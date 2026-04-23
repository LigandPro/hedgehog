"""SYNC synthesizability classifier inference.

This module implements the inference path for the 3D SYNC classifier described
in Liu et al., "SYNC: Measuring and Advancing Synthesizability in
Structure-Based Drug Design" (ICLR 2026), using the public checkpoint from
https://github.com/XYxiyang/SYNC.
"""

from __future__ import annotations

import contextlib
import re
import sys
import types
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from hedgehog.configs.logger import logger

_SYNC_CACHE: dict[str, Any] = {"model": None, "torch": None, "load_failed": False}
_VOCAB_PATH = Path(__file__).resolve().parent / "data" / "sync_smiles_vocab_with_aro.txt"


def _load_torch():
    if _SYNC_CACHE["torch"] is not None:
        return _SYNC_CACHE["torch"]
    try:
        import torch
        from torch import nn
    except Exception as exc:
        logger.warning("SYNC score unavailable: PyTorch could not be imported (%s)", exc)
        _SYNC_CACHE["load_failed"] = True
        return None
    _SYNC_CACHE["torch"] = (torch, nn)
    return _SYNC_CACHE["torch"]


def _read_vocab() -> tuple[list[str], dict[str, int], re.Pattern[str]]:
    vocab = _VOCAB_PATH.read_text(encoding="utf-8").splitlines()
    token_to_id = {token: i for i, token in enumerate(vocab, start=1)}
    pattern = re.compile(
        "(" + "|".join(map(re.escape, sorted(vocab, reverse=True))) + ")"
    )
    return vocab, token_to_id, pattern


_SYNC_VOCAB, _SYNC_TOKEN_TO_ID, _SYNC_TOKEN_PATTERN = _read_vocab()


def _tokenize_atom_symbols(mol: Chem.Mol) -> list[int]:
    atom_string = ""
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            atom_string += f"{atom.GetSymbol()}pi"
        else:
            atom_string += atom.GetSymbol()
    return [_SYNC_TOKEN_TO_ID[token] for token in _SYNC_TOKEN_PATTERN.findall(atom_string)]


def _bond_feature_simple(bond: Chem.rdchem.Bond) -> int:
    bond_type = bond.GetBondType()
    if bond_type == Chem.rdchem.BondType.SINGLE:
        return 1
    if bond_type == Chem.rdchem.BondType.DOUBLE:
        return 2
    if bond_type == Chem.rdchem.BondType.TRIPLE:
        return 3
    if bond_type == Chem.rdchem.BondType.AROMATIC:
        return 4
    return 5


def _embed_smiles(smiles: str, seed: int) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        logger.debug("UFF optimization failed for SYNC input: %s", smiles)
    return Chem.RemoveHs(mol)


def _mol_to_graph_tensors(mol: Chem.Mol, torch):
    atom_tokens = _tokenize_atom_symbols(mol)
    if len(atom_tokens) != mol.GetNumAtoms() or mol.GetNumAtoms() == 0:
        return None

    edge_index: list[list[int]] = []
    edge_attr: list[int] = []
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        feature = _bond_feature_simple(bond)
        edge_index.extend(([begin, end], [end, begin]))
        edge_attr.extend((feature, feature))

    if not edge_index:
        return None

    conformer = mol.GetConformer()
    positions = [
        [
            conformer.GetAtomPosition(i).x,
            conformer.GetAtomPosition(i).y,
            conformer.GetAtomPosition(i).z,
        ]
        for i in range(mol.GetNumAtoms())
    ]
    return {
        "x": torch.tensor(atom_tokens, dtype=torch.long),
        "edge_index": torch.tensor(edge_index, dtype=torch.long).t().contiguous(),
        "edge_attr": torch.tensor(edge_attr, dtype=torch.long),
        "positions": torch.tensor(positions, dtype=torch.float),
    }


def _unsorted_segment_sum(data, segment_ids, num_segments):
    result = data.new_full((num_segments, data.size(1)), 0)
    segment_ids = segment_ids.unsqueeze(-1).expand(-1, data.size(1))
    result.scatter_add_(0, segment_ids, data)
    return result


def _unsorted_segment_mean(data, segment_ids, num_segments):
    result = data.new_full((num_segments, data.size(1)), 0)
    count = data.new_full((num_segments, data.size(1)), 0)
    segment_ids = segment_ids.unsqueeze(-1).expand(-1, data.size(1))
    result.scatter_add_(0, segment_ids, data)
    count.scatter_add_(0, segment_ids, data.new_ones(data.shape))
    return result / count.clamp(min=1)


def _model_classes():
    loaded = _load_torch()
    if loaded is None:
        return None
    torch, nn = loaded

    class E_GCL(nn.Module):
        def __init__(
            self,
            input_nf,
            output_nf,
            hidden_nf,
            edges_in_d=0,
            act_fn=nn.SiLU(),
            residual=True,
            attention=False,
            normalize=False,
            coords_agg="mean",
            tanh=False,
        ):
            super().__init__()
            self.residual = residual
            self.attention = attention
            self.normalize = normalize
            self.coords_agg = coords_agg
            self.tanh = tanh
            self.epsilon = 1e-8
            self.edge_mlp = nn.Sequential(
                nn.Linear(input_nf * 2 + 1 + edges_in_d, hidden_nf),
                act_fn,
                nn.Linear(hidden_nf, hidden_nf),
                act_fn,
            )
            self.node_mlp = nn.Sequential(
                nn.Linear(hidden_nf + input_nf, hidden_nf),
                act_fn,
                nn.Linear(hidden_nf, output_nf),
            )
            layer = nn.Linear(hidden_nf, 1, bias=False)
            torch.nn.init.xavier_uniform_(layer.weight, gain=0.001)
            coord_layers = [nn.Linear(hidden_nf, hidden_nf), act_fn, layer]
            if self.tanh:
                coord_layers.append(nn.Tanh())
            self.coord_mlp = nn.Sequential(*coord_layers)
            if self.attention:
                self.att_mlp = nn.Sequential(nn.Linear(hidden_nf, 1), nn.Sigmoid())

        def edge_model(self, source, target, radial, edge_attr):
            if edge_attr is None:
                out = torch.cat([source, target, radial], dim=1)
            else:
                out = torch.cat([source, target, radial, edge_attr], dim=1)
            out = self.edge_mlp(out)
            if self.attention:
                out = out * self.att_mlp(out)
            return out

        def node_model(self, x, edge_index, edge_attr, node_attr):
            row, _ = edge_index
            agg = _unsorted_segment_sum(edge_attr, row, num_segments=x.size(0))
            if node_attr is not None:
                agg = torch.cat([x, agg, node_attr], dim=1)
            else:
                agg = torch.cat([x, agg], dim=1)
            out = self.node_mlp(agg)
            if self.residual:
                out = x + out
            return out, agg

        def coord_model(self, coord, edge_index, coord_diff, edge_feat):
            row, _ = edge_index
            trans = coord_diff * self.coord_mlp(edge_feat)
            if self.coords_agg == "sum":
                agg = _unsorted_segment_sum(trans, row, num_segments=coord.size(0))
            else:
                agg = _unsorted_segment_mean(trans, row, num_segments=coord.size(0))
            return coord + agg

        def coord2radial(self, edge_index, coord):
            row, col = edge_index
            coord_diff = coord[row] - coord[col]
            radial = torch.sum(coord_diff**2, 1).unsqueeze(1)
            if self.normalize:
                coord_diff = coord_diff / (torch.sqrt(radial).detach() + self.epsilon)
            return radial, coord_diff

        def forward(self, h, edge_index, coord, edge_attr=None, node_attr=None):
            radial, coord_diff = self.coord2radial(edge_index, coord)
            row, col = edge_index
            edge_feat = self.edge_model(h[row], h[col], radial, edge_attr)
            coord = self.coord_model(coord, edge_index, coord_diff, edge_feat)
            h, _ = self.node_model(h, edge_index, edge_feat, node_attr)
            return h, coord, edge_attr

    class EGNN(nn.Module):
        def __init__(
            self,
            in_node_nf,
            hidden_nf,
            out_node_nf,
            in_edge_nf=0,
            device="cpu",
            act_fn=nn.SiLU(),
            n_layers=4,
            residual=True,
            attention=False,
            normalize=False,
            tanh=False,
        ):
            super().__init__()
            self.hidden_nf = hidden_nf
            self.device = device
            self.n_layers = n_layers
            self.embedding_in = nn.Linear(in_node_nf, hidden_nf)
            self.embedding_out = nn.Linear(hidden_nf, out_node_nf)
            for i in range(n_layers):
                self.add_module(
                    f"gcl_{i}",
                    E_GCL(
                        hidden_nf,
                        hidden_nf,
                        hidden_nf,
                        edges_in_d=in_edge_nf,
                        act_fn=act_fn,
                        residual=residual,
                        attention=attention,
                        normalize=normalize,
                        tanh=tanh,
                    ),
                )

        def forward(self, h, x, edges, edge_attr):
            h = self.embedding_in(h)
            for i in range(self.n_layers):
                h, x, _ = self._modules[f"gcl_{i}"](h, edges, x, edge_attr=edge_attr)
            return self.embedding_out(h), x

    class ClassifierHead(nn.Module):
        def __init__(self, d_model: int, num_classes: int, dim_hidden: int):
            super().__init__()
            self.mlp = nn.Sequential(
                nn.Linear(d_model, dim_hidden),
                nn.ReLU(),
                nn.Linear(dim_hidden, dim_hidden),
                nn.ReLU(),
                nn.Linear(dim_hidden, num_classes),
            )

        def forward(self, h):
            return self.mlp(h)

    class Sync3DClassifier(nn.Module):
        def __init__(self):
            super().__init__()
            self.atom_encoder = nn.Embedding(158, 128, padding_idx=0)
            self.edge_encoder = nn.Embedding(6, 128, padding_idx=0)
            self.encoder = EGNN(
                in_node_nf=128,
                hidden_nf=128,
                out_node_nf=128,
                in_edge_nf=128,
                device="cpu",
                n_layers=4,
            )
            self.head = ClassifierHead(d_model=128, num_classes=2, dim_hidden=256)

        def forward(self, batch):
            node = self.atom_encoder(batch["x"]).squeeze()
            edge = self.edge_encoder(batch["edge_attr"]).squeeze()
            node = node.reshape(-1, node.shape[-1])
            node, _ = self.encoder(
                node, batch["positions"], batch["edge_index"], edge
            )
            pooled = node.sum(0, keepdim=True) / node.shape[0]
            return self.head(pooled)

    return torch, E_GCL, EGNN, ClassifierHead, Sync3DClassifier


@contextlib.contextmanager
def _sync_checkpoint_class_aliases():
    classes = _model_classes()
    if classes is None:
        yield
        return
    _, e_gcl, egnn, classifier_head, sync_classifier = classes
    aliases = {
        "synsbdd": types.ModuleType("synsbdd"),
        "synsbdd.models": types.ModuleType("synsbdd.models"),
        "synsbdd.models.syn3d_lin": types.ModuleType("synsbdd.models.syn3d_lin"),
        "synsbdd.models.egnn": types.ModuleType("synsbdd.models.egnn"),
        "synsbdd.models.egnn.egnn": types.ModuleType("synsbdd.models.egnn.egnn"),
        "synsbdd.models.output_head": types.ModuleType("synsbdd.models.output_head"),
    }
    aliases["synsbdd.models.syn3d_lin"].Syn3dlin = sync_classifier
    aliases["synsbdd.models.egnn.egnn"].E_GCL = e_gcl
    aliases["synsbdd.models.egnn.egnn"].EGNN = egnn
    aliases["synsbdd.models.output_head"].ClassifierHead = classifier_head
    previous = {name: sys.modules.get(name) for name in aliases}
    try:
        sys.modules.update(aliases)
        yield
    finally:
        for name, old_module in previous.items():
            if old_module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = old_module


def _project_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _resolve_sync_model_path(config: dict[str, Any] | None) -> Path | None:
    config = config or {}
    custom_path = config.get("sync_model_path")
    if custom_path:
        path = Path(custom_path).expanduser()
        if not path.is_absolute():
            path = _project_root() / path
        if path.exists():
            return path
        logger.warning("SYNC model checkpoint not found: %s", path)
        return None

    path = _project_root() / "modules" / "sync" / "classifier_emb.ckpt"
    if path.exists():
        return path
    if not config.get("sync_auto_install", True):
        logger.warning("SYNC model checkpoint not found: %s", path)
        return None

    try:
        from hedgehog.setup import ensure_sync_model

        return ensure_sync_model(_project_root())
    except Exception as exc:
        logger.warning("SYNC model checkpoint is unavailable (%s)", exc)
        return None


def _load_checkpoint_state_dict(path: Path):
    loaded = _load_torch()
    if loaded is None:
        return None
    torch, _ = loaded
    try:
        payload = torch.load(path, map_location="cpu", weights_only=True)
    except Exception:
        with _sync_checkpoint_class_aliases():
            payload = torch.load(path, map_location="cpu", weights_only=False)

    if isinstance(payload, dict) and "state_dict" in payload:
        return payload["state_dict"]
    if isinstance(payload, dict):
        return payload
    if hasattr(payload, "state_dict"):
        return payload.state_dict()
    raise RuntimeError(f"Unsupported SYNC checkpoint payload: {type(payload)!r}")


def _resolve_device(torch, config: dict[str, Any] | None):
    requested = str((config or {}).get("sync_device", "cpu"))
    if requested.startswith("cuda") and not torch.cuda.is_available():
        logger.warning("SYNC requested %s but CUDA is unavailable; using CPU", requested)
        requested = "cpu"
    return torch.device(requested)


def _load_sync_model(config: dict[str, Any] | None = None):
    if _SYNC_CACHE["model"] is not None:
        return _SYNC_CACHE["model"]
    if _SYNC_CACHE["load_failed"]:
        return False

    classes = _model_classes()
    if classes is None:
        return False
    torch, _, _, _, sync_classifier = classes
    model_path = _resolve_sync_model_path(config)
    if model_path is None:
        _SYNC_CACHE["load_failed"] = True
        return False

    try:
        state_dict = _load_checkpoint_state_dict(model_path)
        model = sync_classifier()
        model.load_state_dict(state_dict, strict=True)
        model.to(_resolve_device(torch, config))
        model.eval()
        _SYNC_CACHE["model"] = model
        logger.info("SYNC model loaded from %s", model_path)
        return model
    except Exception as exc:
        logger.warning("Failed to load SYNC model from %s: %s", model_path, exc)
        _SYNC_CACHE["load_failed"] = True
        return False


def calculate_sync_scores_batch(
    smiles_list: list[str],
    config: dict[str, Any] | None = None,
    progress_cb=None,
) -> list[float]:
    """Calculate SYNC synthesizability probabilities for SMILES."""
    total = len(smiles_list)
    nan_list = [np.nan] * total
    if total == 0:
        return nan_list
    if progress_cb is not None:
        progress_cb(0, total)

    loaded = _load_torch()
    model = _load_sync_model(config)
    if loaded is None or not model:
        if progress_cb is not None:
            progress_cb(total, total)
        return nan_list
    torch, _ = loaded
    seed = int((config or {}).get("sync_conformer_seed", 61453))
    scores: list[float] = []

    with torch.no_grad():
        for index, smiles in enumerate(smiles_list, start=1):
            score = np.nan
            try:
                mol = _embed_smiles(smiles, seed)
                if mol is not None:
                    graph = _mol_to_graph_tensors(mol, torch)
                    if graph is not None:
                        device = next(model.parameters()).device
                        graph = {key: value.to(device) for key, value in graph.items()}
                        logits = model(graph)
                        score = float(torch.softmax(logits, dim=-1)[0, 1].item())
            except Exception as exc:
                logger.debug("SYNC score failed for %s: %s", smiles, exc)
            scores.append(score)
            if progress_cb is not None:
                progress_cb(index, total)
    return scores
