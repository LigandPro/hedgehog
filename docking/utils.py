import os
import json
import uuid
import subprocess
from datetime import datetime
from pathlib import Path

import pandas as pd
import datamol as dm
import shlex

from configs.config_utils import load_config
from logger_config import logger


def _latest_pass_source(base_folder: Path) -> Path | None:
    candidates = [
        base_folder / 'Retrosynthesis' / 'passRetrosynthesisSMILES.csv',
        base_folder / 'StructFilters' / 'passStructFiltersSMILES.csv',
        base_folder / 'Descriptors' / 'passDescriptorsSMILES.csv',
        base_folder / 'sampledMols.csv',
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def _write_smiles_smi_and_csv(df: pd.DataFrame, smi_path: Path | None, csv_path: Path) -> dict:
    smiles_col = next((c for c in df.columns if c.lower() == 'smiles'), None)
    if smiles_col is None:
        raise ValueError('No smiles column in docking input')
    if 'mol_idx' in df.columns:
        names = df['mol_idx'].astype(str)
    elif 'model_name' in df.columns:
        names = df['model_name'].astype(str)
    else:
        names = pd.Series([f'lig_{i}' for i in range(len(df))])

    total = 0
    written = 0
    skipped: list[str] = []
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for smi, name_val, mn, midx in zip(
        df[smiles_col].astype(str),
        names,
        df['model_name'] if 'model_name' in df.columns else pd.Series([None] * len(df)),
        df['mol_idx'] if 'mol_idx' in df.columns else pd.Series([None] * len(df)),
    ):
        total += 1
        try:
            mol = dm.to_mol(smi)
        except Exception:
            mol = None
        if mol is None:
            skipped.append(smi)
            continue
        rows.append({'smiles': smi, 'name': str(name_val), 'model_name': None if pd.isna(mn) else str(mn), 'mol_idx': None if pd.isna(midx) else str(midx)})
        written += 1
    if smi_path is not None:
        with open(smi_path, 'w') as f:
            for r in rows:
                f.write(f"{r['smiles']}\t{r['name']}\n")
    out_df = pd.DataFrame(rows, columns=['smiles','name','model_name','mol_idx'])
    out_df.to_csv(csv_path, index=False)
    if skipped:
        skip_path = csv_path.parent / 'skipped_smiles.txt'
        with open(skip_path, 'w') as sf:
            for s in skipped:
                sf.write(s + '\n')
        logger.warning(f'Some SMILES could not be parsed for docking: {len(skipped)}/{total}. See {skip_path}')
    return { 'smi_path': (str(smi_path) if smi_path is not None else None), 'csv_path': str(csv_path), 'total': total, 'written': written, 'skipped': len(skipped) }


def run_docking(config: dict) -> bool:
    cfg = load_config(config['config_docking'])
    if not cfg.get('run', False):
        logger.info('Docking disabled in config')
        return False

    base_folder = Path(config['folder_to_save']).resolve()
    source = _latest_pass_source(base_folder)
    if source is None:
        logger.warning('No pass*SMILES.csv or sampledMols.csv found for docking input')
        return False

    try:
        df = pd.read_csv(source)
    except Exception as e:
        logger.error(f'Failed to read docking input {source}: {e}')
        return False

    ligands_dir = base_folder / 'Docking'
    ligands_csv = ligands_dir / 'ligands.csv'
    write_stats = _write_smiles_smi_and_csv(df, None, ligands_csv)

    tools_cfg = cfg.get('tools', 'both')
    if isinstance(tools_cfg, str):
        tools_list = [t.strip().lower() for t in tools_cfg.split(',')] if ',' in tools_cfg else [tools_cfg.strip().lower()]
    elif isinstance(tools_cfg, (list, tuple)):
        tools_list = [str(t).strip().lower() for t in tools_cfg]
    else:
        tools_list = ['both']
    if 'both' in tools_list or not tools_list:
        tools_list = ['smina', 'gnina']

    ligands_dir.mkdir(parents=True, exist_ok=True)

    prepared_any = False
    job_id = f"dock_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    scripts_prepared: list[str] = []
    job_ids: dict[str, str] = {}

    if 'smina' in tools_list:
        ini_src = Path(cfg.get('smina_ini', Path(__file__).resolve().parent / 'smina.ini'))
        ini_dst = ligands_dir / 'smina_auto.ini'
        receptor = cfg.get('receptor_pdb')
        lines = []
        inline_ini = cfg.get('smina_ini_content')
        if inline_ini and isinstance(inline_ini, str):
            lines = inline_ini.splitlines()
        elif ini_src.exists():
            try:
                with open(ini_src, 'r') as f:
                    lines = f.read().splitlines()
            except Exception:
                lines = []
        if not lines:
            lines = [
                'output-dir = smina_results',
                'screen-type = smina',
                'metadata-template = {"software": "smina"}',
            ]
        def _set_kv(lines_in: list[str], key: str, value: str) -> list[str]:
            out = []
            found = False
            for ln in lines_in:
                if ln.strip().startswith(f'{key} '):
                    out.append(f'{key} = {value}')
                    found = True
                else:
                    out.append(ln)
            if not found:
                out.append(f'{key} = {value}')
            return out
        if not any(ln.strip().startswith('screen-type') for ln in lines):
            lines.insert(0, 'screen-type = smina')

        lines = _set_kv(lines, 'input-files', f'[{ligands_csv.name}]')
        lines = _set_kv(lines, 'smiles-col', '0')
        lines = _set_kv(lines, 'name-col', '1')
        if receptor:
            lines = _set_kv(lines, 'receptors', f'[{receptor}]')
        with open(ini_dst, 'w') as f:
            f.write('\n'.join(lines) + '\n')

        activate = cfg.get('pyscreener_activate')
        runfile = ligands_dir / 'run_smina.sh'
        with open(runfile, 'w') as f:
            f.write('#!/usr/bin/env bash\nset -eo pipefail\n')
            if activate:
                f.write(f'{activate}\n')
            f.write(f'pyscreener --config {ini_dst.name}\n')
        os.chmod(runfile, 0o755)
        prepared_any = True
        scripts_prepared.append(str(runfile))
        job_ids['smina'] = f"smina_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"

    if 'gnina' in tools_list:
        gnina_bin = cfg.get('gnina_bin', 'gnina')
        receptor = cfg.get('receptor_pdb')
        if not receptor or not Path(receptor).exists():
            logger.error('gnina: receptor_pdb is missing or not found')
            return prepared_any
        cfg_out_dir = cfg.get('gnina_output_dir')
        if cfg_out_dir:
            out_dir_candidate = Path(cfg_out_dir)
            gnina_dir = out_dir_candidate if out_dir_candidate.is_absolute() else (base_folder / out_dir_candidate)
        else:
            gnina_dir = base_folder / 'Docking' / 'gnina_results'
        gnina_dir.mkdir(parents=True, exist_ok=True)
        out_sdf = (gnina_dir / 'gnina_out.sdf').resolve()
        ligands_arg = cfg.get('gnina_ligands')
        ligands_val = str(ligands_arg) if ligands_arg else ligands_csv.name

        needs_conversion = not ligands_val.lower().endswith(('.sdf', '.sdf.gz', '.osd', '.mol2'))
        ligprep_path = cfg.get('ligprep_path')
        ligprep_cmd_line = None
        if needs_conversion:
            if ligprep_path:
                ligands_val_final = 'ligands_prepared.osd'
                ligprep_cmd_line = f"{ligprep_path} -ismi {ligands_val} -osd {ligands_val_final}"
            else:
                ligands_val_final = 'ligands_prepared.sdf'
                sdf_path = ligands_dir / ligands_val_final
                try:
                    import rdkit
                    from rdkit import Chem
                    from rdkit.Chem import AllChem
                    df_ids = pd.read_csv(ligands_csv)
                    smiles_series = df_ids['smiles'] if 'smiles' in df_ids.columns else pd.Series([])
                    name_series = df_ids['name'] if 'name' in df_ids.columns else pd.Series([f'lig_{i}' for i in range(len(smiles_series))])
                    writer = Chem.SDWriter(str(sdf_path))
                    written_local = 0
                    for smi, nm in zip(smiles_series.astype(str), name_series.astype(str)):
                        try:
                            mol = Chem.MolFromSmiles(smi)
                            if mol is None:
                                continue
                            mol = Chem.AddHs(mol)
                            try:
                                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                                AllChem.UFFOptimizeMolecule(mol)
                            except Exception:
                                pass
                            mol.SetProp('_Name', nm)
                            writer.write(mol)
                            written_local += 1
                        except Exception:
                            continue
                    writer.close()
                    if written_local == 0:
                        logger.error('RDKit fallback produced 0 molecules for GNINA SDF.')
                        return prepared_any
                except Exception as e:
                    logger.error(f'RDKit fallback SDF conversion failed: {e}')
                    return prepared_any
        else:
            ligands_val_final = ligands_val

        user_path = cfg.get('gnina_path')
        user_args = cfg.get('gnina_args')
        if user_path and user_args:
            formatted = user_args.format(receptor=receptor, ligands=ligands_val_final, out=str(out_sdf))
            cmd = [str(user_path)] + shlex.split(formatted)
        else:
            cmd = [gnina_bin, '-r', receptor, '-l', ligands_val_final, '-o', str(out_sdf)]
        center = cfg.get('center')
        size = cfg.get('size')
        if isinstance(center, (list, tuple)) and len(center) == 3:
            cmd += ['--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2])]
        if isinstance(size, (list, tuple)) and len(size) == 3:
            cmd += ['--size_x', str(size[0]), '--size_y', str(size[1]), '--size_z', str(size[2])]
        if 'exhaustiveness' in cfg:
            cmd += ['--exhaustiveness', str(cfg['exhaustiveness'])]
        if 'num_modes' in cfg:
            cmd += ['--num_modes', str(cfg['num_modes'])]

        if cfg.get('gnina_autobox_ligand'):
            cmd += ['--autobox_ligand', str(cfg['gnina_autobox_ligand'])]
        if cfg.get('gnina_autobox_add') is not None:
            cmd += ['--autobox_add', str(cfg['gnina_autobox_add'])]
        if cfg.get('gnina_cpu') is not None:
            cmd += ['--cpu', str(cfg['gnina_cpu'])]
        if cfg.get('gnina_seed') is not None:
            cmd += ['--seed', str(cfg['gnina_seed'])]
        if cfg.get('gnina_extra_args'):
            extra = str(cfg['gnina_extra_args']).strip()
            if extra:
                cmd += extra.split()

        runfile = ligands_dir / 'run_gnina.sh'
        gnina_activate = cfg.get('gnina_activate')
        if not gnina_activate:
            env_path = cfg.get('gnina_env_path')
            if env_path:
                conda_sh = cfg.get('conda_sh', os.path.expanduser('~/miniconda3/etc/profile.d/conda.sh'))
                gnina_activate = f'source {conda_sh} && conda activate {env_path}'
        with open(runfile, 'w') as f:
            f.write('#!/usr/bin/env bash\nset -eo pipefail\n')
            if gnina_activate:
                f.write(f'{gnina_activate}\n')
            if ligprep_cmd_line:
                f.write(f'{ligprep_cmd_line}\n')
                f.write('rm -f ligands_raw.smi || true\n')
            f.write(f'mkdir -p "$(dirname "{str(out_sdf)}")"\n')
            f.write(f': > "{str(out_sdf)}"\n')
            f.write(' '.join(cmd) + '\n')
        os.chmod(runfile, 0o755)
        prepared_any = True
        scripts_prepared.append(str(runfile))
        job_ids['gnina'] = f"gnina_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"

    try:
        ligands_dir.mkdir(parents=True, exist_ok=True)
        meta = {
            'job_id': job_id,
            'timestamp': datetime.now().isoformat(timespec='seconds'),
            'source_file': str(source),
            'num_ligands': int(len(df)),
            'receptor_pdb': cfg.get('receptor_pdb'),
            'tools_prepared': [t for t in tools_list if t in ['smina','gnina']],
            'scripts': scripts_prepared,
            'ligands_csv': str(ligands_csv),
            'ligands_counts': write_stats,
            'jobs': {
                k: {
                    'name': k,
                    'job_id': v,
                    'script': str(ligands_dir / f"run_{k}.sh")
                } for k, v in job_ids.items()
            }
        }
        with open(ligands_dir / 'job_meta.json', 'w') as f:
            json.dump(meta, f, indent=2)
        ids_path = ligands_dir / 'job_ids.txt'
        try:
            with open(ids_path, 'w') as jf:
                jf.write(f"overall: {job_id}\n")
                jf.write(f"smina: {job_ids.get('smina','')}\n")
                jf.write(f"gnina: {job_ids.get('gnina','')}\n")
        except Exception:
            pass
        logger.info(f"Docking job id: {job_id}")
    except Exception:
        pass

    final_ok = prepared_any
    try:
        if prepared_any and bool(cfg.get('auto_run', True)):
            background = bool(cfg.get('run_in_background', False))
            run_status = {}

            if 'smina' in tools_list:
                smina_script = ligands_dir / 'run_smina.sh'
                if smina_script.exists():
                    log_path = ligands_dir / 'smina_run.log'
                    if background:
                        with open(log_path, 'ab') as logf:
                            subprocess.Popen(['./run_smina.sh'], stdout=logf, stderr=logf, cwd=str(ligands_dir))
                        run_status['smina'] = 'started_background'
                        if 'smina' in job_ids:
                            run_status['smina_job_id'] = job_ids['smina']
                    else:
                        with open(log_path, 'wb') as logf:
                            subprocess.run(['./run_smina.sh'], check=True, stdout=logf, stderr=logf, cwd=str(ligands_dir))
                        run_status['smina'] = 'completed'
                        if 'smina' in job_ids:
                            run_status['smina_job_id'] = job_ids['smina']
                    smina_out_dir = ligands_dir / 'smina_results'
                    run_status['smina_results_dir'] = str(smina_out_dir) if smina_out_dir.exists() else None

            if 'gnina' in tools_list:
                gnina_script = ligands_dir / 'run_gnina.sh'
                if gnina_script.exists():
                    log_path = ligands_dir / 'gnina_run.log'
                    if background:
                        with open(log_path, 'ab') as logf:
                            subprocess.Popen(['./run_gnina.sh'], stdout=logf, stderr=logf, cwd=str(ligands_dir))
                        run_status['gnina'] = 'started_background'
                        if 'gnina' in job_ids:
                            run_status['gnina_job_id'] = job_ids['gnina']
                    else:
                        with open(log_path, 'wb') as logf:
                            subprocess.run(['./run_gnina.sh'], check=True, stdout=logf, stderr=logf, cwd=str(ligands_dir))
                        run_status['gnina'] = 'completed'
                        if 'gnina' in job_ids:
                            run_status['gnina_job_id'] = job_ids['gnina']
                    cfg_out_dir = cfg.get('gnina_output_dir')
                    if cfg_out_dir:
                        out_dir_candidate = Path(cfg_out_dir)
                        gnina_dir_resolved = out_dir_candidate if out_dir_candidate.is_absolute() else (base_folder / out_dir_candidate)
                    else:
                        gnina_dir_resolved = base_folder / 'Docking' / 'gnina_results'
                    run_status['gnina_output'] = str((gnina_dir_resolved / 'gnina_out.sdf'))
                    run_status['gnina_log'] = str(log_path)

            try:
                meta_path = ligands_dir / 'job_meta.json'
                meta_obj = {}
                if meta_path.exists():
                    with open(meta_path, 'r') as f:
                        meta_obj = json.load(f)
                meta_obj['run_status'] = run_status
                with open(meta_path, 'w') as f:
                    json.dump(meta_obj, f, indent=2)
            except Exception:
                pass
            if not background:
                selected = [t for t in tools_list if t in ['smina','gnina']]
                all_completed = all(run_status.get(t) == 'completed' for t in selected)
                final_ok = all_completed
    except Exception as e:
        logger.error(f'Docking auto-run failed: {e}')
        final_ok = False

    return final_ok
