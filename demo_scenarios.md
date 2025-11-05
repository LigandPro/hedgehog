# HEDGE Folder Management - Usage Scenarios

## Overview
The new hybrid folder management system combines automatic smart logic with explicit control flags to handle results folders intelligently.

---

## Automatic Behavior (No Flags)

### Scenario 1: First Full Run
```bash
uv run hedge run
```
**Result:** Creates `results/test/` (from config)
**Reason:** Folder doesn't exist yet

---

### Scenario 2: Second Full Run
```bash
uv run hedge run
```
**Result:** Creates `results/test1/` (auto-increment)
**Reason:** `results/test/` already contains results from previous run

---

### Scenario 3: Rerun Specific Stage (Typical Use Case)
```bash
# After full run completed, you want to rerun just docking
uv run hedge run --stage docking
```
**Result:** Uses `results/test/` (reuses existing)
**Reason:** Stage rerun without new molecules → finds artifacts from previous stages

**This fixes the original issue!** ✅

---

### Scenario 4: Run Stage with New Molecules
```bash
uv run hedge run --stage descriptors --mols new_molecules.csv
```
**Result:** Creates `results/test2/` (auto-increment)
**Reason:** New input data → fresh folder to avoid mixing results

---

## Explicit Control Flags

### Scenario 5: Force Reuse Existing Folder
```bash
# Force reuse even for full run
uv run hedge run --reuse
```
**Result:** Uses `results/test/` (reuses)
**Use Case:** You want to overwrite/update results intentionally

---

### Scenario 6: Force New Folder
```bash
# Force new folder even for stage rerun
uv run hedge run --stage docking --force-new
```
**Result:** Creates `results/test3/` (new folder)
**Use Case:** Testing different parameters for same stage

---

### Scenario 7: Stage Rerun with New Molecules + Reuse
```bash
uv run hedge run --stage synthesis --mols new.csv --reuse
```
**Result:** Uses `results/test/` (reuses)
**Use Case:** Override auto-increment, explicitly reuse folder

---

## Priority Rules

The system follows this priority:

1. **Explicit `--reuse`** → Always reuse configured folder
2. **Explicit `--force-new`** → Always create incremented folder
3. **Automatic (no flags):**
   - `--stage` WITHOUT `--mols` → **Reuse** (stage rerun scenario)
   - `--stage` WITH `--mols` → **New folder** (new data)
   - Full run (no `--stage`) → **New folder** (protect existing)

---

## Error Handling

### Conflicting Flags
```bash
uv run hedge run --reuse --force-new
```
**Result:** ❌ Error: "Cannot use --reuse and --force-new together"

---

## Summary Table

| Command | Has `--stage`? | Has `--mols`? | Has Flag? | Result |
|---------|---------------|---------------|-----------|--------|
| `uv run hedge run` | ❌ | ❌ | ❌ | New folder (auto-increment) |
| `uv run hedge run --stage docking` | ✅ | ❌ | ❌ | **Reuse folder** ✅ |
| `uv run hedge run --stage desc --mols x.csv` | ✅ | ✅ | ❌ | New folder |
| `uv run hedge run --reuse` | ❌ | ❌ | ✅ | Reuse folder |
| `uv run hedge run --force-new` | ❌ | ❌ | ✅ | New folder |
| `uv run hedge run --stage dock --force-new` | ✅ | ❌ | ✅ | New folder (override) |

---

## Benefits

✅ **Stage reruns work correctly** - finds existing artifacts automatically
✅ **No accidental overwrites** - full runs get fresh folders
✅ **Explicit control when needed** - flags override automatic behavior
✅ **Intuitive defaults** - "do what I mean" behavior
✅ **Flexible for all workflows** - research, production, testing
