# FeAlOY ODS Nanocomposite RVE Pipeline

End-to-end pipeline for generating polycrystal RVE meshes (Neper), preparing
them for MOOSE (sideset/interface splitting), and producing ready-to-run
CZM simulation inputs.

## Pipeline overview

```
generate_odsfe_rve.sh          (neper env)
  -> Neper tessellation (-T) + meshing (-M)
  -> rve_<config>.msh, .tess, .inp, .vtu

make_mesh_prep.py              (neper env)
  -> reads actual mesh bounds from .msh (no hardcoded thresholds)
  -> writes mesh_prep_<config>.i

mesh_prep_<config>.i run via --mesh-only   (moose env)
  -> ParsedGenerateSideset x6 (left/right/front/back/bottom/top)
  -> BreakMeshByBlockGenerator (split_interface=true)
  -> rve_<config>_split.e

extract_czm_boundaries.py      (moose env)
  -> reads split Exodus file, finds BlockX_BlockY interface sidesets
  -> writes rve_<config>_split_interfaces.txt

generate_odsfe_inputs.py       (moose env)
  -> reads interfaces list + mesh bounds
  -> writes config_<config>.i  (ready to run: bulk elasticity, CZM material,
     BCs, postprocessors, executioner)
```

`run_pipeline.sh` chains all of this into one command, switching conda
environments automatically at the right points.

All four generator/prep scripts write their outputs into a single
per-config directory (`rve_<config>_n<n_grains>/`), not the pipeline root.

## Directory layout after a run

```
rve_needles_x_n200/
  rve_needles_x.tess
  rve_needles_x.msh
  rve_needles_x.inp
  rve_needles_x.vtu
  rve_needles_x.stcell
  mesh_prep_needles_x.i
  rve_needles_x_split.e
  rve_needles_x_split_interfaces.txt
  config_needles_x.i          <- final, ready-to-run MOOSE input
  config_needles_x_csv.csv    <- produced when you actually run it
  config_needles_x_out.e
```

## Environment setup

Two separate conda environments are needed -- Neper and MOOSE have
incompatible/unrelated dependency stacks, so keep them apart rather than
trying to merge into one environment.

```bash
conda create -n neper python=3.11
conda activate neper
pip install -r requirements-neper.txt
# + Neper itself (not a pip package -- see https://neper.info/ for install
#   instructions; this pipeline was built against Neper 4.10.2)
conda deactivate

conda create -n moose python=3.11
conda activate moose
pip install -r requirements-moose.txt
# + MOOSE itself and your aragonite-opt executable (built separately, not
#   a pip package)
conda deactivate
```

Set `NEPER_ENV`/`MOOSE_ENV` env vars if your actual conda environment names
differ from `neper`/`moose` (see `run_pipeline.sh` for all overridable
env vars: `ARAGONITE_OPT`, `MPI_RANKS`, `CONDA_BASE`, `NEPER_ENV`,
`MOOSE_ENV`).

**Why mesh-handling libraries live in the `moose` env, not split across
both:** `mesh_prep_<config>.i` runs in `moose` env (it needs the actual
MOOSE binary), and `extract_czm_boundaries.py` + `generate_odsfe_inputs.py`
also run there in the current pipeline -- so `netCDF4` (Exodus reading) and
`meshio` (mesh bounds reading) both need to be available in `moose` env.
`meshio` is also needed in `neper` env, since `generate_odsfe_rve.sh`'s own
`.msh -> .vtu` conversion step and `make_mesh_prep.py` both currently run
there, before the environment switch.

## Usage

### One config, full pipeline, auto-sized domain

```bash
./run_pipeline.sh needles_x 200 1.0 3
```
Arguments: `<config_name> [n_grains=200] [grain_diameq_mm=1.0] [min_grains_across=5]`

Domain size is derived automatically so the domain is at least
`min_grains_across` times a grain's own extent along each axis (computed
per-axis from the config's `aspratio`, not forced into a cube -- see
`generate_odsfe_rve.sh` header comments for the full derivation).

### One config, fixed cross-section ("direct mode")

```bash
./run_pipeline.sh needles_x 200 0.7 5 4.5
```
Adds a 5th argument: `cross_section_mm`. When set, the short (non-elongated)
axes are fixed to exactly this value, and the long axis is derived from
`min_grains_across x` that axis's own grain extent -- **`min_grains_across`
now means "grains along the long axis" directly, not a floor to check
afterward.** `n_grains` is ignored in this mode; actual grain count is
whatever falls out of the resulting box volume, printed at generation time.

Config names: `isotropic`, `needles_x`, `needles_y`, `pancake_x`, `pancake_z`.

### Reproducible tessellation

Neper's tessellation seed is random per run by default (`SEED=$RANDOM`),
so identical arguments can still produce a different mesh -- and,
occasionally, a mesh that hits a 3D-meshing failure that a different seed
would have avoided (see Troubleshooting). Pin it if you need a repeatable
result:
```bash
SEED=12345 ./run_pipeline.sh needles_x 200 1.0 3
```

### Running the final simulation

`config_<config>.i` is a complete MOOSE input once generated -- run it
directly:
```bash
cd rve_needles_x_n200
conda activate moose
mpirun -n 16 ~/projects/aragonite/aragonite-opt -i config_needles_x.i
```

## Troubleshooting -- issues already hit and their fixes

- **`neper -M` fails with `Illegal hash-position` / `BFGS update error2` /
  `Meshing of poly N failed`**: a degenerate/sliver cell in the
  tessellation. `-reg 1` (already in `generate_odsfe_rve.sh`) fixes most
  cases; if it still fails, try a different `SEED` (see above) -- this is
  partly luck-of-the-draw on which random tessellation you get.

- **MOOSE crashes with `ERROR: You may not set a boundary ID of -123`**:
  raw Gmsh-format read issue, seen at very high sideset/nodeset counts
  (tens of thousands). If you're generating a mesh with a very large
  number of grains, this and the next item are the ones to watch for.

- **`Error writing sidesets` (`exodusII_io_helper.C`)**: classic Exodus
  (netCDF-3-based) format has a hard ceiling on total sidesets+nodesets in
  one file. Confirmed hit at ~28,000 sidesets + ~36,000 nodesets. Fix:
  `write_hdf5 = true` in the `[Outputs]/[Exodus]` block (already set in
  `make_mesh_prep.py`'s generated `.i`), which switches to the HDF5-backed
  netCDF-4 Exodus format and removes the ceiling.

- **All nodeset IDs come back as `-1` in `ncdump`**: seen even after
  `write_hdf5 = true`, at very large sideset counts (~20,000+). Root cause
  not fully confirmed -- looked like an uninitialized-value bug in the
  auto-generated sideset-companion-nodeset writer at scale, not a Neper
  tagging issue. The practical fix that's actually been used: keep grain
  counts in the few-hundred range (via `min_grains_across=2` or `3`, or
  fixed cross-section/direct mode) rather than chasing this further into
  MOOSE/libMesh internals.

- **`set -euo pipefail` + `conda activate` -> `CONDA_BUILD: unbound
  variable`**: conda's own activation scripts aren't `nounset`-safe.
  `run_pipeline.sh` already wraps `conda activate` calls in `set +u` / `set
  -u`.

- **Cube-shaped domains blow up grain count for elongated configs**:
  a cube forces the short axes to match whatever the long axis needs.
  Fixed by the anisotropic-box logic in `generate_odsfe_rve.sh` (each axis
  sized independently) -- see script header comments for the full
  before/after numbers.

## Known open items / unverified assumptions

- `extract_czm_boundaries.py`'s exact output filename convention was
  inferred from the calling pattern in earlier pipeline commands, not
  confirmed against its actual source -- verify `rve_<config>_split_
  interfaces.txt` is really what it writes before relying on this
  unattended.
- Whether MOOSE resolves a `.i` file's internal `[Mesh]/file` reference
  relative to the `.i`'s own location or the process's working directory
  was never directly confirmed. `run_pipeline.sh` sidesteps this by `cd`ing
  into the config's own output directory before running MOOSE, which works
  regardless of which convention is correct.
- `write_hdf5` is confirmed to exist and work via the `[Outputs]/[Exodus]`
  block. Whether it's honored by a bare `--mesh-only` invocation
  specifically (rather than a full `Outputs`-block-driven run) was not
  separately confirmed -- `mesh_prep_<config>.i` includes it either way as
  free insurance.
