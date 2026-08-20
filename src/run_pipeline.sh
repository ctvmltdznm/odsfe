#!/usr/bin/env bash
# ============================================================================
# run_pipeline.sh
#
# End-to-end: mesh generation -> mesh_prep.i (auto-built from actual mesh
# bounds, no hardcoded thresholds) -> mesh-only split -> CZM boundary
# extraction -> final MOOSE input generation.
#
# Usage:
#   ./run_pipeline.sh <config_name> [n_grains] [grain_diameq_mm] [min_grains_across] [cross_section_mm]
#   (same arguments as generate_odsfe_rve.sh -- passed straight through)
#
# Example:
#   ./run_pipeline.sh needles_x 200 1.0 3
#   ./run_pipeline.sh needles_x 200 0.7 5 4.5   # direct mode, fixed cross-section
# ============================================================================
set -euo pipefail

CONFIG=${1:?"Usage: $0 <config_name> [n_grains=200] [grain_diameq_mm=1.0] [min_grains_across=5] [cross_section_mm]"}
N_GRAINS=${2:-200}
GRAIN_DIAMEQ=${3:-1.0}
MIN_GRAINS_ACROSS=${4:-5}
CROSS_SECTION_MM=${5:-0}

ARAGONITE_OPT="${ARAGONITE_OPT:-$HOME/projects/aragonite/aragonite-opt}"
MPI_RANKS="${MPI_RANKS:-16}"
NEPER_ENV="${NEPER_ENV:-neper}"
MOOSE_ENV="${MOOSE_ENV:-moose}"

# Conda's `activate` is a shell function injected by `conda init`, not
# available by default in a script's subshell -- source the hook first.
# CONDA_BASE detection: prefer `conda info --base` if conda is already on
# PATH; override via CONDA_BASE env var if that fails in your environment.
CONDA_BASE="${CONDA_BASE:-$(conda info --base 2>/dev/null)}"
if [ -z "$CONDA_BASE" ]; then
    echo "ERROR: could not determine conda base (conda info --base failed)."
    echo "Set CONDA_BASE explicitly, e.g. CONDA_BASE=/home/you/miniconda3"
    exit 1
fi
# conda's own scripts (both the activation hook and per-env activate.d/
# scripts, e.g. activate_zzz_moose-mpi-base.sh) reference variables like
# CONDA_BUILD without a default, which is fine normally but fatal under
# `set -u` above. Relax -u only around conda sourcing/activation, restore
# it immediately after -- keeps the rest of the script's nounset safety net
# intact.
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set -u

# --- Step 1: mesh generation (Neper -> .msh/.inp) --------------------------
echo "=== [1/4] Generating mesh (env: ${NEPER_ENV}) ==="
set +u
conda activate "$NEPER_ENV"
set -u
./generate_odsfe_rve.sh "$CONFIG" "$N_GRAINS" "$GRAIN_DIAMEQ" "$MIN_GRAINS_ACROSS" "$CROSS_SECTION_MM"

OUTDIR="rve_${CONFIG}_n${N_GRAINS}"
MSH_FILE="${OUTDIR}/rve_${CONFIG}.msh"
if [ ! -f "$MSH_FILE" ]; then
    echo "ERROR: expected mesh file not found: ${MSH_FILE}"
    exit 1
fi

# --- Step 2: auto-build mesh_prep_<config>.i from actual mesh bounds -------
# Stays in the neper env -- only needs python3+meshio, which generate_odsfe_
# rve.sh's own Step 3 (vtu conversion) already relies on in that same env.
echo "=== [2/4] Building mesh_prep_${CONFIG}.i from actual mesh bounds ==="
python3 make_mesh_prep.py "$MSH_FILE" "$CONFIG"

# --- Step 3: mesh-only split (sidesets + BreakMeshByBlockGenerator) --------
# cd into OUTDIR: mesh_prep_<config>.i now lives there (Step 2 writes it
# alongside the .msh), and its internal [Mesh]/file reference is a bare
# basename -- whether MOOSE resolves that relative to the .i's own location
# or the process's CWD is unverified, so cd'ing guarantees it works either
# way rather than assuming.
echo "=== [3/4] Running mesh-only split (env: ${MOOSE_ENV}) ==="
set +u
conda activate "$MOOSE_ENV"
set -u
pushd "$OUTDIR" > /dev/null
SPLIT_FILE_NAME="rve_${CONFIG}_split.e"
mpirun -n "$MPI_RANKS" "$ARAGONITE_OPT" \
    -i "mesh_prep_${CONFIG}.i" --mesh-only "$SPLIT_FILE_NAME"

if [ ! -f "$SPLIT_FILE_NAME" ]; then
    echo "ERROR: mesh-only split did not produce ${SPLIT_FILE_NAME} in ${OUTDIR}"
    popd > /dev/null
    exit 1
fi

# --- Step 4: extract CZM boundaries, generate final inputs -----------------
# Still inside OUTDIR -- keeps extract_czm_boundaries.py/generate_odsfe_
# inputs.py's outputs co-located with everything else for this config.
# Stays in moose env -- switch to $NEPER_ENV here instead if these two
# scripts actually need neper-env's python packages rather than moose-env's;
# not verified either way, adjust if step 4 errors on a missing import.
echo "=== [4/4] Extracting CZM boundaries and generating inputs ==="
INTERFACES_FILE_NAME="rve_${CONFIG}_split_interfaces.txt"
python3 "${OLDPWD}/extract_czm_boundaries.py" "$SPLIT_FILE_NAME"
# extract_czm_boundaries.py's actual output filename convention was not
# re-verified here -- confirm INTERFACES_FILE_NAME matches what it actually
# writes before relying on this line; adjust if it writes elsewhere.
python3 "${OLDPWD}/generate_odsfe_inputs.py" "rve_${CONFIG}.msh" "$INTERFACES_FILE_NAME" "$CONFIG"
popd > /dev/null

SPLIT_FILE="${OUTDIR}/${SPLIT_FILE_NAME}"
INTERFACES_FILE="${OUTDIR}/${INTERFACES_FILE_NAME}"

echo "=== Done: ${CONFIG} ==="
echo "  Mesh:       ${MSH_FILE}"
echo "  Split mesh: ${SPLIT_FILE}"
echo "  Interfaces: ${INTERFACES_FILE}"
echo "  Config:     ${OUTDIR}/config_${CONFIG}.i"

