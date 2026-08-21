#!/usr/bin/env bash
# ============================================================================
# generate_odsfe_rve.sh
#
# Generates polycrystal RVE meshes for coarse-grained Fe-10Al-4Cr-4Y2O3
# (FeAlOY) ODS nanocomposite. Grain elongation is GLOBALLY ALIGNED (rolled-
# metal texture, long axes coinciding) using Neper's literal aspratio(rx,ry,rz)
# morphology key. Grain count is the fixed anchor (default 150, coral-
# equivalent); domain size is DERIVED from grain count + absolute grain size,
# not the other way around -- EXCEPT that the derived size is now also
# checked against a minimum-grains-across-the-long-axis floor (see below),
# and the larger of the two wins. For elongated configs this floor usually
# dominates, meaning the actual grain count Neper produces will end up
# higher than N_GRAINS -- expected, and reported by the diagnostic print.
#
# WHY THE FLOOR WAS ADDED: the pure volume-balance domain size is correct
# for isotropic (equiaxed) grains, but says nothing about a grain's actual
# extent along its elongation axis. Checked directly for the aspratio
# values below (all chosen with rx*ry*rz=1, so grain VOLUME stays fixed at
# the equivalent-sphere volume regardless of config): at the previous
# volume-only domain size, needles_x/needles_y's domain was only ~1.07x the
# needle's own length -- i.e. barely wider than one grain, guaranteed to
# span the box end-to-end. Pancakes were ~2.14x -- better, still thin.
# Isotropic was ~4.28x, already close to adequate since equiaxed grains
# aren't elongated in any one direction. Per grain-geometry requirement
# from the group: domain edge must be >= MIN_GRAINS_ACROSS (default 5)
# times a grain's largest extent, so no grain's ends land on two opposite
# domain surfaces simultaneously.
#
# Units: millimeters throughout (domain, grain size). Stress stays in MPa
# downstream in MOOSE (MPa = N/mm^2, consistent with mm lengths).
#
# Outputs:
#   rve_<config>.msh  -- raw Neper mesh
#   rve_<config>.vtu  -- for visualization only (meshio .msh->.vtu is a
#                        reliable path; this is what you already viewed OK)
#   rve_<config>.inp  -- Abaqus format, for MOOSE FileMeshGenerator. Neper
#                        writes this natively -- no meshio Exodus writer
#                        involved, since that path is confirmed broken
#                        (fails in ParaView's vtkIOSSReader / SEACAS-strict
#                        readers). Exodus conversion is dropped from this
#                        script until we find a reliable way to produce it.
#
# Usage:
#   ./generate_odsfe_rve.sh <config_name> [n_grains] [grain_diameq_mm] [min_grains_across] [cross_section_mm]
#   defaults: n_grains=150, grain_diameq_mm=1.0, min_grains_across=5, cross_section_mm=0 (auto)
#
# cross_section_mm (5th arg, optional): if >0, FIXES the short/non-elongated
# axes to this value directly (e.g. the two thin axes of a needle, or the
# one thin axis of a pancake), deriving the remaining long axis/axes from
# the volume-balance target instead of the elongation-floor formula. Use
# this when the auto-computed cross-section looks too small for meaningful
# sampling in the thin direction(s). TRADEOFF: at fixed n_grains, a bigger
# fixed cross-section eats into the long axis's own margin above
# min_grains_across -- the script prints the resulting margin and warns if
# it drops below min_grains_across; increase n_grains to compensate if so.
#
# Example:
#   ./generate_odsfe_rve.sh needles_x                    # 150 grains, 1mm, 5x floor, auto cross-section
#   ./generate_odsfe_rve.sh needles_x 150 0.7             # 150 grains, 0.7mm
#   ./generate_odsfe_rve.sh needles_x 150 1.0 7            # stricter 7x floor
#   ./generate_odsfe_rve.sh needles_x 200 1.0 3 4.5        # fixed 4.5mm cross-section
# ============================================================================
set -euo pipefail

CONFIG=${1:?"Usage: $0 <config_name> [n_grains=150] [grain_diameq_mm=1.0] [min_grains_across=5] [cross_section_mm]"}
N_GRAINS=${2:-200}
GRAIN_DIAMEQ=${3:-1.0}   # mm
MIN_GRAINS_ACROSS=${4:-5}
CROSS_SECTION_MM=${5:-0}   # mm; if >0, FIXES the short (non-elongated) axes
                            # to this value and derives the long axis/axes
                            # from the volume-balance target instead of the
                            # elongation-floor formula. See tradeoff note
                            # printed below -- fixing cross-section bigger,
                            # at the SAME grain count, eats into the long
                            # axis's own margin above MIN_GRAINS_ACROSS.
                            # Bump N_GRAINS if that margin comes out too low.

SIZE_CV=0.31
SEED="${SEED:-$RANDOM}"  # NOTE: random per run unless SEED is preset --
                  # identical parameters will still produce a DIFFERENT
                  # tessellation each invocation otherwise. If a run hits a
                  # meshing failure (e.g. a degenerate/sliver cell), a
                  # rerun with the same arguments may simply avoid it by
                  # luck of the draw. Set SEED=<fixed_value> as an env var
                  # before calling this script for a reproducible
                  # tessellation instead.
MORPHO_TOL="${MORPHO_TOL:-eps<1e-3||iter>=15000}"

case "$CONFIG" in
  isotropic)  ASPRATIO="1,1,1"       ;;
  needles_x)  ASPRATIO="4,0.5,0.5"   ;;
  pancake_z)  ASPRATIO="2,2,0.25"    ;;
  needles_y)  ASPRATIO="0.5,4,0.5"   ;;
  pancake_x)  ASPRATIO="0.25,2,2"    ;;
  *)
    echo "Unknown config '${CONFIG}'. Options: isotropic needles_x pancake_z needles_y pancake_x"
    exit 1
    ;;
esac

# --- Derive domain edges: anisotropic box, not a cube ----------------------
# Each axis gets its OWN size, satisfying its own elongation-floor
# requirement using that axis's actual grain extent -- not forced to match
# whatever the LONGEST axis needs, which is what a cube domain did and is
# why needles_x blew up to 11,602 grains at min_across=5 (the short axes
# were forced 8x larger than they needed to be, just to match the long
# axis's requirement).
#
# Method: compute each axis's own floor (min_grains_across * that axis's
# actual grain extent). If the resulting box's volume is below the
# volume-balance target for N_GRAINS, scale ALL THREE axes up by the same
# factor (preserves aspect ratio exactly, never violates the floor --
# scaling up only adds margin). This typically brings grain count back
# down near N_GRAINS for elongated configs, instead of the 50-100x
# overshoot the cube approach produced.
DOMAIN_SIZES=$(python3 -c "
import math
n = ${N_GRAINS}
d = ${GRAIN_DIAMEQ}
rx, ry, rz = ${ASPRATIO}
min_across = ${MIN_GRAINS_ACROSS}
cross_section = ${CROSS_SECTION_MM}

prod = rx * ry * rz
k = (d / 2.0) / prod**(1/3)
comps = (rx, ry, rz)
extents = [2.0 * k * r for r in comps]  # actual grain extent per axis
vol_per_grain = (math.pi / 6.0) * d**3
vol_target = n * vol_per_grain

max_comp = max(comps)
is_long = [abs(c - max_comp) < 1e-9 for c in comps]
n_long = sum(is_long)
n_short = 3 - n_long

if cross_section > 0 and 0 < n_long < 3:
    # DIRECT mode: cross_section and min_across together fully define the
    # box -- no volume-balance scaling, no margin to check afterward,
    # because min_across is now the literal definer of long-axis length,
    # not a floor checked against a volume-derived value. n_grains (n) is
    # NOT used here at all; grain count is a reported RESULT, not a target.
    L = [(min_across * extents[i]) if is_long[i] else cross_section for i in range(3)]
    binding = f'direct: cross-section={cross_section:g}mm fixed, {min_across:g}x grain-lengths along long axis (n_grains arg ignored in this mode)'
elif cross_section > 0:
    # cross_section given but config has no long/short distinction
    # (isotropic, or n_long==3) -- direct mode has no clear meaning here,
    # fall back to the volume-balance/elongation-floor logic below.
    L_floor = [min_across * e for e in extents]
    vol_floor = L_floor[0] * L_floor[1] * L_floor[2]
    if vol_floor < vol_target:
        scale = (vol_target / vol_floor) ** (1/3)
        L = [l * scale for l in L_floor]
        binding = 'volume balance (cross_section ignored -- no long/short axis distinction for this config)'
    else:
        L = L_floor
        binding = 'elongation floor (cross_section ignored -- no long/short axis distinction for this config)'
else:
    # Existing behavior: each axis satisfies its own elongation-floor
    # requirement, scaled up uniformly if needed to hit the volume target.
    L_floor = [min_across * e for e in extents]
    vol_floor = L_floor[0] * L_floor[1] * L_floor[2]
    if vol_floor < vol_target:
        scale = (vol_target / vol_floor) ** (1/3)
        L = [l * scale for l in L_floor]
        binding = 'volume balance (scaled up uniformly, aspect ratio preserved)'
    else:
        L = L_floor
        binding = 'elongation floor (already exceeds volume-balance target)'

direct_mode = 1 if (cross_section > 0 and 0 < n_long < 3) else 0
vol_final = L[0] * L[1] * L[2]
est_grains = vol_final / vol_per_grain
margins = [L[i] / extents[i] for i in range(3)]
min_margin = min(margins)
# In direct mode the long-axis margin is exact by construction (always
# passes trivially) and the short-axis figure is just whatever cross_section
# produces -- not a constraint the user asked to enforce, so don't warn.
margin_warning = 0 if direct_mode else (1 if min_margin < min_across else 0)
print(f'{L[0]:.4f} {L[1]:.4f} {L[2]:.4f} {extents[0]:.4f} {extents[1]:.4f} {extents[2]:.4f} '
      f'{est_grains:.0f} {min_margin:.2f} {margin_warning} {direct_mode} {binding}')
")
read -r LX LY LZ EXT_X EXT_Y EXT_Z EST_GRAINS MIN_MARGIN MARGIN_WARNING DIRECT_MODE BINDING <<< "$DOMAIN_SIZES"

echo "=== Config: ${CONFIG} (aspratio ${ASPRATIO}) ==="
echo "=== Grain extents per axis: x=${EXT_X}mm y=${EXT_Y}mm z=${EXT_Z}mm ==="
echo "=== Domain (anisotropic box): x=${LX}mm y=${LY}mm z=${LZ}mm ==="
echo "=== Bound by: ${BINDING} ==="
if [ "$DIRECT_MODE" = "1" ]; then
    echo "=== Estimated grain count: ${EST_GRAINS} (n_grains arg ignored in direct mode) ==="
else
    echo "=== Estimated grain count: ${EST_GRAINS} (target was ${N_GRAINS}) ==="
fi
echo "=== Smallest per-axis margin (domain/grain-extent): ${MIN_MARGIN}x (requested floor: ${MIN_GRAINS_ACROSS}x) ==="
if [ "$MARGIN_WARNING" = "1" ]; then
    echo "=== NOTE: fixed cross-section pushed a per-axis margin BELOW the ${MIN_GRAINS_ACROSS}x floor. ==="
    echo "===       Increase n_grains if this margin is too tight for your needs.        ==="
fi

OUTDIR="rve_${CONFIG}_n${N_GRAINS}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# --- Step 1: tessellation ---------------------------------------------------
# -n from_morpho + absolute diameq: grain size stays literal (mm), count
#   falls out from domain volume / mean grain volume -- for elongated
#   configs where the elongation floor dominates, this means the actual
#   grain count will be HIGHER than N_GRAINS (larger domain, same grain
#   size), which is expected and reported below, not a bug.
# Periodicity is OFF: it was producing fragmented/wrapped mesh output for
# elongated grains (needles and pancakes split into disconnected pieces
# across the periodic boundary) and complicates simple face-based Dirichlet
# BCs anyway. Non-periodic tessellation clipped to a plain box is guaranteed
# space-filling with no such fragmentation. With the elongation floor above,
# grains touching/crossing the domain boundary should now be rare rather
# than the dominant failure mode -- still possible for individual outliers
# given the lognormal size distribution's tail, but no longer systematic.
neper -T \
  -n from_morpho \
  -id "${SEED}" \
  -morpho "diameq:lognormal(${GRAIN_DIAMEQ},${SIZE_CV}),aspratio(${ASPRATIO})" \
  -domain "cube(${LX},${LY},${LZ})" \
  -morphooptistop "${MORPHO_TOL}" \
  -reg 1 \
  -statcell diameq,vol \
  -o "rve_${CONFIG}"

N_ACTUAL=$(wc -l < "rve_${CONFIG}.stcell" 2>/dev/null || echo "?")
echo "=== Grain count Neper actually generated: ${N_ACTUAL} (target was ${N_GRAINS}) ==="

# --- Step 2: mesh, native msh + Abaqus inp ----------------------------------
neper -M "rve_${CONFIG}.tess" \
  -order 1 \
  -elttype tet \
  -format msh,inp \
  -o "rve_${CONFIG}"

# --- Step 3: vtu for visualization only (proven-reliable meshio path) ------
python3 -c "
import meshio
m = meshio.read('rve_${CONFIG}.msh')
meshio.write('rve_${CONFIG}.vtu', m)
print('Wrote rve_${CONFIG}.vtu for viewing in ParaView')
"

echo "=== Done ==="
echo "View:  ${OUTDIR}/rve_${CONFIG}.vtu"
echo "MOOSE: ${OUTDIR}/rve_${CONFIG}.inp  (point FileMeshGenerator at this -- untested, try it and report back)"
