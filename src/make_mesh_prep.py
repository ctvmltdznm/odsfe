#!/usr/bin/env python3
"""
make_mesh_prep.py

Generates mesh_prep_<config>.i from the ACTUAL bounds of a given .msh file
-- no hardcoded thresholds, ever. This replaces hand-editing mesh_prep.i's
left/right/front/back/bottom/top combinatorial_geometry values per config,
which is exactly what caused the earlier 4.713mm-vs-20mm mismatch.

write_hdf5=true is included as cheap defensive insurance (Exodus classic
format has a hard ceiling on total sidesets+nodesets that we hit directly
at ~20,000+ boundaries -- see project notes). It costs nothing at normal
scale and avoids re-discovering that failure mode if grain count ever
creeps up again.

Deliberately does NOT include any nodeset-deletion/cleanup logic. That
approach was considered and rejected: BoundaryDeletionGenerator's exact
scope (sideset-only vs. sideset+nodeset pair) was never verified, and at
the grain counts we're now targeting (via min_grains_across=2 or 3) the
underlying ID-collapse problem shouldn't occur in the first place. If it
resurfaces, diagnose fresh rather than assuming this script's old
workaround still applies -- the mechanism was never fully confirmed.

Usage:
    python3 make_mesh_prep.py <msh_file> <config_name> [--tol-frac 0.001]

Example:
    python3 make_mesh_prep.py rve_needles_x_n200/rve_needles_x.msh needles_x
    # writes mesh_prep_needles_x.i
"""
import argparse
from pathlib import Path
import meshio

parser = argparse.ArgumentParser()
parser.add_argument("msh_file")
parser.add_argument("config")
parser.add_argument("--tol-frac", type=float, default=0.001,
                     help="Sideset threshold as a fraction of the domain's "
                          "x-extent (default 0.001 = 0.1%%, matching the "
                          "convention already used in earlier mesh_prep.i "
                          "files).")
parser.add_argument("--outdir", type=str, default=None,
                     help="default: same directory as msh_file")
args = parser.parse_args()

outdir = Path(args.outdir) if args.outdir else Path(args.msh_file).parent
outdir.mkdir(parents=True, exist_ok=True)

print(f"Reading bounds from {args.msh_file} ...")
m = meshio.read(args.msh_file)
pts = m.points
bounds = [(float(pts[:, i].min()), float(pts[:, i].max())) for i in range(3)]
extents = [hi - lo for lo, hi in bounds]
tol = args.tol_frac * extents[0]

print("Domain bounds (mm):")
for axis, (lo, hi) in zip("xyz", bounds):
    print(f"  {axis}: {lo:.6f} to {hi:.6f}  (extent {hi - lo:.6f})")
print(f"Sideset tolerance: {tol:.6f} mm ({args.tol_frac*100:.3g}% of x-extent)")

# (low_name, high_name, axis_letter, axis_index)
faces = [("left", "right", "x", 0), ("front", "back", "y", 1), ("bottom", "top", "z", 2)]

out_path = outdir / f"mesh_prep_{args.config}.i"
with open(out_path, "w") as f:
    f.write(f"""[Mesh]
  [file]
    type = FileMeshGenerator
    file = {Path(args.msh_file).name}
  []
""")
    prev = "file"
    for lo_name, hi_name, axis, idx in faces:
        lo, hi = bounds[idx]
        f.write(f"""  [{lo_name}]
    type = ParsedGenerateSideset
    input = {prev}
    combinatorial_geometry = '{axis} < {lo + tol:.6f}'
    new_sideset_name = '{lo_name}'
  []
  [{hi_name}]
    type = ParsedGenerateSideset
    input = {lo_name}
    combinatorial_geometry = '{axis} > {hi - tol:.6f}'
    new_sideset_name = '{hi_name}'
  []
""")
        prev = hi_name
    f.write(f"""  [break]
    type = BreakMeshByBlockGenerator
    input = {prev}
    split_interface = true
  []
[]

[Outputs]
  [out]
    type = Exodus
    write_hdf5 = true
  []
[]
""")

print(f"Wrote {out_path}")

