#!/usr/bin/env python3
"""
extract_czm_boundaries.py

After BreakMeshByBlockGenerator (split_interface=true) runs, MOOSE creates one
sideset per grain-grain contact, named 'BlockX_BlockY'. At 150 grains there
can be hundreds of these -- too many to hand-type into `boundary = '...'`
like the 7-grain coral example does.

This script does two things depending on what you give it:

  1. Given a raw Neper .msh file: converts it to Exodus via meshio (MOOSE's
     FileMeshGenerator has had trouble reading Neper's .msh directly -- this
     sidesteps that). It does NOT contain interface sidesets yet, since those
     only get created when BreakMeshByBlockGenerator actually runs inside
     MOOSE. You still need to run mesh-only after this before interfaces
     exist to extract.

  2. Given an already-split Exodus file (the output of running your [Mesh]
     block, including BreakMeshByBlockGenerator, through `--mesh-only`):
     extracts the grain-grain interface sideset names and writes them as a
     ready-to-paste boundary string.

Full workflow:
  1. python3 extract_czm_boundaries.py rve_ar4.msh
       -> writes rve_ar4.e
  2. Point [Mesh]/[file]/file at rve_ar4.e in a minimal mesh_prep.i containing
     your full [Mesh] block (ParsedGenerateSideset x6 + BreakMeshByBlockGenerator).
  3. moose-opt -i mesh_prep.i --mesh-only production_split.e
  4. python3 extract_czm_boundaries.py production_split.e
       -> writes production_split_interfaces.txt, ready to paste into
          `boundary = '...'`

Requires netCDF4 and meshio (both already in your toolchain).
"""
import sys
import os
import re
import netCDF4


def convert_msh_to_exodus(msh_path):
    import meshio
    exo_path = os.path.splitext(msh_path)[0] + ".e"
    m = meshio.read(msh_path)
    meshio.write(exo_path, m, file_format="exodus")
    print(f"# Converted {msh_path} -> {exo_path}")
    print(f"# Element blocks found: {sorted(set(c.type for c in m.cells))}")
    print("#")
    print("# This is the raw (pre-split) mesh -- no interface sidesets exist")
    print("# yet. Point FileMeshGenerator at this .e file, run your [Mesh]")
    print("# block (through BreakMeshByBlockGenerator) with --mesh-only, then")
    print("# rerun this script on THAT output to get the interface list.")
    return exo_path


def get_sideset_names(exodus_path):
    ds = netCDF4.Dataset(exodus_path, "r")
    names = []
    if "ss_names" in ds.variables:
        raw = ds.variables["ss_names"][:]
        for row in raw:
            # exodus stores fixed-width char arrays; strip null padding
            s = b"".join(c for c in row.tobytes().split(b"\x00")[:1]).decode(
                "utf-8", errors="ignore"
            )
            if not s:
                # fallback: join all bytes and strip nulls if the above missed it
                s = row.tobytes().split(b"\x00")[0].decode("utf-8", errors="ignore")
            names.append(s)
    ds.close()
    return names


def get_sideset_elem_side_pairs(ds, sideset_index_1based):
    """Read (element_id, side_id) pairs for one sideset by its 1-based index
    in ss_prop1/ss_names ordering."""
    elem_var = f"elem_ss{sideset_index_1based}"
    side_var = f"side_ss{sideset_index_1based}"
    if elem_var not in ds.variables or side_var not in ds.variables:
        return set()
    elems = ds.variables[elem_var][:]
    sides = ds.variables[side_var][:]
    return set(zip(elems.tolist(), sides.tolist()))


def check_boundary_overlap(exodus_path, names, interfaces, domain_faces):
    """Check whether any interface sideset shares (element, side) pairs with
    ANY other sideset that isn't itself a BlockX_BlockY interface. A genuine
    interior BlockX_BlockY interface should NEVER overlap with a true
    exterior domain-face sideset -- if it does, that's a direct, concrete
    cause of MOOSE's 'Element ... missing a neighbor ... has interface
    kernel(s) defined on the boundary' error, rather than a guess.

    IMPORTANT: this checks by POSITIONAL index (1-based, matching ss order
    in the file), not by name. ~1000+ sidesets in a typical run come back
    with an empty/blank name (likely companion primary/secondary sides that
    BreakMeshByBlockGenerator creates internally per interface) -- a
    name-keyed dict silently collapses all of those onto one key and skips
    checking the rest. Positional indexing checks every sideset in the file,
    named or not.
    """
    ds = netCDF4.Dataset(exodus_path, "r")

    interface_set = set(interfaces)
    domain_face_set = set(domain_faces) - {""}   # real domain names only

    # Build (index, name, category) for every sideset by position.
    all_pairs = {}   # idx (1-based) -> set of (elem, side)
    idx_category = {}  # idx -> 'interface' | 'domain' | 'anonymous'
    for i, name in enumerate(names):
        idx = i + 1
        pairs = get_sideset_elem_side_pairs(ds, idx)
        all_pairs[idx] = pairs
        if name in interface_set:
            idx_category[idx] = 'interface'
        elif name in domain_face_set:
            idx_category[idx] = 'domain'
        else:
            idx_category[idx] = 'anonymous'

    interface_idxs = [i for i, c in idx_category.items() if c == 'interface']
    non_interface_idxs = [i for i, c in idx_category.items() if c != 'interface']

    problems = []
    for iidx in interface_idxs:
        ipairs = all_pairs[iidx]
        if not ipairs:
            continue
        for oidx in non_interface_idxs:
            opairs = all_pairs[oidx]
            if not opairs:
                continue
            overlap = ipairs & opairs
            if overlap:
                problems.append((names[iidx - 1], oidx, idx_category[oidx], overlap))

    ds.close()
    return problems


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <mesh_file.msh | split_mesh.e>")
        sys.exit(1)

    path = sys.argv[1]

    if path.endswith(".msh"):
        convert_msh_to_exodus(path)
        return

    names = get_sideset_names(path)

    interface_pattern = re.compile(r"^Block\d+_Block\d+$")
    interfaces = sorted(
        [n for n in names if interface_pattern.match(n)],
        key=lambda s: [int(x) for x in re.findall(r"\d+", s)],
    )
    named_non_interface = [n for n in names if n not in interfaces and n != ""]
    n_anonymous = sum(1 for n in names if n == "")

    if not interfaces:
        print(f"# No 'BlockX_BlockY' sidesets found in {path}.")
        print("# This looks like the pre-split mesh -- BreakMeshByBlockGenerator")
        print("# hasn't run yet. Run your full [Mesh] block (through")
        print("# BreakMeshByBlockGenerator) with --mesh-only first, then rerun")
        print("# this script on that output.")
        print(f"# Sidesets present: {names}")
        return

    print(f"# Found {len(names)} total sidesets in {path}")
    print(f"# Named domain-face sidesets (not CZM interfaces): {named_non_interface}")
    print(f"# Anonymous/blank-named sidesets: {n_anonymous} "
          f"(likely BreakMeshByBlockGenerator's internal companion/secondary "
          f"sides -- checked positionally below, not skipped)")
    print(f"# Grain-grain interface sidesets: {len(interfaces)}")

    # Check for the specific failure mode that caused
    # "Element ... missing a neighbor ... has interface kernel(s)":
    # an interface sideset sharing (elem, side) pairs with ANY other
    # sideset (named or anonymous) -- checked by position, not name, since
    # name-based lookup silently missed ~1000+ anonymous sidesets before.
    print()
    print("# Checking for interface/other-sideset overlap "
          "(root cause of the 'missing neighbor' crash, if present)...")
    problems = check_boundary_overlap(path, names, interfaces, named_non_interface)
    if problems:
        print(f"# FOUND {len(problems)} OVERLAPPING (interface, other-sideset) PAIR(S):")
        for iname, oidx, ocat, overlap in problems[:20]:
            print(f"#   {iname} overlaps sideset #{oidx} ({ocat}) on "
                  f"{len(overlap)} face(s): {sorted(overlap)[:5]}"
                  f"{'...' if len(overlap) > 5 else ''}")
        if len(problems) > 20:
            print(f"#   ... and {len(problems) - 20} more")
        print("# These interfaces share faces with another sideset (domain")
        print("# boundary or an anonymous companion side) -- this is the")
        print("# concrete, checkable cause of the crash, not a guess.")
    else:
        print("# No overlap found across ALL sidesets, named and anonymous.")
        print("# The threshold/overlap hypothesis is now ruled out with")
        print("# actual data -- the crash has a different cause (worth")
        print("# looking at a genuine mesh degeneracy, e.g. a sliver/")
        print("# non-manifold element near a triple junction).")

    print()
    print("# --- Paste this into `boundary = '...'` for both the CohesiveZone")
    print("#     physics block and the HomogenizedExponentialCZM material ---")
    print(" ".join(interfaces))

    # Also write to a file since 150-grain lists can be long / unwieldy in a terminal
    out_path = path.rsplit(".", 1)[0] + "_interfaces.txt"
    with open(out_path, "w") as f:
        f.write(" ".join(interfaces))
    print(f"\n# Also written to: {out_path}")


if __name__ == "__main__":
    main()
