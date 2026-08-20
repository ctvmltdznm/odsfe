#!/usr/bin/env python3
"""
generate_odsfe_inputs.py
=========================
Produces a ready-to-run MOOSE .i file for a FeAlOY needle/pancake RVE,
following the same structural conventions as fill_and_generate_configs_v2.py
(coral pipeline): displacement-controlled loading, HomogenizedExponentialCZM
via [Physics][SolidMechanics][QuasiStatic][CohesiveZone], matching
postprocessor set, proven executioner block.

CHANGE LOG (this version):
  - fill_method: symmetric_isotropic -> symmetric_isotropic_E_nu. The old
    value silently misread C_ijkl as Lame's lambda/mu instead of (E, nu) --
    see the elasticity tensor bug writeup in the two-grain calibration
    thread. This affected every prior RVE run's bulk material.
  - delta_0_normal/delta_0_tangent, eta, quality_std_dev, damage_viscosity:
    updated to the values calibrated in the two-grain stiffness/eta sweeps
    (see czm_model_description.md), not the original placeholders.
  - Loading rate: now a FIXED constant (0.000942 mm per time unit),
    independent of domain size, matching the rate validated across the
    isotropic RVE runs -- NOT recomputed per-config from dx*0.01/end_time
    as before. Overridable via --rate if a config genuinely needs
    something different, but the default is this shared, already-tested
    value, not a per-config re-derivation.
  - end_time default: 5.0 -> 200.0, matching the slow-rate/long-end_time
    combination that got the isotropic runs through genuine peak+softening
    instead of stalling early.
  - Output files (config_<name>.i) now write into the SAME DIRECTORY as
    the input mesh file, not the current working directory -- avoids
    scattering config_*.i files across the pipeline root when running many
    configs. If you want a different output location entirely, pass
    --outdir explicitly.

Usage:
    python3 generate_odsfe_inputs.py <mesh_file> <interfaces_txt> <config_name> [options]

    mesh_file       - the mesh MOOSE will read (verify this loads FIRST)
    interfaces_txt  - output of extract_czm_boundaries.py (space-separated
                       BlockX_BlockY list)
    config_name     - label used in output filename / comments

Options:
    --end-time <float>   default 200.0
    --rate <float>       mm per time unit, default 0.000942 (fixed,
                          domain-size-independent -- see CHANGE LOG above)
    --outdir <path>      default: same directory as mesh_file

Example:
    python3 generate_odsfe_inputs.py rve_needles_x_n200/rve_needles_x.msh \\
        rve_needles_x_n200/rve_needles_x_split_interfaces.txt needles_x
"""
import argparse
import sys
import meshio
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("mesh_file")
parser.add_argument("interfaces_txt")
parser.add_argument("config_name")
parser.add_argument("--end-time", type=float, default=200.0)
parser.add_argument("--rate", type=float, default=0.000942,
                     help="mm per time unit -- fixed, NOT scaled by domain "
                          "size. See CHANGE LOG in module docstring.")
parser.add_argument("--outdir", type=str, default=None,
                     help="default: same directory as mesh_file")
args = parser.parse_args()

MESH_FILE      = args.mesh_file
INTERFACES_TXT = args.interfaces_txt
CONFIG_NAME    = args.config_name
END_TIME       = args.end_time
RATE           = args.rate
OUTDIR         = Path(args.outdir) if args.outdir else Path(MESH_FILE).parent

# -- Read interface boundary list ---------------------------------------------
interfaces = Path(INTERFACES_TXT).read_text().strip()
n_interfaces = len(interfaces.split())
print(f"Loaded {n_interfaces} interfaces from {INTERFACES_TXT}")

# -- Mesh bounds (mm, from the actual mesh -- no separate JSON step needed) --
m = meshio.read(MESH_FILE)
pts = m.points
xmin, ymin, zmin = pts.min(axis=0)
xmax, ymax, zmax = pts.max(axis=0)
dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
tol_x, tol_y, tol_z = dx * 0.001, dy * 0.001, dz * 0.001
print(f"Bounds: X[{xmin:.4f},{xmax:.4f}] Y[{ymin:.4f},{ymax:.4f}] Z[{zmin:.4f},{zmax:.4f}] (mm)")

# -- CZM parameters ------------------------------------------------------------
# Calibrated via the two-grain stiffness-matching and eta sweeps (see
# czm_model_description.md for the full derivation). delta_0_normal/
# delta_0_tangent are tied by elasticity (E/G ratio), not independently
# fit -- see the model description doc, section 1.
CZM = {
    'normal_strength':          20.0,
    'shear_strength_s':         12.0,
    'shear_strength_t':         12.0,
    'delta_0_normal':           2.8425461e-05,
    'delta_0_tangent':          7.1601097e-05,
    'mu':                       1.0,
    'eta':                      0.3,
    'quality_std_dev':          0.15,
    'spatial_quality_std_dev':  0.15,
    'damage_viscosity':         5,
    'spatial_random_seed':      12345,
}

def fmt(v):
    return f"{v:g}" if isinstance(v, float) else f"{v}"

def czm_mat_block():
    lines = ['  [czm_grain_boundary]', '    type = HomogenizedExponentialCZM',
              f"    boundary = '{interfaces}'"]
    lines += [f'    {k:<24} = {fmt(v)}' for k, v in CZM.items()]
    lines.append('  []')
    return '\n'.join(lines)

def sidesets_block():
    faces = [
        ('left',   f"x < {xmin+tol_x:.6f}", 'file'),
        ('right',  f"x > {xmax-tol_x:.6f}", 'left'),
        ('front',  f"y < {ymin+tol_y:.6f}", 'right'),
        ('back',   f"y > {ymax-tol_y:.6f}", 'front'),
        ('bottom', f"z < {zmin+tol_z:.6f}", 'back'),
        ('top',    f"z > {zmax-tol_z:.6f}", 'bottom'),
    ]
    lines = []
    for name, expr, inp in faces:
        lines.append(f"  [{name}]\n    type = ParsedGenerateSideset\n"
                      f"    input = {inp}\n"
                      f"    combinatorial_geometry = '{expr}'\n"
                      f"    new_sideset_name = '{name}'\n  []")
    return '\n'.join(lines)

# -- BCs: displacement-controlled, LOADING ALONG X ----------------------------
# Right face (x=xmax) pulled in x; left face (x=xmin) fixed in x.
# fix_front_y/fix_bottom_z pin the lateral DOFs on faces OTHER than the
# fixed/loaded x-faces (front=y_min, bottom=z_min), preventing rigid-body
# translation in y/z without over-constraining the loaded face itself.
strain_rate_fn = f"{RATE:g} * t"

content = f"""# {'='*70}
# FeAlOY RVE -- config: {CONFIG_NAME}
# {'='*70}
# Mesh: {Path(MESH_FILE).name}  ({n_interfaces} grain-grain interfaces)
# Bulk: isotropic linear elastic, E=105000 MPa, nu=0.29 (see [Materials])
# CZM:  see czm_model_description.md for the full model/calibration writeup.
# Loading: uniaxial, along X, displacement-controlled, FIXED rate
#   ({RATE:g} mm/time-unit, independent of domain size -- see module
#   docstring CHANGE LOG).
#   "Applied stress at fracture" = peak of the stress_xx_avg postprocessor
#   over the run (NOT a prescribed load -- force control is unstable past
#   the softening peak).
#   Interface-normal stress time history = czm_max_normal_traction
#   (SideExtremeValue over all interfaces, per-timestep max).
#
# VERIFY boundary names with a --mesh-only run before trusting this file.
# {'='*70}

[Mesh]
  [file]
    type = FileMeshGenerator
    file = '{Path(MESH_FILE).name}'
    use_for_exodus_restart = true
  []
{sidesets_block()}
  [break]
    type = BreakMeshByBlockGenerator
    input = top
    split_interface = true
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [bulk]
        strain = SMALL
        incremental = true
        add_variables = true
        volumetric_locking_correction = true
        generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz
                           strain_xx strain_yy strain_zz vonmises_stress'
      []
      [CohesiveZone]
        [grain_boundary]
          boundary = '{interfaces}'
          strain = SMALL
          generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction
                             jump_x jump_y jump_z normal_jump tangent_jump'
        []
      []
    []
  []
[]

[AuxVariables]
  [damage]          order = CONSTANT  family = MONOMIAL []
  [delta_eff_czm]   order = CONSTANT  family = MONOMIAL []
[]

[AuxKernels]
  [czm_damage_aux]
    type = MaterialRealAux
    variable = damage
    property = damage
    execute_on = TIMESTEP_END
    check_boundary_restricted = false
    boundary = '{interfaces}'
  []
  [delta_eff_czm_aux]
    type = MaterialRealAux
    variable = delta_eff_czm
    property = delta_eff
    execute_on = TIMESTEP_END
    check_boundary_restricted = false
    boundary = '{interfaces}'
  []
[]

[Materials]
  [elasticity]
    type = ComputeElasticityTensor
    fill_method = symmetric_isotropic_E_nu
    C_ijkl = '105000 0.29'   # E = 105000 MPa, nu = 0.29
  []
  [stress]
    type = ComputeLinearElasticStress
  []

{czm_mat_block()}
[]

[BCs]
  # Uniaxial loading along X, displacement-controlled.
  [fix_left_x]
    type = DirichletBC  variable = disp_x  boundary = left   value = 0
  []
  [fix_front_y]
    type = DirichletBC  variable = disp_y  boundary = front  value = 0
  []
  [fix_bottom_z]
    type = DirichletBC  variable = disp_z  boundary = bottom value = 0
  []
  [load_right_x]
    type = FunctionDirichletBC  variable = disp_x  boundary = right
    function = '{strain_rate_fn}'
  []
[]

[Postprocessors]
  # -- Bulk / nominal applied stress -------------------------------------------
  # Peak of stress_xx_avg over the run = "applied stress at fracture".
  [stress_xx_avg]  type = ElementAverageValue  variable = stress_xx []
  [stress_xx_max]  type = ElementExtremeValue  variable = stress_xx []
  [strain_xx_avg]  type = ElementAverageValue  variable = strain_xx []
  [vonmises_avg]   type = ElementAverageValue  variable = vonmises_stress []

  # -- Interface (grain-boundary-normal) stress --------------------------------
  [czm_avg_normal_traction]
    type = SideAverageValue  variable = normal_traction
    boundary = '{interfaces}'
  []
  [czm_max_normal_traction]
    type = SideExtremeValue  variable = normal_traction
    boundary = '{interfaces}'
  []
  [czm_avg_tangent_traction]
    type = SideAverageValue  variable = tangent_traction
    boundary = '{interfaces}'
  []
  [czm_avg_normal_jump]
    type = SideAverageValue  variable = normal_jump
    boundary = '{interfaces}'
  []
  [czm_avg_tangent_jump]
    type = SideAverageValue  variable = tangent_jump
    boundary = '{interfaces}'
  []
  [czm_avg_damage]
    type = SideAverageValue  variable = damage
    boundary = '{interfaces}'
  []
  [czm_max_damage]
    type = SideExtremeValue  variable = damage
    boundary = '{interfaces}'
  []
  [czm_max_delta_eff]
    type = SideExtremeValue  variable = delta_eff_czm
    boundary = '{interfaces}'
  []
  # -- Spatial fields are in Exodus -- open in ParaView for per-interface detail
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-pc_type -pc_gamg_agg_nsmooths -ksp_type -ksp_gmres_restart -ksp_max_it -ksp_rtol'
  petsc_options_value = 'gamg     2                     gmres     300                200        1e-4'

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-2
  nl_max_its = 50
  l_max_its  = 200

  line_search = 'l2'

  dtmax    = 0.001
  dtmin    = 1e-6
  end_time = {END_TIME:g}

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.002
    optimal_iterations = 8
    iteration_window   = 2
    growth_factor       = 1.25
    cutback_factor       = 0.25
  []
[]

[Outputs]
  [csv]
    type = CSV
  []
  [out]
    type = Exodus
    time_step_interval = 10
  []
  print_linear_residuals = false
  perf_graph = false
  [console]
    type = Console
    verbose = false
  []
[]
"""

OUTDIR.mkdir(parents=True, exist_ok=True)
out_path = OUTDIR / f"config_{CONFIG_NAME}.i"
out_path.write_text(content)
print(f"Written: {out_path}")
print(f"\nVERIFY boundary names before trusting this file:")
print(f"  moose-opt -i {out_path} --mesh-only")

