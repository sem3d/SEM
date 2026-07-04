# SEM2D anisotropic-from-file — test cases

Run on the cluster after building SEM2D. Generate the test material with
`python3 gen_test_cstar2d.py`, rename the chosen file to the name in
`material.spec` (e.g. `Cstar.h5`), and run.

The three gates below are listed in increasing risk: do them in order.

--------------------------------------------------------------------------------
## Generating the mesh

Each test needs `mesh4spec.NNNN.h5` + `material.input`. Two ways:

1. **mesher2D (on the fly)** — preferred. Reads `mat.dat` + `mater.in`, writes
   `mesh4spec.NNNN.h5` and `material.input`. `mesh.input` feeds the stdin menu
   (line 1 = NPROCS, line 2 = `1` for on-the-fly):

       .../build/MESH2D/mesher2D < mesh.input > outputmesh.log

   `mat.dat`/`mater.in`/`mesh.input` are provided in each test dir for the
   `[0,0.05] x [0,0.03]` 10x6 grid. (Overwrites the static `material.input`
   with identical content.)

2. **gen_test_mesh2d.py** — fallback that writes the HDF5 mesh directly, then
   partition with `mesh2dc Nproc mesh_input.h5 mesh4spec`. Use the bundled
   `material.input` as-is.

For these isotropic-limit gates use `has_pml=0` (no PML). To add PML, set the
three PML lines of `mat.dat` (lateral W/E are always on when `has_pml`>0):

    2              # has_pml = number of PML element layers (npml)
    0 1            # pml_top pml_bottom (1 = add PML on that side; here bottom only)
    5 2 10. 0. 0.  # ngllPML  npow  Apow  omegac  kc  (PML element NGLL + attenuation)
    4              # mesh_type

mesher2D extends the grid by `npml` element layers on each active side and emits
the derived `P` materials + the PML descriptor block
(`Filtering npow Apow Px Left Pz Down omegac kc`) into `material.input`, as the
SEM2D solver (`Domain.F90 read_material_file`) expects.

--------------------------------------------------------------------------------
## Test 1 — ELASTIC isotropic-limit   [should pass first]

Goal: an aniso-from-file run whose Cij IS isotropic must match a plain isotropic
run. Validates the reader + the Cij->Acoeff kernel + the rho-from-Cstar mass.

material.spec:
    material 0 {
      domain = solid;
      deftype = Hooke_Aniso;
      spacedef = file;
      filename0 = "Cstar.h5";
    };

Use Cstar_elastic_iso.h5 (rho=2000, vp=3000, vs=1700) -> rename to Cstar.h5.
Reference run: the SAME mesh/source with a plain material.input isotropic medium
(vp=3000, vs=1700, rho=2000), no material.spec.

PASS: traces from the two runs match (to round-off). Look for log lines:
  [aniso] subdomain 0 ...
  [aniso] elastic Cstar tensor read & interpolated to GLL ...
  [aniso] CG Acoeff + mass rebuilt from Cstar on N element(s).

--------------------------------------------------------------------------------
## Test 2 — FLUID isotropic-limit   [pure fluid]

Goal: a fluid-aniso run with iRho=(1/rho)I, iKappa=1/kappa must match an isotropic
acoustic run. Validates the potential domain (mass 1/kappa, stiffness rho^-1) + the
fluid source.

material.spec:
    material 0 {
      domain = fluid;
      deftype = Fluid_Aniso;
      spacedef = file;
      filename0 = "Cstar.h5";
    };

Use Cstar_fluid_iso.h5 (rho=1000, vp=1500) -> rename to Cstar.h5.

input.spec source (pressure pulse in the fluid):
    source {
      coords = 0.025  0.015 ;      # inside the domain, NOT on an element boundary
      type = fluidpulse;
      func = ricker;
      ... (tau/freq as usual) ...
    };
    # receiver(s) via the usual capteurs/stations mechanism

PASS: a clean acoustic pulse propagates at vp=1500 (check arrival time at a
receiver at known distance). Log:
  [aniso] fluid-aniso material read & interpolated to GLL ...
  [aniso] fluid-aniso mass(1/kappa) + stiffness built on N element(s).
RECEIVER NOTE: for the fluid, the recorded Veloc column 0 = VelPhi = -pressure.
  So pressure(t) = -(column 0). Column 1 stays ~0.

--------------------------------------------------------------------------------
## Test 3 — FLUID anisotropic density   [pure fluid]

Use Cstar_fluid_aniso.h5 (rho_x=1000, rho_z=1600) -> rename to Cstar.h5.
Same material.spec / source as Test 2.

CHECK: horizontal vs vertical wavefront speeds differ. With v_i = sqrt(kappa/rho_i):
  vx = sqrt(kappa/1000) = 1500 m/s ;  vz = sqrt(kappa/1600) ~ 1186 m/s.
Measure arrival times at receivers placed at equal distance along x and along z;
the ratio should be ~ vz/vx = sqrt(1000/1600) = 0.79.

--------------------------------------------------------------------------------
## Test 4 — SOLID-FLUID interface   [EXPERIMENTAL / NOT verified]

Solid<->fluid-aniso coupling is implemented but UNVERIFIED — it also requires
the interface-DOF separation (see solid_fluid_coupling_2d.F90.scaffold header and
the plan). Only attempt after Tests 1-3 pass.

Setup: a mesh with a solid subdomain (material 0) over a fluid-aniso subdomain
(material 1), sharing a flat horizontal interface.
    material 0 { domain = solid; deftype = Hooke_Aniso; spacedef = file; filename0 = "Cstar_solid.h5"; };
    material 1 { domain = fluid; deftype = Fluid_Aniso; spacedef = file; filename0 = "Cstar_fluid.h5"; };
CHECK (when it works): reflection/transmission coefficients at the interface vs the
analytic acoustic-elastic values; total energy conserved.

--------------------------------------------------------------------------------
## mat.dat (mesh) reminder
Keep the bounding box consistent with gen_test_cstar2d.py (default [0,0.05]x[0,0.03]).
NO trailing dots in numbers (e.g. "0.05" not "0.05."). Put the source coords strictly
inside the domain and not exactly on an element boundary.

--------------------------------------------------------------------------------
## Input decks provided (everything except the .h5 and the mesh)

  t1_elastic_iso/   input.spec, material.input (solid), material.spec (Hooke_Aniso), stations.txt
  t2_fluid/         input.spec, material.input (fluid), material.spec (Fluid_Aniso), stations.txt
                    -> reuse for fluid_iso AND fluid_aniso, just swap Cstar.h5.

To run a case:
  1. python3 gen_test_cstar2d.py        (in this dir)
  2. cp the right file to the test dir as Cstar.h5
       t1: Cstar_elastic_iso.h5 ; t2 iso: Cstar_fluid_iso.h5 ; t2 aniso: Cstar_fluid_aniso.h5
  3. generate the mesh "mesh4spec" with the 2D mesher for the SAME box [0,0.05]x[0,0.03]
     (NOT provided here -- needs your 2D mesher; keep the box consistent with gen_test_cstar2d.py).
  4. run SEM2D in the test dir.

For the T1 isotropic-limit REFERENCE run: same dir but WITHOUT material.spec (rename it away)
-> SEM2D uses the plain isotropic material.input; traces must match the material.spec run.

Notes / caveats:
- material.input columns: type Vp Vs Rho NGLLx <ignored> NGLLz Dt Qp Qs. The Dt field is a
  placeholder (1e-8); the actual dt is set from time_scheme%courant. ngll in input.spec and the
  NGLL columns should be consistent.
- source type=fluidpulse (type 3) is the fluid source. type=pressure (type 7) now STOPs with an
  error in 2D (not implemented) -- use fluidpulse.
- source coords are mid-domain; if they land exactly on an element boundary, nudge them slightly.
- stations.txt: one "x z" per line; receivers 1&2 are at equal distance (0.015) along x and z
  for the test-3 speed-ratio check.

--------------------------------------------------------------------------------
## 2D mesh: use the automatic mesher MESH2D (input.create provided)

The 2D on-the-fly mesher is MESH2D (SEM/MESH2D, program create_model_2D in
create_2D.f90). It reads "input.create" (the 2D analogue of the 3D mat.dat) and
builds a regular n_elem_x x n_elem_z mesh. An input.create is provided in each test
dir (domain [0,0.05]x[0,0.03], 100x60 elements, ngll=5).

  file_out  = mesh4spec            -> matches mesh_file in input.spec
  mater_out = material_mesher.input -> a THROWAWAY (single isotropic material the
              mesher writes). Do NOT use it for the run: use the provided
              material.input, which sets the correct type ("S" solid / "F" fluid).
              (For the fluid test the mesher would write a solid-type line.)

Steps: run the MESH2D mesher with input.create to produce mesh4spec (HDF5 mesh via
main_mesh2d if your build partitions it), then run SEM2D with the provided
input.spec / material.input / material.spec / Cstar.h5.

--------------------------------------------------------------------------------
## MESH — definitive workflow (supersedes the input.create notes above)

MESH2D builds only `mesh2dc` (a PARTITIONER): it reads an existing quad mesh and
splits it; it does NOT create geometry. The on-the-fly creator (create_2D.f90 /
input.create) is NOT compiled by the build, so use this instead:

  1. python3 gen_test_mesh2d.py            -> mesh_input.h5  (/Nodes, /Sem2D/Quad4, /Sem2D/Mat)
  2. mesh2dc <Nproc> mesh_input.h5 mesh4spec   (argv, NOT stdin!) -> mesh4spec.0000.h5 ...
       e.g. SLURM:
       srun /home/${SLURM_JOB_USER}/WIP/SEM/build/MESH2D/mesh2dc 4 mesh_input.h5 mesh4spec > outputmesh.log
  3. run SEM2D (mesh_file = "mesh4spec" in input.spec) with material.input/material.spec/Cstar.h5.

The SAME mesh_input.h5 serves all tests (single material, index 0): solid vs fluid is
set in material.input (material 0 = "S" or "F"), not in the mesh.
NOTE: mesh2dc takes command-line ARGUMENTS, not stdin -- the 3D-style "mesher < mesh.input"
does NOT apply in 2D. (input.create / create_model_2D would only be needed if you build that tool.)
