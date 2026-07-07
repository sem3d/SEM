# Solid-fluid interface, anisotropic-from-file (both sides read a Cstar.h5).
# Material indices match the mesher output (mater.in): 0 = solid (S), 1 = fluid (F).
# The Cstar fields are isotropic-equivalent constants, so this run should match a
# plain isotropic S/F run (t5_solid_fluid_iso) to round-off.
material 0 {
  domain = solid;
  deftype = Hooke_Aniso;
  spacedef = file;
  filename0 = "Cstar_solid.h5";
};
material 1 {
  domain = fluid;
  deftype = Fluid_Aniso;
  spacedef = file;
  filename0 = "Cstar_fluid.h5";
};
