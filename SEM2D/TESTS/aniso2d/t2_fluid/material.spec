# Anisotropic-density fluid (reads Cstar.h5 -> IDensTensor2d + invKappa2d).
#   isotropic-limit test : Cstar.h5 = Cstar_fluid_iso.h5   (must match isotropic acoustic vp=1500)
#   anisotropy test      : Cstar.h5 = Cstar_fluid_aniso.h5 (vx=1500, vz~1186; ratio ~0.79)
material 0 {
  domain = fluid;
  deftype = Fluid_Aniso;
  spacedef = file;
  filename0 = "Cstar.h5";
};
