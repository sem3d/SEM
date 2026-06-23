# Anisotropic-from-file elastic material (reads Cstar.h5 -> Cij2d).
# For the isotropic-limit test, Cstar.h5 = Cstar_elastic_iso.h5 (Cij of an
# isotropic medium vp=3000 vs=1700 rho=2000); result must match the plain
# isotropic run (this directory WITHOUT material.spec).
material 0 {
  domain = solid;
  deftype = Hooke_Aniso;
  spacedef = file;
  filename0 = "Cstar.h5";
};
