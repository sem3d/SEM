
# Definition materiaux

material 0 {
domain = solid;
deftype = Kappa_Mu_Rho;
spacedef = file;
filename0 = "mat/h5/Mat_0_Kappa.h5";
filename1 = "mat/h5/Mat_0_Mu.h5";
filename2 = "mat/h5/Mat_0_Density.h5";
};

material 3 { copy = 0; };
material 4 { copy = 0; };
material 5 { copy = 0; };
material 6 { copy = 0; };
material 7 { copy = 0; };
material 8 { copy = 0; };
material 9 { copy = 0; };
material 10 { copy = 0; };

material 1 {
domain = solid;
deftype = Kappa_Mu_Rho;
spacedef = file;
filename0 = "mat/h5/Mat_1_Kappa.h5";
filename1 = "mat/h5/Mat_1_Mu.h5";
filename2 = "mat/h5/Mat_1_Density.h5";
};

material 11 { copy = 1; };
material 12 { copy = 1; };
material 13 { copy = 1; };
material 14 { copy = 1; };
material 15 { copy = 1; };
material 16 { copy = 1; };
material 17 { copy = 1; };
material 18 { copy = 1; };


