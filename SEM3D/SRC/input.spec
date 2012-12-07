.false.                # Acceleration scheme for Newmark
.true.                 # Velocity scheme for Newmark
800                    # Duration of the simulation
0.5                    # alpha (Newmark parameter)
0.5                    # beta (Newmark parameter)
1                      # gamma (Newmark parameter)
mesh4spec.             # input mesh file
CUB  T                 # Model (homo, prem, 3D_Berkeley or CUB) and anisotropy (T or F)
elastic_n4.hetero      # material properties file
.true.                 # save traces
.false.                # save snapshots
stations               # station location file
100                    # Time interval to skip snapshots
.true.                 # introduce a source
1                      # number of sources
29.46 65.89 6351000    # coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
2                      # Type (1 Impulse, 2 Moment Tensor)
1                      # Direction (0 x, 1 y, 2 z) (only for Impulse)
1e+18 1e+18 -2e+18     # M11 M22 M33 (only for Moment Tensor)
0e+18 0e+18 0e+18      # M12 M13 M23 (only for Moment Tensor)
3                      # Function (1 Gaussian, 2 Ricker, 3 TF(Heaviside))
300                    # tau
0.025                  # source main frequency (only for Ricker)
0.001 0.009 0.012 0.02 # source frequency band f1,f2,f3,f4 (only for TF(Heaviside))
0                      # number of solids for attenuation (0 if no attenuation)
1000  50               # attenuation period band
1                      # model period 
