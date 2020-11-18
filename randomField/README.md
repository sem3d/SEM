# randomFields

## INSTALL

Follow the guidelines in `INSTALL.md`

## USING INSIDE SEM

We suppose that we're running a test case in a folder called TEST. The location of this folder is arbitrary.

Inside TEST you will generate your mesh and then run the properties generation.
This means that you call the executable randomfield.exe from your TEST folder (RF_main_input have the informations
to the random library). This will generate the properties h5 files needed to run the simulation. They will be
located in the TEST/mat/h5 folder. 

Once properties are generated you SEM simulation.


## INPUTS

### RF_main_input

This file indicates that we are using the library with SEM and the folder of materials.
It REMAINS UNCHANGED for every SEM simulation.
Its contents are:

$application 2
$folder "./mat"

Inside input.spec:
-----------------

material {
        type = random;
};


Inside mater.in (automatic meshing) or material.input:
-----------------------------------------------------
The Domains having the random properties are identified by the leter R.
On the line after the main declaration of random subdomain there should be the statistical parameters.

Example of mater.in file of 2 random subdomains (1st and 3rd)  and 1 constant subdomain (2nd).

3
R   3535.533   2500.000   2800.000   5   0.0  0.0
1
1   20.00   20.00  20.00  2  0.5  -1
1   20.00   20.00  20.00  2  0.5  -1
1   20.00   20.00  20.00  2  0.5  -1
S   3535.533   2500.000   2800.000   5   0.0  0.0
R   3535.533   2500.000   2800.000   5   0.0  0.0
1
1   10.00   10.00  10.00  2  0.15  0
1   20.00   20.00  20.00  2  0.25  0
1   30.00   40.00  50.00  2  0.0   3

Same example inside a material.input file (the statistical parameters go all togethersequentially, after)
Obs: the comment lines (#) ARE NEEDED. In the case with PMLs the PML properties come before the Random properties

3
R   3535.533   2500.000   2800.000   5   0.0 0.0
S   3535.533   2500.000   2800.000   5   0.0 0.0
R   3535.533   2500.000   2800.000   5   0.0 0.0
# Random properties
# Parametrization Choice (0 for Kappa, 1 for Lambda)
# Rho            : corrMod, corrL_x, corrL_y, corrL_z, margiF, CV, seedStart
# Kappa or Lambda: corrMod, corrL_x, corrL_y, corrL_z, margiF, CV, seedStart
# Mu             : corrMod, corrL_x, corrL_y, corrL_z, margiF, CV, seedStart
1
1   20.00   20.00  20.00  2  0.5  -1
1   20.00   20.00  20.00  2  0.5  -1
1   20.00   20.00  20.00  2  0.5  -1
1
1   10.00   10.00  10.00  2  0.15  0
1   20.00   20.00  20.00  2  0.25  0
1   30.00   40.00  50.00  2  0.0   3

For each "R" subdomain the 4 folowing lines define the statistical and parameterization inputs
In the first line the integer informs the Parameterization Choice - PC (0 for Kappa, 1 for Lambda)
The second, third and forth lines are the statistical inputs for each property. These properties
will be:

If PC = 0
1st. Density
2nd. Kappa
3rd. Mu

If PC = 1 
1st. Density
2nd. Lambda
3rd. Mu 


For each of these properties the statistical inputs are, from left to right
Correlation Model: 1 for gaussian
Correlation Length in X
Correlation Length in Y
Correlation Length in Z
First order marginal (1 for gaussian, 2 for lognormal)
Coefficient of variation
Random seed. If it's smaller than 0 the seed will be determined by the processors clock.

Graphically it means:

0   ---------------------------------------- Parameterization Choice (PC = 0 for Kappa and PC = 1 for Lambda)
1   10.00   10.00  10.00  2  0.15  0 -- Density Statistical parameters
1   20.00   20.00  20.00  2  0.25  0 -- Kappa Statistical parameters (It would be the parameters for Lambda if PC = 1)
1   30.00   40.00  50.00  2  0.0   3 -- Mu Statistical parameters
|    |       |      |     |  |     |
|    |       |      |     |  |     |_Random Seed
|    |       |      |     |  |_Coefficient of Variation
|    |       |      |     |_First Order Marginal
|    |       |      |_Z Correlation Length
|    |       |_Y Correlation Legth
|    |_X Correlation Length
|_First Order Marginal

FAQ
===

-Can I make only one of the properties constant?
Yes, putting 0.0 to the Coefficient of Variation of this parameter.

-What does it mean the Random Seed?
Is the input for the generator of pseudo-random number generator. It defines the random sequence.

-How can I make 2 properties have a 100% correlation?
Putting the same seed AND THE SAME CORRELATION LENGTHS for these properties.

-Why do we have the "Parametrization Choice"?
In most of the cases you want to parameterize your simulation with Kappa and Mu (CP = 0), cause they are independent eigen values
of your behaviour matrix and you are sure that your problem is physically coherent. Sometimes, when thinking in terms of
V_p and V_s, working with Lambda (putting it as a constant) can lead to regimes where there is no mode conversion and this
is specially interesting when you want to validate random media propagation.
