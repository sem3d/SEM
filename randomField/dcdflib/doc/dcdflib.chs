










                                    DCDFLIB

            Library of Fortran Routines for Cumulative Distribution
                 Functions, Inverses, and Other Parameters

                               (February, 1994)






                    Summary Documentation of Each Routine








                            Compiled and Written by:

                                 Barry W. Brown
                                  James Lovato
                                  Kathy Russell








                     Department of Biomathematics, Box 237
                     The University of Texas, M.D. Anderson Cancer Center
                     1515 Holcombe Boulevard
                     Houston, TX      77030


 This work was supported by grant CA-16672 from the National Cancer Institute.


WHICH and STATUS are integer; all other arguments are Double Precision.

-------------------- DISTRIBUTION                 

WHICH   PARAMETERS     INPUT RANGE       SEARCH RANGE       REQUIREMENTS


-------------------- Beta 

      SUBROUTINE CDFBET( WHICH, P, Q, X, Y, A, B, STATUS, BOUND )

  1     P and Q        [0,1];[0,1]       -----------        SUM to 1.0
  2     X and Y        [0,1];[0,1]       [0,1],[0,1]        SUM to 1.0
  3     A              (0,infinity)      [1E-300,1E300]
  4     B              (0,infinity)      [1E-300,1E300]

-------------------- Binomial

      SUBROUTINE CDFBIN ( WHICH, P, Q, S, XN, PR, OMPR, STATUS, BOUND )

  1     P and Q        [0,1];[0,1]       -----------        SUM to 1.0
  2     S              [0,XN]            [0,XN]
  3     XN             (0,infinity)      [1E-300,1E300]
  4     PR and OMPR    [0,1];[0,1]       [0,1];[0,1]        SUM to 1.0

-------------------- Chi-square

      SUBROUTINE CDFCHI( WHICH, P, Q, X, DF, STATUS, BOUND )

  1     P and Q        [0,1],(0,1]       -----------        SUM to 1.0
  2     X              [0,infinity]      [0,1E300]
  3     DF             (0,infinity)      [1E-300,1E300]

-------------------- Noncentral Chi-square

      SUBROUTINE CDFCHN( WHICH, P, Q, X, DF, PNONC, STATUS, BOUND )

  1     P and Q        [0,1-1E-16],none  -----------
  2     X              [0,infinity]      [0,1E300]
  3     DF             (0,infinity)      [1E-300,1E300]
  4     PNONC          [0,infinity)      [0,1E4]

NOTE: We do not yet have a method to calculation the Noncentral Chi-Square
distribution acurately near 1;  therefore, Q is not used by CDFCHN.  There
are no input requirements of Q, and when WHICH is 1, Q is returned as 1-P.

-------------------- F


      SUBROUTINE CDFF( WHICH, P, Q, F, DFN, DFD, STATUS, BOUND )

  1     P and Q        [0,1],(0,1]       -----------        SUM to 1.0
  2     F              [0,infinity)      [0,1E300]
  3     DFN            (0,infinity)      [1E-300,1E300]
  4     DFD            (0,infinity)      [1E-300,1E300]

-------------------- Noncentral F

      SUBROUTINE CDFFNC( WHICH, P, Q, F, DFN, DFD, PNONC, STATUS, BOUND )

  1     P and Q        [0,1-1E-16],none  -----------
  2     F              [0,infinity)      [0,1E300]
  3     DFN            (0,infinity)      [1E-300,1E300]
  4     DFD            (0,infinity)      [1E-300,1E300]
  5     PNONC          [0,infinity)      [0,1E4]

NOTE: We do not yet have a method to calculation the Noncentral F
distribution acurately near 1;  therefore, Q is not used by CDFF. 
There are no input requirements of Q, and when WHICH is 1, Q is returned
as 1-P.

-------------------- Gamma

      SUBROUTINE CDFGAM( WHICH, P, Q, X, SHAPE, SCALE, STATUS, BOUND )

  1     P and Q        [0,1],(0,1]       -----------        SUM to 1.0
  2     X              [0,infinity)      [0,1E300]
  3     SHAPE          (0,infinity)      [1E-300,1E300]
  4     SCALE          (0,infinity)      [1E-300,1E300]

-------------------- Negative Binomial

      SUBROUTINE CDFNBN ( WHICH, P, Q, S, XN, PR, OMPR, STATUS, BOUND )


  1     P and Q        [0,1];(0,1]       -----------        SUM to 1.0
  2     S              [0,infinity)      [0,1E300]
  3     XN             [0,infinity)      [0,1E300]
  4     PR and OMPR    [0,1];[0,1]       [0,1];[0,1]        SUM to 1.0

-------------------- Normal

      SUBROUTINE CDFNOR( WHICH, P, Q, X, MEAN, SD, STATUS, BOUND )

  1     P and Q        (0,1],(0,1]       -----------        SUM to 1.0
  2     X              (-inf.,inf.)      -----------
  3     MEAN           (-inf.,inf.)      -----------
  4     SD             (0,infinity)      -----------

-------------------- Poisson

      SUBROUTINE CDFPOI( WHICH, P, Q, S, XLAM, STATUS, BOUND )

  1     P and Q        [0,1],(0,1]       -----------        SUM to 1.0
  2     S              [0,infinity)      [0,1E300]
  3     XLAM           [0,infinity)      [0,1E300]

-------------------- Student's t

      SUBROUTINE CDFT( WHICH, P, Q, T, DF, STATUS, BOUND )

  1     P and Q        (0,1],(0,1]       -----------        SUM to 1.0
  2     T              (-inf.,inf.)      [-1E300,1E300]
  3     DF             (0,infinity)      [1E-300,1E10]



