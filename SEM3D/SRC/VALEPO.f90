      SUBROUTINE VALEPO(N,X,Y,DY,D2Y)                                           
!*************************************************************                  
!   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N                   
!   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT                       
!   N  = DEGREE OF THE POLYNOMIAL                                               
!   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                            
!   Y  = VALUE OF THE POLYNOMIAL IN X                                           
!   DY = VALUE OF THE FIRST DERIVATIVE IN X                                     
!   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                                    
!*************************************************************                  
                                                                                
         Y   = 1.D0                                                             
         DY  = 0.D0                                                             
         D2Y = 0.D0                                                             
      IF (N .EQ. 0) RETURN                                                      
                                                                                
         Y   = X                                                                
         DY  = 1.D0                                                             
         D2Y = 0.D0                                                             
      IF(N .EQ. 1) RETURN                                                       
                                                                                
         YP   = 1.D0                                                            
         DYP  = 0.D0                                                            
         D2YP = 0.D0                                                            
      DO 1 I=2,N                                                                
         C1 = DFLOAT(I)                                                         
         C2 = 2.D0*C1-1.D0                                                      
         C4 = C1-1.D0                                                           
         YM = Y                                                                 
         Y  = (C2*X*Y-C4*YP)/C1                                                 
         YP = YM                                                                
         DYM  = DY                                                              
         DY   = (C2*X*DY-C4*DYP+C2*YP)/C1                                       
         DYP  = DYM                                                             
         D2YM = D2Y                                                             
         D2Y  = (C2*X*D2Y-C4*D2YP+2.D0*C2*DYP)/C1                               
         D2YP = D2YM                                                            
1     CONTINUE                                                                  
                                                                                
      RETURN                                                                    
      END                                                                       
!                                                                               
