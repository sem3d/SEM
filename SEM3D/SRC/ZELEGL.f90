      SUBROUTINE ZELEGL(N,ET,VN)                                                
!********************************************************************           
!   COMPUTES THE NODES RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA           
!   N  = ORDER OF THE FORMULA                                                   
!   ET = VECTOR OF THE NODES, ET(I), I=0,N                                      
!   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N           
!********************************************************************           
      DIMENSION ET(0:*), VN(0:*)                                                
      IF (N .EQ. 0) RETURN                                                      
                                                                                
         N2 = (N-1)/2                                                           
         SN = DFLOAT(2*N-4*N2-3)                                                
         ET(0) = -1.D0                                                          
         ET(N) = 1.D0                                                           
         VN(0) = SN                                                             
         VN(N) = 1.D0                                                           
      IF (N .EQ. 1) RETURN                                                      
                                                                                
         ET(N2+1) = 0.D0                                                        
         X = 0.D0                                                               
      CALL VALEPO(N,X,Y,DY,D2Y)                                                 
         VN(N2+1) = Y                                                           
      IF(N .EQ. 2) RETURN                                                       
                                                                                
         PI = 3.14159265358979323846D0                                          
         C  = PI/DFLOAT(N)                                                      
      DO 1 I=1,N2                                                               
         ETX = DCOS(C*DFLOAT(I))                                                
      DO 2 IT=1,8                                                               
      CALL VALEPO(N,ETX,Y,DY,D2Y)                                               
         ETX = ETX-DY/D2Y                                                       
2     CONTINUE                                                                  
         ET(I) = -ETX                                                           
         ET(N-I) = ETX                                                          
         VN(I) = Y*SN                                                           
         VN(N-I) = Y                                                            
1     CONTINUE                                                                  
                                                                                
      RETURN                                                                    
      END                                                                       
!                                                                               
