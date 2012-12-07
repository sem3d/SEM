      SUBROUTINE WELEGL(N,ET,VN,WT)                                             
!**********************************************************************         
!   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA         
!   N  = ORDER OF THE FORMULA                                                   
!   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N                               
!   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N           
!   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N                                    
!**********************************************************************         
      DIMENSION ET(0:*), VN(0:*), WT(0:*)                                       
      IF (N .EQ. 0) RETURN                                                      
                                                                                
          N2 = (N-1)/2                                                          
          DN = DFLOAT(N)                                                        
          C  = 2.D0/(DN*(DN+1.D0))                                              
      DO 1 I=0,N2                                                               
          X = ET(I)                                                             
          Y = VN(I)                                                             
          WTX = C/(Y*Y)                                                         
          WT(I) = WTX                                                           
          WT(N-I) = WTX                                                         
1     CONTINUE                                                                  
                                                                                
      IF(N-1 .EQ. 2*N2) RETURN                                                  
          X = 0.D0                                                              
          Y = VN(N2+1)                                                          
          WT(N2+1) = C/(Y*Y)                                                    
                                                                                
      RETURN                                                                    
      END                                                                       
!                                                                               
