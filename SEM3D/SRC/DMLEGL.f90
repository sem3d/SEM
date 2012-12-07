      SUBROUTINE DMLEGL(N,NM,ET,VN,DMA)                                         
!***********************************************************************        
!  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE                
!  LEGENDRE GAUSS-LOBATTO NODES                                                 
!  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX                      
!  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT        
!  ET  = VECTOR OF THE NODES, ET(I), I=0,N                                      
!  VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N           
!  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N                              
!***********************************************************************        
      DIMENSION ET(0:*), VN(0:*), DMA(0:NM,0:*)                                 
          DMA(0,0) = 0.D0                                                       
      IF (N .EQ. 0) RETURN                                                      
                                                                                
      DO 1 I=0,N                                                                
          VI = VN(I)                                                            
          EI = ET(I)                                                            
      DO 2 J=0,N                                                                
      IF (I .NE. J) THEN                                                        
          VJ = VN(J)                                                            
          EJ = ET(J)                                                            
          DMA(I,J) = VI/(VJ*(EI-EJ))                                            
      ELSE                                                                      
          DMA(I,I) = 0.D0                                                       
      ENDIF                                                                     
2     CONTINUE                                                                  
1     CONTINUE                                                                  
                                                                                
          DN = DFLOAT(N)                                                        
          C  = .25D0*DN*(DN+1.D0)                                               
          DMA(0,0) = -C                                                         
          DMA(N,N) = C                                                          
                                                                                
      RETURN                                                                    
      END                                                                       
!                                                                               
