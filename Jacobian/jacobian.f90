      MODULE JACOBIAN 
       USE PARAM 
       IMPLICIT NONE

       REAL(DP) :: eps = 1D-6
       
      CONTAINS
!------------------------------------------------------------------------
         FUNCTION dfdxdy(f,xi,yi)
          !USE EQ2
          USE PARAM
          IMPLICIT NONE
          REAL(DP), INTENT(IN) :: xi,yi
          REAL(DP)             :: dfdxdy
          REAL(DP), EXTERNAL   :: f       ! function 
            
           dfdxdy = (f(xi,yi+eps)-f(xi,yi))/eps 
          
          RETURN 
          END FUNCTION
!------------------------------------------------------------------------
          FUNCTION dfdx(f,xi)
           !USE EQ1
           USE PARAM
           IMPLICIT NONE
           REAL(DP), INTENT(IN) :: xi
           REAL(DP)             :: dfdx
           REAL(DP), EXTERNAL   :: f
            
            dfdx = (f(xi+eps)-f(xi))/eps 
           
          END FUNCTION 
!------------------------------------------------------------------------
      END MODULE JACOBIAN 

