      MODULE PARAM
        INTEGER, PARAMETER :: DP = KIND(1.D0)
      END MODULE
!------------------------------------------------------------------------      
      MODULE eq4   ! EQUAZIONE DI RUNGE
       !USE PARAM
       IMPLICIT NONE  
       
       REAL(8),PARAMETER  :: y0 = 1./26.0
       REAL(8),PARAMETER  :: t0 = -5.D0
       REAL(8),PARAMETER  :: tF = 5.D0
       REAL(8),PARAMETER  :: dt = 0.01
       INTEGER, PARAMETER :: Ns = CEILING((tF-t0)/dt)
       CHARACTER(LEN=3)   :: eqn = '004'
       !:INTEGER            :: i, stat
!------------------------------------------------------------------------      
      CONTAINS
!------------------------------------------------------------------------      
         FUNCTION NUMSOL(a,b)
         IMPLICIT NONE
          REAL(8)              :: NumSol
          REAL(8) , INTENT(IN) :: a,b  

             NumSol = -2.*a*b*b ! FUNZIONE DI RUNGE
          
          RETURN 
         END FUNCTION   
!------------------------------------------------------------------------
         SUBROUTINE realSol()
          IMPLICIT NONE
          REAL(8), PARAMETER      :: delta = 0.08
          INTEGER , PARAMETER     :: Nr = CEILING((tF-t0)/delta) 
          REAL(8) , DIMENSION(Nr) :: x,y
          INTEGER                 :: i, stat

          x(1) = t0
          y(1) = y0

          x(2:Nr) = (/ ((x(1) + delta*(i-1)) , i=2,Nr  ) /)
          y(2:Nr) = (/ ( 1.0/(1.0+(x(i)**2. )) , i=2,Nr )  /)

       OPEN(UNIT=10,FILE='realSol_004.out',STATUS='UNKNOWN',IOSTAT=stat)
          IF(STAT .NE. 0) THEN 
             PRINT*, 'Error opening file realSol_004.out' 
             STOP
          ELSE 
             WRITE(10,FMT='(2F10.5)') ( x(i),y(i), i=1,Nr  ) 
          END IF   
       CLOSE(10)   
      END SUBROUTINE
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      END MODULE     
