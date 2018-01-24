      MODULE eq3
       
       IMPLICIT NONE
        REAL(8),PARAMETER :: t0  = -2.D0
        REAL(8),PARAMETER :: tF  = 2.D0
        REAL(8),PARAMETER :: y0  = exp(2.D0)
        REAL(8),PARAMETER :: dt  = 0.05
        INTEGER ,PARAMETER :: Ns  = CEILING((TF-T0)/dt)
        CHARACTER(LEN=3)   :: eqn = '003'


      CONTAINS
!-----------------------------------------------------------------------
      FUNCTION NumSol(a,b) 
       !USE PARAM
       IMPLICIT NONE 
       REAL(8)    :: NumSol
       REAL(8), INTENT(IN) :: a,b
       
         NumSol = a*b 

       RETURN 
       END FUNCTION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE realSol()
       IMPLICIT NONE
       REAL(8), PARAMETER     :: delta=0.08
       INTEGER, PARAMETER     :: Nr = CEILING((tF-t0)/delta)     
       REAL(8), DIMENSION(Nr) :: t,y
       INTEGER :: i, stat
       t(1)   = t0
       y(1)   = y0
       t(2:Nr) = (/ ( (t0 + delta*(i-1)) , i=2,Nr )  /)
       y(2:Nr) = (/ ( (exp((t(i)**2.d0)/2.d0) ), i=2,Nr )/)
       !DO i=1,N
       ! y(i) =  
       !END DO 
      OPEN(UNIT=10,FILE='realSol_003.out',STATUS='UNKNOWN',IOSTAT=stat) 
            IF(stat .NE. 0) THEN 
                  Print*, 'Error opening File!'
                  CALL EXIT(1)
            ELSE    
                   WRITE(10,FMT='(2F10.5)') ( t(i), y(i) , i=1,Nr )
            END IF
      CLOSE(10)
      RETURN 
      END SUBROUTINE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      END MODULE 
