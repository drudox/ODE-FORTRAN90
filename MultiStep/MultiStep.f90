      MODULE PARAM 
!module containi
       IMPLICIT NONE
        INTEGER , PARAMETER :: SP = KIND(1.0)
        INTEGER , PARAMETER :: DP = KIND(1.D0)
      END MODULE
!-----------------------------------------------------------------------
      MODULE MultiStep 
        USE PARAM
        IMPLICIT NONE 
!-----------------------------------------------------------------------
!  @brief  This module contains the subroutine for solving ODE using 
!          general multi-step method. like the mid point, Leap-Frog and 
!          so on ...
!   
!
!  @ author Marco Ghiani 
!  @ date   Glasgow Dec 2017
!
!-----------------------------------------------------------------------


      CONTAINS 
!------------------------------------------------------------------------      
         SUBROUTINE LEAPFROG(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------
!     This subroutine perform the 2 order accurancy solution of a given
!     Differential Problem (ODE) using a multi step scheme :
!           >>> LEAP FROG Method <<< 
!
!------------------------------------------------------------------------
          IMPLICIT NONE 
          REAL(DP), DIMENSION(Ns),INTENT(INOUT) :: t,u
          REAL(DP), EXTERNAL                    :: f
          REAL(DP), INTENT(IN)                  :: t0,tf,y0,dt
          INTEGER , INTENT(IN)                  :: Ns
          CHARACTER(LEN=*) , INTENT(IN)         :: eqn 
          INTEGER                    :: i, stat
          CHARACTER(LEN=4)           :: f_ext   = '.out'  
          CHARACTER(LEN=9)           :: f_name1 = 'LeapFrog_'
          CHARACTER(LEN=17)          :: path    = 'result/multistep/'
          CHARACTER(LEN=33)          :: filename  
          REAL(DP)                   :: k1 , k2
!------------------------------------------------------------------------
          
          filename = path//f_name1//eqn//f_ext  
            
!--------- Initial Value             
          t(1)  = t0
          t(Ns) = tf
          u(1)  = y0

!---------- initiation time vector 
          
          t(2:Ns) = (/ (t0 + dt*(i-1), i=2,Ns)  /)
          
          WRITE (6,'(A)') "Running Leap Frog Solver ... "
          
!     --- ''Start-up`` the solver 
          
          k1 = f(t(1), u(1) )
          k2 = f(t(1)+dt/2. , u(1) + k1* dt/2.) ! RK 2nd Order solution
            
          u(2) = u(1) + dt*k2 



          DO i=2,Ns
            u(i+1) = u(i-1) + 2.*dt * f(t(i),u(i)) 
          END DO   
          
          WRITE (6,'(A)') "... Done "
          
          OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
            IF(stat .NE. 0)THEN 
              PRINT*, 'Error opening File: ' , filename, ': (EXIT 1)'
              CALL EXIT(1)    
            ELSE 
              WRITE(10,FMT='(2F10.5)') (t(i),u(i), i=1,Ns)
            END IF
          CLOSE(10)
          RETURN   
         END SUBROUTINE   
!------------------------------------------------------------------------
      END MODULE    
