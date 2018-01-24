      MODULE PARAM 
        IMPLICIT NONE
        INTEGER , PARAMETER :: SP = KIND(1.0)
        INTEGER , PARAMETER :: DP = KIND(1.D0)
      END MODULE
!------------------------------------------------------------------------
      MODULE Euler 
        USE PARAM
        !USE eq
        IMPLICIT NONE 

      CONTAINS 
!------------------------------------------------------------------------      
         SUBROUTINE BWDEULER(t,u,f,df,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------
!   this subroutine perform the numerical solution of a ODE using   
!   implicit (backward Euler) method . 
!   For solving the non linear algebric equations the Newton-Raphson 
!   method is used.       
! 
!------------------------------------------------------------------------
           IMPLICIT NONE
           
           REAL(DP), DIMENSION(Ns), INTENT(INOUT) :: t,u
           REAL(DP), EXTERNAL              :: f,df
           REAL(DP), PARAMETER             :: toll=1D-12
           REAL(DP), INTENT(IN)            :: t0,tf,y0,dt
           INTEGER , INTENT(IN)            :: Ns
           CHARACTER(LEN=*) , INTENT(IN)   :: eqn 
           REAL(DP)                        :: uOLD, uNEW, error
           INTEGER                         :: i,  stat
           CHARACTER(LEN=9)                :: f_name = 'bwdEuler_'
           CHARACTER(LEN=4)                :: f_ext = '.out'
           CHARACTER(LEN=13)               :: path = 'result/euler/'
           CHARACTER(LEN=30)               :: filename

!----------------------------------------------------------------------- 
           
          ! IF( .NOT. PRESENT(h)) THEN 
          !   h = CEILING((tf-t0)/h)
          ! END IF  

           filename = path//f_name//eqn//f_ext
           
!--------- INITIAL CONDITION             
           
           t(1) = t0
           t(Ns) = tf
           u(1) = y0
           
           t(2:Ns) = (/( t0 + dt*(i-1), i=2,Ns  )/)
!------------
           WRITE (6,'(A)') "Running Backward Euler Solver ..."
           DO I=1,Ns
             uOld = u(i) + f(t(i+1),u(i))
             error=1.D0
             DO WHILE(error .GT. toll)
               
               uNew  = uOld - (uOld - (u(i) + dt* f(t(i+1),uOld) )) /  &
     &                       (1- dt*df(f,t(i+1),uOld)) 
               
               error = ABS(uNew-uOld)
               uOld = uNew   
             END DO
             u(i+1) = uNew
           END DO       
           WRITE (6,'(A)') "... Done "

!------------           
           
           OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
            IF(stat .NE. 0) THEN 
              PRINT*, 'Error opening file ', filename , ': (EXIT 1)'
              CALL EXIT(1)
            ELSE 
                WRITE(10,'(2f10.5)') (t(i),u(i) , i=1,Ns)
            END IF
            
            CLOSE(10)
           RETURN   
         END SUBROUTINE 

      SUBROUTINE FWDEULER(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------
!     @brief FwdEuler Subroutine :
!      Compute the solution of a non stiff ODE using explicit Euler (FWD)
!      
!      The result will be write into a file ( the variable filename 
!      contains the name of the file )
!------------------------------------------------------------------------
         IMPLICIT NONE 
         REAL(DP), DIMENSION(Ns),INTENT(INOUT) :: t,u
         REAL(DP), EXTERNAL                    :: f
         REAL(DP), INTENT(IN)                  :: t0,tf,y0,dt
         INTEGER , INTENT(IN)                  :: Ns
         CHARACTER(LEN=*) , INTENT(IN)         :: eqn 
         INTEGER                               :: i, stat
         CHARACTER(LEN=4)                      :: f_ext = '.out'
         CHARACTER(LEN=9)                      :: fname = 'fwdEuler_' 
         CHARACTER(LEN=13)                     :: path  = 'result/euler/'
         CHARACTER(LEN=29)                     :: filename  
!------------------------------------------------------------------------
          
          filename = path//fname//eqn//f_ext  
            
!--------- Initial Value             
          t(1)  = t0
          t(Ns) = tf
          u(1)  = y0

!---------- VECTOR U CONSTRUCTION 
          
          t(2:Ns) = (/ (t0 + dt*(i-1), i=2,Ns)  /)
          
          WRITE (6,'(A)') "Running Fordward Euler Solver ... "

          DO i=1,Ns
            u(i+1) = u(i) + dt* f(t(i),u(i)) 
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
!------------------------------------------------------------------------
      END MODULE    
