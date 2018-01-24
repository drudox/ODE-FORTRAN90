      MODULE PARAM 
       IMPLICIT NONE
       INTEGER , PARAMETER :: DP = KIND(1.D0)
       INTEGER , PARAMETER :: SP = KIND(1.0)
      END MODULE

      MODULE RUNGEKUTTA 
!------------------------------------------------------------------------ 
! @brief 
! Modules containg the commons RK methods (non stiff ODE)
! RK2 (modified - Euler / Heun)
! RK4th (classic Runge Kutta 4th order)
! RK-fehlberg (4-5 order)
! RK-Merson
!------------------------------------------------------------------------
       USE PARAM 
       !USE EQ
       IMPLICIT NONE
       
      CONTAINS
!------------------------------------------------------------------------
        SUBROUTINE RungeKutta4th(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
         IMPLICIT NONE 
         REAL(DP), DIMENSION(Ns),INTENT(INOUT) :: t, u
         REAL(DP), EXTERNAL                    :: f   
         REAL(DP), INTENT(IN)                  :: t0,tf,y0,dt
         INTEGER , INTENT(IN)                  :: Ns
         CHARACTER(LEN=*) , INTENT(IN)         :: eqn 
         CHARACTER(LEN=4)                      :: f_name = 'RK4_' 
         CHARACTER(LEN=4)                      :: f_ext  = '.out'
         CHARACTER(LEN=30)                     :: filename
         CHARACTER(LEN=18)                     :: path = 'result/rungekutta/'
         REAL(DP)                              :: k1, k2, k3, k4   
         INTEGER                               :: i,  stat
!-----------------------------------------------------------------------
         filename = path//f_name//eqn//f_ext
         
         t(1)  = t0
         t(Ns) = tF
         u(1)  = y0
            

          WRITE(6,FMT='(A)') "Running Runge Kutta 4th order (std RK) "
                    
         t(2:Ns) = (/( (t0 + dt*(i-1)) , i=2,Ns )/)   
         
        DO i=1,Ns                  ! computing RK coef. and y solution  
            k1  = f(t(i)      , u(i))
            k2  = f(t(i)+dt/2 , u(i)+dt/2*k1 )
            k3  = f(t(i)+dt/2 , u(i)+dt/2*k2 )
            k4  = f(t(i)+dt   , u(i)+dt*k3   )
         u(i+1) = u(i) + dt/6.*(k1 + 2*k2 + 2*k3 + k4)
        END DO
         
         OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
            IF(stat .NE. 0) THEN 
             PRINT*, 'Error opening file ', filename ,': (EXIT 1)'
            ELSE 
              WRITE(10,FMT='(2F10.5)') (t(i), u(i) , i=1,Ns)
            END IF
          CLOSE(10)   
          
          WRITE(6,FMT='(A)') "... Done"
         RETURN 
         END SUBROUTINE 
 
!-----------------------------------------------------------------------      
      SUBROUTINE HEUN(t,u,f,t0,tf,y0,dt,Ns,eqn)
!-----------------------------------------------------------------------         
!        
!
!-----------------------------------------------------------------------         
         IMPLICIT NONE 
         REAL(DP), DIMENSION(Ns)        :: t,u
         REAL(DP), EXTERNAL             :: f   
         REAL(DP), INTENT(IN)           :: t0,tf,y0,dt
         INTEGER , INTENT(IN)           :: Ns
         CHARACTER(LEN=*) , INTENT(IN)  :: eqn 
         CHARACTER(LEN=5)               :: f_name = 'Heun_' 
         CHARACTER(LEN=4)               :: f_ext  = '.out'
         CHARACTER(LEN=30)              :: filename
         CHARACTER(LEN=18)              :: path = 'result/rungekutta/'
         REAL(DP)                       :: k1, k2  
         INTEGER                        :: i,  stat
            
         filename = path//f_name//eqn//f_ext
         
         t(1)  = t0
         t(Ns) = tF
         u(1)  = y0
      
         t(2:Ns) = (/( (t0 + dt*(i-1)) , i=2,Ns )/)   
        
         WRITE(6,FMT='(A)') "Running Heun Solver (RungeKutta 2th order) "
 
         DO i=1,Ns
              k1 = f(t(i),u(i))
              k2 = f(t(i)+dt , u(i) + dt*k1)
          u(i+1) = u(i) + dt/2*(k1+k2) 
         END DO
       
         OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
            IF(stat .NE. 0) THEN 
             PRINT*, 'Error opening file ', filename ,': (EXIT 1)'
            ELSE 
              WRITE(10,FMT='(2F10.5)') (t(i), u(i) , i=1,Ns)
            END IF
          CLOSE(10)   
          
          WRITE(6,FMT='(A)') "... Done"

         RETURN 
        END SUBROUTINE     
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
        SUBROUTINE ModifiedEuler(t,u,f,t0,tf,y0,dt,Ns,eqn)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
         IMPLICIT NONE 
         REAL(DP), DIMENSION(Ns),INTENT(INOUT) :: u,t
         REAL(DP), EXTERNAL                    :: f   
         REAL(DP), INTENT(IN)                  :: t0,tf,y0,dt
         INTEGER , INTENT(IN)                  :: Ns
         CHARACTER(LEN=*) , INTENT(IN)         :: eqn 
         CHARACTER(LEN=5)                      :: f_name = 'MEul_' 
         CHARACTER(LEN=18)                     :: path = 'result/rungekutta/'
         CHARACTER(LEN=4)                      :: f_ext  = '.out'
         CHARACTER(LEN=30)                     :: filename
         REAL(DP)                              :: k1, k2  
         INTEGER                               :: i,  stat

         filename = path//f_name//eqn//f_ext
         
         t(1)  = t0
         t(Ns) = tF
         u(1)  = y0
      
         t(2:Ns) = (/( (t0 + dt*(i-1)) , i=2,Ns )/)   
         
         WRITE(6,FMT='(A)') "Running Modified Euler Solver (RungeKutta 2th order) "
         
         DO i=1,Ns
          k1 = f(t(i),u(i))
          k2 = f(t(i)+dt/2 , u(i)+dt/2*k1)
          u(i+1) = u(i) + dt*k2
         END DO
        
          OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
            IF(stat .NE. 0) THEN 
             PRINT*, 'Error opening file ', filename ,': (EXIT 1)'
            ELSE 
              WRITE(10,FMT='(2F10.5)') (t(i), u(i) , i=1,Ns)
            END IF
          CLOSE(10)   

          WRITE(6,FMT='(A)') "... Done"

         RETURN 
        END SUBROUTINE 

!------------------------------------------------------------------------         
      SUBROUTINE RungeKuttaMerson(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------         
!
!
!
!------------------------------------------------------------------------         
        IMPlICIT NONE
        REAL(DP), DIMENSION(Ns), INTENT(INOUT) :: t,u
        REAL(DP), EXTERNAL                     :: f
        REAL(DP), INTENT(IN)                   :: t0,tf,y0,dt
        INTEGER , INTENT(IN)                   :: Ns
        CHARACTER(LEN=*) , INTENT(IN)          :: eqn 
        REAL(DP)                               :: k1,k2,k3,k4,k5
        INTEGER                                :: stat, i
        CHARACTER(LEN=9)                       :: f_name = 'RKMerson_'
        CHARACTER(LEN=4)                       :: f_ext = '.out'
        CHARACTER(LEN=18)                      :: path = 'result/rungekutta/'
        CHARACTER(LEN=35)                      :: filename 
      !  REAL(DP)                               :: 
        
        filename = path//f_name//eqn//f_ext

        t(1) = t0 
        t(Ns)= tf
        u(1) = y0

        t(2:Ns) = (/( (t0 + dt*(i-1)), i=2,Ns) /)
         
        WRITE(6,FMT='(A)') "Running RungeKutta - Merson 5th order "
        
        DO i=1,Ns
         k1 = dt * f( t(i)       , u(i) )
         k2 = dt * f( t(i)+dt/3. , u(i)+k1/3. )
         k3 = dt * f( t(i)+dt/3. , u(i)+1./6.*(k1+k2) )
         k4 = dt * f( t(i)+dt/2. , u(i)+1./8.*(k1+3.*k3) )
         k5 = dt * f( t(i)+dt    , u(i)+1./2.*(k1-3.*k3+4*k4))

         u(i+1) = u(i) + 1./6. * (k1+ 4*k4 + k5)
        END DO 



          OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
            IF(stat .NE. 0) THEN 
             PRINT*, 'Error opening file ', filename ,': (EXIT 1)'
            ELSE 
              WRITE(10,FMT='(2F10.5)') (t(i), u(i) , i=1,Ns)
            END IF
          CLOSE(10)   

          WRITE(6,FMT='(A)') "... Done"

         RETURN 
        END SUBROUTINE    
!------------------------------------------------------------------------         
      SUBROUTINE RungeKuttaFehlberg45(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------         
!    Adaptative step solver - Runge Kutta Fehlberg 4th to 5th order         
!    evalutate the truncation error each step (using 4th order compared
!    with 5th order solution ) and increase or decrease the step size 
!
!------------------------------------------------------------------------         
       IMPLICIT NONE 
       REAL(DP), DIMENSION(Ns*4),INTENT(INOUT) :: u,t
       REAL(DP)                              :: u1, u2
       REAL(DP), EXTERNAL                    :: f   
       REAL(DP), INTENT(IN)                  :: t0,tf,y0,dt
       REAL(DP)                              :: h                     ! h = dt
       INTEGER , INTENT(IN)                  :: Ns
       CHARACTER(LEN=*) , INTENT(IN)         :: eqn 
       CHARACTER(LEN=11)                     :: f_name = 'RKFehlberg_' 
       CHARACTER(LEN=18)                     :: path = 'result/rungekutta/'
       CHARACTER(LEN=4)                      :: f_ext  = '.out'
       CHARACTER(LEN=38)                     :: filename
       REAL(DP)                              :: k1, k2, k3, k4, k5, k6 
       INTEGER                               :: i,  stat, Nstep =1 
       REAL(DP)                              :: R , delta , s 
       REAL(DP), PARAMETER                   :: toll=1.0e-12 , eps=1.0d-3


        filename = path//f_name//eqn//f_ext

        t(1) = t0 
        !t(Ns)= tf
        u(1) = y0
        u1   = y0
        u2   = y0
        i = 1
        h=dt

       ! t(2:Ns) = (/( (t0 + dt*(i-1)), i=2,Ns) /)
         
        WRITE(6,FMT='(A)') "Running RungeKutta - Fehlberg (adaptative) 4-5th order "
            
       ! Nstep = 1
       !  print*, Nstep
        DO WHILE( t(i) .LT. tf )
         k1 = h* f(t(i)          , u(i)                               )
         k2 = h* f(t(i)+h/4     ,  u(i)+ k1/4                         )
         k3 = h* f(t(i)+3./8.*h ,  u(i)+ 3./32.*k1 + 9./32.*k2        )
         k4 = h* f(t(i)+12./13.*h, u(i)+ 1932./2197.*k1- 7200./2197.*k2&
    &                                  + 7296./2197.*k3                )
         k5 = h* f(t(i)+ h      ,  u(i)+ 439./216*k1 - 8.*k2           &
    &                                  + 3680./513.* k3 - 845./4104.*k4)
         k6 = h* f(t(i)+h/2     ,  u(i)- 8./27.*k1+ 2.*k2-3544./2565*k3&
    &                                  + 1859./4104.*k4 - 11./40.*k5   )
        
        ! 4th order solution
         u1 = u(i)+( 25./216.*k1 + 1408./2565.*k3 +2197./4101*k4 -1./5*k5)

        ! 5th order solution
         u2 = u(i)+(16./135.*k1+6656./12825.*k3 + 28561./56430.*k4   &
    &             - 9./50.*k5 + 2./55.*k6  )

          R = ABS(u1-u2)/h
          s = 0.8409*((eps/R)**(1./4.))
          
        !  print*, Nstep , i

          IF(R .LT. eps ) THEN
            t(i+1) = t(i)+h
            u(i+1) = u1
            h      = s * h 
            i      = i + 1
            Nstep  = Nstep + 1
          ELSE 
            h = s*h
          END IF 
            
        END DO
 
 
         OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
            IF(stat .NE. 0) THEN 
             PRINT*, 'Error opening file ', filename ,': (EXIT 1)'
            ELSE 
              WRITE(10,FMT='(2F10.5)') (t(i), u(i) , i=1,Nstep)
            END IF
          CLOSE(10)   

          WRITE(6,FMT='(A)') "... Done"


        
         RETURN 
        END SUBROUTINE    
!------------------------------------------------------------------------         


        SUBROUTINE CrankNicholson(t,u,f,t0,tf,y0,dt,Ns,eqn)
!-----------------------------------------------------------------------
!   Perform the solution of a given IVProblem using:
!   PREDICTOR (fwd Euler) - implicit CORRECTOR scheme (CRANK-NICHOLSON)
!   
!-----------------------------------------------------------------------
         IMPLICIT NONE 
         REAL(DP), DIMENSION(Ns),INTENT(INOUT) :: u,t
         REAL(DP), DIMENSION(NS)               :: up  
         REAL(DP), EXTERNAL                    :: f   
         REAL(DP), INTENT(IN)                  :: t0,tf,y0,dt
         INTEGER , INTENT(IN)                  :: Ns
         CHARACTER(LEN=*) , INTENT(IN)         :: eqn 
         CHARACTER(LEN=10)                     :: f_name = 'CrankNich_' 
         CHARACTER(LEN=18)                     :: path = 'result/rungekutta/'
         CHARACTER(LEN=4)                      :: f_ext  = '.out'
         CHARACTER(LEN=35)                     :: filename
         REAL(DP)                              :: k1, k2  
         INTEGER                               :: i,  stat

         filename = path//f_name//eqn//f_ext
         
         !---------    set initial Value
         t(1)  = t0
         t(Ns) = tF
         u(1)  = y0
         up(1) = y0

         t(2:Ns) = (/( (t0 + dt*(i-1)) , i=2,Ns )/)   
         
         WRITE(6,FMT='(A)') "Running Crank-Nicolson (Predictor - Corrector) Solver (RK 2th order) "
         
            

         DO i=1,Ns
          
          u(i+1) = u(i) + dt*f(t(i),u(i))   ! PREDICTOR: FORWARD EULER 

          u(i+1) = u(i) + dt/2. * (f(t(i),u(i)) + f(t(i)+dt , u(i+1) ) ) ! CORRECTOR: CRANK NICHOLSON

         END DO
        
          OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
            IF(stat .NE. 0) THEN 
             PRINT*, 'Error opening file ', filename ,': (EXIT 1)'
            ELSE 
              WRITE(10,FMT='(2F10.5)') (t(i), u(i) , i=1,Ns)
            END IF
          CLOSE(10)   

          WRITE(6,FMT='(A)') "... Done"

         RETURN 

        END SUBROUTINE 
!------------------------------------------------------------------------         


!------------------------------------------------------------------------         
      END MODULE

