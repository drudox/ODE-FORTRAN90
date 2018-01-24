      MODULE PARAM
       INTEGER, PARAMETER :: DP = KIND(1.D0)
      END MODULE 
!------------------------------------------------------------------------            
      MODULE ADAMS_BASHFORTH
       USE PARAM 
       IMPLICIT NONE 
             
      CONTAINS
       
       SUBROUTINE AdamsBashforth2nd(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------            
!------------------------------------------------------------------------
!  @brief SUBROUTINE that compute 2nd orders Adams Bashforth 
!                    multi step solution of ODE 
!
!  @author Marco Ghiani 
!  @date Dec 2017 
!  
!  @ input t  : vector of time
!          u  : numerical result of solver
!          t0 : initial time
!          tf : finale time 
!          y0 : initial value
!          dt : step size
!          Ns : Number of step ( tf-t0/dt )
!          eqn: string to identify the output file (related to the eq.)        
!------------------------------------------------------------------------
!------------------------------------------------------------------------            
         IMPLICIT NONE
         REAL(DP) , DIMENSION(Ns),INTENT(INOUT) :: t,u
         REAL(DP) , EXTERNAL                    :: f
         REAL(DP), INTENT(IN)                   :: t0,tf,y0,dt
         INTEGER , INTENT(IN)                   :: Ns
         CHARACTER(LEN=*) , INTENT(IN)          :: eqn 
         CHARACTER(LEN=4)                       :: f_name = 'AB2_'
         CHARACTER(LEN=4)                       :: f_ext  = '.out'
         CHARACTER(LEN=13)                      :: path = 'result/adams/'
         CHARACTER(LEN=24)                      :: filename
         REAL(DP)                               :: k1,k2
         INTEGER                                :: i,stat
!------------------------------------------------------------------------            
         
         filename = path//f_name//eqn//f_ext
            
         t(1)  = t0
         t(Ns) = tF

         u(1) =  y0 
         
         t(2:Ns) = (/( (t0+ dt*(i-1)), i=2,Ns )/)    
         
         WRITE(6,FMT='(A)') "Running Adams Bashforth 2 step - 2th order"
         ! compute first point using RK 2nd order
         k1 = f(t(1),u(1)) 
         k2 = f(t(1)+dt, u(1)+k1*dt)   
         
         u(2) = u(1) + dt/2.*(k1+k2)   
            
         DO i=2,Ns
             u(i+1) = u(i) + dt/2.*(3.*f(t(i),u(i)) - f(t(i-1),u(i-1)))      
         END DO   

            
         WRITE(6,FMT='(A)') "... Done"
                  
                  
         OPEN(UNIT=10, FILE=filename, STATUS='UNKNOWN', IOSTAT=stat)
            IF(stat .NE. 0) THEN 
             PRINT*, 'Error opening file ', filename ,': (EXIT -1)'
             CALL EXIT(1)
            ELSE 
              WRITE(10,FMT='(2F20.5)') (t(i), u(i) , i=1,Ns)
            END IF
          CLOSE(10)       


        RETURN     
       END SUBROUTINE     
!-----------------------------------------------------------------------      
      
      SUBROUTINE AdamsBashforth3th(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------            
!------------------------------------------------------------------------
!  @brief SUBROUTINE that compute 3th orders Adams Bashforth 
!                    multi step solution of ODE ! 
!
!  @author Marco Ghiani 
!  @date Dec 2017 
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------            
       IMPLICIT NONE
!-------------------
       REAL(DP), DIMENSION(Ns), INTENT(INOUT) :: t,u
       REAL(DP), EXTERNAL                     :: f
       REAL(DP), INTENT(IN)                   :: t0,tf,y0,dt
       INTEGER , INTENT(IN)                   :: Ns
       CHARACTER(LEN=*) , INTENT(IN)          :: eqn 
       INTEGER                                :: I , stat
       CHARACTER(LEN=4)                       :: f_name = 'AB3_'
       CHARACTER(LEN=4)                       :: f_ext  = '.out'
       CHARACTER(LEN=13)                      :: path = 'result/adams/'
       CHARACTER(LEN=24)                      :: filename
       REAL(DP)                               :: k1,k2
!-------------------       
       
       filename = path//f_name//eqn//f_ext
      
       t(1)  = t0
       t(Ns) = tF
      
       u(1)  = y0

       t(2:Ns) = (/( (t0 + dt*(i-1)), i=2,Ns  )/)

       k1  = f(t(1),u(1))
       k2  = f(t(1)+dt, u(1)+k1*dt)
       WRITE(6,FMT='(A)') "Running Adams Bashforth 3 step - 3th order"

 ! compute first point 
       u(2) = u(1) + dt/2.*(k1+k2)

       k1  = f(t(2),u(2))
       k2  = f(t(2)+dt, u(2)+k1*dt)
       ! compute second point           
       u(3) = u(2) + dt/2.*(k1+k2)

       DO i=3,Ns
         u(i+1) = u(i) + dt/12.*(23.*f(t(i),u(i)) - 16*f(t(i-1),u(i-1))   &
     &    + 5*(f(t(i-2),u(i-2) ) ) )       
       END DO 
       
       WRITE(6,FMT='(A)') "... Done"
       
       OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
         IF(STAT .NE. 0.) THEN
            PRINT*, 'Error opening file, ' , filename , ': program terminate with EXIT 1'
            CALL EXIT(1)
         ELSE
            WRITE(10,FMT='(2F10.5)') (t(i),u(i), i=1,Ns)   
         END IF 
         CLOSE(10)
       RETURN 
      END SUBROUTINE AdamsBashforth3th


      SUBROUTINE AdamsBashforth4th(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------
!  @brief SUBROUTINE that compute 4th orders Adams Bashforth 
!                    multi step solution of ODE ! 
!
!  @author Marco Ghiani 
!  @date Dec 2017 
!
!------------------------------------------------------------------------
       IMPLICIT NONE
!-------------------
       REAL(DP), DIMENSION(Ns), INTENT(INOUT) :: t,u
       REAL(DP), EXTERNAL                     :: f
       REAL(DP), INTENT(IN)                   :: t0,tf,y0,dt
       INTEGER , INTENT(IN)                   :: Ns
       CHARACTER(LEN=*) , INTENT(IN)          :: eqn 
       INTEGER                                :: I , stat
       CHARACTER(LEN=4)                       :: f_name = 'AB4_'
       CHARACTER(LEN=4)                       :: f_ext  = '.out'
       CHARACTER(LEN=13)                      :: path   = 'result/adams/'
       CHARACTER(LEN=24)                      :: filename
       REAL(DP)                               :: k1,k2,k3,k4
!-------------------       
       
       filename = path//f_name//eqn//f_ext
      
       t(1)  = t0
       t(Ns) = tF
      
       u(1)  = y0
      
       !> compute the  vector      
       t(2:Ns) = (/( (t0 + dt*(i-1)), i=2,Ns  )/)

        ! compute first points using Runge Kutta 4th Order 
       DO i =1,3  
         k1      = f( t(i)       , u(i)          )
         k2      = f( t(i)+dt/2. , u(i)+dt/2.*k1 )
         k3      = f( t(i)+dt/2. , u(i)+dt/2.*k2 )
         k4      = f( t(i)+ dt   , u(i)+dt   *k3 )
         
         u(i+1)  = u(i) + dt/6.* (k1 + 2*k2+ 2*k3 + k4)  
       END DO  
           
       WRITE(6,FMT='(A)') "Running Adams Bashforth 4 step - 4th order"

       ! compute first point 

       DO i=4,Ns
         u(i+1) = u(i) + dt/24. *(55.*f(t(i),u(i))  &
                                 -59.*f(t(i-1),u(i-1))   &
     &                           +37.*(f(t(i-2),u(i-2) ) ) &
                                  -9.*(f(t(i-3),u(i-3) ) )  )       
       END DO 
       
       WRITE(6,FMT='(A)') "... Done"

       OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
         IF(STAT .NE. 0.) THEN
            PRINT*, 'Error opening file, ' , filename , ': program terminate with EXIT 1'
            CALL EXIT(1)
         ELSE
            WRITE(10,FMT='(2F25.3)') (t(i),u(i), i=1,Ns)   
         END IF 
         CLOSE(10)
       RETURN 
      END SUBROUTINE AdamsBashforth4th
!------------------------------------------------------------------------
!
      SUBROUTINE AdamsBashforth5th(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------
!  @brief SUBROUTINE that compute 5th orders Adams Bashforth 
!                    multi step solution of ODE  
!
!  @author Marco Ghiani 
!  @date Dec 2017 
!
!------------------------------------------------------------------------
       IMPLICIT NONE
!-------------------
       REAL(DP), DIMENSION(Ns), INTENT(INOUT) :: t,u
       REAL(DP), EXTERNAL                     :: f
       REAL(DP), INTENT(IN)                   :: t0,tf,y0,dt
       INTEGER , INTENT(IN)                   :: Ns
       CHARACTER(LEN=*) , INTENT(IN)          :: eqn 
       INTEGER                                :: I , stat
       CHARACTER(LEN=4)                       :: f_name = 'AB5_'
       CHARACTER(LEN=4)                       :: f_ext  = '.out'
       CHARACTER(LEN=13)                      :: path = 'result/adams/'
       CHARACTER(LEN=24)                      :: filename
       REAL(DP)                               :: k1,k2,k3,k4
!-------------------       
       
       filename = path//f_name//eqn//f_ext
      
       t(1)  = t0
       t(Ns) = tF
      
       u(1)  = y0
      

       !> compute the  vector      
       t(2:Ns) = (/( (t0 + dt*(i-1)), i=2,Ns  )/)
      

        ! compute the first 4 points using Runge Kutta 4th Order 
       DO i =1,4  
         k1      = f( t(i)       , u(i)          )
         k2      = f( t(i)+dt/2. , u(i)+dt/2.*k1 )
         k3      = f( t(i)+dt/2. , u(i)+dt/2.*k2 )
         k4      = f( t(i)+ dt   , u(i)+dt   *k3 )
         
         u(i+1)  = u(i) + dt/6.* (k1 + 2*k2+ 2*k3 + k4)  
       END DO  
           
       WRITE(6,FMT='(A)') "Running Adams Bashforth 5 step - 5th order"

       ! compute first point 

       DO i=5,Ns
         u(i+1) = u(i) + dt  *(1901.0/720.0  *f(t(i),u(i))          &
     &                         -1387.0/360.0 *f(t(i-1),u(i-1))      &
     &                         +109.0/30.0   *(f(t(i-2),u(i-2) ) )  &
     &                         -637.0/360.0  *(f(t(i-3),u(i-3) ) )  &
     &                         +251.0/720.0  *(f(t(i-4),u(i-4) )))       
       END DO 
       
       WRITE(6,FMT='(A)') "... Done"

       OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
         IF(STAT .NE. 0.) THEN
            PRINT*, 'Error opening file, ' , filename , ': program terminate with EXIT 1'
            CALL EXIT(1)
         ELSE
            WRITE(10,FMT='(2F10.5)') (t(i),u(i), i=1,Ns)   
         END IF 
         CLOSE(10)
       RETURN 
      END SUBROUTINE AdamsBashforth5th

!------------------------------------------------------------------------            
      END MODULE
