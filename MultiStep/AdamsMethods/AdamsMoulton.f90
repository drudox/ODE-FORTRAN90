      MODULE PARAM
       IMPLICIT NONE
       INTEGER, PARAMETER :: DP = KIND(1.D0)
       INTEGER, PARAMETER :: SP = KIND(1.0)
      END MODULE 
!------------------------------------------------------------------------      
      MODULE adams_moulton
       USE PARAM
       IMPLICIT NONE
      
      CONTAINS 
!------------------------------------------------------------------------
        SUBROUTINE ADAMSMOULTON2nd(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
         IMPLICIT NONE

         REAL(DP), DIMENSION(Ns), INTENT(INOUT) :: t,u 
         REAL(DP), EXTERNAL                     :: f
         REAL(DP), INTENT(IN)                   :: t0,tf,y0,dt
         INTEGER , INTENT(IN)                   :: Ns
         CHARACTER(LEN=*) , INTENT(IN)          :: eqn 
         REAL(DP)                :: fpred , fold, up , uOld, uNew 
         REAL(DP)                :: error , fnew , k1,k2
         REAL(DP), PARAMETER     :: toll = 1D-15
         INTEGER                 :: i, j , stat
         CHARACTER(LEN=13)      :: path = 'result/adams/'
         CHARACTER(LEN=4) :: f_name   = 'AM2_'
         CHARACTER(LEN=4) :: f_ext    = '.out'
         CHARACTER(LEN=24):: filename 

!-----------         
         filename = path//f_name//eqn//f_ext
         
         t(1)  = t0
         t(Ns) = tF
         
         u(1)  = y0

         t(2:Ns) = (/ ((t0 + dt*(i-1)), i=2,Ns)  /)

      ! compute first point 
         k1    = f(t(1) ,u(1))
         k2    = f(t(1) , u(1)+dt*k1)
         
         u(2)  = u(1) + dt/2.*(k1+k2)  

!>      Predictor  ADAMS BASHFORTH 2ND
      
         DO I=2,Ns

           up = u(i) + dt/2.*(3.*f(t(i),u(i)) - f(t(i-1),u(i-1))) 
 
           fpred  = f(t(i+1),up) 

            
!>      CORRECTOR ADAMS MOULTON 2ND
           error = 1.D0
           fold = fpred
           j=0
           DO WHILE(error .GT. toll)
            uOld = u(i) + dt/2.*( 1.*fOld + f(t(i),u(i)) )
            
            fOld = f(t(i+1),uOld)
            
            uNew = u(i) + dt/2.* ( 1.*fOld + f(t(i),u(i) ) ) 
            
            fNew = f(t(i+1),uNew)

            error = abs(uNew-uOld)
            fOld = fNew
            j = j+1  
           ! print*, 'iter AM2: ',j , 'error ' , error
           END DO 
           u(i+1) = uNew
         END DO

         OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
            IF( STAT .NE. 0) THEN 
               PRINT*, 'Error opening File ', filename, ' program terminate with EXIT 1'
               CALL EXIT(1)   
            ELSE 
               WRITE(10,FMT='(2F10.5)') (t(i),u(i), i=1,Ns)
            END IF
         CLOSE(10)
      
      END SUBROUTINE
!------------------------------------------------------------------------      

      SUBROUTINE ADAMSMOLULTON3th(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------      
!
!
!------------------------------------------------------------------------      

       IMPLICIT NONE
       REAL(DP),DIMENSION(Ns),INTENT(INOUT) :: t, u 
       REAL(DP),EXTERNAL                    :: f
       REAL(DP), INTENT(IN)                   :: t0,tf,y0,dt
       INTEGER , INTENT(IN)                   :: Ns
       CHARACTER(LEN=*) , INTENT(IN)          :: eqn 
       REAL(DP)               :: error, up ,uOld , uNew, toll=1D-15
       REAL(DP)               :: fpred , fNew ,fOld
       REAL(DP)               :: k1,k2 !,k3,k4
       INTEGER                :: i,stat
       CHARACTER(LEN=13)      :: path = 'result/adams/'
       CHARACTER(LEN=4)       :: f_name= 'AM3_'
       CHARACTER(LEN=4)       :: f_ext = '.out'
       CHARACTER(LEN=24)      :: filename
!----------       
       filename = path//f_name//eqn//f_ext

       t(1) = t0
       t(Ns)= tF
      
       u(1) = y0
       
       ! compute first points using Runge Kutta 2nd Order accuracy
       DO i =1,2  
         k1      = f(t(i)       , u(i)          )
         k2      = f(t(i)+dt    , u(i)+dt   *k1 )
         u(i+1)  = u(i) + dt/2.* (k1 + k2)  
       END DO  
           

       t(2:Ns) = (/( (t(1) + dt*(i-1)), i=2,Ns )/)
!
!>          PREDICTOR ADAMS BASHFORTH  3th order            
       DO i=3,Ns
        up = u(i) + dt/12.*(23. * f(t(i),u(i))      &
     &                    - 16. * f(t(i-1),u(i-1))  &
     &                     + 5. * f(t(i-2),u(i-2))   )   

        fpred = f(t(i+1),up)
        error = 1.D0
        fold = fpred
        !j=0
!>           CORRECTOR ADAMS MOULTON 3th ORDER
        DO WHILE( error .GT. toll )
         
         uOld = u(i) + dt/12. * ( 5.* fold              &
     &                          + 8.* f( t(i),u(i)   )  & 
     &                          - 1.* f(t(i-1),u(i-1)) )
         
         fold = f(t(i+1),uOld) 

         uNew = u(i) + dt/12. *( 5 * fold              &  
     &                         + 8 * f ( t(i),u(i) )   & 
     &                         - 1.* f(t(i-1),u(i-1))  )
         
         fnew = f(t(i+1),uNew )
         
         error = abs(uNew-uOld)
         
         fold = fnew

        END DO 
        u(i+1) = unew
       END DO
      
      OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
         IF (stat .NE. 0) THEN 
          PRINT*, 'Error opening file ' , filename, ': (EXIT 1)' 
         ELSE 
          WRITE(10,FMT='(2F10.5)') ( t(i),u(i), i=1,Ns)
         END IF 
      RETURN 
      END SUBROUTINE ADAMSMOLULTON3th

!------------------------------------------------------------------------      
      SUBROUTINE ADAMSMOLULTON4th(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------      
!
!
!------------------------------------------------------------------------      
       IMPLICIT NONE
       REAL(DP),DIMENSION(Ns),INTENT(INOUT) :: t, u 
       REAL(DP),EXTERNAL                    :: f
       REAL(DP), INTENT(IN)                 :: t0,tf,y0,dt
       INTEGER , INTENT(IN)                 :: Ns
       CHARACTER(LEN=*) , INTENT(IN)        :: eqn 
       REAL(DP)               :: error, up ,uOld , uNew, toll=1D-15
       REAL(DP)               :: fpred , fNew ,fOld
       REAL(DP)               :: k1,k2,k3,k4
       INTEGER                :: i,stat
       CHARACTER(LEN=13)      :: path = 'result/adams/'
       CHARACTER(LEN=4)       :: f_name= 'AM4_'
       CHARACTER(LEN=4)       :: f_ext = '.out'
       CHARACTER(LEN=24)      :: filename
!----------       
       filename = path//f_name//eqn//f_ext

       t(1) = t0
       t(Ns)= tF
      
       u(1) = y0
       
       ! compute first points using Runge Kutta 4th Order 
       DO i =1,3  
         k1      = f( t(i)       , u(i)          )
         k2      = f( t(i)+dt/2. , u(i)+dt/2.*k1 )
         k3      = f( t(i)+dt/2. , u(i)+dt/2.*k2 )
         k4      = f( t(i)+ dt   , u(i)+dt   *k3 )
         
         u(i+1)  = u(i) + dt/6.* (k1 + 2*k2+ 2*k3 + k4)  
       END DO  
           

       t(2:Ns) = (/( (t(1) + dt*(i-1)), i=2,Ns )/)    ! time vector 
      
      WRITE(6,'(A)') "Running Adams Bashforth (Predictor) , Adams Moulton 4th order (Corrector)"

!
!>          PREDICTOR ADAMS BASHFORTH  4th order            
       DO i=4,Ns
        up = u(i) + dt/24.*(55. * f(t(i),u(i))      &
     &                    - 59. * f(t(i-1),u(i-1))  &
     &                    + 37. * f(t(i-2),u(i-2))  & 
     &                    +  9. * f(t(i-3),u(i-3))  )   

        fpred = f(t(i+1),up)
        error = 1.D0
        fold = fpred
        !j=0
!>           CORRECTOR ADAMS MOULTON 4th ORDER
        DO WHILE( error .GT. toll )
         
         uOld = u(i) + dt/24. * ( 9.* fold              &
     &                          +19.* f( t(i),u(i)   )  & 
     &                          - 5.* f(t(i-1),u(i-1))  &
     &                          + 1. *f(t(i-2),u(i-2))  )
         
         fold = f(t(i+1),uOld) 

         uNew = u(i) + dt/24. *(  9.* fold              &
     &                          +19.* f( t(i),u(i)   )  & 
     &                          - 5.* f(t(i-1),u(i-1))  &
     &                          + 1. *f(t(i-2),u(i-2))  )


         fnew = f(t(i+1),uNew )
         
         error = abs(uNew-uOld)
         
         fold = fnew

        END DO 
        u(i+1) = unew
       END DO
      WRITE(6,'(A)') "... Done"

      OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
         IF (stat .NE. 0) THEN 
          PRINT*, 'Error opening file ' , filename, ': (EXIT 1)' 
         ELSE 
          WRITE(10,FMT='(2F10.5)') (t(i),u(i), i=1,Ns)
         END IF 
      RETURN 
      END SUBROUTINE ADAMSMOLULTON4th
!------------------------------------------------------------------------      
      SUBROUTINE ADAMSMOLULTON5th(t,u,f,t0,tf,y0,dt,Ns,eqn)
!------------------------------------------------------------------------      
!
!
!------------------------------------------------------------------------      
       IMPLICIT NONE
       REAL(DP),DIMENSION(Ns),INTENT(INOUT) :: t, u 
       REAL(DP),EXTERNAL                    :: f
       REAL(DP), INTENT(IN)                   :: t0,tf,y0,dt
       INTEGER , INTENT(IN)                   :: Ns
       CHARACTER(LEN=*) , INTENT(IN)          :: eqn 
       REAL(DP)               :: error, up ,uOld , uNew, toll=1D-13
       REAL(DP)               :: fpred , fNew ,fOld
       REAL(DP)               :: k1,k2,k3,k4,k5
       INTEGER                :: i,stat
       CHARACTER(LEN=4)       :: f_name= 'AM5_'
       CHARACTER(LEN=4)       :: f_ext = '.out'
       CHARACTER(LEN=13)      :: path = 'result/adams/'
       CHARACTER(LEN=24)      :: filename
!----------       
       
       
       filename = path//f_name//eqn//f_ext
       

       t(1) = t0
       t(Ns)= tF
      
       u(1) = y0
       
       ! compute first 4 points using Runge-Kutta-Merson 5th Order (start-up the method )
       DO i =1,4  
         k1      = dt*f( t(i)       , u(i)                        )
         k2      = dt*f( t(i)+dt/3. , u(i)+k1/3.                  )
         k3      = dt*f( t(i)+dt/3. , u(i)+1./6.*(k1+k2)          )
         k4      = dt*f( t(i)+dt/2. , u(i)+1./8.*(k1+3.*k3)       )
         k5      = dt*f( t(i)+dt    , u(i)+1./2.*(k1-3.*k3+4.*k4) )
         u(i+1)  = u(i) + 1./6.* (k1 + 4*k4 + k5)  
       END DO  
           
      ! DO i =1,4  
      !   k1      = f( t(i)       , u(i)                   )
      !   k2      = f( t(i)+dt/2. , u(i)+dt/2.*k1          )
      !   k3      = f( t(i)+dt/2. , u(i)+dt/2.*k2          )
      !   k4      = f( t(i)+dt    , u(i)+dt*k3             )
      !   
      !   u(i+1)  = u(i) + dt/6. (k1 +2.*k2 + 2.*k3+  k4 )  
      ! END DO  
     !  DO i =1,4  
     !    k1      = f( t(i)       , u(i)          )
     !    k2      = f( t(i)+dt/2. , u(i)+dt/2.*k1 )
     !    k3      = f( t(i)+dt/2. , u(i)+dt/2.*k2 )
     !    k4      = f( t(i)+ dt   , u(i)+dt   *k3 )
     !    
     !    u(i+1)  = u(i) + dt/6.* (k1 + 2*k2+ 2*k3 + k4)  
     !  END DO  
 


       t(2:Ns) = (/( (t(1) + dt*(i-1)), i=2,Ns )/)    ! time vector 
      
      WRITE(6,'(A)') "Running Adams Bashforth (Predictor) , Adams Moulton 5th order (Corrector)"

!
!>          PREDICTOR ADAMS BASHFORTH  5th order            
       DO i=5,Ns
        up = u(i) + dt  *(1901.0/720.0  *f(t(i),u(i))          &
     &                   -1387.0/360.0 *f(t(i-1),u(i-1))       &
     &                    +109.0/30.0   *(f(t(i-2),u(i-2) ) )  &
     &                    -637.0/360.0  *(f(t(i-3),u(i-3) ) )  &
     &                    +251.0/720.0  *(f(t(i-4),u(i-4) )))     


        fpred = f(t(i+1),up)
        error = 1.D0
        fold = fpred
        !j=0
!>           CORRECTOR ADAMS MOULTON 5th ORDER
        DO WHILE( error .GT. toll )
         
         uOld = u(i) + dt/720. * ( 251. * fold              &
     &                           + 646. * f( t(i),u(i)   )  & 
     &                           - 264. * f( t(i-1),u(i-1))  &
     &                           + 106. * f( t(i-2),u(i-2))  &
     &                           -  19. * f( t(i-3),u(i-3)) )
         
         fold = f(t(i+1),uOld) 

         uNew = u(i) + dt/720. * ( 251. * fold              &
     &                           + 646. * f( t(i),u(i)   )  & 
     &                           - 264. * f( t(i-1),u(i-1))  &
     &                           + 106. * f( t(i-2),u(i-2))  &
     &                           -  19. * f( t(i-3),u(i-3)) )


         fnew = f(t(i+1),uNew )
         
         error = abs(uNew-uOld)
         !print*, error
         fold = fnew

        END DO 
        u(i+1) = unew
       END DO
      WRITE(6,'(A)') "... Done"

      OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',IOSTAT=stat)
         IF (stat .NE. 0) THEN 
          PRINT*, 'Error opening file ' , filename, ': (EXIT 1)' 
         ELSE 
          WRITE(10,FMT='(2F10.5)') (t(i),u(i), i=1,Ns)
         END IF 
      RETURN 
      END SUBROUTINE ADAMSMOLULTON5th


!------------------------------------------------------------------------      
      END MODULE adams_moulton





