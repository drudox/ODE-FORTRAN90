      PROGRAM main 
       USE PARAM
       USE EQ2
       USE EULER 
       USE JACOBIAN
       USE MULTISTEP
       USE RUNGEKUTTA
       USE ADAMS_BASHFORTH
       USE ADAMS_MOULTON

       
       IMPLICIT NONE

       REAL(DP),DIMENSION(Ns) :: t1,u1      
       REAL(DP),DIMENSION(Ns*4) :: t2,u2      
       

       CALL BWDEULER(t1,u1,NumSol,dfdxdy,t0,tf,y0,dt,Ns,eqn)     
       CALL FWDEULER(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)     
       
       CALL LEAPFROG(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)     
       
       CALL RungeKutta4th(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL HEUN(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL ModifiedEuler(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL CRANKNICHOLSON(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)     
     !  CALL RungeKuttaFehlberg45(t2,u2,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL RungeKuttaMerson(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)

       
       CALL AdamsBashforth2nd(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL AdamsBashforth3th(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL AdamsBashforth4th(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL AdamsBashforth5th(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL AdamsMoulton2nd(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL ADAMSMOLULTON3th(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL ADAMSMOLULTON4th(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       CALL ADAMSMOLULTON5th(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)
       
       CALL realSOL()
       


       STOP
      END PROGRAM
