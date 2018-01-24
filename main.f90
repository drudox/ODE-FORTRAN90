      PROGRAM main 
       USE PARAM
       USE EQ
       USE EULER 
       USE JACOBIAN
       IMPLICIT NONE

       REAL(DP),DIMENSION(Ns) :: t1,u1      
       

       CALL BWDEULER(t1,u1,NumSol,dfdxdy,t0,tf,y0,dt,Ns,eqn)     
       CALL FWDEULER(t1,u1,NumSol,t0,tf,y0,dt,Ns,eqn)     
       CALL realSOL()
       
       !print*,Ns


       STOP
      END PROGRAM
