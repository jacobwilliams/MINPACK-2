      subroutine siacfj(m,n,x,fvec,fjac,ldfjac,task,nint,sigma)
      character*(*) task
      integer m,n,ldfjac,nint
      real sigma
      real x(n),fvec(m),fjac(ldfjac,n)
c     **********
c
c     Subroutine siacfj
c
c     This subroutine computes the function and Jacobian matrix of the
c     Isomerization of Alpha-Pinene - Collocation Formulation problem.
c     The problem formulation is based on the work of I. Tjoa and L.
c     Biegler and arises in the determination of reaction coefficients
c     in the thermal isomerization of alpha-pinene.
c
c     The alpha-pinene problem is modeled by the system initial value
c     problems
c
c                    y1'  =  -(T1 + T2)*y1
c                    y2'  =  T1*y1
c                    y3'  =  T2*y1 - (T3 + T4)*y3 + T5*y5
c                    y4'  =  T3*y3
c                    y5'  =  T4*y3 - T5*y5
c
c     with initial data
c
c                    y1(0) = 100,
c                    y2(0) = y3(0) = y4(0) = y5(0) = 0.
c
c     and unknown coefficients T1 ... T5.
c
c     This is a data fitting parameter estimation problem with
c     differential equation constraints. The method of collocation
c     is used to discretize the constraint equations which become
c     part of the optimization problem.  The functions and Jacobians
c     for both the data residual equations and the constraint
c     equations are computed.
c
c     The subroutine statement is:
c
c       subroutine siacfj(m,n,x,fvec,fjac,ldfjac,task,nint,sigma)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions.
c            For the alpha-pinene problem, m = 25*nint + 40.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables.
c            For the alpha-pinene problem n = 25*nint + 5.
c         On exit n is unchanged.
c
c       x is a real array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a real array of dimension m.
c         On entry fvec need not be specified.
c         On exit fvec contains the function evaluated at x if
c            task = 'F' or 'FJ'.
c
c       fjac is a real array of dimension (ldfjac,n).
c         On entry fjac need not be specified.
c         On exit fjac contains the Jacobian matrix evaluated at x if
c            task = 'J' or 'FJ'.
c
c       ldfjac is an integer variable.
c          On entry ldfjac is the leading dimension of fjac.
c          On exit ldfjac is unchanged.
c
c       task is a character variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c            'F'      Evaluate the function at x.
c            'J'      Evaluate the Jacobian matrix at x.
c            'FJ'     Evaluate the function and the Jacobian at x.
c            'XS'     Set x to the standard starting point xs.
c
c         On exit task is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the
c            k-stage collocation.
c         On exit nint is unchanged.
c
c       sigma is a real variable.
c         On entry sigma is the penalty constant for the constraint
c            equations.
c         On exit sigma is unchanged
c
c     Subroutines called:
c
c     MINPACK Supplied ....... siarfj, siaofj.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer s,cpts,ntimes
      parameter (s=5,cpts=4,ntimes=8)

      integer i,j,ideg,npi,neqn
      integer deg(s)

      data (deg(i),i=1,s)/1,1,1,1,1/

c     Check input arguments for errors.

      if (nint.le.0) task = 'ERROR: NINT MUST BE .GT. 0'
      if (m.ne.25*nint+40) task = 'ERROR: M MUST .EQ. 25*NINT + 40'
      if (n.ne.25*nint+5) task = 'ERROR: N MUST .EQ. 25*NINT + 5'
      if (task(1:5).eq.'ERROR') return

      if (task .eq. 'XS') then
         call siaofj(m,n,x,fvec,fjac,ldfjac,task,nint)
         return
      endif

c     Initialize.

      ideg = 0
      do 10 i = 1,s
         ideg = ideg + deg(i)
   10 continue
      npi = s*cpts + ideg
      neqn = nint*npi

c     Evaluate the data residual equations.

      call siarfj(m,n,x,fvec,fjac,ldfjac,task,nint)

c     Move the information down.

      if (task.eq.'F' .or. task.eq.'FJ') then
         do 20 i = 1,s*ntimes
            fvec(neqn+i) = fvec(i)
   20    continue
      end if

      if (task.eq.'J' .or. task.eq.'FJ') then
         do 40 j = 1,n
            do 30 i = 1,s*ntimes
               fjac(neqn+i,j) = fjac(i,j)
   30       continue
   40    continue
      end if

c     Evaluate the collocation constraint equations.

      call siaofj(m,n,x,fvec,fjac,ldfjac,task,nint)

c     Penalize the collocation constraints by sigma.

      if (task.eq.'F' .or. task.eq.'FJ') then
         do 50 i = 1,neqn
            fvec(i) = sigma*fvec(i)
   50    continue
      end if

      if (task.eq.'J' .or. task.eq.'FJ') then
         do 70 j = 1,n
            do 60 i = 1,neqn
               fjac(i,j) = sigma*fjac(i,j)
   60       continue
   70    continue
      end if

      return

      end
