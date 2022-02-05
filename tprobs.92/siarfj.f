      subroutine siarfj(m,n,x,fvec,fjac,ldfjac,task,nint)
      character*(*) task
      integer n,m,ldfjac,nint
      real x(n),fvec(m),fjac(ldfjac,n)
c     **********
c
c     Subroutine siarfj
c
c     This subroutine computes the function and Jacobian matrix of the
c     data residual equations for the Isomerization of Alpha-Pinene -
c     Collocation Formulation problem.  The problem formulation is based
c     on the work of I. Tjoa and L. Biegler and arises in determining
c     reaction coefficients in thermal isomerization of alpha-pinene.
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
c     part of the optimization problem.
c
c     The subroutine statement is:
c
c       subroutine siacfj(m,n,x,fvec,fjac,ldfjac,task,nint,sigma)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions.
c            For the alpha-pinene problem, m = 40.
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
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer s,cpts,maxdeg,ntimes,dim
      parameter (s=5,cpts=4,maxdeg=1,ntimes=8,dim=maxdeg+cpts-1)
      real len,zero,one
      parameter (len=4.0E4,zero=0.0,one=1.0)

      integer i,j,k,mm,e,ideg,npi,neqn
      integer deg(s),var(s),eqn(s),sumdeg(0:s),intt(ntimes)
      real h,hm,tijh,nf
      real obs(ntimes,s),dt(ntimes,2),dw(maxdeg+1,cpts+maxdeg,s),
     +     tnfhk(ntimes,0:dim,0:dim,0:maxdeg),w(maxdeg+1,s)

      data (deg(i), i = 1, s)/1,1,1,1,1/
      data (dt(i,1), i = 1, ntimes)
     +     /1230.0,3060.0,4920.0,7800.0,10680.0,15030.0,22620.0,36420.0/
      data (obs(i,1),i=1,ntimes)
     +     /88.35,76.40,65.10,50.40,37.50,25.90,14.00,4.50/
      data (obs(i,2),i=1,ntimes)/7.3,15.6,23.1,32.9,42.7,49.1,57.4,63.1/
      data (obs(i,3),i=1,ntimes)/2.3,4.5,5.3,6.0,6.0,5.9,5.1,3.8/
      data (obs(i,4),i=1,ntimes)/0.4,0.7,1.1,1.5,1.9,2.2,2.6,2.9/
      data (obs(i,5),i=1,ntimes)/1.75,2.8,5.8,9.3,12.0,17.0,21.0,25.7/

c     Check input arguments for errors.

      if (nint .le. 0) task = 'ERROR: NINT MUST BE .GT. 0'
      if (m .ne. 40) task = 'ERROR: M MUST .EQ. 40'
      if (n .ne. 25*nint+5) task = 'ERROR: N MUST .EQ. 25*NINT + 5'
      if (task(1:5).eq.'ERROR') return

c     Initialization.

      h = len/real(nint)

c     Determine in which time interval each observation lies and
c     what the value of the independent variable is relative to
c     that interval.

      do 10 i = 1,ntimes
         k = int(dt(i,1)/h)
         dt(i,2) = (dt(i,1)-real(k)*h)/h
         intt(i) = k + 1
   10 continue

c     Store all possible combinations of zeta, h, and factorial.

      hm = one
      do 50 mm = 0,maxdeg
         do 40 i = 1,ntimes
            tijh = hm
            do 30 j = 0,dim
               nf = one
               do 20 k = 0,dim
                  tnfhk(i,j,k,mm) = tijh/nf
                  nf = nf*real(k+1)
   20          continue
               tijh = tijh*dt(i,2)
   30       continue
   40    continue
         hm = hm*h
   50 continue

      ideg = 0
      sumdeg(0) = 0
      do 60 i = 1,s
         ideg = ideg + deg(i)
         sumdeg(i) = sumdeg(i-1) + deg(i)
   60 continue
      npi = s*cpts + ideg
      neqn = nint*npi

      if (task.eq.'F' .or. task.eq.'FJ') then
         do 70 i = 1, ntimes*s
            fvec(i) = zero
   70    continue
      end if

      if (task.eq.'J' .or. task.eq.'FJ') then
         do 90 j = 1, neqn + s
            do 80 i = 1, ntimes*s
               fjac(i,j) = zero
   80       continue
   90    continue
      end if

c     Set up the data residual equations.

      do 170 k = 1, ntimes
         i = intt(k)
         do 130 e = 1, s
            var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
c           eqn(e) = neqn + (k - 1)*s + e
            eqn(e) = (k-1)*s + e
            do 120 mm = 1, deg(e) + 1
               w(mm,e) = zero
               do 100 j = mm,deg(e)
                  w(mm,e) = w(mm,e) + tnfhk(k,j-mm,j-mm,j-mm)*
     +            x(var(e)+j)
                  dw(mm,j,e) = tnfhk(k,j-mm,j-mm,j-mm)
  100          continue
               do 110 j = 1,cpts
                  w(mm,e) = w(mm,e) + tnfhk(k,deg(e)+j-mm,deg(e)+j-mm,
     +            deg(e)-mm+1)*x(var(e)+deg(e)+j)
                  dw(mm,deg(e)+j,e) = tnfhk(k,deg(e)+j-mm,deg(e)+j-mm,
     +            deg(e)-mm+1)
  110          continue
  120       continue
  130    continue

         if (task.eq.'F' .or. task.eq.'FJ') then
            do 140 e = 1,s
               fvec(eqn(e)) = w(1,e) - obs(k,e)
  140       continue
         end if

         if (task.eq.'J' .or. task.eq.'FJ') then
            do 160 e = 1,s
               do 150 j = 1,cpts + deg(e)
                  fjac(eqn(e),var(e)+j) = dw(1,j,e)
  150          continue
  160       continue
         end if

  170 continue

      return

      end
