      subroutine sodcfg(nx,ny,x,f,fgrad,task,lambda)
      character*(*) task
      integer nx,ny
      real f,lambda
      real x(nx*ny),fgrad(nx*ny)
c     **********
c
c     Subroutine sodcfg
c
c     This subroutine computes the function and gradient of the Optimal
c     Design with Composite Materials problem formulated by J. Goodman,
c     R. Kohn, and L. Reyna. The problem arises in the determination of
c     the placement of two elastic materials in a rod with square cross-
c     section so as to obtain maximal torsional rigidity.  The nonlinear
c     functional that arises depends on a Lagrange multiplier lambda in
c     (0,1). A particularly interesting value is lambda = .008. A finite
c     element triangulation of the unit square is considered. The domain
c     can be altered by changing the appropriate parameter statements.
c
c     The subroutine statement is:
c
c       subroutine sodcfg(nx,ny,x,f,fgrad,task,lambda)
c
c     where
c
c       nx is an integer variable.
c         On entry nx is the number of grid points in the first
c            coordinate direction.
c         On exit nx is unchanged.
c
c       ny is an integer variable.
c         On entry ny is the number of grid points in the second
c            coordinate direction.
c         On exit ny is unchanged.
c
c       x is a real array of dimension nx*ny.
c         On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
c            x is set according to task.
c
c       f is a real variable.
c         On entry f need not be specified.
c         On exit f is set to the function evaluated at x if task = 'F'
c            or 'FG'.
c
c       fgrad is a real array of dimension nx*ny.
c         On entry fgrad need not be specified.
c         On exit fgrad contains the gradient evaluated at x if
c            task = 'G' or 'FG'.
c
c       task is a character variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c             'F'     Evaluate the function at x.
c             'G'     Evaluate the gradient at x.
c             'FG'    Evaluate the function and the gradient at x.
c             'XS'    Set x to the standard starting point xs.
c
c         On exit task is unchanged.
c
c       lambda is a real variable.
c         On entry lambda is the Lagrange multiplier.
c         On exit lambda is unchanged.
c
c     Subprograms called
c
c       MINPACK-supplied   ...   sodcps
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      real b,t,l,r,mu1,mu2,zero,p5,two
      parameter (b=0.0,t=1.0,l=0.0,r=1.0)
      parameter (mu1=1.0,mu2=2.0)
      parameter (zero=0.0,p5=0.5,two=2.0)

      integer i,j,k
      real dvdx,dvdy,vl,vr,vb,vt,v,area,gradv,hx,hy,hxhy,dpsi,dpsip,t1,
     +     t2,temp

c     Initialization.

      hx = (r-l)/real(nx+1)
      hy = (t-b)/real(ny+1)
      hxhy = hx*hy
      area = p5*hxhy

c     Compute the break points.

      t1 = sqrt(two*lambda*mu1/mu2)
      t2 = sqrt(two*lambda*mu2/mu1)

c     Compute the standard starting point if task = 'XS'.

      if (task.eq.'XS') then
         do 20 j = 1,ny
            temp = real(min(j,ny-j+1))*hy
            do 10 i = 1,nx
               k = nx* (j-1) + i
               x(k) = - (min(real(min(i,nx-i+1))*hx,temp))**2
   10       continue
   20    continue

         return

      end if

c     Evaluate the function if task = 'F', the gradient if task = 'G',
c     or both if task = 'FG'.

      if (task.eq.'F' .or. task.eq.'FG') f = zero
      if (task.eq.'G' .or. task.eq.'FG') then
         do 30 k = 1,nx*ny
            fgrad(k) = zero
   30    continue
      end if

c     Computation of the function and the gradient over the lower
c     triangular elements.

      do 50 j = 0,ny
         do 40 i = 0,nx
            k = nx* (j-1) + i
            v = zero
            vr = zero
            vt = zero
            if (j.ge.1 .and. i.ge.1) v = x(k)
            if (i.lt.nx .and. j.gt.0) vr = x(k+1)
            if (i.gt.0 .and. j.lt.ny) vt = x(k+nx)
            dvdx = (vr-v)/hx
            dvdy = (vt-v)/hy
            gradv = dvdx**2 + dvdy**2
            if (task.eq.'F' .or. task.eq.'FG') then
               call sodcps(gradv,mu1,mu2,t1,t2,dpsi,'PSI',lambda)
               f = f + dpsi
            end if

            if (task.eq.'G' .or. task.eq.'FG') then
               call sodcps(gradv,mu1,mu2,t1,t2,dpsip,'PSIP',lambda)
               if (i.ge.1 .and. j.ge.1) fgrad(k) = fgrad(k) -
     +         two* (dvdx/hx+dvdy/hy)*dpsip
               if (i.lt.nx .and. j.gt.0) fgrad(k+1) = fgrad(k+1) +
     +         two* (dvdx/hx)*dpsip
               if (i.gt.0 .and. j.lt.ny) fgrad(k+nx) = fgrad(k+nx) +
     +         two* (dvdy/hy)*dpsip
            end if

   40    continue
   50 continue

c     Computation of the function and the gradient over the upper
c     triangular elements.

      do 70 j = 1,ny + 1
         do 60 i = 1,nx + 1
            k = nx* (j-1) + i
            vb = zero
            vl = zero
            v = zero
            if (i.le.nx .and. j.gt.1) vb = x(k-nx)
            if (i.gt.1 .and. j.le.ny) vl = x(k-1)
            if (i.le.nx .and. j.le.ny) v = x(k)
            dvdx = (v-vl)/hx
            dvdy = (v-vb)/hy
            gradv = dvdx**2 + dvdy**2
            if (task.eq.'F' .or. task.eq.'FG') then
               call sodcps(gradv,mu1,mu2,t1,t2,dpsi,'PSI',lambda)
               f = f + dpsi
            end if

            if (task.eq.'G' .or. task.eq.'FG') then
               call sodcps(gradv,mu1,mu2,t1,t2,dpsip,'PSIP',lambda)
               if (i.le.nx .and. j.gt.1) fgrad(k-nx) = fgrad(k-nx) -
     +         two* (dvdy/hy)*dpsip
               if (i.gt.1 .and. j.le.ny) fgrad(k-1) = fgrad(k-1) -
     +         two* (dvdx/hx)*dpsip
               if (i.le.nx .and. j.le.ny) fgrad(k) = fgrad(k) +
     +         two* (dvdx/hx+dvdy/hy)*dpsip
            end if

   60    continue
   70 continue

c     Scale the function.

      if (task.eq.'F' .or. task.eq.'FG') f = area*f

c     Integrate v over the domain.

      if (task.eq.'F' .or. task.eq.'FG') then
         temp = zero
         do 80 k = 1,nx*ny
            temp = temp + x(k)
   80    continue
         f = f + hxhy*temp
      end if

      if (task.eq.'G' .or. task.eq.'FG') then
         do 90 k = 1,nx*ny
            fgrad(k) = area*fgrad(k) + hxhy
   90    continue
      end if

      return

      end

      subroutine sodcps(t,mu1,mu2,t1,t2,result,task,lambda)
      character*(*) task
      real t,mu1,mu2,t1,t2,result,lambda
c     **********
c
c     This subroutine computes the function psi(t) and the scaled
c     functions psi'(t)/t and psi''(t)/t from the Optimal Design
c     with Composite Materials problem.
c
c     The subroutine statement is:
c
c       subroutine sodcps(t,mu1,mu2,t1,t2,result,task,lambda)
c
c     where
c
c       t is a real variable.
c         On entry t is the variable t
c         On exit t is unchanged
c
c       mu1 is a real variable.
c         On entry mu1 is the reciprocal shear modulus of material 1.
c         On exit mu1 is unchanged.
c
c       mu2 is a real variable.
c         On entry mu2 is the reciprocal shear modulus of material 2.
c         On exit mu2 is unchanged.
c
c       t1 is a real variable.
c         On entry t1 is the first breakpoint.
c         On exit t1 is unchanged.
c
c       t2 is a real variable.
c         On entry t2 is the second breakpoint.
c         On exit t2 is unchanged.
c
c       result is a real variable.
c         On entry result need not be specified.
c         On exit result is set according to task.
c
c       task is a character variable.
c         On entry task specifies the action of the subroutine:
c
c             task        action
c             ----        ------
c            'PSI'        evaluate the function psi(t).
c            'PSIP'       evaluate the scaled function psi'(t)/t.
c            'PSIPP'      evaluate the scaled function psi''(t)/t.
c
c        On exit task is unchanged.
c
c       lambda is a real variable
c         On entry lambda is the Lagrange multiplier.
c         On exit lambda is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      real zero,p25,p5
      parameter (zero=0.0,p25=0.25,p5=0.5)

      real c1,c2,sqrtt

c     Compute the constants of integration.

      c1 = -lambda*mu1
      c2 = lambda* (mu2-mu1)

      sqrtt = sqrt(t)

      if (task.eq.'PSI') then
         if (sqrtt.le.t1) then
            result = p5*mu2*t

         else if (sqrtt.gt.t1 .and. sqrtt.lt.t2) then
            result = mu2*t1*sqrtt + c1

         else if (sqrtt.ge.t2) then
            result = p5*mu1*t + c2
         end if

      end if

      if (task.eq.'PSIP') then
         if (sqrtt.le.t1) then
            result = p5*mu2

         else if (sqrtt.gt.t1 .and. sqrtt.lt.t2) then
            result = p5*mu2*t1/sqrtt

         else if (sqrtt.ge.t2) then
            result = p5*mu1
         end if

      end if

      if (task.eq.'PSIPP') then
         if (sqrtt.le.t1) then
            result = zero

         else if (sqrtt.gt.t1 .and. sqrtt.lt.t2) then
            result = -p25*mu2*t1/ (sqrtt*t)

         else if (sqrtt.ge.t2) then
            result = zero
         end if

      end if

      return

      end
