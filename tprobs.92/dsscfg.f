      subroutine dsscfg(nx,ny,x,f,fgrad,task,lambda)
      character*(*) task
      integer nx,ny
      double precision f,lambda
      double precision x(nx*ny),fgrad(nx*ny)
c     **********
c
c     Subroutine dsscfg
c
c     This subroutine computes the function and gradient of the
c     Steady State Combustion problem formulated as energy functional. 
c     This problem arises in the study of steady state solid fuel 
c     ignition. A finite element approximation to the energy functional 
c     is obtained by using piecewise linear basis functions over a 
c     triangulation of the unit square.
c
c     The subroutine statement is:
c
c       subroutine dsscfg(nx,ny,x,f,fgrad,task,lambda)
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
c       x is a double precision array of dimension nx*ny.
c         On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
c            x is set according to task.
c
c       f is a double precision variable.
c         On entry f need not be specified.
c         On exit f is set to the function evaluated at x if task = 'F'
c            or 'FG'.
c
c       fgrad is a double precision array of dimension nx*ny.
c         On entry fgrad need not be specified.
c         On exit fgrad contains the gradient evaluated at x if 
c            task = 'G' or 'FG'.
c
c       task is a character variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c            'F'      Evaluate the function at x.
c            'G'      Evaluate the gradient at x.
c            'FG'     Evaluate the function and the gradient at x.
c            'XS'     Set x to the standard starting point xs.
c            'XL'     Set x to the lower bound xl.
c            'XU'     Set x to the upper bound xu.
c
c         On exit task is unchanged.
c
c       lambda is a double precision variable.
c         On entry lambda is a nonnegative Frank-Kamenetski parameter.
c         On exit lambda is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      double precision zero,p5,one,three
      parameter(zero=0.0d0,p5=0.5d0,one=1.0d0,three=3.0d0)

      integer i,j,k
      double precision fexp,fquad,dvdx,dvdy,hx,hy,v,vb,vl,vr,vt,
     +       expv,expvr,expvt,expvb,expvl,area,temp,temp1

c     Initialization.

      hx = one/dble(nx + 1)
      hy = one/dble(ny + 1)
      area = p5*hx*hy

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then

         temp1 = lambda/(lambda + one)
         do 20 j = 1, ny
            temp = dble(min(j,ny-j+1))*hy
            do 10 i = 1, nx
               k = nx*(j - 1) + i
               x(k) = temp1*sqrt(min(dble(min(i,nx-i+1))*hx,temp))
   10       continue
   20    continue

         return

      endif 

c     Compute the function if task = 'F', the gradient if task = 'G', or
c     both if task = 'FG'.

      if (task .eq. 'F' .or. task .eq. 'FG') then
         fquad = zero
         fexp = zero
      endif
      if (task .eq. 'G' .or. task .eq. 'FG') then
         do 30 k = 1, nx*ny
            fgrad(k) = zero
   30    continue
      endif

c     Computation of the function and the gradient over the lower 
c     triangular elements.  The trapezoidal rule is used to estimate 
c     the integral of the exponential term.

      do 50 j = 0, ny
         do 40 i = 0, nx
            k = nx*(j - 1) + i
            v = zero
            vr = zero
            vt = zero
            if (i .ne. 0 .and. j .ne. 0) v = x(k)
            if (i .ne. nx .and. j .ne. 0) vr = x(k+1)
            if (i .ne. 0 .and. j .ne. ny) vt = x(k+nx)
            dvdx = (vr - v)/hx
            dvdy = (vt - v)/hy
            expv = exp(v)
            expvr = exp(vr)
            expvt = exp(vt)
            if (task .eq. 'F' .or. task .eq. 'FG') then
               fquad = fquad + dvdx**2 + dvdy**2
               fexp = fexp - lambda*(expv + expvr + expvt)/three
            endif
            if (task .eq. 'G' .or. task .eq. 'FG') then
               if (i .ne. 0 .and. j .ne. 0) fgrad(k) = fgrad(k) - 
     +                             dvdx/hx - dvdy/hy - lambda*expv/three
               if (i .ne. nx .and. j .ne. 0) fgrad(k+1) = fgrad(k+1) + 
     +                                      dvdx/hx - lambda*expvr/three
               if (i .ne. 0 .and. j .ne. ny) fgrad(k+nx) = fgrad(k+nx) + 
     +                                      dvdy/hy - lambda*expvt/three 
            endif
   40    continue
   50 continue

c     Computation of the function and the gradient over the upper 
c     triangular elements.  The trapezoidal rule is used to estimate
c     the integral of the exponential term.

      do 70 j = 1, ny + 1
         do 60 i = 1, nx + 1
            k = nx*(j - 1) + i
            vb = zero
            vl = zero
            v = zero
            if (i .ne. nx + 1 .and. j .ne. 1) vb = x(k-nx)
            if (i .ne. 1 .and. j .ne. ny + 1) vl = x(k-1)
            if (i .ne. nx + 1 .and. j .ne. ny + 1) v = x(k)
            dvdx = (v - vl)/hx
            dvdy = (v - vb)/hy
            expvb = exp(vb)
            expvl = exp(vl)
            expv = exp(v)
            if (task .eq. 'F' .or. task .eq. 'FG') then
               fquad = fquad + dvdx**2 + dvdy**2
               fexp = fexp - lambda*(expvb + expvl + expv)/three
            endif
            if (task .eq. 'G' .or. task .eq. 'FG') then
               if (i .ne. nx + 1 .and. j .ne. 1) fgrad(k-nx) = 
     +                        fgrad(k-nx) - dvdy/hy - lambda*expvb/three
               if (i .ne. 1 .and. j .ne. ny + 1) fgrad(k-1) = 
     +                         fgrad(k-1) - dvdx/hx - lambda*expvl/three
               if (i .ne. nx + 1 .and. j .ne. ny + 1) fgrad(k) = 
     +                  fgrad(k) + dvdx/hx + dvdy/hy - lambda*expv/three
            endif
   60    continue
   70 continue

c     Scale the result.

      if (task .eq. 'F' .or. task .eq. 'FG') then
         f = area*(p5*fquad + fexp)
      endif
      if (task .eq. 'G' .or. task .eq. 'FG') then
         do 80 k = 1, nx*ny
            fgrad(k) = area*fgrad(k)
   80    continue
      endif

      return

      end
