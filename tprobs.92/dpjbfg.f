      subroutine dpjbfg(nx,ny,x,f,fgrad,task,ecc,b)
      character*(*) task
      integer nx,ny
      double precision f,ecc,b
      double precision x(nx*ny),fgrad(nx*ny)
c     **********
c
c     Subroutine dpjbfg
c
c     This subroutine computes the function and gradient of the 
c     Pressure Distribution in a Journal Bearing problem as formulated 
c     by G. Cimatti. This problem arises in the determination of the
c     pressure distribution in a thin film of lubricant between two 
c     circular cylinders with eccentricity ecc. A finite element 
c     approximation to the journal bearing problem is obtained by using 
c     piecewise linear basis functions over a triangulation of the 
c     domain D = (0,2*pi)X(0,2*b).
c
c     The subroutine statement is:
c
c       subroutine dpjbfg(nx,ny,x,f,fgrad,task,ecc,b)
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
c             'F'     Evaluate the function at x.
c             'G'     Evaluate the gradient at x.
c             'FG'    Evaluate the function and the gradient at x.
c             'XS'    Set x to the standard starting point xs.
c             'XL'    Set x to the lower bound xl.
c
c         On exit task is unchanged.
c
c       ecc is a double precision variable
c         On entry ecc is the eccentricity in (0,1).
c         On exit ecc is unchanged
c
c       b is a double precision variable
c         On entry b defines the domain as D = (0,2*pi) X (0,2*b).
c         On exit b is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      double precision zero,p5,one,two,four,six
      parameter(zero=0.0d0,p5=0.5d0,one=1.0d0,two=2.0d0,four=4.0d0,
     +          six=6.0d0)

      integer i,j,k
      double precision flin,fquad,trule,v,vb,vl,vr,vt,xi,pi,hx,hy,hxhy,
     +       temp,dvdx,dvdy,ehxhy

      double precision p

      p(xi) = (1 + ecc*cos(xi))**3

c     Initialization. 

      pi = four*atan(one)
      hx = two*pi/dble(nx + 1)
      hy = two*b/dble(ny + 1)
      hxhy = hx*hy
      ehxhy = ecc*hxhy

c     Compute the lower bound xl for x if task = 'XL'.

      if (task .eq. 'XL') then
         do 10 k = 1, nx*ny
            x(k) = zero
   10    continue

         return

      endif
     
c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         do 30 i = 1, nx
            temp = max(sin(dble(i)*hx),zero)
            do 20 j = 1, ny  
                k = nx*(j - 1) + i
                x(k) = temp
   20       continue
   30    continue

         return

      endif
  
c     Compute the function if task = 'F', the gradient if task = 'G', or
c     both if task = 'FG'.

      if (task .eq. 'F' .or. task .eq. 'FG') then
         fquad = zero
         flin = zero
      endif
      if (task .eq. 'G' .or. task .eq. 'FG') then
         do 40 k = 1, nx*ny
            fgrad(k) = zero
   40    continue
      endif

c     Computation of the quadratic part of the function and 
c     corresponding components of the gradient over the 
c     lower triangular elements.
      
      do 60 i = 0, nx
         xi = dble(i)*hx
         trule = hxhy*(p(xi) + p(xi+hx) + p(xi))/six
         do 50 j = 0, ny
            k = nx*(j - 1) + i
            v = zero
            vr = zero
            vt = zero
            if (i .ne. 0 .and. j .ne. 0) v = x(k)
            if (i .ne. nx .and. j .ne. 0) vr = x(k+1)
            if (i .ne. 0 .and. j .ne. ny) vt = x(k+nx)
            dvdx = (vr - v)/hx
            dvdy = (vt - v)/hy
            if (task .eq. 'F' .or. task .eq. 'FG') fquad = fquad + 
     +                                         trule*(dvdx**2 + dvdy**2)
            if (task .eq. 'G' .or. task .eq. 'FG') then
               if (i .ne. 0 .and. j .ne. 0) fgrad(k) = fgrad(k) - 
     +                                         trule*(dvdx/hx + dvdy/hy)
               if (i .ne. nx .and. j .ne. 0) fgrad(k+1) = fgrad(k+1) + 
     +                                                     trule*dvdx/hx
               if (i .ne. 0 .and. j .ne. ny) fgrad(k+nx) = fgrad(k+nx) + 
     +                                                     trule*dvdy/hy
            endif
   50    continue
   60 continue

c     Computation of the quadratic part of the function and
c     corresponding components of the gradient over the upper
c     triangular elements.

      do 80 i = 1, nx + 1
         xi = dble(i)*hx
         trule = hxhy*(p(xi) + p(xi-hx) + p(xi))/six
         do 70 j = 1, ny + 1
            k = nx*(j - 1) + i
            vb = zero
            vl = zero
            v = zero
            if (i .ne. nx + 1 .and. j .ne. 1) vb = x(k-nx)
            if (i .ne. 1 .and. j .ne. ny + 1) vl = x(k-1)
            if (i .ne. nx + 1 .and. j .ne. ny + 1) v = x(k)
            dvdx = (v - vl)/hx
            dvdy = (v - vb)/hy
            if (task .eq. 'F' .or. task .eq. 'FG') fquad = fquad + 
     +                                         trule*(dvdx**2 + dvdy**2)
            if (task .eq. 'G' .or. task .eq. 'FG') then
               if (i .le. nx .and. j .gt. 1) fgrad(k-nx) = fgrad(k-nx) - 
     +                                                     trule*dvdy/hy
               if (i .gt. 1 .and. j .le. ny) fgrad(k-1) = fgrad(k-1) - 
     +                                                     trule*dvdx/hx
               if (i .le. nx .and. j .le. ny) fgrad(k) = fgrad(k) + 
     +                                         trule*(dvdx/hx + dvdy/hy)
            endif
   70    continue
   80 continue

c     Computation of the linear part of the function and 
c     corresponding components of the gradient.

      if (task .eq. 'F' .or. task .eq. 'FG') then
         do 100 i = 1, nx
            temp = sin(dble(i)*hx)
            do 90 j = 1, ny
               k = nx*(j - 1) + i
               flin = flin + temp*x(k)
   90       continue
  100    continue
         f = p5*fquad - ehxhy*flin
      endif

      if (task .eq. 'G' .or. task .eq. 'FG') then
         do 120 i = 1, nx
            temp = ehxhy*sin(dble(i)*hx)
            do 110 j = 1, ny
               k = nx*(j - 1) + i
               fgrad(k) = fgrad(k) - temp
  110       continue
  120    continue
      endif

      return

      end

      subroutine dpjbds(nx,ny,w,ecc,b)
      integer nx,ny
      double precision ecc,b
      double precision w(nx*ny)
c     **********
c
c     This subroutine computes a scaling for the Pressure Distribution 
c     in a Journal Bearing problem as the square roots of the diagonal 
c     elements of the Hessian matrix.
c
c     The subroutine statement is:
c
c       subroutine dpjbds(nx,ny,w,ecc,b)
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
c       w is a double precision array of dimension nx*ny.
c         On entry w need not be specified.
c         On exit w contains the square roots of the diagonal elements
c            of the Hessian matrix.
c
c       ecc is a double precision variable.
c         On entry ecc is the eccentricity in (0,1).
c         On exit ecc is unchanged.
c
c       b is a double precision variable.
c         On entry b defines the domain as D = (0,2*pi) X (0,2*b).
c         On exit b is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      double precision one,two,three,four,six
      parameter(one=1.0d0,two=2.0d0,three=3.0d0,four=4.0d0,six=6.0d0)

      integer i,j,k
      double precision xi,pl,pc,pr,rulel,ruler,ruleu,ruled,hx,hy,hxhy,pi

      double precision p

      p(xi) = (1 + ecc*cos(xi))**3

      pi = four*atan(one)
      hx = two*pi/dble(nx + 1)
      hy = two*b/dble(ny + 1)
      hxhy = hx*hy

c     Computation of the scaling.

      do 20 i = 1, nx
         xi = dble(i)*hx
         do 10 j = 1, ny
            k = nx*(j - 1) + i
            pl = p(xi-hx)
            pc = p(xi)
            pr = p(xi+hx)
            ruler = hxhy*(three*pc + three*pr)/six
            rulel = hxhy*(three*pc + three*pl)/six
            ruleu = hxhy*(four*pc + pr + pl)/six
            ruled = ruleu
            w(k) = sqrt((ruler/hx)*(one/hx) + (ruleu/hy)*(one/hy) +
     +                        (rulel/hx)*(one/hx) + (ruled/hy)*(one/hy))
   10    continue
   20 continue

      return

      end
