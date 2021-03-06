      subroutine smsahs(nx,ny,x,s,y,bottom,top,left,right)
      integer nx,ny
      real x(nx*ny),s(nx*ny),y(nx*ny),bottom(nx+2),top(nx+2),left(ny+2),
     +     right(ny+2)
c     **********
c
c     Subroutine smsahs
c
c     This subroutine computes the product
c
c                            H(x)*s = y
c
c     where H(x) is the Hessian for the Minimal Surface problem
c     exaluted at x.
c
c     The subroutine statement is:
c
c       subroutine smsahs(nx,ny,x,s,y,bottom,top,left,right)
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
c         On entry x specifies the vector x.
c         On exit x is unchanged.
c
c       s is a real array of dimension nx*ny.
c         On entry s contains the vector s.
c         On exit s is unchanged.
c
c       y is a real array of dimension nx*ny.
c         On entry out need not be specified.
c         On exit y contains H(x)*s.
c
c       bottom is a real array of dimension nx + 2.
c         On entry bottom must contain appropriate boundary data.
c         On exit bottom is unchanged.
c
c       top is a real work array of dimension nx + 2.
c         On entry top must contain appropriate boundary data.
c         On exit top is unchanged.
c
c       left is a real work array of dimension ny + 2.
c         On entry left must contain appropriate boundary data.
c         On exit left is unchanged.
c
c       right is a real work array of dimension ny + 2.
c         On entry right must contain appropriate boundary data.
c         On exit right is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      real zero,p5,one,two
      parameter (zero=0.0,p5=0.5,one=1.0,two=2.0)
      real b,t,l,r
      parameter (b=-.50,t=.50,l=-.50,r=.50)

      integer i,j,k
      real hx,hy,vl,vr,vb,vt,v,dvdx,dvdy,dvdxhx,dvdyhy,area,tu,fu,fu3,
     +     tl,fl,fl3,z,zr,zt,zl,zb,dzdx,dzdxhx,dzdy,dzdyhy

      hx = (r-l)/real(nx+1)
      hy = (t-b)/real(ny+1)
      area = p5*hx*hy

      do 10 k = 1,nx*ny
         y(k) = zero
   10 continue

c     Computation of H(x)*s over the lower triangular elements.

      do 30 j = 0,ny
         do 20 i = 0,nx
            k = nx* (j-1) + i
            if (i.ne.0 .and. j.ne.0) then
               v = x(k)
               z = s(k)

            else
               if (j.eq.0) v = bottom(i+1)
               if (i.eq.0) v = left(j+1)
               z = zero
            end if

            if (i.ne.nx .and. j.ne.0) then
               vr = x(k+1)
               zr = s(k+1)

            else
               if (i.eq.nx) vr = right(j+1)
               if (j.eq.0) vr = bottom(i+2)
               zr = zero
            end if

            if (i.ne.0 .and. j.ne.ny) then
               vt = x(k+nx)
               zt = s(k+nx)

            else
               if (i.eq.0) vt = left(j+2)
               if (j.eq.ny) vt = top(i+1)
               zt = zero
            end if

            dvdx = (vr-v)/hx
            dvdy = (vt-v)/hy
            dzdx = (zr-z)/hx
            dzdy = (zt-z)/hy
            dvdxhx = dvdx/hx
            dvdyhy = dvdy/hy
            dzdxhx = dzdx/hx
            dzdyhy = dzdy/hy
            tl = one + dvdx**2 + dvdy**2
            fl = sqrt(tl)
            fl3 = fl*tl
            if (i.ne.0 .and. j.ne.0) y(k) = y(k) +
     +      (dvdx*dzdx+dvdy*dzdy)* (dvdxhx+dvdyhy)/fl3 -
     +      (dzdxhx+dzdyhy)/fl
            if (i.ne.nx .and. j.ne.0) y(k+1) = y(k+1) -
     +      (dvdx*dzdx+dvdy*dzdy)*dvdxhx/fl3 + dzdxhx/fl
            if (i.ne.0 .and. j.ne.ny) y(k+nx) = y(k+nx) -
     +      (dvdx*dzdx+dvdy*dzdy)*dvdyhy/fl3 + dzdyhy/fl
   20    continue
   30 continue

c     Computation of H(x)*s over the upper triangular elements.

      do 50 j = 1,ny + 1
         do 40 i = 1,nx + 1
            k = nx* (j-1) + i
            if (i.ne.nx+1 .and. j.ne.1) then
               vb = x(k-nx)
               zb = s(k-nx)

            else
               if (j.eq.1) vb = bottom(i+1)
               if (i.eq.nx+1) vb = right(j)
               zb = zero
            end if

            if (i.ne.1 .and. j.ne.ny+1) then
               vl = x(k-1)
               zl = s(k-1)

            else
               if (j.eq.ny+1) vl = top(i)
               if (i.eq.1) vl = left(j+1)
               zl = zero
            end if

            if (i.ne.nx+1 .and. j.ne.ny+1) then
               v = x(k)
               z = s(k)

            else
               if (i.eq.nx+1) v = right(j+1)
               if (j.eq.ny+1) v = top(i+1)
               z = zero
            end if

            dvdx = (v-vl)/hx
            dvdy = (v-vb)/hy
            dzdx = (z-zl)/hx
            dzdy = (z-zb)/hy
            dvdxhx = dvdx/hx
            dvdyhy = dvdy/hy
            dzdxhx = dzdx/hx
            dzdyhy = dzdy/hy
            tu = one + dvdx**2 + dvdy**2
            fu = sqrt(tu)
            fu3 = fu*tu
            if (i.ne.nx+1 .and. j.ne.ny+1) y(k) = y(k) -
     +      (dvdx*dzdx+dvdy*dzdy)* (dvdxhx+dvdyhy)/fu3 +
     +      (dzdxhx+dzdyhy)/fu
            if (i.ne.1 .and. j.ne.ny+1) y(k-1) = y(k-1) +
     +      (dvdx*dzdx+dvdy*dzdy)*dvdxhx/fu3 - dzdxhx/fu
            if (i.ne.nx+1 .and. j.ne.1) y(k-nx) = y(k-nx) +
     +      (dvdx*dzdx+dvdy*dzdy)*dvdyhy/fu3 - dzdyhy/fu
   40    continue
   50 continue

c     Scale the result by the area.

      do 60 k = 1,nx*ny
         y(k) = area*y(k)
   60 continue

      return

      end
