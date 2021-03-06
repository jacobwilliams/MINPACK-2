      subroutine spjbhs(nx,ny,s,y,ecc,b)
      integer nx,ny
      real ecc,b
      real s(nx*ny),y(nx*ny)
c     **********
c
c     Subroutine spjbhs
c
c     This subroutine computes the product
c
c                            H*s = y
c
c     where H is the Hessian for the Pressure Distribution in a
c     Journal Bearing problem.
c
c     The subroutine statement is:
c
c       subroutine spjbhs(nx,ny,s,y,ecc,b)
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
c       s is a real array of dimension nx*ny.
c         On entry s contains the vector s.
c         On exit s is unchanged.
c
c       y is a real array of dimension nx*ny.
c         On entry out need not be specified.
c         On exit y contains H*s.
c
c       ecc is a real variable
c         On entry ecc is the eccentricity in (0,1).
c         On exit ecc is unchanged
c
c       b is a real variable
c         On entry b defines the domain as D = (0,2*pi) X (0,2*b).
c         On exit b is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      real zero,one,two,four,six
      parameter (zero=0.0,one=1.0,two=2.0,four=4.0,six=6.0)

      integer i,j,k
      real hx,hxhx,hy,hyhy,vl,vr,vb,vt,v,pi,hxhy,trule,xi
      real p

      p(xi) = (one+ecc*cos(xi))**3

      pi = four*atan(one)
      hx = two*pi/real(nx+1)
      hy = two*b/real(ny+1)
      hxhy = hx*hy
      hxhx = one/ (hx*hx)
      hyhy = one/ (hy*hy)

      do 10 k = 1,nx*ny
         y(k) = zero
   10 continue

c     Computation of H*s over the lower triangular elements.

      do 30 i = 0,nx
         xi = real(i)*hx
         trule = hxhy* (p(xi)+p(xi+hx)+p(xi))/six
         do 20 j = 0,ny
            k = nx* (j-1) + i
            v = zero
            vr = zero
            vt = zero
            if (i.ne.0 .and. j.ne.0) v = s(k)
            if (i.ne.nx .and. j.ne.0) then
               vr = s(k+1)
               y(k+1) = y(k+1) + trule*hxhx* (vr-v)
            end if

            if (i.ne.0 .and. j.ne.ny) then
               vt = s(k+nx)
               y(k+nx) = y(k+nx) + trule*hyhy* (vt-v)
            end if

            if (i.ne.0 .and. j.ne.0) y(k) = y(k) +
     +      trule* (hxhx* (v-vr)+hyhy* (v-vt))
   20    continue
   30 continue

c     Computation of H*s over the upper triangular elements.

      do 50 i = 1,nx + 1
         xi = real(i)*hx
         trule = hxhy* (p(xi)+p(xi-hx)+p(xi))/six
         do 40 j = 1,ny + 1
            k = nx* (j-1) + i
            vb = zero
            vl = zero
            v = zero
            if (i.ne.nx+1 .and. j.ne.ny+1) v = s(k)
            if (i.ne.nx+1 .and. j.ne.1) then
               vb = s(k-nx)
               y(k-nx) = y(k-nx) + trule*hyhy* (vb-v)
            end if

            if (i.ne.1 .and. j.ne.ny+1) then
               vl = s(k-1)
               y(k-1) = y(k-1) + trule*hxhx* (vl-v)
            end if

            if (i.ne.nx+1 .and. j.ne.ny+1) y(k) = y(k) +
     +      trule* (hyhy* (v-vb)+hxhx* (v-vl))
   40    continue
   50 continue

      return

      end
