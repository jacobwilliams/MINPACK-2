      subroutine ssschs(nx,ny,x,s,y,lambda)
      integer nx, ny
      real lambda
      real x(nx*ny), s(nx*ny), y(nx*ny)
c     **********
c
c     Subroutine ssschs
c
c     This subroutine computes the product f''(x)*s = y, where
c     f''(x) is the Hessian matrix of the steady state combustion
c     problem exaluated at x.
c
c     The subroutine statement is
c
c       subroutine ssschs(nx,ny,x,s,y,lambda)
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
c         On exit y contains f''(x)*s.
c
c       lambda is a real variable.
c         On entry lambda is the multiplier from the Bratu problem.
c         On exit lambda is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      real one, p5, three, zero
      parameter (zero=0.0,p5=0.5,one=1.0,three=3.0)

      integer i, j, k
      real area, hx, hxhx, hy, hyhy, v, vb, vl, vr, vt

      hx = one/real(nx+1)
      hy = one/real(ny+1)
      hxhx = one/(hx*hx)
      hyhy = one/(hy*hy)
      area = p5*hx*hy

      do 10 k = 1, nx*ny
         y(k) = zero
   10 continue

c     Computation of f''(x)*s over the lower triangular elements.

      do 30 j = 0, ny
         do 20 i = 0, nx
            k = nx*(j-1) + i
            v = zero
            vr = zero
            vt = zero
            if (i .ne. 0 .and. j .ne. 0) v = s(k)
            if (i .ne. nx .and. j .ne. 0) then
               vr = s(k+1)
               y(k+1) = y(k+1) + hxhx*(vr-v) -
     +                  lambda*exp(x(k+1))*vr/three
            end if
            if (i .ne. 0 .and. j .ne. ny) then
               vt = s(k+nx)
               y(k+nx) = y(k+nx) + hyhy*(vt-v) -
     +                   lambda*exp(x(k+nx))*vt/three
            end if
            if (i .ne. 0 .and. j .ne. 0) y(k) = y(k) + hxhx*(v-vr) +
     +          hyhy*(v-vt) - lambda*exp(x(k))*v/three
   20    continue
   30 continue

c     Computation of f''(x)*s over the upper triangular elements.

      do 50 j = 1, ny + 1
         do 40 i = 1, nx + 1
            k = nx*(j-1) + i
            vb = zero
            vl = zero
            v = zero
            if (i .ne. nx+1 .and. j .ne. ny+1) v = s(k)
            if (i .ne. nx+1 .and. j .ne. 1) then
               vb = s(k-nx)
               y(k-nx) = y(k-nx) + hyhy*(vb-v) -
     +                   lambda*exp(x(k-nx))*vb/three
            end if
            if (i .ne. 1 .and. j .ne. ny+1) then
               vl = s(k-1)
               y(k-1) = y(k-1) + hxhx*(vl-v) -
     +                  lambda*exp(x(k-1))*vl/three
            end if
            if (i .ne. nx+1 .and. j .ne. ny+1) y(k) = y(k) +
     +          hxhx*(v-vl) + hyhy*(v-vb) - lambda*exp(x(k))*v/three
   40    continue
   50 continue

      do 60 k = 1, nx*ny
         y(k) = area*y(k)
   60 continue

      end
