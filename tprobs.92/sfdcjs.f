      subroutine sfdcjs(nx,ny,x,s,y,r)
      integer nx,ny
      real r
      real x(nx*ny),s(nx*ny),y(nx*ny)
c     *******
c
c     Subroutine sfdcjs
c
c     This subroutine computes the product
c
c                            J(x)*s = y
c
c     where J is the Jacobian for the flow in a channel problem
c     at the point x.
c
c     The subroutine statement is:
c
c        subroutine sfdcjs(nx,ny,x,s,y,r)
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
c         On exit y contains the product J*s.
c
c       r is a real variable.
c         On entry r is the Reynolds number.
c         On exit r is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick
c
c     **********
      real zero,one,two,three,four
      parameter (zero=0.0,one=1.0,two=2.0,three=3.0,four=4.0)

      integer i,j,k,n
      real pbb,pbl,pb,pbr,pll,pl,p,pr,prr,ptl,pt,ptr,ptt,dpdy,dpdx,
     +     pblap,pllap,prlap,ptlap,hy,hx,hy2,hx2

      n = nx*ny
      hx = one/real(nx+1)
      hy = one/real(ny+1)
      hy2 = hy*hy
      hx2 = hx*hx

      do 20 i = 1,ny
         do 10 j = 1,nx
            k = (i-1)*nx + j
            if (i.eq.1 .or. j.eq.1) then
               pbl = zero

            else
               pbl = x(k-nx-1)
            end if

            if (i.eq.1) then
               pb = zero
               pbb = x(k)

            else if (i.eq.2) then
               pb = x(k-nx)
               pbb = zero

            else
               pb = x(k-nx)
               pbb = x(k-2*nx)
            end if

            if (i.eq.1 .or. j.eq.nx) then
               pbr = zero

            else
               pbr = x(k-nx+1)
            end if

            if (j.eq.1) then
               pl = zero
               pll = x(k)

            else if (j.eq.2) then
               pl = x(k-1)
               pll = zero

            else
               pl = x(k-1)
               pll = x(k-2)
            end if

            p = x(k)
            if (j.eq.nx-1) then
               pr = x(k+1)
               prr = zero

            else if (j.eq.nx) then
               pr = zero
               prr = x(k)

            else
               pr = x(k+1)
               prr = x(k+2)
            end if

            if (i.eq.ny .or. j.eq.1) then
               ptl = zero

            else
               ptl = x(k+nx-1)
            end if

            if (i.eq.ny-1) then
               pt = x(k+nx)
               ptt = zero

            else if (i.eq.ny) then
               pt = zero
               ptt = x(k) + two*hy

            else
               pt = x(k+nx)
               ptt = x(k+2*nx)
            end if

            if (i.eq.ny .or. j.eq.nx) then
               ptr = zero

            else
               ptr = x(k+nx+1)
            end if

            dpdy = (pt-pb)/ (two*hy)
            dpdx = (pr-pl)/ (two*hx)

            pblap = (pbr-two*pb+pbl)/hx2 + (p-two*pb+pbb)/hy2
            pllap = (p-two*pl+pll)/hx2 + (ptl-two*pl+pbl)/hy2
            prlap = (prr-two*pr+p)/hx2 + (ptr-two*pr+pbr)/hy2
            ptlap = (ptr-two*pt+ptl)/hx2 + (ptt-two*pt+p)/hy2

            y(k) = zero
            if (i.gt.2) y(k) = y(k) + (one/hy2/hy2-
     +      r*dpdx/hy2/ (two*hy))*s(k-2*nx)
            if (j.gt.2) y(k) = y(k) + (one/hx2/hx2+
     +      r*dpdy/hx2/ (two*hx))*s(k-2)
            if (j.lt.nx-1) y(k) = y(k) +
     +      (one/hx2/hx2-r*dpdy/hx2/ (two*hx))*s(k+2)
            if (i.lt.ny-1) y(k) = y(k) +
     +      (one/hy2/hy2+r*dpdx/hy2/ (two*hy))*s(k+2*nx)
            if (i.ne.1 .and. j.ne.1) y(k) = y(k) +
     +      (two/hy2/hx2+r* (dpdy/hy2/ (two*hx)-dpdx/hx2/ (two*hy)))*
     +      s(k-nx-1)
            if (i.ne.1 .and. j.ne.nx) y(k) = y(k) +
     +      (two/hy2/hx2-r* (dpdy/hy2/ (two*hx)+dpdx/hx2/ (two*hy)))*
     +      s(k-nx+1)
            if (i.ne.ny .and. j.ne.1) y(k) = y(k) +
     +      (two/hy2/hx2+r* (dpdy/hy2/ (two*hx)+dpdx/hx2/ (two*hy)))*
     +      s(k+nx-1)
            if (i.ne.ny .and. j.ne.nx) y(k) = y(k) +
     +      (two/hy2/hx2-r* (dpdy/hy2/ (two*hx)-dpdx/hx2/ (two*hy)))*
     +      s(k+nx+1)
            if (i.ne.1) y(k) = y(k) + (-four* (one/hy2/hx2+one/hy2/hy2)+
     +      r* ((prlap-pllap)/ (two*hx)/ (two*hy)+dpdx* (one/hx2+
     +      one/hy2)/hy))*s(k-nx)
            if (j.ne.1) y(k) = y(k) + (-four* (one/hy2/hx2+one/hx2/hx2)-
     +      r* ((ptlap-pblap)/ (two*hx)/ (two*hy)+dpdy* (one/hx2+
     +      one/hy2)/hx))*s(k-1)
            if (j.ne.nx) y(k) = y(k) + (-four*
     +      (one/hy2/hx2+one/hx2/hx2)+r* ((ptlap-
     +      pblap)/ (two*hx)/ (two*hy)+dpdy* (one/hx2+one/hy2)/hx))*
     +      s(k+1)
            if (i.ne.ny) y(k) = y(k) + (-four*
     +      (one/hy2/hx2+one/hy2/hy2)-r* ((prlap-
     +      pllap)/ (two*hx)/ (two*hy)+dpdx* (one/hx2+one/hy2)/hy))*
     +      s(k+nx)

            y(k) = y(k) + (two* (four/hx2/hy2+three/hx2/hx2+
     +      three/hy2/hy2))*s(k)
            if (i.eq.1) y(k) = y(k) + (one/hy2/hy2-
     +      r* (dpdx/hy2/ (two*hy)))*s(k)
            if (j.eq.1) y(k) = y(k) + (one/hx2/hx2+
     +      r* (dpdy/hx2/ (two*hx)))*s(k)
            if (j.eq.nx) y(k) = y(k) + (one/hx2/hx2-
     +      r* (dpdy/hx2/ (two*hx)))*s(k)
            if (i.eq.ny) y(k) = y(k) + (one/hy2/hy2+
     +      r* (dpdx/hy2/ (two*hy)))*s(k)
   10    continue
   20 continue

c     Scale the result.

      do 30 k = 1,n
         y(k) = hx2*hy2*y(k)
   30 continue

      return

      end
