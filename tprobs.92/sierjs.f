      subroutine sierjs(n,x,s,y,nint)
      integer n,nint
      real x(n),s(n),y(n)
c     **********
c
c     Subroutine sierjs
c
c     This subroutine computes the product
c
c                            J(x)*s = y
c
c     where J(x) is the Jacobian matrix of the Inverse Elastic
c     Rod problem.
c
c     The subroutine statement is:
c
c       subroutine sierjs(n,x,s,y,nint)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c            For the Inverse Elastic Rod problem n = 15*nint+3.
c         On exit n is unchanged.
c
c       x is a real array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       s is a real array of dimension n.
c         On entry s contains the vector s.
c         On exit s is unchanged.
c
c       y is a real array of dimension n.
c         On entry out need not be specified.
c         On exit y contains J*s.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the
c            k-stage collocation.
c         On exit nint is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota
c     Brett M. Averick, R.S. Maier, G.L. Xue
c
c     **********
      integer sigma,cpts,maxdeg,dim
      parameter (sigma=3,cpts=4,maxdeg=1,dim=maxdeg+cpts-1)
      real len,zero,one
      parameter (len=1.0,zero=0.0,one=1.0)

      integer i,j,k,l,m,e,ideg,npi
      integer deg(sigma),var(sigma),eqn(sigma),sumdeg(0:sigma)
      real rho(cpts),dw(maxdeg+1,cpts+maxdeg,sigma),
     +     rhnfhk(cpts,0:dim,0:dim,0:maxdeg),w(maxdeg+1,sigma)
      real h,hm,rhoijh,nf

      data (rho(i),i = 1, cpts)
     +     /0.694318413734436035E-1,0.330009490251541138,
     +      0.669990539550781250,0.930568158626556396/
      data (deg(i), i = 1, sigma)/1,1,1/

c     Check input arguments for errors.

      if (nint.le.0) print *,'ERROR: NINT MUST BE .GT. 0'
      if (n.ne.15*nint+3) print *,'ERROR: N MUST .EQ. 15*NINT + 3'

      h = len/real(nint)

      hm = one
      do 40 m = 0,maxdeg
         do 30 i = 1,cpts
            rhoijh = hm
            do 20 j = 0,dim
               nf = one
               do 10 k = 0,dim
                  rhnfhk(i,j,k,m) = rhoijh/nf
                  nf = nf*real(k+1)
   10          continue
               rhoijh = rhoijh*rho(i)
   20       continue
   30    continue
         hm = hm*h
   40 continue

      ideg = 0
      sumdeg(0) = 0
      do 50 i = 1,sigma
         ideg = ideg + deg(i)
         sumdeg(i) = sumdeg(i-1) + deg(i)
   50 continue
      npi = sigma*cpts + ideg

      do 60 i = 1,n
         y(i) = zero
   60 continue

      do 90 k = 1,sigma
         do 80 j = 1,maxdeg + 1
            w(j,k) = zero
            do 70 l = 1,cpts + maxdeg
               dw(j,l,k) = zero
   70       continue
   80    continue
   90 continue

      do 100 i = 1,sigma
         y(i) = y(i) + s((i-1)*cpts+sumdeg(i-1)+1)
  100 continue

      do 170 i = 1,nint
         do 160 k = 1,cpts
            do 140 e = 1,sigma
               var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
               eqn(e) = sigma + (i-1)*npi + (e-1)*cpts
               do 130 m = 1,deg(e) + 1
                  w(m,e) = zero
                  do 110 j = m,deg(e)
                     w(m,e) = w(m,e) + x(var(e)+j)*rhnfhk(k,j-m,j-m,j-m)
                     dw(m,j,e) = rhnfhk(k,j-m,j-m,j-m)
  110             continue
                  do 120 j = 1,cpts
                     w(m,e) = w(m,e) + x(var(e)+deg(e)+j)*
     +               rhnfhk(k,deg(e)+j-m,deg(e)+j-m,deg(e)-m+1)
                     dw(m,deg(e)+j,e) = rhnfhk(k,deg(e)+j-m,deg(e)+j-m,
     +               deg(e)-m+1)
  120             continue
  130          continue
  140       continue

            do 150 j = 1,cpts + deg(1)
               y(eqn(1)+k) = y(eqn(1)+k) + s(var(1)+j)*dw(2,j,1)
               y(eqn(3)+k) = y(eqn(3)+k) - s(var(1)+j)*x(n-2)*dw(1,j,1)
               y(eqn(2)+k) = y(eqn(2)+k) + s(var(2)+j)*dw(2,j,2)
               y(eqn(3)+k) = y(eqn(3)+k) + s(var(2)+j)*x(n-1)*dw(1,j,2)
               y(eqn(1)+k) = y(eqn(1)+k) +
     +         s(var(3)+j)*sin(w(1,3))*dw(1,j,3)
               y(eqn(2)+k) = y(eqn(2)+k) -
     +         s(var(3)+j)*cos(w(1,3))*dw(1,j,3)
               y(eqn(3)+k) = y(eqn(3)+k) + s(var(3)+j)*dw(2,j,3)
  150       continue
            y(eqn(3)+k) = y(eqn(3)+k) - s(n-2)*w(1,1)
            y(eqn(3)+k) = y(eqn(3)+k) + s(n-1)*w(1,2)
            y(eqn(3)+k) = y(eqn(3)+k) - s(n)*one
  160    continue
  170 continue

      do 250 i = 1,nint
         do 210 e = 1,sigma
            var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
            eqn(e) = sigma + (i-1)*npi + sigma*cpts + sumdeg(e-1)
            do 200 m = 1,deg(e)
               w(m,e) = zero
               do 180 j = m,deg(e)
                  w(m,e) = w(m,e) + rhnfhk(1,0,j-m,j-m)*x(var(e)+j)
                  dw(m,j,e) = rhnfhk(1,0,j-m,j-m)
  180          continue
               do 190 j = 1,cpts
                  w(m,e) = w(m,e) + x(var(e)+deg(e)+j)*
     +            rhnfhk(1,0,deg(e)+j-m,deg(e)-m+1)
                  dw(m,deg(e)+j,e) = rhnfhk(1,0,deg(e)+j-m,deg(e)-m+1)
  190          continue
  200       continue
  210    continue
         if (i.eq.nint) go to 260

         do 240 e = 1,sigma
            do 230 m = 1,deg(e)
               y(eqn(e)+m) = y(eqn(e)+m) + s(var(e)+npi+m)
               do 220 j = 1,cpts + deg(e)
                  y(eqn(e)+m) = y(eqn(e)+m) - s(var(e)+j)*dw(m,j,e)
  220          continue
  230       continue
  240    continue
  250 continue

  260 continue

      var(1) = n - npi - 3
      var(2) = var(1) + 5
      var(3) = var(2) + 5
      do 270 j = 1,5
         y(n-2) = y(n-2) + s(var(1)+j)*dw(1,j,1)
         y(n-1) = y(n-1) + s(var(2)+j)*dw(1,j,2)
         y(n) = y(n) + s(var(3)+j)*dw(1,j,3)
  270 continue

      return

      end
