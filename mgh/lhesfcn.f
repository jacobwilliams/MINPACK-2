      subroutine hesfcn(n,x,h,ldh,nprob)
      integer n,ldh,nprob
      double precision x(n),h(ldh)
c     **********
c
c     subroutine hesfcn
c
c     This subroutine defines the hessian matrices of eighteen
c     nonlinear unconstrained minimization problems. The problem
c     dimensions are as described in the prologue comments of objfcn.
c     The upper triangle of the (symmetric) hessian matrix is
c     computed columnwise and stored as a one-dimensional array.
c
c     The subroutine statement is
c
c       subroutine hesfcn(n,x,h,ldh,nprob)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c       h is an output array of length (n*(n+1))/2 which contains 
c         the hessian matrix of the nprob objective function 
c         evaluated at x.
c
c       ldh is a positive integer input variable not less than
c         (n*(n+1))/2 which specifies the dimension of the array h.
c
c       nprob is a positive integer input variable which defines the
c         number of the problem. nprob must not exceed 18.
c
c     Subprograms called
c
c       FORTRAN-supplied ... abs,atan,cos,dble,exp,log,sign,sin,sqrt
c
c     Argonne National Laboratory. MINPACK Project. march 1980.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
c
c     **********
      integer i,iev,j,k
      integer jj,jjp1,jjp2,jjp3,jn,jn1,j1j1,kj,kn,kn1,nn,ntr,n1n,n1n1
      equivalence (kj,jj,nn,n1n1)
      double precision ap,arg,cp0001,cp1,cp25,cp5,c1p5,c2p25,
     *                 c2p625,c3p5,c19p8,c25,c29,c100,c200,c10000,d1,
     *                 d2,eight,fifty,five,four,one,r,s1,s2,s3,t,t1,
     *                 t2,t3,ten,th,three,tpi,twenty,two,zero
      double precision d3,r1,r2,r3,u1,u2,v,v1,v2
      double precision fvec(50),fvec1(50),y(15)
      double precision six,xnine,twelve,c120,c200p2,c202,c220p2,c360,
     *                 c400,c1200
      data six,xnine,twelve,c120,c200p2,c202,c220p2,c360,c400,c1200
     *     /6.0d0,9.0d0,1.2d1,1.2d2,2.002d2,2.02d2,2.202d2,3.6d2,
     *      4.0d2,1.2d3/
      data zero,one,two,three,four,five,eight,ten,twenty,fifty
     *     /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,8.0d0,1.0d1,2.0d1,
     *      5.0d1/
      data cp0001,cp1,cp25,cp5,c1p5,c2p25,c2p625,c3p5,c19p8,c25,c29,
     *     c100,c200,c10000
     *     /1.0d-4,1.0d-1,2.5d-1,5.0d-1,1.5d0,2.25d0,2.625d0,3.5d0,
     *      1.98d1,2.5d1,2.9d1,1.0d2,2.0d2,1.0d4/
      data ap /1.0d-5/
      data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),
     *     y(12),y(13),y(14),y(15)
     *     /9.0d-4,4.4d-3,1.75d-2,5.4d-2,1.295d-1,2.42d-1,3.521d-1,
     *      3.989d-1,3.521d-1,2.42d-1,1.295d-1,5.4d-2,1.75d-2,4.4d-3,
     *      9.0d-4/
c
c     Hessian routine selector.
c
      go to (10,20,60,100,110,170,210,290,330,380,390,450,490,580,620,
     *       660,670,680), nprob
c
c     Helical valley function.
c
   10 continue
      tpi = eight*atan(one)
      th = sign(cp25,x(2))
      if (x(1) .gt. zero) th = atan(x(2)/x(1))/tpi
      if (x(1) .lt. zero) th = atan(x(2)/x(1))/tpi + cp5
      arg = x(1)**2 + x(2)**2
      r = sqrt(arg)
      t = x(3) - ten*th
      s1 = ten*t/(tpi*arg)
      t1 = ten/tpi
      t2 = t1/arg
      t3 = (x(1)/r - t1*t2*x(1) - two*x(2)*s1)/arg
      h(1) = c200
     *         *(one - x(2)/arg*(x(2)/r - t1*t2*x(2) + two*x(1)*s1))
      h(2) = c200*(s1 + x(2)*t3)
      h(3) = c200*(one - x(1)*t3)
      h(4) = c200*t2*x(2)
      h(5) = -c200*t2*x(1)
      h(6) = c202
      return
c
c     Biggs exp6 function.
c
   20 continue
      do 40 kj = 1, 21
         h(kj) = zero
   40    continue
      do 50 i = 1, 13
         d1 = dble(i)/ten
         d2 = exp(-d1) - five*exp(-ten*d1) + three*exp(-four*d1)
         s1 = exp(-d1*x(1))
         s2 = exp(-d1*x(2))
         s3 = exp(-d1*x(5))
         t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2
         th = d1*t
         r1 = d1*s1
         r2 = d1*s2
         r3 = d1*s3
         h(1) = h(1) + r1*(th + x(3)*r1)
         h(2) = h(2) - r1*r2
         h(3) = h(3) - r2*(th - x(4)*r2)
         h(4) = h(4) - s1*(th + x(3)*r1)
         h(6) = h(6) + s1**2
         h(7) = h(7) + r1*s2
         h(8) = h(8) + s2*(th - x(4)*r2)
         h(9) = h(9) - s1*s2
         h(10) = h(10) + s2**2
         h(11) = h(11) + r1*r3
         h(12) = h(12) - r2*r3
         h(15) = h(15) + r3*(th + x(6)*r3)
         h(16) = h(16) - r1*s3
         h(17) = h(17) + r2*s3
         h(18) = h(18) + s1*s3
         h(19) = h(19) - s2*s3
         h(20) = h(20) - s3*(th + x(6)*r3)
         h(21) = h(21) + s3**2
   50    continue
      h(1) = two*x(3)*h(1)
      h(2) = two*x(3)*x(4)*h(2)
      h(3) = two*x(4)*h(3)
      h(4) = two*h(4)
      h(5) = two*x(4)*h(7)
      h(6) = two*h(6)
      h(7) = two*x(3)*h(7)
      h(8) = two*h(8)
      h(9) = two*h(9)
      h(10) = two*h(10)
      h(11) = two*x(3)*x(6)*h(11)
      h(12) = two*x(4)*x(6)*h(12)
      h(13) = two*x(6)*h(16)
      h(14) = two*x(6)*h(17)
      h(15) = two*x(6)*h(15)
      h(16) = two*x(3)*h(16)
      h(17) = two*x(4)*h(17)
      h(18) = two*h(18)
      h(19) = two*h(19)
      h(20) = two*h(20)
      h(21) = two*h(21)
      return
c
c     Gaussian function.
c
   60 continue
      do 80 kj = 1, 6
         h(kj) = zero
   80    continue
      do 90 i = 1, 15
         d1 = cp5*dble(i-1)
         d2 = c3p5 - d1 - x(3)
         arg = -cp5*x(2)*d2**2
         r = exp(arg)
         t = x(1)*r - y(i)
         s1 = r*t
         s2 = d2*s1
         t1 = s2 + d2*x(1)*r**2
         t2 = d2*t1
         h(1) = h(1) + r**2
         h(2) = h(2) - t2
         h(3) = h(3) + d2**2*t2
         h(4) = h(4) + t1
         h(5) = h(5) + two*s2 - d2*x(2)*t2
         h(6) = h(6) + x(2)*t2 - s1
   90    continue
      h(1) = two*h(1)
      h(2) = h(2)
      h(3) = cp5*x(1)*h(3)
      h(4) = two*x(2)*h(4)
      h(5) = x(1)*h(5)
      h(6) = two*x(1)*x(2)*h(6)
      return
c
c     Powell badly scaled function.
c
  100 continue
      t1 = c10000*x(1)*x(2) - one
      s1 = exp(-x(1))
      s2 = exp(-x(2))
      t2 = s1 + s2 - one - cp0001
      h(1) = two*((c10000*x(2))**2 + s1*(s1 + t2))
      h(2) = two*(c10000*(one + two*t1) + s1*s2)
      h(3) = two*((c10000*x(1))**2 + s2*(s2 + t2))
      return
c
c     Box 3-dimensional function.
c
  110 continue
      do 130 kj = 1, 6
         h(kj) = zero
  130    continue
      do 140 i = 1, 10
         d1 = dble(i)
         d2 = d1/ten
         s1 = exp(-d2*x(1))
         s2 = exp(-d2*x(2))
         s3 = exp(-d2) - exp(-d1)
         t = s1 - s2 - s3*x(3)
         th = d2*t
         r1 = d2*s1
         r2 = d2*s2
         h(1) = h(1) + r1*(th + r1)
         h(2) = h(2) - r1*r2
         h(3) = h(3) - r2*(th - r2)
         h(4) = h(4) + r1*s3
         h(5) = h(5) - r2*s3
         h(6) = h(6) + s3**2
  140    continue
      do 160 kj = 1, 6
         h(kj) = two*h(kj)
  160    continue
      return
c
c     Variably dimensioned function.
c
  170 continue
      t1 = zero
      do 180 j = 1, n
         t1 = t1 + dble(j)*(x(j) - one)
  180    continue
c     t = t1*(one + two*t1**2)
      t2 = two + twelve*t1**2
      kj = 0
      do 200 j = 1, n
         do 190 k = 1, j
            kj = kj + 1
            h(kj) = dble(k*j)*t2
  190       continue
         h(jj) = h(jj) + two
  200    continue
      return
c
c     Watson function.
c
  210 continue
      ntr = (n*(n + 1))/2
      do 230 kj = 1, ntr
         h(kj) = zero
  230    continue
      do 280 i = 1, 29
         d1 = dble(i)/c29
         s1 = zero
         d2 = one
         do 240 j = 2, n
            s1 = s1 + dble(j-1)*d2*x(j)
            d2 = d1*d2
  240       continue
         s2 = zero
         d2 = one
         do 250 j = 1, n
            s2 = s2 + d2*x(j)
            d2 = d1*d2
  250       continue
         t = s1 - s2**2 - one
         s3 = two*d1*s2
         d2 = two/d1
         th = two*d1**2*t
         kj = 0
         do 270 j = 1, n
            v = dble(j-1) - s3
            d3 = one/d1
            do 260 k = 1, j
               kj = kj + 1
               h(kj) = h(kj) + d2*d3*(v*(dble(k-1) - s3) - th)
               d3 = d1*d3
  260          continue
            d2 = d1*d2
  270       continue
  280    continue
      t1 = x(2) - x(1)**2 - one
      h(1) = h(1) + eight*x(1)**2 + two - four*t1
      h(2) = h(2) - four*x(1)
      h(3) = h(3) + two
      return
c
c     Penalty function I.
c
  290 continue
      t1 = -cp25
      do 300 j = 1, n
         t1 = t1 + x(j)**2
  300    continue
      d1 = two*ap
      th = four*t1
      kj = 0
      do 320 j = 1, n
         t2 = eight*x(j)
         do 310 k = 1, j
            kj = kj + 1
            h(kj) = x(k)*t2
  310       continue
         h(jj) = h(jj) + d1 + th
  320    continue
      return
c
c     Penalty function II.
c
  330 continue
      t1 = -one
      do 340 j = 1, n
         t1 = t1 + dble(n-j+1)*x(j)**2
  340    continue
      d1 = exp(cp1)
      d2 = one
      th = four*t1
      kj = 0
      do 370 j = 1, n
         t2 = eight*dble(n-j+1)*x(j)
         do 350 k = 1, j
            kj = kj + 1
            h(kj) = dble(n-k+1)*x(k)*t2
  350       continue
         h(jj) = h(jj) + dble(n-j+1)*th
         s1 = exp(x(j)/ten)
         if (j .gt. 1) then
            s3 = s1 + s2 - d2*(d1 + one)
            h(jj) = h(jj) + ap*s1*(s3 + three*s1 - one/d1)/fifty
            h(jj-1) = h(jj-1) + ap*s1*s2/fifty
            h(j1j1) = h(j1j1) + ap*s2*(s2 + s3)/fifty
            end if
         s2 = s1
         d2 = d1*d2
         j1j1 = jj
  370    continue
      h(1) = h(1) + two
      return
c
c     Brown badly scaled function.
c
  380 continue
c     t1 = x(1) - c1pd6
c     t2 = x(2) - c2pdm6
      t3 = x(1)*x(2) - two
      h(1) = two*(one + x(2)**2)
      h(2) = four*(one + t3)
      h(3) = two*(one + x(1)**2)
      return
c
c     Brown and Dennis function.
c
  390 continue
      do 410 kj = 1, 10
         h(kj) = zero
  410    continue
      do 420 i = 1, 20
         d1 = dble(i)/five
         d2 = sin(d1)
         t1 = x(1) + d1*x(2) - exp(d1)
         t2 = x(3) + d2*x(4) - cos(d1)
         t = t1**2 + t2**2
c        s1 = t1*t
c        s2 = t2*t
         s3 = two*t1*t2
         r1 = t + two*t1**2
         r2 = t + two*t2**2
         h(1) = h(1) + r1
         h(2) = h(2) + d1*r1
         h(3) = h(3) + d1**2*r1
         h(4) = h(4) + s3
         h(5) = h(5) + d1*s3
         h(6) = h(6) + r2
         h(7) = h(7) + d2*s3
         h(8) = h(8) + d1*d2*s3
         h(9) = h(9) + d2*r2
         h(10) = h(10) + d2**2*r2
  420    continue
      do 440 kj = 1, 10
         h(kj) = four*h(kj)
  440    continue
      return
c
c     Gulf research and development function.
c
  450 continue
      do 470 kj = 1, 6
         h(kj) = zero
  470    continue
      d1 = two/three
      do 480 i = 1, 99
         arg = dble(i)/c100
         r = (-fifty*log(arg))**d1 + c25 - x(2)
         t1 = abs(r)**x(3)/x(1)
         t2 = exp(-t1)
         t = t2 - arg
         s1 = t1*t2*t
         s2 = t1*(s1 + t2*(t1*t2 - t))
         r1 = log(abs(r))
         r2 = r1*s2
         h(1) = h(1) + s2 - s1
         h(2) = h(2) + s2/r
         h(3) = h(3) + (s1 + x(3)*s2)/r**2
         h(4) = h(4) - r2
         h(5) = h(5) + (s1 - x(3)*r2)/r
         h(6) = h(6) + r1*r2
  480    continue
      h(1) = two*h(1)/x(1)**2
      h(2) = two*x(3)*h(2)/x(1)
      h(3) = two*x(3)*h(3)
      h(4) = two*h(4)/x(1)
      h(5) = two*h(5)
      h(6) = two*h(6)
      return
c
c     Trigonometric function.
c
  490 continue
      u2 = cos(x(n))
      s1 = u2
      if (n .gt. 1) then
         u1 = cos(x(n-1))
         s1 = s1 + u1
         ntr = ((n - 2)*(n - 1))/2
         jn1 = ntr
         do 500 j = 1, n-2
            jn1 = jn1 + 1
            h(jn1) = cos(x(j))
            s1 = s1 + h(jn1)
  500       continue
         end if
      v2 = sin(x(n))
      s2 = dble(2*n) - v2 - s1 - dble(n)*u2
      r2 = dble(2*n)*v2 - u2
      kj = 0
      if (n .gt. 1) then
         v1 = sin(x(n-1))
         s2 = s2 + dble(2*n-1) - v1 - s1 - dble(n-1)*u1
         r1 = dble(2*n-1)*v1 - u1
         jn1 = ntr
         do 520 j = 1, n-2
            jn = jn1 + n
            jn1 = jn1 + 1
            h(jn) = sin(x(j))
            t = dble(n+j) - h(jn) - s1 - dble(j)*h(jn1)
            s2 = s2 + t
  520       continue
         jn1 = ntr
         do 540 j = 1, n-2
            jn = jn1 + n
            jn1 = jn1 + 1
            v = dble(j)*h(jn1) + h(jn)
            t = dble(n+j) - s1 - v
            t1 = dble(n+j)*h(jn) - h(jn1)
            kn1 = ntr
            do 530 k = 1, j
               kj = kj + 1
               kn = kn1 + n
               kn1 = kn1 + 1
               th = dble(k)*h(kn) - h(kn1)
               h(kj) = two*(h(kn)*t1 + h(jn)*th)
  530          continue
            h(jj) = h(jj) + two*(h(jn1)*s2 + v*t + th**2)
  540       continue
         do 550 k = 1, n-2
            kn = kj + n
            kj = kj + 1
            kn1 = kj
            th = dble(k)*h(kn) - h(kn1)
            h(kn1) = two*(h(kn)*r1 + v1*th)
            h(kn) = two*(h(kn)*r2 + v2*th)
  550       continue
         v = dble(n-1)*u1 + v1
         t = dble(2*n-1) - s1 - v
         th = dble(n-1)*v1 - u1
         n1n = kj + n
         kj = kj + 1
         h(n1n1) = two*(v1*(r1 + th) + u1*s2 + v*t + th**2)
         h(n1n) = two*(v1*r2 + v2*th)
         end if
      v = dble(n)*u2 + v2
      t = dble(2*n) - s1 - v
      th = dble(n)*v2 - u2
      kj = kj + n
      h(nn) = two*(v2*(r2 + th) + u2*s2 + v*t + th**2)
      return
c
c     Extended Rosenbrock function.
c
  580 continue
      ntr = (n*(n + 1))/2
      do 600 kj = 1, ntr
         h(kj) = zero
  600    continue
      jjp1 = -1
      do 610 j = 1, n, 2
c        t1 = one - x(j)
         jj = jjp1 + j + 1
         jjp1 = jj + j
         h(jj) = c1200*x(j)**2 - c400*x(j+1) + two
         h(jjp1) = -c400*x(j)
         h(jjp1+1) = c200
  610    continue
      return
c
c     Extended Powell function.
c
  620 continue
      ntr = (n*(n + 1))/2
      do 640 kj = 1, ntr
         h(kj) = zero
  640    continue
      jjp3 = -3
      do 650 j = 1, n, 4
c        t = x(j) + ten*x(j+1)
c        t1 = x(j+2) - x(j+3)
c        s1 = five*t1
         t2 = x(j+1) - two*x(j+2)
c        s2 = four*t2**3
         t3 = x(j) - x(j+3)
c        s3 = twenty*t3**3
         r2 = twelve*t2**2
         r3 = c120*t3**2
         jj = jjp3 + j + 3
         jjp1 = jj + j
         jjp2 = jjp1 + j + 1
         jjp3 = jjp2 + j + 2
         h(jj) = two + r3
         h(jjp1) = twenty
         h(jjp1+1) = c200 + r2
         h(jjp2+1) = -two*r2
         h(jjp2+2) = ten + four*r2
         h(jjp3) = -r3
         h(jjp3+2) = -ten
         h(jjp3+3) = ten + r3
  650    continue
      return
c
c     Beale function.
c
  660 continue
      s1 = one - x(2)
      t1 = c1p5 - x(1)*s1
      s2 = one - x(2)**2
      t2 = c2p25 - x(1)*s2
      s3 = one - x(2)**3
      t3 = c2p625 - x(1)*s3
      h(1) = two*(s1**2 + s2**2 + s3**2)
      h(2) = two
     *         *(t1 + x(2)*(two*t2 + three*x(2)*t3)
     *           - x(1)*(s1 + x(2)*(two*s2 + three*x(2)*s3)))
      h(3) = two*x(1)
     *         *(x(1) + two*t2
     *           + x(2)*(six*t3 + x(1)*x(2)*(four + xnine*x(2)**2)))
      return
c
c     Wood function.
c
  670 continue
      s1 = x(2) - x(1)**2
c     s2 = one - x(1)
c     s3 = x(2) - one
      t1 = x(4) - x(3)**2
c     t2 = one - x(3)
c     t3 = x(4) - one
      h(1) = c400*(two*x(1)**2 - s1) + two
      h(2) = -c400*x(1)
      h(3) = c220p2
      h(4) = zero
      h(5) = zero
      h(6) = c360*(two*x(3)**2 - t1) + two
      h(7) = zero
      h(8) = c19p8
      h(9) = -c360*x(3)
      h(10) = c200p2
      return
c
c     Chebyquad function.
c
  680 continue
      do 690 i = 1, n
         fvec(i) = zero
  690    continue
      do 710 j = 1, n
         t1 = one
         t2 = two*x(j) - one
         t = two*t2
         do 700 i = 1, n
            fvec(i) = fvec(i) + t2
            th = t*t2 - t1
            t1 = t2
            t2 = th
  700       continue
  710    continue
      d1 = one/dble(n)
      iev = -1
      do 720 i = 1, n
         fvec(i) = d1*fvec(i)
         if (iev .gt. 0) fvec(i) = fvec(i) + one/(dble(i)**2 - one)
         iev = -iev
  720    continue
      kj = 0
      do 770 j = 1, n
         do 730 k = 1, j
            kj = kj + 1
            h(kj) = zero
  730       continue
         t1 = one
         t2 = two*x(j) - one
         t = two*t2
         s1 = zero
         s2 = two
         r1 = zero
         r2 = zero
         do 740 i = 1, n
            h(jj) = h(jj) + fvec(i)*r2
            th = eight*s2 + t*r2 - r1
            r1 = r2
            r2 = th
            fvec1(i) = d1*s2
            th = four*t2 + t*s2 - s1
            s1 = s2
            s2 = th
            th = t*t2 - t1
            t1 = t2
            t2 = th
  740       continue
         kj = kj - j
         do 760 k = 1, j
            kj = kj + 1
            v1 = one
            v2 = two*x(k) - one
            v = two*v2
            u1 = zero
            u2 = two
            do 750 i = 1, n
               h(kj) = h(kj) + fvec1(i)*u2
               th = four*v2 + v*u2 - u1
               u1 = u2
               u2 = th
               th = v*v2 - v1
               v1 = v2
               v2 = th
  750          continue
  760       continue
  770    continue
      d2 = two*d1
      ntr = (n*(n + 1))/2
      do 790 kj = 1, ntr
         h(kj) = d2*h(kj)
  790    continue
  800 continue
      return
c
c     Last card of subroutine hesfcn.
c
      end
      subroutine expand(n,alin,amat,lda)
      integer n,lda
      double precision alin(*),amat(lda,n)
c     **********
c
c     This subroutine expands the upper triangle of a symmetric 
c     matrix A stored as a linear array to full two-dimensional form.
c
c     The subroutine statement is
c
c       subroutine expand(n,alin,amat,lda)
c
c     where
c
c       n is a positive integer input variable set to the order of A.
c
c       alin is an input array of length (n*(n+1))/2 which contains
c         the upper triangle of A stored linearly by columns.
c
c       amat is an output n by n array which contains the expanded
c         full matrix A. If desired, amat can be equivalent to alin.
c
c       lda is a positive integer input variable not less than n
c         which specifies the leading dimension of the array amat.
c
c     Argonne National Laboratory. MINPACK Project. February 1984.
c     Burton S. Garbow
c
c     **********
      integer j,k,kj
c
c     Copy the linear array into the upper triangle of the matrix.
c
      kj = (n*(n+1))/2
      do 20 j = n, 1, -1
         do 10 k = j, 1, -1
            amat(k,j) = alin(kj)
            kj = kj - 1
   10       continue
   20    continue
c
c     Reflect the upper triangle to fill the square matrix.
c
      do 40 j = 1, n
         do 30 k = j+1, n
            amat(k,j) = amat(j,k)
   30       continue
   40    continue
      return
c
c     Last card of subroutine expand.
c
      end
