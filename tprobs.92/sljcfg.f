      subroutine sljcfg(n,x,f,fgrad,task,ndim,natoms)
      character*(*) task
      integer n,natoms,ndim
      real f
      real x(n),fgrad(n)
c     **********
c
c     Subroutine sljcfg
c
c     This subroutine computes the function and gradient of the Leonard-
c     Jones Clusters (Molecular Conformation)  Problem.  This problem
c     arises in the study of low-energy states of proteins and in the
c     study of cluster statics.  This subroutine can be used for both
c     the two and three-dimensional cases.
c
c     The subroutine statement is:
c
c       subroutine sljcfg(n,x,f,fgrad,task,natoms,ndim)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c            For the 2-dimensional problem n = 2*natoms.
c            For the 3-dimensional problem n = 3*natoms.
c         On exit n is unchanged.
c
c       x is a real array of dimension n.
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
c       fgrad is a real array of dimension n.
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
c       natoms is an integer variable.
c         On entry natoms is the number of atoms in the cluster.
c         On exit natoms is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, R.S. Maier and G.L. Xue
c
c     **********
      real zero,p5,one,two,three,six
      parameter (zero=0.0,p5=0.5,one=1.0,two=2.0,three=3.0,six=6.0)

      integer i,j,k,isqrtn,icrtn,ileft,il,jl,ctr
      real rij,temp,xx,yy,zz

c     Check Input for Errors.

      if (ndim.eq.2 .and. n.ne.2*natoms)
     +   task = 'ERROR: N MUST .EQ. 2*NATOMS'
      if (ndim.eq.3 .and. n.ne.3*natoms)
     +   task = 'ERROR: N MUST .EQ. 3*NATOMS'
      if (ndim.ne.2 .and. ndim.ne.3) task =
     +   'ERROR: NDIM MUST .EQ. 2 OR 3'
      if (task(1:5).eq.'ERROR') return

c     Compute the standard starting point if task = 'XS'

      if (task.eq.'XS') then
         if (ndim.eq.2) then
            isqrtn = int(sqrt(real(natoms)))
            ileft = natoms - isqrtn*isqrtn
            xx = zero
            yy = zero
            do 20 j = 1,isqrtn + ileft/isqrtn + 1
               do 10 i = 1,min(isqrtn,natoms- (j-1)*isqrtn)
                  ctr = (j-1)*isqrtn + i
                  x(2*ctr-1) = xx
                  x(2*ctr) = yy
                  xx = xx + one
   10          continue
               yy = yy + one
               xx = zero
   20       continue

         else if (ndim.eq.3) then
            icrtn = int((natoms+p5)** (one/three))
            ileft = natoms - icrtn*icrtn*icrtn
            xx = zero
            yy = zero
            zz = zero
            do 50 k = 1,icrtn + ileft/ (icrtn*icrtn) + 1
               jl = min(icrtn, (natoms- (k-1)*icrtn*icrtn)/icrtn+1)
               do 40 j = 1,jl
                  il = min(icrtn,natoms- (k-1)*icrtn*icrtn- (j-1)*icrtn)
                  do 30 i = 1,il
                     ctr = (k-1)*icrtn*icrtn + (j-1)*icrtn + i
                     x(3*ctr-2) = xx
                     x(3*ctr-1) = yy
                     x(3*ctr) = zz
                     xx = xx + one
   30             continue
                  yy = yy + one
                  xx = zero
   40          continue
               yy = zero
               zz = zz + one
   50       continue
         end if

         return

      end if

c     Evaluate the function if task = 'F', the gradient if task = 'G',
c     or both if task = 'FG'.

      if (task.eq.'F' .or. task.eq.'FG') then
         f = zero
         if (ndim.eq.2) then
            do 70 j = 2,natoms
               do 60 i = 1,j - 1
                  xx = x(2*j-1) - x(2*i-1)
                  yy = x(2*j) - x(2*i)
                  rij = xx**2 + yy**2
                  temp = one/rij/rij/rij
                  f = f + temp* (temp-two)
   60          continue
   70       continue

         else if (ndim.eq.3) then
            do 90 j = 2,natoms
               do 80 i = 1,j - 1
                  xx = x(3*j-2) - x(3*i-2)
                  yy = x(3*j-1) - x(3*i-1)
                  zz = x(3*j) - x(3*i)
                  rij = xx*xx + yy*yy + zz*zz
                  temp = one/rij/rij/rij
                  f = f + temp* (temp-two)
   80          continue
   90       continue
         end if

         if (task.eq.'F') return

      end if

c     Compute the gradient

      if (task.eq.'G' .or. task.eq.'FG') then
         do 100 i = 1,n
            fgrad(i) = zero
  100    continue
         if (ndim.eq.2) then
            do 120 j = 2,natoms
               do 110 i = 1,j - 1
                  xx = x(2*j-1) - x(2*i-1)
                  yy = x(2*j) - x(2*i)
                  rij = xx**2 + yy**2
                  temp = one/rij/rij/rij
                  temp = six*temp* (two*temp-two)/rij
                  fgrad(2*j-1) = fgrad(2*j-1) - xx*temp
                  fgrad(2*j) = fgrad(2*j) - yy*temp
                  fgrad(2*i-1) = fgrad(2*i-1) + xx*temp
                  fgrad(2*i) = fgrad(2*i) + yy*temp
  110          continue
  120       continue

         else if (ndim.eq.3) then
            do 140 j = 2,natoms
               do 130 i = 1,j - 1
                  xx = x(3*j-2) - x(3*i-2)
                  yy = x(3*j-1) - x(3*i-1)
                  zz = x(3*j) - x(3*i)
                  rij = xx*xx + yy*yy + zz*zz
                  temp = one/rij/rij/rij
                  temp = six*temp* (two*temp-two)/rij
                  fgrad(3*j-2) = fgrad(3*j-2) - xx*temp
                  fgrad(3*j-1) = fgrad(3*j-1) - yy*temp
                  fgrad(3*j) = fgrad(3*j) - zz*temp
                  fgrad(3*i-2) = fgrad(3*i-2) + xx*temp
                  fgrad(3*i-1) = fgrad(3*i-1) + yy*temp
                  fgrad(3*i) = fgrad(3*i) + zz*temp
  130          continue
  140       continue
         end if

         return

      end if

      end
