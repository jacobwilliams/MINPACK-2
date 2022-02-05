      subroutine scprfj(n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer n,ldfjac
      real x(n),fvec(n),fjac(ldfjac,n)
c     **********
c
c     Subroutine scprfj
c
c     This subroutine computes the function and the Jacobian matrix
c     of the Combustion of Propane - Reduced Formulation - problem of
c     K. Meintjes and A. Morgan. This problem describes the combustion
c     of propane in air using element variables to eliminate square
c     roots in the function evaluations.
c
c     The subroutine statement is:
c
c       subroutine scprfj(n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c            For the reduced combustion of propane problem n = 5.
c         On exit n is unchanged.
c
c       x is a real array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a real array of dimension n.
c         On entry fvec need not be specified.
c         On exit fvec contains the function evaluated at x if
c            task = 'F' or 'FJ'.
c
c       fjac is a real array of dimension (ldfjac,n).
c         On entry fjac need not be specified.
c         On exit fjac contains the Jacobian matrix evaluated at x if
c            task = 'J' or 'FJ'.
c
c       ldfjac is an integer variable.
c          On entry ldfjac is the leading dimension of fjac.
c          On exit ldfjac is unchanged.
c
c       task is a character variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c             'F'     Evaluate the function at x.
c             'J'     Evaluate the Jacobian matrix at x.
c             'FJ'    Evaluate the function and the Jacobian at x.
c             'XS'    Set x to the standard starting point xs.
c             'XL'    Set x to the lower bound xl.
c
c         On exit task is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      real p,rr
      parameter (p=40.0,rr=10.0)
      real zero,one,two,three,four,eight
      parameter (zero=0.0,one=1.0,two=2.0,three=3.0,four=4.0,eight=8.0)

      integer i
      real r(10),k(10)
      real sqrtp

      data k /0.0,0.0,0.0,0.0,1.930E-1,2.597E-3,3.448E-3,1.799E-5,
     +        2.155E-4,3.846E-5/

c     Check input arguments for errors.

      if (n.ne.5) then
         task = 'ERROR: N MUST .EQ. 5'
         return

      end if

c     Initialization.

      sqrtp = sqrt(p)
      r(5) = k(5)
      r(6) = k(6)/sqrtp
      r(7) = k(7)/sqrtp
      r(8) = k(8)/p
      r(9) = k(9)/sqrtp
      r(10) = k(10)/p

c     Compute a lower bound for x if task = 'XL'.

      if (task.eq.'XL') then
         do 10 i = 1,n
            x(i) = zero
   10    continue

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task.eq.'XS') then
         x(1) = 5.0E-3
         x(2) = 5.0E-3
         x(3) = 5.0E-2
         x(4) = 5.0E-1
         x(5) = 5.0E-2

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task.eq.'F' .or. task.eq.'FJ') then
         fvec(1) = x(1)*x(2) + x(1) - three*x(5)
         fvec(2) = two*x(1)*x(2) + x(1) + two*r(10)*x(2)**2 +
     +   x(2)*x(3)**2 + r(7)*x(2)*x(3) + r(9)*x(2)*x(4) + r(8)*x(2) -
     +   rr*x(5)
         fvec(3) = two*x(2)*x(3)**2 + r(7)*x(2)*x(3) +
     +   two*r(5)*x(3)**2 + r(6)*x(3) - eight*x(5)
         fvec(4) = r(9)*x(2)*x(4) + two*x(4)**2 - four*rr*x(5)
         fvec(5) = x(1)*x(2) + x(1) + r(10)*x(2)**2 + x(2)*x(3)**2 +
     +   r(7)*x(2)*x(3) + r(9)*x(2)*x(4) + r(8)*x(2) + r(5)*x(3)**2 +
     +   r(6)*x(3) + x(4)**2 - one

         if (task.eq.'F') return

      end if

      if (task.eq.'J' .or. task.eq.'FJ') then
         fjac(1,1) = x(2) + one
         fjac(2,1) = two*x(2) + one
         fjac(3,1) = zero
         fjac(4,1) = zero
         fjac(5,1) = x(2) + one

         fjac(1,2) = x(1)
         fjac(2,2) = two*x(1) + four*r(10)*x(2) + x(3)**2 + r(7)*x(3) +
     +   r(9)*x(4) + r(8)
         fjac(3,2) = two*x(3)**2 + r(7)*x(3)
         fjac(4,2) = r(9)*x(4)
         fjac(5,2) = x(1) + two*r(10)*x(2) + x(3)**2 + r(7)*x(3) +
     +   r(9)*x(4) + r(8)

         fjac(1,3) = zero
         fjac(2,3) = two*x(2)*x(3) + r(7)*x(2)
         fjac(3,3) = four*x(2)*x(3) + r(7)*x(2) + four*r(5)*x(3) + r(6)
         fjac(4,3) = zero
         fjac(5,3) = two*x(2)*x(3) + r(7)*x(2) + two*r(5)*x(3) + r(6)

         fjac(1,4) = zero
         fjac(2,4) = r(9)*x(2)
         fjac(3,4) = zero
         fjac(4,4) = r(9)*x(2) + four*x(4)
         fjac(5,4) = r(9)*x(2) + two*x(4)

         fjac(1,5) = -three
         fjac(2,5) = -rr
         fjac(3,5) = -eight
         fjac(4,5) = -four*rr
         fjac(5,5) = zero

         return

      end if

      end
