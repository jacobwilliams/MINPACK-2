      subroutine sctsfj(m,n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer m,n,ldfjac
      real x(n),fvec(m),fjac(ldfjac,n)
c     **********
c
c     Subroutine sctsfj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the Coating Thickness Standardization problem communicated by
c     Janet Rogers of the National Institute of Standards and
c     Technology. This is a multiple response data fitting problem that
c     arises from the need to nondestructively determine any
c     nonuniformity in the lead-tin coating on samples of standard
c     reference materials. This problem has data which has been supplied
c     by Susannah Schiller of the National Institute of Standards and
c     Technology.
c
c     NOTE: This subroutine assumes that the data is stored in the file
c           'scts.dat'.
c
c     The subroutine statement is:
c
c       subroutine sctsfj(m,n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions.
c            For the coating thickness standardization problem m = 252.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables.
c            For the coating thickness standardization problem n = 134.
c         On exit n is unchanged.
c
c       x is a real array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a real array of dimension m.
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
c
c         On exit task is unchanged.
c
c     MINPACK-2 Project. October 1991.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer nread
      parameter (nread=99)

      integer mdiv4,mdiv2
      parameter (mdiv4=63,mdiv2=126)
      real scale1,scale2,zero,one,hund
      parameter (scale1=4.08,scale2=0.417,zero=0.0,one=1.0,hund=1.0E2)

      integer i,j
      logical first
      real indvar(mdiv4,2),y(mdiv2)

      save first,indvar,y

      data first/.true./

c     Check input arguments for errors.

      if (m.ne.252) task = 'ERROR: M MUST .EQ. 252'
      if (n.ne.134) task = 'ERROR: N MUST .EQ. 134'
      if (task(1:5).eq.'ERROR') return

c     Initialization on the first call to the subroutine.

      if (first) then
         first = .false.
         open (nread,file='scts.dat',status='old')
         do 10 i = 1,mdiv4
            read (nread,FMT=*) indvar(i,1),indvar(i,2),y(i),y(i+mdiv4)
   10    continue
         close (nread)
         do 20 i = 1,mdiv4
            indvar(i,1) = indvar(i,1)/hund
            y(i) = y(i)/hund
            y(i+mdiv4) = y(i+mdiv4)/hund
   20    continue
      end if

c     Compute the standard starting point if task = 'XS'.

      if (task.eq.'XS') then
         x(1) = -8.0
         x(2) = 13.0
         x(3) = 1.2
         x(4) = 0.2
         x(5) = 0.1
         x(6) = 6.0
         x(7) = 5.5
         x(8) = -5.2
         do 30 i = 1,mdiv2
            x(8+i) = zero
   30    continue

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task.eq.'F' .or. task.eq.'FJ') then
         do 40 i = 1,mdiv4
            fvec(i) = x(1) + x(2)* (indvar(i,1)+x(8+i)) +
     +      x(3)* (indvar(i,2)+x(8+i+mdiv4)) +
     +      x(4)* (indvar(i,1)+x(8+i))* (indvar(i,2)+x(8+i+mdiv4)) -
     +      y(i)
            fvec(i+mdiv4) = x(5) + x(6)* (indvar(i,1)+x(8+i)) +
     +      x(7)* (indvar(i,2)+x(8+i+mdiv4)) +
     +      x(8)* (indvar(i,1)+x(8+i))* (indvar(i,2)+x(8+i+mdiv4)) -
     +      y(i+mdiv4)
            fvec(i+2*mdiv4) = scale1*x(8+i)
            fvec(i+3*mdiv4) = scale2*x(8+i+mdiv4)
   40    continue

         if (task.eq.'F') return

      end if

      if (task.eq.'J' .or. task.eq.'FJ') then
         do 60 j = 1,8 + 2*mdiv4
            do 50 i = 1,m
               fjac(i,j) = zero
   50       continue
   60    continue

         do 70 i = 1,mdiv4
            fjac(i,1) = one
            fjac(i,2) = indvar(i,1) + x(8+i)
            fjac(i,3) = indvar(i,2) + x(8+i+mdiv4)
            fjac(i,4) = (indvar(i,1)+x(8+i))* (indvar(i,2)+x(8+i+mdiv4))
            fjac(i,8+i) = x(2) + x(4)* (indvar(i,2)+x(8+i+mdiv4))
            fjac(i,8+i+mdiv4) = x(3) + x(4)* (indvar(i,1)+x(8+i))
            fjac(i+mdiv4,5) = one
            fjac(i+mdiv4,6) = indvar(i,1) + x(8+i)
            fjac(i+mdiv4,7) = indvar(i,2) + x(8+i+mdiv4)
            fjac(i+mdiv4,8) = (indvar(i,1)+x(8+i))*
     +      (indvar(i,2)+x(8+i+mdiv4))
            fjac(i+mdiv4,8+i) = x(6) + x(8)* (indvar(i,2)+x(8+i+mdiv4))
            fjac(i+mdiv4,8+i+mdiv4) = x(7) + x(8)* (indvar(i,1)+x(8+i))
            fjac(i+2*mdiv4,8+i) = scale1
            fjac(i+3*mdiv4,8+i+mdiv4) = scale2
   70    continue

         return

      end if

      end
