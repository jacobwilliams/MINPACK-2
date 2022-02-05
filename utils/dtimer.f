      subroutine dtimer(ttime)
      double precision ttime
c     *********
c
c     Subroutine dtimer
c
c     This subroutine is used to determine user time. In a typical
c     application, the user time for a code segment requires calls
c     to subroutine dtimer to determine the initial and final time.
c
c     The subroutine statement is
c
c       subroutine dtimer(ttime)
c
c     where
c
c       ttime is a double precision variable.
c         On input ttime need not be specified.
c         On output ttime specifies the user time.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and Jorge J. More'.
c
c     **********

c     Declarations for machines with etime.

      real tarray(2)
      real etime

c     Declarations for machines with mclock (ibm).

c     integer mclock

c     For machines with etime (sun).
c     The first element of the array tarray specifies user time

      ttime = etime(tarray)
      ttime = dble(tarray(1))

c     For machines with mclock (ibm).

c     ttime = dble(mclock())/100.0d0

      end
