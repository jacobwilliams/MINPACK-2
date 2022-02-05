      subroutine dsfisp(nx,ny,nnz,indrow,indcol)
      integer nx,ny,nnz
      integer indrow(*),indcol(*)
c     **********
c
c     Subroutine dsfisp
c
c     This subroutine defines the sparsity pattern of the Jacobian 
c     matrix for the Solid Fuel Ignition  Problem. Given the number 
c     of grid points nx and ny dbrssp provides column and row indices 
c     for the nonzeros in the Jacobian. The number of nonzeros, nnz, 
c     is also provided.
c
c     The subroutine statement is:
c
c       subroutine dsfisp(nx,ny,nnz,indrow,indcol)
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
c       nnz is an integer variable.
c         On entry nnz need not be specified
c         On exit nnz is set to the number of nonzeros in the 
c            sparsity structure. 
c
c       indrow is an integer array of dimension at least nnz.
c         On entry indrow need not be specified.
c         On exit indrow contains the row indices of the nonzeros
c            in the Jacobian matrix.
c
c       indcol is an integer array of dimension at least nnz.
c         On entry indcol need not be specified.
c         On exit indcol contains the column indices of the nonzeros
c            in the Jacobian matrix.
c
c     MINPACK-2 Project. March 1992.
c     Argonne National Laboratory and University of Minnesota.
c     R.S. Maier
c
c     **********
      integer i,j,k

c     Compute the sparsity structure.

      nnz = 0
      do 20 j = 1, ny
         do 10 i = 1, nx
            k = (j - 1)*nx + i
            nnz = nnz + 1
            indrow(nnz) = k
            indcol(nnz) = k
            if (i .ne. 1) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k - 1
            endif
            if (i .ne. nx) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k + 1
            endif
            if (j .ne. 1) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k - nx
            endif
            if (j .ne. ny) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k + nx
            endif
   10    continue
   20 continue

      return

      end
