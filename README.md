
This directory contains MINPACK-2 software for 
the solution of systems of nonlinear equations, 
nonlinear least squares problems, and minimization problems.

The initial software development work concentrated on the
software for the solution of medium-scale problems on vector
and shared-memory architectures. A Levenberg-Marquardt method 
for nonlinear least squares and a trust region Newton method 
for minimization problems are currently available.
The BLAS and LAPACK routines have been used wherever possible
to enhance performance on a variety of architectures.
Current work centers on software for large-scale problems.
A variable-storage variable metric method
and a trust region Newton method are currently available.
Future work will focus on software for distributed memory architectures.

The codes use a reverse communication interface for flexibility and 
ease of use. The codes are accompanied by auxiliary subroutines.
For example, subroutines for computing the covariance
matrix in least squares problems, for estimating
the noise in the function, for estimating the difference
parameter, and for estimating sparse Jacobian and Hessian matrices.

Software will be added to this directory as it becomes available.
The MINPACK-2 test problem collection is available in 
single and double precision; all other software is in double precision.
Single precision versions of the software can be generated easily,
but are only available on request.

The README file in each subdirectory contains additional
information on the status of the software.

Please report any comments and errors to: 

	Jorge J. More' 
        Mathematics and Computer Science Division
        Argonne National Laboratory
        9700 South Cass Avenue
        Argonne, IL  60439
        email: more@mcs.anl.gov

Last modification: September 1, 1997

Organization
============

   csrch
   -----

   This directory has the files for the line search subroutine.

      dcsrch.f
      dcstep.f
      csrch.tar.gz

   vmlm
   -----

   This directory has the files for the limited memory 
   variable metric method. 

      dvmlm.f
      dlmmv.f
      vmlm.tar.gz

   tprobs
   -------

   This directory contains the MINPACK-2 test problem collection.
```
   dficfj.f dficjs.f dficsp.f --- flow in a channel
   dsfdfj.f dsfdjs.f dsfdsp.f --- swirling flow between disks
   dierfj.f dierjs.f diersp.f --- incompressible elastic rods
   dsfifj.f dsfijs.f dsfisp.f --- solid fule ignition
   dfdcfj.f dfdcjs.f dfdcsp.f --- flow in a driven cavity
   dhhdfj.f                   --- human heart dipole
   dcpffj.f                   --- combustion of propane: full formulation
   dcprfj.f                   --- combustion of propane: reduced formulation

   diacfj.f                   --- Isomerization of alpha-pinene: collocation
   diadfj.f                   --- Isomerization of alpha-pinene: direct
   diaofj.f                   --- Isomerization of alpha-pinene: constraints
   diarfj.f                   --- Isomerization of alpha-pinene: residuals
   dctsfj.f                   --- Coating thickness standardization
   dedffj.f                   --- Exponential data fitting
   dgdffj.f                   --- Gaussian data fitting
   datrfj.f                   --- Analysis of thermistor resistance
   daerfj.f                   --- Analysis on an enzyme reaction
   dchqfj.f                   --- Chebychev quadrature

   deptfg.f depths.f deptsp.f --- elastic-plastic torsion
   dpjbfg.f dpjbhs.f dpjbsp.f --- pressure distribution in a journal bearing
            dpjbds.f
   dmsafg.f dmsahs.f dmsasp.f --- minimal surfaces
   dmsabc.f
   dodcfg.f dodchs.f dodcsp.f --- optimal design with composites
   dodcps.f
   dljcfg.f                   --- Lennard-Jones clusters
   dgl1fg.f dgl1hs.f dgl1sp.f --- 1-d Ginzburg-Landau
   dsscfg.f dsschs.f dsscsp.f --- steady-state combustion
   dgl2fg.f dgl2hs.f dgl2sp.f --- 2-d Ginzburg-Landau
            dgl2co.f
```
   sockets
   -------

   This directory contains the sockets files for
   the MINPACK-2 test problem collection.
```
   fic.m --- flow in a channel
   sfd.m --- swirling flow between disks
   ier.m --- incompressible elastic rods
   sfi.m --- solid fule ignition
   fdc.m --- flow in a driven cavity
   ept.m --- elastic-plastic torsion
   pjb.m --- pressure distribution in a journal bearing
   msa.m --- minimal surfaces
   odc.m --- optimal design with composites
   gl1.m --- inhomogeneous superconductors: 1-d Ginzburg-Landau
   ssc.m --- steady-state combustion
   gl2.m --- inhomogeneous superconductors: 2-d Ginzburg-Landau
```
   utils
   -----

   This directory contains utilty programs needed by MINPACK-2.
```
   dnrm2.f  --- computes the l2 norm
   dpmeps.f --- determines the machine precision
   dtimer.f --- computes the current time
   surn01.f --- computes a (single precision) uniform random number in (0,1)
```
   lapack
   ------

   This directory contains the LAPACK files needed by MINPACK-2.
```
   dgeqrf.f dgeqr2.f --- Computes the QR factorization
   dpotf2.f dpotrf.f --- Computes the Cholesky factorization
   dormqr.f dorm2r.f --- Aplies the matrix Q
   dorg2r.f dorgqr.f --- Generates the matrix Q
   dlarf.f  dlarfb.f
   dlarfg.f dlarft.f --- Applies a reflector
   dlamc1.f dlamc2.f
   dlamc3.f dlamc4.f
   dlamc5.f          --- Computes machine constants
   ilaenv.f          --- Returns problem-dependent parametrs
   dlapy2.f          --- Computation of sqrt(x**2+y**2)
   lsame.f           --- Compares strings
   xerbla.f          --- xerbla.f
```
   blas
   ----

   This directory contains the blas files needed by MINPACK-2.
```
   Level 1: dasum.f daxpy.f dcopy.f ddot.f dscal.f dswap.f idamax.f

   Level 2: dger.f dtrsv.f xerbla.f

   Level 3: dgemm.f dsyrk.f dtrsm.f
```
   tprobs.92
   ---------

   This directory contains the 1992 version of the
   MINPACK-2 test problem collection.

   mgh
   ---

   This directory has the More', Garbow, Hillstrom MINPACK-1 test problems.
```
   vecfcn.f --- nonlinear equations problems
   vecjac.f --- nonlinear equations Jacobain matrices
   ssqfcn.f --- nonlinear least squares problems
   ssqjac.f --- nonlinear least squares Jacobian matrices
   objfcn.f --- unconstrained minimization problems
   grdfcn.f --- unconstrained minimization gradients
   hesfcn.f --- unconstrained minimization Hessians
  lhesfcn.f --- unconstrained minimization Hessians (compressed form)
   neq.data --- nonlinear equations sample data
   ssq.data --- nonlinear least squares sample data
    um.data --- unconstrained minimization sample data
```

