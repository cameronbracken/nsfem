!program slaptest
!    implicit none
!    integer,parameter::n=5,nelt=11
!    integer::i
!    integer,dimension(nelt)::ia,ja
!    double precision,dimension(nelt)::a
!    double precision,dimension(n)::b
!
!    
!    a  = (/51d0,12d0,11d0,33d0,15d0,53d0,55d0,22d0,35d0,44d0,21d0/) 
!    ia = (/5   ,1   ,1   ,3   ,1   ,5   ,5   ,2   ,3   ,4   ,2/) 
!    ja = (/1   ,2   ,1   ,3   ,5   ,3   ,5   ,2   ,5   ,4   ,1/) 
!    b  = (/1d0,2d0,3d0,4d0,5d0/) 
!
!    call slapit(a,ia,ja,b,nelt,n)
!    do i=1,n
!        write(*,*)b(i)
!    end do
!
!
!end program slaptest


subroutine slapit(a,ia,ja,b,nelt,n,ierr,quiet)

!!****************************************************************************80
!!
!!slapit is a modification of DLAP_DGMRES_PRB written by John Burkhardt as an
!!interface to the DLAP GMRES solver. DLAP is the ouble precision implementation
!!of SLAP the Sparse Linear Algebra Package. 
!!
!!This routine provides an easy calling interface to DGEMRS which has a large
!!and nasty call statement.  Basically all of the things that are not needed for
!!a normal solution to Ax=b are hard coded into the routine.
!!
!!~Inputs~
!!
!!nelt  -[integer] the number of nonzero entries in A
!!n     -[integer] the number of rows or columns in A
!!ia    -[integer (nelt)] the row numbers
!!ja    -[integer (nelt)] the column numbers
!!A     -[double (nelt)] the nonzero matrix entries
!!b     -[double (n)] the rhs vector
!!
!!~Outputs~
!!
!!b     -[double (n)] the solution vector replaces the rhs vector on output
!!
!!
!!  Createded by John Burkhardt:  24 July 2007
!!
!!  Modified by Cameron Bracken:  08 March 2008
!!
  implicit none

  integer :: maxl
  integer:: n,nelt
  integer::lrgw,ligw

  double precision::err,tol
  integer,dimension(nelt)::ia,ja
  double precision,dimension(nelt)::a
  double precision,dimension(n):: b,x
  
  integer::ierr
  integer,dimension(:),allocatable::igwk
  integer::isym,iter,itmax,itol,iunit,i
  integer,dimension(1)::iwork
  external matvec_triad
  external msolve_identity
  double precision,dimension(:),allocatable:: rgwk
  double precision,dimension(1):: rwork
  double precision,dimension(n):: sb,sx
  logical::quiet
  
  maxl=100!int(n/4)
  
  ligw=131 + 16*N
  allocate(igwk(ligw))

  lrgw=1+n*(maxl+6)+maxl*(maxl+3)
  !write(*,*)'Allocated',lrgw,'for RGWK'
  allocate(rgwk(lrgw))

  isym = 0
  itol = 0
  tol = 0d0
  itmax = maxl*(10+1)
  iunit = 0
  sb(1:n) = 1d0
  sx(1:n) = 1d0

  igwk(1) = maxl
  igwk(2) = maxl
  igwk(3) = 0
  igwk(4) = 0
  igwk(5) = 10


    !Zero out X, because it will be used as an initial guess.
  x=0d0

  call dgmres ( n, b, x, nelt, ia, ja, a, isym, matvec_triad, &
    msolve_identity, itol, tol, itmax, iter, err, ierr, iunit, sb, &
    sx, rgwk, lrgw, igwk, ligw, rwork, iwork )

    !overwrite rhs vector with solution vector
  b=x

    if(.not.quiet)then
      write(*,*)
      write(*,*)'Number of iterations      ', iter
      write(*,*)'Minimum length for RGWK   ', igwk(6)
      write(*,*)'Error estimate            ', err
      write(*,*)'Error code is', ierr,'should be 0'
    end if

return
end

subroutine matvec_triad ( n, x, y, nelt, ia, ja, a, isym )

!*****************************************************************************80
!
!! MATVEC_TRIAD computes A*X for a matrix A stored in SLAP Triad form.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vectors.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) Y(N), the product A * X.
!
!    Input, integer NELT, the number of nonzero entries in A.
!
!    Input, integer IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), the data
!    structure storing the sparse matrix.
!
!    Input, integer ISYM, is 0 if all nonzero entries of the matrix
!    are stored, and 1 if only the diagonal and upper or lower triangle
!    are stored.
!
  implicit none

  integer::n,nelt,isym,k,i

  double precision,dimension(nelt):: a
  integer,dimension(nelt)::ia,ja
  double precision,dimension(n):: x,y
  y=0d0


  do k = 1, nelt
    y(ia(k)) = y(ia(k)) + a(k) * x(ja(k))
  end do

  return
end

subroutine msolve_identity ( n, r, z, nelt, ia, ja, a, isym, rwork, iwork )

!*****************************************************************************80
!
!! MSOLVE_IDENTITY applies the identity matrix preconditioner.
!
!  Discussion:
!
!    Most SLAP solver routines require a preconditioner routine
!    that can solve M * Z = R.  If no preconditioning is required,
!    then you can simply behave as though the preconditioning matrix
!    M was the identity matrix.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vectors.
!
!    Input, real ( kind = 8 ) R(N), the right hand side.
!
!    Output, real ( kind = 8 ) Z(N), the solution of M * Z = R.
!
!    Input, integer NELT, the number of nonzero entries in A.
!
!    Input, integer IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), the data
!    structure storing the sparse matrix.
!
!    Input, integer ISYM, is 0 if all nonzero entries of the matrix
!    are stored, and 1 if only the diagonal and upper or lower triangle
!    are stored.
!
!    Input, real ( kind = 8 ) RWORK(*), a real array that
!    can be used to pass information to the preconditioner.
!
!    Input, integer IWORK(*), an integer array that
!    can be used to pass information to the preconditioner.
!
  implicit none

  integer::n,nelt,isym

  double precision,dimension(nelt)::a
  integer,dimension(nelt)::ia,ja
  integer::iwork(*)
  double precision::rwork(*)
  double precision,dimension(n)::r,z
  write(*,*)'msolve'

  z(1:n) = r(1:n)

  return
end

