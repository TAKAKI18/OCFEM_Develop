program test_marquardt_pardiso
  use, intrinsic :: ieee_arithmetic
  use sparse_triplet_builder, only: dp
  use marquardt_pardiso_solver, only: solve_marquardt_pardiso
  implicit none
  integer :: ir(5),ic(9),k,i,status,nnz_h,info,ipiv(4)
  real(dp) :: a(9),f(4),rhs(4),d(4),jmat(4,4),h(4,4),expected(4,1),x(4),relres,lambda
  ir=[1,3,6,8,10]
  ic=[1,2,1,2,3,3,4,1,4]
  a=[2.0_dp,-1.0_dp,0.5_dp,3.0_dp,1.0_dp,4.0_dp,0.25_dp,-0.5_dp,2.0_dp]
  jmat=0.0_dp
  do i=1,4
    do k=ir(i),ir(i+1)-1
      jmat(i,ic(k))=a(k)
    end do
  end do
  f=[0.25_dp,-0.5_dp,1.0_dp,0.125_dp]
  rhs=-matmul(transpose(jmat),f)
  d=1.0_dp+sum(jmat*jmat,dim=1)
  do k=0,2
    lambda=0.125_dp*k
    h=matmul(transpose(jmat),jmat)
    do i=1,4
      h(i,i)=h(i,i)+lambda*d(i)
    end do
    expected(:,1)=rhs
    call dgesv(4,1,h,4,ipiv,expected,4,info)
    if(info/=0) error stop 'DENSE_REFERENCE_FAILED'
    call solve_marquardt_pardiso(4,ir,ic,a,rhs,d,lambda,1.0d-11,x,relres,status,nnz_h)
    if(status/=0.or.relres>1.0d-11) error stop 'PARDISO_PRECISION_FAILED'
    if(maxval(abs(x-expected(:,1)))>1.0d-11) error stop 'PARDISO_REFERENCE_MISMATCH'
    if(nnz_h<=4) error stop 'NORMAL_MATRIX_WAS_DIAGONAL_ONLY'
  end do
  print *, 'PASS: full CSR normal matrix, lambda=0/positive, dense LU reference and repeated cleanup'
  rhs=0.0_dp
  call solve_marquardt_pardiso(4,ir,ic,a,rhs,d,1.0_dp,1.0d-11,x,relres,status,nnz_h)
  if(status/=0.or.any(x/=0.0_dp)) error stop 'ZERO_RHS_FAILED'
  print *, 'PASS: zero right hand side'
  rhs=1.0_dp
  a(1)=ieee_value(0.0_dp,ieee_positive_inf)
  call solve_marquardt_pardiso(4,ir,ic,a,rhs,d,1.0_dp,1.0d-11,x,relres,status,nnz_h)
  if(status/=-1002) error stop 'NONFINITE_MATRIX_NOT_REJECTED'
  a(1)=2.0_dp
  call solve_marquardt_pardiso(4,ir,ic,a,rhs,d,-1.0_dp,1.0d-11,x,relres,status,nnz_h)
  if(status/=-1001) error stop 'NEGATIVE_DAMPING_NOT_REJECTED'
  a=0.0_dp
  call solve_marquardt_pardiso(4,ir,ic,a,rhs,d,0.0_dp,1.0d-11,x,relres,status,nnz_h)
  if(status/=-1004) error stop 'EMPTY_NORMAL_ROW_NOT_REJECTED'
  print *, 'PASS: nonfinite matrix, negative damping and empty normal rows rejected'
  call parallel_case()
  print *, 'MARQUARDT_PARDISO_TESTS_PASSED'
contains
  subroutine parallel_case()
    integer,parameter :: n=1200
    integer :: rows(n+1),cols(2*n),i,k,next,exit_status,nnz_h
    real(dp) :: vals(2*n),exact(n),b(n),ax(n),diag(n),d(n),answer(n),relative
    do i=1,n
      rows(i)=2*i-1
      cols(2*i-1)=i
      next=mod(i,n)+1
      cols(2*i)=next
      vals(2*i-1)=2.0_dp+0.01_dp*mod(i,7)
      vals(2*i)=-0.25_dp
      exact(i)=sin(0.13_dp*i)
    end do
    rows(n+1)=2*n+1
    ax=0.0_dp
    b=0.0_dp
    diag=0.0_dp
    do i=1,n
      do k=rows(i),rows(i+1)-1
        ax(i)=ax(i)+vals(k)*exact(cols(k))
        diag(cols(k))=diag(cols(k))+vals(k)**2
      end do
    end do
    do i=1,n
      do k=rows(i),rows(i+1)-1
        b(cols(k))=b(cols(k))+vals(k)*ax(i)
      end do
    end do
    d=1.0_dp+diag
    b=b+0.125_dp*d*exact
    call solve_marquardt_pardiso(n,rows,cols,vals,b,d,0.125_dp,1.0d-11,answer,relative,exit_status,nnz_h)
    if(exit_status/=0.or.relative>1.0d-11.or.maxval(abs(answer-exact))>1.0d-10) &
      error stop 'PARALLEL_NORMAL_PARDISO_FAILED'
    print *, 'PASS: 1200-row off-diagonal normal matrix / 8 MKL threads',nnz_h,relative
  end subroutine
end program
