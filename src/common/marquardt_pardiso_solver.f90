module marquardt_pardiso_solver
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use sparse_triplet_builder, only: dp, csr_transpose, csr_matmul_ozaki_diag_upper, &
    oz_add_product_slice6_module, oz_add_slice6_module
  implicit none
  private
  public :: solve_marquardt_pardiso, solve_marquardt_pardiso_augmented, &
    solve_marquardt_pardiso_augmented_cached, release_marquardt_pardiso_cache

  type :: augmented_pardiso_cache_type
    integer :: n=0, naug=0, nnz_j=0, nnz_augmented=0
    integer(int64) :: pt(64)=0_int64
    integer :: iparm(64)=0
    integer, allocatable :: k_ir(:), k_ic(:), perm(:)
    real(dp), allocatable :: k_a(:)
    logical :: analyzed=.false., factored=.false.
    real(dp) :: factored_lambda=huge(1.0_dp)
  end type

  type(augmented_pardiso_cache_type), save :: augmented_cache
contains
  subroutine solve_marquardt_pardiso(n, ir, ic, a, rhs, d, lambda, tol, x, relres, status, nnz_normal)
    integer, intent(in) :: n, ir(:), ic(:)
    real(dp), intent(in) :: a(:), rhs(:), d(:), lambda, tol
    real(dp), intent(out) :: x(:), relres
    integer, intent(out) :: status, nnz_normal
    integer, allocatable :: jt_ir(:), jt_ic(:), h_ir(:), h_ic(:), perm(:)
    real(dp), allocatable :: jt_a(:), h_a(:), b(:), jx(:), hx(:), residual(:), correction(:)
    integer(int64), parameter :: max_normal_nnz=600000000_int64
    integer(int64), parameter :: default_max_memory_kb=891289600_int64
    integer(int64) :: pt(64), symbolic_peak_kb, numerical_peak_kb, max_memory_kb
    integer, parameter :: max_refine=8
    integer :: iparm(64), i, s, nnzt, nnz, error, release_error, phase, refine_iter
    integer :: env_status, env_length, read_status, build_status
    real(dp) :: rhsnorm, slices(6), max_memory_gib
    character(len=64) :: env_value
    logical :: symbolic_only

    x=0.0_dp
    relres=huge(1.0_dp)
    status=-1001
    nnz_normal=0
    if (n<1.or.size(ir)/=n+1.or.size(rhs)/=n.or.size(d)/=n.or.size(x)/=n) return
    if (.not.ieee_is_finite(lambda).or..not.ieee_is_finite(tol)) return
    if (lambda<0.0_dp.or.tol<=0.0_dp.or.ir(1)/=1) return
    nnz=ir(n+1)-1
    if (nnz<0.or.nnz>size(a).or.nnz>size(ic)) return
    do i=1,n
      if (ir(i+1)<ir(i).or.ir(i)<1.or.ir(i+1)>nnz+1) return
    end do
    if (any(ic(1:nnz)<1).or.any(ic(1:nnz)>n)) return
    status=-1002
    if (.not.all(ieee_is_finite(a(1:nnz)))) return
    if (.not.all(ieee_is_finite(rhs)).or..not.all(ieee_is_finite(d))) return
    if (any(d<=0.0_dp)) return

    call csr_transpose(n,n,ir,ic,a,jt_ir,jt_ic,jt_a,nnzt)

!   Store only the exact upper triangle of J^T*J+lambda*D.  The product
!   still uses the six-slice compensated accumulation used elsewhere.
    call csr_matmul_ozaki_diag_upper(n,jt_ir,jt_ic,jt_a,ir,ic,a,d,lambda, &
      max_normal_nnz,h_ir,h_ic,h_a,nnz_normal,build_status)
    if (build_status/=0) then
      status=build_status
      return
    end if
    if (.not.all(ieee_is_finite(h_a))) return
    status=-1004
    if (any(h_ir(2:n+1)==h_ir(1:n))) return
    allocate(perm(n),b(n),jx(n),hx(n),residual(n),correction(n))
    perm=0
    b=rhs
    pt=0_int64
    iparm=0
    iparm(1)=1
    iparm(2)=2
    iparm(8)=10
    iparm(10)=0
    iparm(11)=0
    iparm(13)=0
    iparm(27)=1
    iparm(35)=0
!   Keep the factorization in core.  This HPC run reserves enough RAM for
!   the phase-11 estimate and avoids the 2 GiB default OOC core limit.
    iparm(60)=0

    max_memory_kb=default_max_memory_kb
    env_value=''
    call get_environment_variable('OCFEM_PARDISO_MAX_MEMORY_GIB',env_value, &
      length=env_length,status=env_status)
    if (env_status==0.and.env_length>0) then
      read(env_value(1:min(env_length,len(env_value))),*,iostat=read_status) max_memory_gib
      if (read_status==0.and.max_memory_gib>0.0_dp) &
        max_memory_kb=int(max_memory_gib*1024.0_dp*1024.0_dp,int64)
    end if
    symbolic_only=.false.
    env_value=''
    call get_environment_variable('OCFEM_PARDISO_SYMBOLIC_ONLY',env_value, &
      length=env_length,status=env_status)
    if (env_status==0.and.env_length>0) then
      symbolic_only=trim(adjustl(env_value(1:min(env_length,len(env_value)))))=='1'
    end if

!   Analyze first.  PARDISO mtype=2 consumes one triangle of a real
!   symmetric positive-definite matrix and reports memory in iparm(15:17).
    phase=11
    call pardiso(pt,1,1,2,phase,n,h_a,h_ir,h_ic,perm,1,iparm,0,b,x,error)
    if (error/=0) then
      phase=-1
      call pardiso(pt,1,1,2,phase,n,h_a,h_ir,h_ic,perm,1,iparm,0,b,x,release_error)
      status=error
      return
    end if
    symbolic_peak_kb=int(iparm(15),int64)
    numerical_peak_kb=int(iparm(16),int64)+int(iparm(17),int64)
    write(*,*) 'PARDISO_PHASE11 N/NNZ_UPPER=',n,nnz_normal
    write(*,*) 'PARDISO_PHASE11 SYMBOLIC_KB/PERMANENT_KB/FACTORS_KB/PEAK_KB/LIMIT_KB=', &
      iparm(15),iparm(16),iparm(17),max(symbolic_peak_kb,numerical_peak_kb),max_memory_kb

    if (symbolic_only) then
      phase=-1
      call pardiso(pt,1,1,2,phase,n,h_a,h_ir,h_ic,perm,1,iparm,0,b,x,release_error)
      status=-2003
      return
    end if
    if (max(symbolic_peak_kb,numerical_peak_kb)>max_memory_kb) then
      phase=-1
      call pardiso(pt,1,1,2,phase,n,h_a,h_ir,h_ic,perm,1,iparm,0,b,x,release_error)
      status=-2002
      return
    end if

    phase=22
    call pardiso(pt,1,1,2,phase,n,h_a,h_ir,h_ic,perm,1,iparm,0,b,x,error)
    if (error==0) then
      phase=33
      call pardiso(pt,1,1,2,phase,n,h_a,h_ir,h_ic,perm,1,iparm,0,b,x,error)
    end if

!   Verify and refine with the true matrix action.  PARDISO's own residual
!   is based on the stored normal matrix; this check also catches errors in
!   J^T*J assembly and keeps the Newton update tied to the requested tol.
    if (error==0) then
      rhsnorm=stable_norm(rhs)
      do refine_iter=0,max_refine
        if (.not.all(ieee_is_finite(x))) then
          error=-1002
          exit
        end if
        call checked_matvec(n,ir,ic,a,x,jx,error)
        if (error/=0) exit
        call checked_matvec(n,jt_ir,jt_ic,jt_a,jx,hx,error)
        if (error/=0) exit
        do i=1,n
          if (.not.ieee_is_finite(lambda*d(i))) then
            error=-1002
            exit
          end if
          if (.not.ieee_is_finite((lambda*d(i))*x(i))) then
            error=-1002
            exit
          end if
          slices=0.0_dp
          call oz_add_slice6_module(slices,rhs(i))
          call oz_add_slice6_module(slices,-hx(i))
          call oz_add_product_slice6_module(slices,-lambda*d(i),x(i))
          residual(i)=0.0_dp
          do s=6,1,-1
            residual(i)=residual(i)+slices(s)
          end do
        end do
        if (error/=0) exit
        if (.not.all(ieee_is_finite(residual))) then
          error=-1002
          exit
        end if
        if (rhsnorm==0.0_dp) then
          relres=stable_norm(residual)
        else
          relres=stable_norm(residual)/rhsnorm
        end if
        if (.not.ieee_is_finite(relres)) then
          error=-1002
          exit
        end if
        write(*,*) 'PARDISO_TRUE_RESIDUAL_REFINEMENT ITER/RELRES/TOL=',refine_iter,relres,tol
        if (relres<=tol.or.refine_iter==max_refine) exit
        b=residual
        correction=0.0_dp
        phase=33
        call pardiso(pt,1,1,2,phase,n,h_a,h_ir,h_ic,perm,1,iparm,0,b,correction,error)
        if (error/=0) exit
        x=x+correction
      end do
    end if

    phase=-1
    call pardiso(pt,1,1,2,phase,n,h_a,h_ir,h_ic,perm,1,iparm,0,b,x,release_error)
    status=error
    if (status==0.and.release_error/=0) status=release_error
    if (status/=0) return
    status=0
    if (relres>tol) status=1
  end subroutine

  subroutine solve_marquardt_pardiso_augmented_cached(n,ir,ic,a,f,d,lambda,tol, &
      x,relres,status,nnz_augmented,rebuild_factor)
    integer, intent(in) :: n,ir(:),ic(:)
    real(dp), intent(in) :: a(:),f(:),d(:),lambda,tol
    real(dp), intent(out) :: x(:),relres
    integer, intent(out) :: status,nnz_augmented
    logical, intent(in) :: rebuild_factor
    integer :: nnz,naug,i,k,pos,error,phase,release_error
    integer :: env_status,env_length,read_status
    integer(int64), parameter :: default_max_memory_kb=891289600_int64
    integer(int64) :: symbolic_peak_kb,numerical_peak_kb,max_memory_kb
    real(dp), allocatable :: b(:),solution(:)
    real(dp) :: max_memory_gib
    character(len=64) :: env_value
    logical :: pattern_matches,symbolic_only

    x=0.0_dp
    relres=huge(1.0_dp)
    status=-1001
    nnz_augmented=0
    if (n<1.or.size(ir)/=n+1.or.size(f)/=n.or.size(d)/=n.or.size(x)/=n) return
    if (.not.ieee_is_finite(lambda).or..not.ieee_is_finite(tol)) return
    if (lambda<0.0_dp.or.tol<=0.0_dp.or.ir(1)/=1) return
    nnz=ir(n+1)-1
    if (nnz<0.or.nnz>size(a).or.nnz>size(ic)) return
    if (.not.all(ieee_is_finite(a(1:nnz)))) return
    if (.not.all(ieee_is_finite(f)).or..not.all(ieee_is_finite(d))) return
    if (any(d<=0.0_dp)) return

    naug=2*n
    nnz_augmented=nnz+2*n
    pattern_matches=augmented_cache%analyzed.and.augmented_cache%n==n.and. &
      augmented_cache%nnz_j==nnz.and.augmented_cache%nnz_augmented==nnz_augmented
    if (pattern_matches) then
      pos=1
      do i=1,n
        if (augmented_cache%k_ir(i)/=pos.or.augmented_cache%k_ic(pos)/=i) then
          pattern_matches=.false.
          exit
        end if
        pos=pos+1
        do k=ir(i),ir(i+1)-1
          if (augmented_cache%k_ic(pos)/=n+ic(k)) then
            pattern_matches=.false.
            exit
          end if
          pos=pos+1
        end do
        if (.not.pattern_matches) exit
      end do
    end if
    if (.not.pattern_matches) then
      call release_marquardt_pardiso_cache()
      augmented_cache%n=n
      augmented_cache%naug=naug
      augmented_cache%nnz_j=nnz
      augmented_cache%nnz_augmented=nnz_augmented
      allocate(augmented_cache%k_ir(naug+1),augmented_cache%k_ic(nnz_augmented), &
        augmented_cache%k_a(nnz_augmented),augmented_cache%perm(naug))
      pos=1
      augmented_cache%k_ir(1)=1
      do i=1,n
        augmented_cache%k_ic(pos)=i
        pos=pos+1
        do k=ir(i),ir(i+1)-1
          augmented_cache%k_ic(pos)=n+ic(k)
          pos=pos+1
        end do
        augmented_cache%k_ir(i+1)=pos
      end do
      do i=1,n
        augmented_cache%k_ic(pos)=n+i
        pos=pos+1
        augmented_cache%k_ir(n+i+1)=pos
      end do
      if (pos-1/=nnz_augmented) then
        status=-1004
        call release_marquardt_pardiso_cache()
        return
      end if
      augmented_cache%perm=0
      augmented_cache%pt=0_int64
      augmented_cache%iparm=0
      augmented_cache%iparm(1)=1
      augmented_cache%iparm(2)=2
      augmented_cache%iparm(8)=10
      augmented_cache%iparm(11)=1
      augmented_cache%iparm(13)=1
      augmented_cache%iparm(27)=1
      augmented_cache%iparm(35)=0
      augmented_cache%iparm(60)=0
    end if

    allocate(b(naug),solution(naug))
    b=0.0_dp
    b(1:n)=-f
    solution=0.0_dp

    if (.not.augmented_cache%analyzed) then
      call update_augmented_values(n,ir,a,d,lambda)
      phase=11
      call pardiso(augmented_cache%pt,1,1,-2,phase,naug,augmented_cache%k_a, &
        augmented_cache%k_ir,augmented_cache%k_ic,augmented_cache%perm,1, &
        augmented_cache%iparm,0,b,solution,error)
      if (error/=0) then
        status=error
        call release_marquardt_pardiso_cache()
        return
      end if
      augmented_cache%analyzed=.true.
      symbolic_peak_kb=int(augmented_cache%iparm(15),int64)
      numerical_peak_kb=int(augmented_cache%iparm(16),int64)+ &
        int(augmented_cache%iparm(17),int64)
      max_memory_kb=default_max_memory_kb
      env_value=''
      call get_environment_variable('OCFEM_PARDISO_MAX_MEMORY_GIB',env_value, &
        length=env_length,status=env_status)
      if (env_status==0.and.env_length>0) then
        read(env_value(1:min(env_length,len(env_value))),*,iostat=read_status) max_memory_gib
        if (read_status==0.and.max_memory_gib>0.0_dp) &
          max_memory_kb=int(max_memory_gib*1024.0_dp*1024.0_dp,int64)
      end if
      write(*,*) 'PARDISO_AUG_CACHE_PHASE11 N2/NNZ_UPPER=',naug,nnz_augmented
      write(*,*) 'PARDISO_AUG_CACHE_MEMORY PEAK_KB/LIMIT_KB=', &
        max(symbolic_peak_kb,numerical_peak_kb),max_memory_kb
      if (max(symbolic_peak_kb,numerical_peak_kb)>max_memory_kb) then
        status=-2002
        call release_marquardt_pardiso_cache()
        return
      end if
      symbolic_only=.false.
      env_value=''
      call get_environment_variable('OCFEM_PARDISO_SYMBOLIC_ONLY',env_value, &
        length=env_length,status=env_status)
      if (env_status==0.and.env_length>0) &
        symbolic_only=trim(adjustl(env_value(1:min(env_length,len(env_value)))))=='1'
      if (symbolic_only) then
        status=-2003
        return
      end if
    end if

    if (rebuild_factor.or..not.augmented_cache%factored) then
      call update_augmented_values(n,ir,a,d,lambda)
      phase=22
      call pardiso(augmented_cache%pt,1,1,-2,phase,naug,augmented_cache%k_a, &
        augmented_cache%k_ir,augmented_cache%k_ic,augmented_cache%perm,1, &
        augmented_cache%iparm,0,b,solution,error)
      if (error/=0) then
        augmented_cache%factored=.false.
        status=error
        return
      end if
      augmented_cache%factored=.true.
      augmented_cache%factored_lambda=lambda
      write(*,*) 'PARDISO_AUG_CACHE_FACTOR_REBUILT LAMBDA=',lambda
    else
      write(*,*) 'PARDISO_AUG_CACHE_FACTOR_REUSED LAMBDA_FACTORED/CURRENT=', &
        augmented_cache%factored_lambda,lambda
    end if

    phase=33
    call pardiso(augmented_cache%pt,1,1,-2,phase,naug,augmented_cache%k_a, &
      augmented_cache%k_ir,augmented_cache%k_ic,augmented_cache%perm,1, &
      augmented_cache%iparm,0,b,solution,error)
    if (error/=0) then
      status=error
      return
    end if
    x=solution(n+1:naug)
    if (.not.all(ieee_is_finite(x))) then
      status=-1002
      return
    end if
    call normal_equation_relres(n,ir,ic,a,f,d,lambda,x,relres,error)
    if (error/=0) then
      status=error
      return
    end if
    status=0
    if (relres>tol) status=1
  end subroutine

  subroutine update_augmented_values(n,ir,a,d,lambda)
    integer, intent(in) :: n,ir(:)
    real(dp), intent(in) :: a(:),d(:),lambda
    integer :: i,k,pos
    pos=1
    do i=1,n
      augmented_cache%k_a(pos)=1.0_dp
      pos=pos+1
      do k=ir(i),ir(i+1)-1
        augmented_cache%k_a(pos)=a(k)
        pos=pos+1
      end do
    end do
    do i=1,n
      augmented_cache%k_a(pos)=-lambda*d(i)
      pos=pos+1
    end do
  end subroutine

  subroutine normal_equation_relres(n,ir,ic,a,f,d,lambda,x,relres,error)
    integer, intent(in) :: n,ir(:),ic(:)
    real(dp), intent(in) :: a(:),f(:),d(:),lambda,x(:)
    real(dp), intent(out) :: relres
    integer, intent(out) :: error
    integer, allocatable :: jt_ir(:),jt_ic(:)
    real(dp), allocatable :: jt_a(:),rhs(:),jx(:),hx(:),residual(:)
    integer :: nnzt,i,s
    real(dp) :: rhsnorm,slices(6)
    allocate(rhs(n),jx(n),hx(n),residual(n))
    call csr_transpose(n,n,ir,ic,a,jt_ir,jt_ic,jt_a,nnzt)
    call checked_matvec(n,jt_ir,jt_ic,jt_a,f,rhs,error)
    if (error/=0) return
    rhs=-rhs
    call checked_matvec(n,ir,ic,a,x,jx,error)
    if (error/=0) return
    call checked_matvec(n,jt_ir,jt_ic,jt_a,jx,hx,error)
    if (error/=0) return
    do i=1,n
      slices=0.0_dp
      call oz_add_slice6_module(slices,rhs(i))
      call oz_add_slice6_module(slices,-hx(i))
      call oz_add_product_slice6_module(slices,-lambda*d(i),x(i))
      residual(i)=0.0_dp
      do s=6,1,-1
        residual(i)=residual(i)+slices(s)
      end do
    end do
    if (.not.all(ieee_is_finite(residual))) then
      error=-1002
      return
    end if
    rhsnorm=stable_norm(rhs)
    if (rhsnorm==0.0_dp) then
      relres=stable_norm(residual)
    else
      relres=stable_norm(residual)/rhsnorm
    end if
    if (.not.ieee_is_finite(relres)) then
      error=-1002
      return
    end if
    error=0
  end subroutine

  subroutine release_marquardt_pardiso_cache()
    real(dp), allocatable :: b(:),x(:)
    integer :: phase,error
    if (augmented_cache%analyzed.and.allocated(augmented_cache%k_a)) then
      allocate(b(augmented_cache%naug),x(augmented_cache%naug))
      b=0.0_dp
      x=0.0_dp
      phase=-1
      call pardiso(augmented_cache%pt,1,1,-2,phase,augmented_cache%naug, &
        augmented_cache%k_a,augmented_cache%k_ir,augmented_cache%k_ic, &
        augmented_cache%perm,1,augmented_cache%iparm,0,b,x,error)
    end if
    if (allocated(augmented_cache%k_ir)) deallocate(augmented_cache%k_ir)
    if (allocated(augmented_cache%k_ic)) deallocate(augmented_cache%k_ic)
    if (allocated(augmented_cache%k_a)) deallocate(augmented_cache%k_a)
    if (allocated(augmented_cache%perm)) deallocate(augmented_cache%perm)
    augmented_cache%n=0
    augmented_cache%naug=0
    augmented_cache%nnz_j=0
    augmented_cache%nnz_augmented=0
    augmented_cache%pt=0_int64
    augmented_cache%iparm=0
    augmented_cache%analyzed=.false.
    augmented_cache%factored=.false.
    augmented_cache%factored_lambda=huge(1.0_dp)
  end subroutine

  subroutine solve_marquardt_pardiso_augmented(n, ir, ic, a, f, d, lambda, tol, x, relres, status, nnz_augmented)
    integer, intent(in) :: n, ir(:), ic(:)
    real(dp), intent(in) :: a(:), f(:), d(:), lambda, tol
    real(dp), intent(out) :: x(:), relres
    integer, intent(out) :: status, nnz_augmented
    integer, allocatable :: jt_ir(:), jt_ic(:), k_ir(:), k_ic(:), perm(:)
    real(dp), allocatable :: jt_a(:), k_a(:), b(:), solution(:), jx(:), hx(:), rhs(:), residual(:)
    integer(int64), parameter :: default_max_memory_kb=891289600_int64
    integer(int64) :: pt(64), symbolic_peak_kb, numerical_peak_kb, max_memory_kb
    integer :: iparm(64), i, k, s, pos, nnz, nnzt, naug, error, release_error, phase
    integer :: env_status, env_length, read_status
    real(dp) :: rhsnorm, slices(6), max_memory_gib
    character(len=64) :: env_value
    logical :: symbolic_only

    x=0.0_dp
    relres=huge(1.0_dp)
    status=-1001
    nnz_augmented=0
    if (n<1.or.size(ir)/=n+1.or.size(f)/=n.or.size(d)/=n.or.size(x)/=n) return
    if (.not.ieee_is_finite(lambda).or..not.ieee_is_finite(tol)) return
    if (lambda<0.0_dp.or.tol<=0.0_dp.or.ir(1)/=1) return
    nnz=ir(n+1)-1
    if (nnz<0.or.nnz>size(a).or.nnz>size(ic)) return
    do i=1,n
      if (ir(i+1)<ir(i).or.ir(i)<1.or.ir(i+1)>nnz+1) return
    end do
    if (any(ic(1:nnz)<1).or.any(ic(1:nnz)>n)) return
    status=-1002
    if (.not.all(ieee_is_finite(a(1:nnz)))) return
    if (.not.all(ieee_is_finite(f)).or..not.all(ieee_is_finite(d))) return
    if (any(d<=0.0_dp)) return

!   The symmetric indefinite augmented system
!
!       [ I    J_s       ] [q] = [-f_s]
!       [ J_s^T -lambda D] [z]   [  0 ]
!
!   has the Schur complement
!       (J_s^T J_s + lambda D) z = -J_s^T f_s.
!   Store only its upper triangle for PARDISO MTYPE=-2.  Unlike the
!   normal equation, this needs only NNZ(J_s)+2*N stored entries.
    naug=2*n
    nnz_augmented=nnz+2*n
    allocate(k_ir(naug+1),k_ic(nnz_augmented),k_a(nnz_augmented))
    pos=1
    k_ir(1)=1
    do i=1,n
      k_ic(pos)=i
      k_a(pos)=1.0_dp
      pos=pos+1
      do k=ir(i),ir(i+1)-1
        k_ic(pos)=n+ic(k)
        k_a(pos)=a(k)
        pos=pos+1
      end do
      k_ir(i+1)=pos
    end do
    do i=1,n
      k_ic(pos)=n+i
      k_a(pos)=-lambda*d(i)
      pos=pos+1
      k_ir(n+i+1)=pos
    end do
    if (pos-1/=nnz_augmented) then
      status=-1004
      return
    end if
    if (.not.all(ieee_is_finite(k_a))) return

    allocate(perm(naug),b(naug),solution(naug),jx(n),hx(n),rhs(n),residual(n))
    perm=0
    b=0.0_dp
    b(1:n)=-f
    solution=0.0_dp
    pt=0_int64
    iparm=0
    iparm(1)=1
    iparm(2)=2
    iparm(8)=10
!   Symmetric-indefinite pivoting benefits from nonsymmetric scaling and
!   maximum-weight matching; both are explicitly enabled with user iparm.
    iparm(11)=1
    iparm(13)=1
    iparm(27)=1
    iparm(35)=0
!   Use the same in-core policy for the augmented-system solver.
    iparm(60)=0

    max_memory_kb=default_max_memory_kb
    env_value=''
    call get_environment_variable('OCFEM_PARDISO_MAX_MEMORY_GIB',env_value, &
      length=env_length,status=env_status)
    if (env_status==0.and.env_length>0) then
      read(env_value(1:min(env_length,len(env_value))),*,iostat=read_status) max_memory_gib
      if (read_status==0.and.max_memory_gib>0.0_dp) &
        max_memory_kb=int(max_memory_gib*1024.0_dp*1024.0_dp,int64)
    end if
    symbolic_only=.false.
    env_value=''
    call get_environment_variable('OCFEM_PARDISO_SYMBOLIC_ONLY',env_value, &
      length=env_length,status=env_status)
    if (env_status==0.and.env_length>0) then
      symbolic_only=trim(adjustl(env_value(1:min(env_length,len(env_value)))))=='1'
    end if

    phase=11
    call pardiso(pt,1,1,-2,phase,naug,k_a,k_ir,k_ic,perm,1,iparm,0,b,solution,error)
    if (error/=0) then
      phase=-1
      call pardiso(pt,1,1,-2,phase,naug,k_a,k_ir,k_ic,perm,1,iparm,0,b,solution,release_error)
      status=error
      return
    end if
    symbolic_peak_kb=int(iparm(15),int64)
    numerical_peak_kb=int(iparm(16),int64)+int(iparm(17),int64)
    write(*,*) 'PARDISO_AUG_PHASE11 N2/NNZ_UPPER=',naug,nnz_augmented
    write(*,*) 'PARDISO_AUG_PHASE11 SYMBOLIC_KB/PERMANENT_KB/FACTORS_KB/PEAK_KB/LIMIT_KB=', &
      iparm(15),iparm(16),iparm(17),max(symbolic_peak_kb,numerical_peak_kb),max_memory_kb

    if (symbolic_only) then
      phase=-1
      call pardiso(pt,1,1,-2,phase,naug,k_a,k_ir,k_ic,perm,1,iparm,0,b,solution,release_error)
      status=-2003
      return
    end if
    if (max(symbolic_peak_kb,numerical_peak_kb)>max_memory_kb) then
      phase=-1
      call pardiso(pt,1,1,-2,phase,naug,k_a,k_ir,k_ic,perm,1,iparm,0,b,solution,release_error)
      status=-2002
      return
    end if

    phase=22
    call pardiso(pt,1,1,-2,phase,naug,k_a,k_ir,k_ic,perm,1,iparm,0,b,solution,error)
    if (error==0) then
      phase=33
      call pardiso(pt,1,1,-2,phase,naug,k_a,k_ir,k_ic,perm,1,iparm,0,b,solution,error)
    end if
    phase=-1
    call pardiso(pt,1,1,-2,phase,naug,k_a,k_ir,k_ic,perm,1,iparm,0,b,solution,release_error)
    status=error
    if (status==0.and.release_error/=0) status=release_error
    if (status/=0) return
    status=-1002
    x=solution(n+1:naug)
    if (.not.all(ieee_is_finite(x))) return

!   Report the true normal-equation residual so RELRES has exactly the
!   same meaning as in the MTYPE=2 implementation.
    call csr_transpose(n,n,ir,ic,a,jt_ir,jt_ic,jt_a,nnzt)
    call checked_matvec(n,jt_ir,jt_ic,jt_a,f,rhs,error)
    if (error/=0) return
    rhs=-rhs
    call checked_matvec(n,ir,ic,a,x,jx,error)
    if (error/=0) return
    call checked_matvec(n,jt_ir,jt_ic,jt_a,jx,hx,error)
    if (error/=0) return
    do i=1,n
      slices=0.0_dp
      call oz_add_slice6_module(slices,rhs(i))
      call oz_add_slice6_module(slices,-hx(i))
      call oz_add_product_slice6_module(slices,-lambda*d(i),x(i))
      residual(i)=0.0_dp
      do s=6,1,-1
        residual(i)=residual(i)+slices(s)
      end do
    end do
    if (.not.all(ieee_is_finite(residual))) return
    rhsnorm=stable_norm(rhs)
    if (rhsnorm==0.0_dp) then
      relres=stable_norm(residual)
    else
      relres=stable_norm(residual)/rhsnorm
    end if
    if (.not.ieee_is_finite(relres)) return
    status=0
    if (relres>tol) status=1
  end subroutine

  subroutine checked_matvec(n,ir,ic,a,x,y,error)
    integer, intent(in) :: n,ir(:),ic(:)
    real(dp), intent(in) :: a(:),x(:)
    real(dp), intent(out) :: y(:)
    integer, intent(out) :: error
    integer :: i,k,s
    real(dp) :: slices(6)
    error=-1002
    do i=1,n
      slices=0.0_dp
      do k=ir(i),ir(i+1)-1
        if (.not.ieee_is_finite(a(k)*x(ic(k)))) return
        call oz_add_product_slice6_module(slices,a(k),x(ic(k)))
        if (.not.all(ieee_is_finite(slices))) return
      end do
      y(i)=0.0_dp
      do s=6,1,-1
        y(i)=y(i)+slices(s)
      end do
      if (.not.ieee_is_finite(y(i))) return
    end do
    error=0
  end subroutine

  real(dp) function stable_norm(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: scale
    scale=maxval(abs(x))
    value=0.0_dp
    if (scale>0.0_dp) value=scale*sqrt(sum((x/scale)**2))
  end function
end module marquardt_pardiso_solver
