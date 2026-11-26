module sparse_triplet_builder
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)

  type :: triplet_matrix
    integer :: nrow = 0
    integer :: ncol = 0
    integer :: nnz = 0
    integer :: capacity = 0
    integer, allocatable :: row(:)
    integer, allocatable :: col(:)
    real(dp), allocatable :: val(:)
  end type triplet_matrix

contains

  subroutine triplet_init(mat, nrow, ncol, initial_capacity)
    type(triplet_matrix), intent(inout) :: mat
    integer, intent(in) :: nrow, ncol, initial_capacity
    integer :: cap

    cap = max(1, initial_capacity)
    if (allocated(mat%row)) deallocate(mat%row)
    if (allocated(mat%col)) deallocate(mat%col)
    if (allocated(mat%val)) deallocate(mat%val)
    mat%nrow = nrow
    mat%ncol = ncol
    mat%nnz = 0
    mat%capacity = cap
    allocate(mat%row(cap), mat%col(cap), mat%val(cap))
  end subroutine triplet_init

  subroutine triplet_add(mat, irow, icol, value)
    type(triplet_matrix), intent(inout) :: mat
    integer, intent(in) :: irow, icol
    real(dp), intent(in) :: value

    if (value == 0.0_dp) return
    if (irow < 1 .or. irow > mat%nrow .or. icol < 1 .or. icol > mat%ncol) then
      write(*,*) 'triplet_add: index out of range', irow, icol, mat%nrow, mat%ncol
      stop
    end if
    if (mat%nnz >= mat%capacity) call triplet_grow(mat)

    mat%nnz = mat%nnz + 1
    mat%row(mat%nnz) = irow
    mat%col(mat%nnz) = icol
    mat%val(mat%nnz) = value
  end subroutine triplet_add

  subroutine triplet_set(mat, irow, icol, value)
    type(triplet_matrix), intent(inout) :: mat
    integer, intent(in) :: irow, icol
    real(dp), intent(in) :: value

    if (irow < 1 .or. irow > mat%nrow .or. icol < 1 .or. icol > mat%ncol) then
      write(*,*) 'triplet_set: index out of range', irow, icol, mat%nrow, mat%ncol
      stop
    end if
    call triplet_add(mat, irow, icol, value)
  end subroutine triplet_set

  subroutine triplet_add_to(mat, irow, icol, value)
    type(triplet_matrix), intent(inout) :: mat
    integer, intent(in) :: irow, icol
    real(dp), intent(in) :: value

    if (value == 0.0_dp) return
    if (irow < 1 .or. irow > mat%nrow .or. icol < 1 .or. icol > mat%ncol) then
      write(*,*) 'triplet_add_to: index out of range', irow, icol, mat%nrow, mat%ncol
      stop
    end if
    call triplet_add(mat, irow, icol, value)
  end subroutine triplet_add_to

  real(dp) function triplet_get(mat, irow, icol)
    type(triplet_matrix), intent(in) :: mat
    integer, intent(in) :: irow, icol
    integer :: k

    triplet_get = 0.0_dp
    if (irow < 1 .or. irow > mat%nrow .or. icol < 1 .or. icol > mat%ncol) then
      write(*,*) 'triplet_get: index out of range', irow, icol, mat%nrow, mat%ncol
      stop
    end if
    do k = 1, mat%nnz
      if (mat%row(k) == irow .and. mat%col(k) == icol) then
        triplet_get = triplet_get + mat%val(k)
      end if
    end do
  end function triplet_get

  subroutine triplet_clear(mat)
    type(triplet_matrix), intent(inout) :: mat
    mat%nnz = 0
  end subroutine triplet_clear

  subroutine triplet_append(mat, src)
    type(triplet_matrix), intent(inout) :: mat
    type(triplet_matrix), intent(in) :: src
    integer :: first, last

    if (src%nnz <= 0) return
    if (mat%nrow /= src%nrow .or. mat%ncol /= src%ncol) then
      write(*,*) 'triplet_append: matrix shape mismatch', &
                 mat%nrow, mat%ncol, src%nrow, src%ncol
      stop
    end if

    do while (mat%capacity < mat%nnz + src%nnz)
      call triplet_grow(mat)
    end do
    first = mat%nnz + 1
    last = mat%nnz + src%nnz
    mat%row(first:last) = src%row(1:src%nnz)
    mat%col(first:last) = src%col(1:src%nnz)
    mat%val(first:last) = src%val(1:src%nnz)
    mat%nnz = last
  end subroutine triplet_append

  subroutine triplet_grow(mat)
    type(triplet_matrix), intent(inout) :: mat
    integer :: new_capacity
    integer, allocatable :: new_row(:), new_col(:)
    real(dp), allocatable :: new_val(:)

    new_capacity = max(2, mat%capacity * 2)
    allocate(new_row(new_capacity), new_col(new_capacity), new_val(new_capacity))
    if (mat%nnz > 0) then
      new_row(1:mat%nnz) = mat%row(1:mat%nnz)
      new_col(1:mat%nnz) = mat%col(1:mat%nnz)
      new_val(1:mat%nnz) = mat%val(1:mat%nnz)
    end if
    call move_alloc(new_row, mat%row)
    call move_alloc(new_col, mat%col)
    call move_alloc(new_val, mat%val)
    mat%capacity = new_capacity
  end subroutine triplet_grow

  subroutine triplet_to_csr(mat, ir, ic, a, nnz_out)
    type(triplet_matrix), intent(in) :: mat
    integer, allocatable, intent(out) :: ir(:), ic(:)
    real(dp), allocatable, intent(out) :: a(:)
    integer, intent(out) :: nnz_out

    integer :: k, k2, pos, i, s, start_pos, end_pos
    integer, allocatable :: order(:), next_pos(:), row_ptr(:), row_nnz(:)
    integer, allocatable :: sorted_row(:), sorted_col(:)
    real(dp), allocatable :: sorted_val(:)
    real(dp) :: slices(6)

    allocate(row_ptr(mat%nrow + 1))
    row_ptr = 0
    do k = 1, mat%nnz
      row_ptr(mat%row(k) + 1) = row_ptr(mat%row(k) + 1) + 1
    end do
    row_ptr(1) = 1
    do i = 1, mat%nrow
      row_ptr(i + 1) = row_ptr(i + 1) + row_ptr(i)
    end do

    allocate(order(mat%nnz), next_pos(mat%nrow))
    next_pos = row_ptr(1:mat%nrow)
    do k = 1, mat%nnz
      pos = next_pos(mat%row(k))
      order(pos) = k
      next_pos(mat%row(k)) = next_pos(mat%row(k)) + 1
    end do

    allocate(sorted_row(mat%nnz), sorted_col(mat%nnz), sorted_val(mat%nnz))
    !$omp parallel do private(k) if(mat%nnz > 10000)
    do k = 1, mat%nnz
      sorted_row(k) = mat%row(order(k))
      sorted_col(k) = mat%col(order(k))
      sorted_val(k) = mat%val(order(k))
    end do
    !$omp end parallel do

    call sort_rows_by_column(mat%nrow, row_ptr, sorted_col, sorted_val)

    allocate(row_nnz(mat%nrow), ir(mat%nrow + 1))
    row_nnz = 0
    !$omp parallel do private(i,k,start_pos,end_pos) if(mat%nrow > 1000)
    do i = 1, mat%nrow
      start_pos = row_ptr(i)
      end_pos = row_ptr(i + 1) - 1
      do k = start_pos, end_pos
        if (k == start_pos) then
          row_nnz(i) = row_nnz(i) + 1
        else if (sorted_col(k) /= sorted_col(k - 1)) then
          row_nnz(i) = row_nnz(i) + 1
        end if
      end do
    end do
    !$omp end parallel do

    ir(1) = 1
    do i = 1, mat%nrow
      ir(i + 1) = ir(i) + row_nnz(i)
    end do
    nnz_out = ir(mat%nrow + 1) - 1
    allocate(ic(nnz_out), a(nnz_out))

    !$omp parallel do private(i,k,k2,s,pos,start_pos,end_pos,slices) if(mat%nrow > 1000)
    do i = 1, mat%nrow
      start_pos = row_ptr(i)
      end_pos = row_ptr(i + 1) - 1
      pos = ir(i) - 1
      k = start_pos
      do while (k <= end_pos)
        pos = pos + 1
        ic(pos) = sorted_col(k)
        slices = 0.0_dp
        k2 = k
        do while (k2 <= end_pos .and. sorted_col(k2) == sorted_col(k))
          call oz_add_slice6_module(slices, sorted_val(k2))
          k2 = k2 + 1
        end do
        a(pos) = 0.0_dp
        do s = 6, 1, -1
          a(pos) = a(pos) + slices(s)
        end do
        k = k2
      end do
    end do
    !$omp end parallel do
  end subroutine triplet_to_csr

  subroutine sort_rows_by_column(nrow, ir, col, val)
    integer, intent(in) :: nrow
    integer, intent(in) :: ir(:)
    integer, intent(inout) :: col(:)
    real(dp), intent(inout) :: val(:)
    integer :: i, left, right

    !$omp parallel do private(i,left,right) if(nrow > 1000)
    do i = 1, nrow
      left = ir(i)
      right = ir(i + 1) - 1
      if (right > left) call insertion_sort_pair(col, val, left, right)
    end do
    !$omp end parallel do
  end subroutine sort_rows_by_column

  subroutine insertion_sort_pair(col, val, left, right)
    integer, intent(inout) :: col(:)
    real(dp), intent(inout) :: val(:)
    integer, intent(in) :: left, right
    integer :: i, j, c
    real(dp) :: v

    do i = left + 1, right
      c = col(i)
      v = val(i)
      j = i - 1
      do while (j >= left)
        if (col(j) <= c) exit
        col(j + 1) = col(j)
        val(j + 1) = val(j)
        j = j - 1
      end do
      col(j + 1) = c
      val(j + 1) = v
    end do
  end subroutine insertion_sort_pair

  subroutine csr_transpose(nrow, ncol, ir, ic, a, irt, ict, at, nnzt)
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: ir(nrow + 1), ic(:)
    real(dp), intent(in) :: a(:)
    integer, allocatable, intent(out) :: irt(:), ict(:)
    real(dp), allocatable, intent(out) :: at(:)
    integer, intent(out) :: nnzt
    integer :: i, k, col, pos
    integer, allocatable :: next_pos(:)

    nnzt = ir(nrow + 1) - 1
    allocate(irt(ncol + 1), ict(nnzt), at(nnzt), next_pos(ncol))
    irt = 0
    do k = 1, nnzt
      irt(ic(k) + 1) = irt(ic(k) + 1) + 1
    end do
    irt(1) = 1
    do i = 1, ncol
      irt(i + 1) = irt(i + 1) + irt(i)
    end do
    next_pos = irt(1:ncol)
    do i = 1, nrow
      do k = ir(i), ir(i + 1) - 1
        col = ic(k)
        pos = next_pos(col)
        ict(pos) = i
        at(pos) = a(k)
        next_pos(col) = next_pos(col) + 1
      end do
    end do
    deallocate(next_pos)
  end subroutine csr_transpose

  subroutine csr_matmul_ozaki_diag(n, ira, ica, aval, irb, icb, bval, diag_add, lambda, irc, icc, cval, nnzc)
    integer, intent(in) :: n
    integer, intent(in) :: ira(n + 1), ica(:), irb(n + 1), icb(:)
    real(dp), intent(in) :: aval(:), bval(:), diag_add(n), lambda
    integer, allocatable, intent(out) :: irc(:), icc(:)
    real(dp), allocatable, intent(out) :: cval(:)
    integer, intent(out) :: nnzc

    integer, allocatable :: marker(:), active(:), row_count(:)
    real(dp), allocatable :: work(:,:)
    integer :: i, ka, kb, krow, jcol, nactive, idx, s, pos
    real(dp) :: aik, cij

    allocate(row_count(n))
    row_count = 0

!   Pass 1: count the final nonzeros per row.  Each C row is independent;
!   thread-local marker/work arrays avoid shared triplet insertion.
    !$omp parallel private(marker,active,work,i,ka,kb,krow,jcol,nactive,idx,s,aik,cij) if(n > 1000)
    allocate(marker(n), active(min(n, 1024)), work(6, min(n, 1024)))
    marker = 0
    !$omp do schedule(dynamic)
    do i = 1, n
      nactive = 0
      do ka = ira(i), ira(i + 1) - 1
        krow = ica(ka)
        aik = aval(ka)
        do kb = irb(krow), irb(krow + 1) - 1
          jcol = icb(kb)
          if (marker(jcol) == 0) then
            nactive = nactive + 1
            if (nactive > size(active)) call ensure_active_capacity(active, work, nactive, n)
            active(nactive) = jcol
            marker(jcol) = nactive
            work(:, nactive) = 0.0_dp
          end if
          call oz_add_product_slice6_module(work(:, marker(jcol)), aik, bval(kb))
        end do
      end do

      if (lambda /= 0.0_dp) then
        if (marker(i) == 0) then
          nactive = nactive + 1
          if (nactive > size(active)) call ensure_active_capacity(active, work, nactive, n)
          active(nactive) = i
          marker(i) = nactive
          work(:, nactive) = 0.0_dp
        end if
        call oz_add_product_slice6_module(work(:, marker(i)), lambda, diag_add(i))
      end if

      row_count(i) = 0
      do idx = 1, nactive
        jcol = active(idx)
        cij = 0.0_dp
        do s = 6, 1, -1
          cij = cij + work(s, idx)
        end do
        if (cij /= 0.0_dp) row_count(i) = row_count(i) + 1
        marker(jcol) = 0
      end do
    end do
    !$omp end do
    deallocate(marker, active, work)
    !$omp end parallel

    allocate(irc(n + 1))
    irc(1) = 1
    do i = 1, n
      irc(i + 1) = irc(i) + row_count(i)
    end do
    nnzc = irc(n + 1) - 1
    allocate(icc(nnzc), cval(nnzc))

!   Pass 2: write each row directly into its CSR segment.
    !$omp parallel private(marker,active,work,i,ka,kb,krow,jcol,nactive,idx,s,pos,aik,cij) if(n > 1000)
    allocate(marker(n), active(min(n, 1024)), work(6, min(n, 1024)))
    marker = 0
    !$omp do schedule(dynamic)
    do i = 1, n
      nactive = 0
      do ka = ira(i), ira(i + 1) - 1
        krow = ica(ka)
        aik = aval(ka)
        do kb = irb(krow), irb(krow + 1) - 1
          jcol = icb(kb)
          if (marker(jcol) == 0) then
            nactive = nactive + 1
            if (nactive > size(active)) call ensure_active_capacity(active, work, nactive, n)
            active(nactive) = jcol
            marker(jcol) = nactive
            work(:, nactive) = 0.0_dp
          end if
          call oz_add_product_slice6_module(work(:, marker(jcol)), aik, bval(kb))
        end do
      end do

      if (lambda /= 0.0_dp) then
        if (marker(i) == 0) then
          nactive = nactive + 1
          if (nactive > size(active)) call ensure_active_capacity(active, work, nactive, n)
          active(nactive) = i
          marker(i) = nactive
          work(:, nactive) = 0.0_dp
        end if
        call oz_add_product_slice6_module(work(:, marker(i)), lambda, diag_add(i))
      end if

      pos = irc(i)
      do idx = 1, nactive
        jcol = active(idx)
        cij = 0.0_dp
        do s = 6, 1, -1
          cij = cij + work(s, idx)
        end do
        if (cij /= 0.0_dp) then
          icc(pos) = jcol
          cval(pos) = cij
          pos = pos + 1
        end if
        marker(jcol) = 0
      end do
    end do
    !$omp end do
    deallocate(marker, active, work)
    !$omp end parallel

    call sort_rows_by_column(n, irc, icc, cval)
    deallocate(row_count)
  end subroutine csr_matmul_ozaki_diag

  subroutine csr_matmul_ozaki_diag_upper(n, ira, ica, aval, irb, icb, bval, &
      diag_add, lambda, max_nnz, irc, icc, cval, nnzc, status)
    use, intrinsic :: iso_fortran_env, only: int64
    integer, intent(in) :: n
    integer, intent(in) :: ira(n + 1), ica(:), irb(n + 1), icb(:)
    real(dp), intent(in) :: aval(:), bval(:), diag_add(n), lambda
    integer(int64), intent(in) :: max_nnz
    integer, allocatable, intent(out) :: irc(:), icc(:)
    real(dp), allocatable, intent(out) :: cval(:)
    integer, intent(out) :: nnzc, status

    integer, allocatable :: marker(:), active(:), row_count(:)
    real(dp), allocatable :: work(:,:)
    integer :: i, ka, kb, krow, jcol, nactive, idx, s, pos
    integer(int64) :: total_nnz
    real(dp) :: aik, cij

    status = -1001
    nnzc = 0
    allocate(row_count(n))
    row_count = 0

!   Count only the upper triangle of A*B.  For B=J and A=J^T this is
!   the exact symmetric normal matrix required by PARDISO mtype=2.
    !$omp parallel private(marker,active,work,i,ka,kb,krow,jcol,nactive,idx,s,aik,cij) if(n > 1000)
    allocate(marker(n), active(min(n, 1024)), work(6, min(n, 1024)))
    marker = 0
    !$omp do schedule(dynamic)
    do i = 1, n
      nactive = 0
      do ka = ira(i), ira(i + 1) - 1
        krow = ica(ka)
        aik = aval(ka)
        do kb = irb(krow), irb(krow + 1) - 1
          jcol = icb(kb)
          if (jcol < i) cycle
          if (marker(jcol) == 0) then
            nactive = nactive + 1
            if (nactive > size(active)) call ensure_active_capacity(active, work, nactive, n)
            active(nactive) = jcol
            marker(jcol) = nactive
            work(:, nactive) = 0.0_dp
          end if
          call oz_add_product_slice6_module(work(:, marker(jcol)), aik, bval(kb))
        end do
      end do

      if (lambda /= 0.0_dp) then
        if (marker(i) == 0) then
          nactive = nactive + 1
          if (nactive > size(active)) call ensure_active_capacity(active, work, nactive, n)
          active(nactive) = i
          marker(i) = nactive
          work(:, nactive) = 0.0_dp
        end if
        call oz_add_product_slice6_module(work(:, marker(i)), lambda, diag_add(i))
      end if

      do idx = 1, nactive
        jcol = active(idx)
        cij = 0.0_dp
        do s = 6, 1, -1
          cij = cij + work(s, idx)
        end do
        if (cij /= 0.0_dp) row_count(i) = row_count(i) + 1
        marker(jcol) = 0
      end do
    end do
    !$omp end do
    deallocate(marker, active, work)
    !$omp end parallel

    total_nnz = 0_int64
    do i = 1, n
      total_nnz = total_nnz + int(row_count(i), int64)
    end do
    if (total_nnz > max_nnz .or. total_nnz >= int(huge(1), int64)) then
      status = -2001
      deallocate(row_count)
      return
    end if

    nnzc = int(total_nnz)
    allocate(irc(n + 1), icc(nnzc), cval(nnzc))
    irc(1) = 1
    do i = 1, n
      irc(i + 1) = irc(i) + row_count(i)
    end do

!   Recompute the values and write each independent row directly to CSR.
    !$omp parallel private(marker,active,work,i,ka,kb,krow,jcol,nactive,idx,s,pos,aik,cij) if(n > 1000)
    allocate(marker(n), active(min(n, 1024)), work(6, min(n, 1024)))
    marker = 0
    !$omp do schedule(dynamic)
    do i = 1, n
      nactive = 0
      do ka = ira(i), ira(i + 1) - 1
        krow = ica(ka)
        aik = aval(ka)
        do kb = irb(krow), irb(krow + 1) - 1
          jcol = icb(kb)
          if (jcol < i) cycle
          if (marker(jcol) == 0) then
            nactive = nactive + 1
            if (nactive > size(active)) call ensure_active_capacity(active, work, nactive, n)
            active(nactive) = jcol
            marker(jcol) = nactive
            work(:, nactive) = 0.0_dp
          end if
          call oz_add_product_slice6_module(work(:, marker(jcol)), aik, bval(kb))
        end do
      end do

      if (lambda /= 0.0_dp) then
        if (marker(i) == 0) then
          nactive = nactive + 1
          if (nactive > size(active)) call ensure_active_capacity(active, work, nactive, n)
          active(nactive) = i
          marker(i) = nactive
          work(:, nactive) = 0.0_dp
        end if
        call oz_add_product_slice6_module(work(:, marker(i)), lambda, diag_add(i))
      end if

      pos = irc(i)
      do idx = 1, nactive
        jcol = active(idx)
        cij = 0.0_dp
        do s = 6, 1, -1
          cij = cij + work(s, idx)
        end do
        if (cij /= 0.0_dp) then
          icc(pos) = jcol
          cval(pos) = cij
          pos = pos + 1
        end if
        marker(jcol) = 0
      end do
    end do
    !$omp end do
    deallocate(marker, active, work)
    !$omp end parallel

    call sort_rows_by_column(n, irc, icc, cval)
    deallocate(row_count)
    status = 0
  end subroutine csr_matmul_ozaki_diag_upper

  subroutine csr_matmul_ozaki_diag_values(n, ira, ica, aval, irb, icb, bval, diag_add, lambda, irc, icc, cval)
    integer, intent(in) :: n
    integer, intent(in) :: ira(n + 1), ica(:), irb(n + 1), icb(:), irc(n + 1), icc(:)
    real(dp), intent(in) :: aval(:), bval(:), diag_add(n), lambda
    real(dp), intent(inout) :: cval(:)
    integer, allocatable :: marker(:), active(:)
    real(dp), allocatable :: work(:,:)
    integer :: i, ka, kb, kc, krow, jcol, local_col, nrow_nnz, s, pos
    real(dp) :: aik, cij

!   Reuse the CSR pattern built by csr_matmul_ozaki_diag and refresh only
!   values.  This is much cheaper than rebuilding row structure each step.
    cval = 0.0_dp
    !$omp parallel private(marker,active,work,i,ka,kb,kc,krow,jcol,local_col,nrow_nnz,s,pos,aik,cij) if(n > 1000)
    allocate(marker(n), active(min(n, 1024)), work(6, min(n, 1024)))
    marker = 0
    !$omp do schedule(dynamic)
    do i = 1, n
      nrow_nnz = irc(i + 1) - irc(i)
      if (nrow_nnz > size(active)) call ensure_active_capacity(active, work, nrow_nnz, n)
      do pos = 1, nrow_nnz
        kc = irc(i) + pos - 1
        marker(icc(kc)) = pos
        work(:, pos) = 0.0_dp
      end do

      do ka = ira(i), ira(i + 1) - 1
        krow = ica(ka)
        aik = aval(ka)
        do kb = irb(krow), irb(krow + 1) - 1
          jcol = icb(kb)
          local_col = marker(jcol)
          if (local_col > 0) call oz_add_product_slice6_module(work(:, local_col), aik, bval(kb))
        end do
      end do

      if (lambda /= 0.0_dp) then
        local_col = marker(i)
        if (local_col > 0) call oz_add_product_slice6_module(work(:, local_col), lambda, diag_add(i))
      end if

      do pos = 1, nrow_nnz
        kc = irc(i) + pos - 1
        cij = 0.0_dp
        do s = 6, 1, -1
          cij = cij + work(s, pos)
        end do
        cval(kc) = cij
        marker(icc(kc)) = 0
      end do
    end do
    !$omp end do
    deallocate(marker, active, work)
    !$omp end parallel
  end subroutine csr_matmul_ozaki_diag_values

  subroutine drop_csr_ilut_pattern(n, ir, ic, a, drop_rel, max_per_row, iro, ico, ao, nnzo)
    integer, intent(in) :: n, max_per_row
    integer, intent(in) :: ir(n + 1), ic(:)
    real(dp), intent(in) :: a(:), drop_rel
    integer, allocatable, intent(out) :: iro(:), ico(:)
    real(dp), allocatable, intent(out) :: ao(:)
    integer, intent(out) :: nnzo

    integer :: i, k, nrow, ncand, nkeep, max_row, pos, diag_pos
    integer, allocatable :: row_count(:), cand_col(:)
    real(dp), allocatable :: cand_val(:), cand_abs(:)
    real(dp) :: row_max, drop_abs

    max_row = 1
    do i = 1, n
      max_row = max(max_row, ir(i + 1) - ir(i))
    end do
    allocate(row_count(n), cand_col(max_row), cand_val(max_row), cand_abs(max_row))

    do i = 1, n
      row_max = 0.0_dp
      do k = ir(i), ir(i + 1) - 1
        row_max = max(row_max, abs(a(k)))
      end do
      drop_abs = drop_rel * row_max
      ncand = 0
      do k = ir(i), ir(i + 1) - 1
        if (ic(k) == i .or. abs(a(k)) >= drop_abs) then
          ncand = ncand + 1
          cand_col(ncand) = ic(k)
          cand_val(ncand) = a(k)
          cand_abs(ncand) = abs(a(k))
        end if
      end do
      if (ncand == 0) then
        row_count(i) = 1
      else
        call sort_drop_candidates(ncand, cand_col, cand_val, cand_abs)
        row_count(i) = min(ncand, max(1, max_per_row))
      end if
    end do

    allocate(iro(n + 1))
    iro(1) = 1
    do i = 1, n
      iro(i + 1) = iro(i) + row_count(i)
    end do
    nnzo = iro(n + 1) - 1
    allocate(ico(nnzo), ao(nnzo))

    do i = 1, n
      row_max = 0.0_dp
      do k = ir(i), ir(i + 1) - 1
        row_max = max(row_max, abs(a(k)))
      end do
      drop_abs = drop_rel * row_max
      ncand = 0
      diag_pos = 0
      do k = ir(i), ir(i + 1) - 1
        if (ic(k) == i .or. abs(a(k)) >= drop_abs) then
          ncand = ncand + 1
          cand_col(ncand) = ic(k)
          cand_val(ncand) = a(k)
          cand_abs(ncand) = abs(a(k))
          if (ic(k) == i) diag_pos = ncand
        end if
      end do
      if (ncand == 0) then
        ico(iro(i)) = i
        ao(iro(i)) = 1.0_dp
      else
        call sort_drop_candidates(ncand, cand_col, cand_val, cand_abs)
        nkeep = min(ncand, max(1, max_per_row))
        diag_pos = 0
        do k = 1, ncand
          if (cand_col(k) == i) diag_pos = k
        end do
        if (diag_pos > nkeep) then
          cand_col(nkeep) = cand_col(diag_pos)
          cand_val(nkeep) = cand_val(diag_pos)
          cand_abs(nkeep) = cand_abs(diag_pos)
        end if
        call sort_kept_by_column(nkeep, cand_col, cand_val)
        pos = iro(i)
        do k = 1, nkeep
          ico(pos) = cand_col(k)
          ao(pos) = cand_val(k)
          pos = pos + 1
        end do
      end if
    end do

    deallocate(row_count, cand_col, cand_val, cand_abs)
  end subroutine drop_csr_ilut_pattern

  subroutine sort_drop_candidates(ncand, col, val, absv)
    integer, intent(in) :: ncand
    integer, intent(inout) :: col(:)
    real(dp), intent(inout) :: val(:), absv(:)
    integer :: i, j, c
    real(dp) :: v, av

    do i = 2, ncand
      c = col(i)
      v = val(i)
      av = absv(i)
      j = i - 1
      do while (j >= 1)
        if (absv(j) >= av) exit
        col(j + 1) = col(j)
        val(j + 1) = val(j)
        absv(j + 1) = absv(j)
        j = j - 1
      end do
      col(j + 1) = c
      val(j + 1) = v
      absv(j + 1) = av
    end do
  end subroutine sort_drop_candidates

  subroutine sort_kept_by_column(nkeep, col, val)
    integer, intent(in) :: nkeep
    integer, intent(inout) :: col(:)
    real(dp), intent(inout) :: val(:)
    integer :: i, j, c
    real(dp) :: v

    do i = 2, nkeep
      c = col(i)
      v = val(i)
      j = i - 1
      do while (j >= 1)
        if (col(j) <= c) exit
        col(j + 1) = col(j)
        val(j + 1) = val(j)
        j = j - 1
      end do
      col(j + 1) = c
      val(j + 1) = v
    end do
  end subroutine sort_kept_by_column

  subroutine ensure_active_capacity(active, work, needed, nmax)
    integer, allocatable, intent(inout) :: active(:)
    real(dp), allocatable, intent(inout) :: work(:,:)
    integer, intent(in) :: needed, nmax
    integer :: old_cap, new_cap
    integer, allocatable :: new_active(:)
    real(dp), allocatable :: new_work(:,:)

    old_cap = size(active)
    if (needed <= old_cap) return
    new_cap = min(nmax, max(needed, old_cap * 2))
    allocate(new_active(new_cap), new_work(6, new_cap))
    if (old_cap > 0) then
      new_active(1:old_cap) = active(1:old_cap)
      new_work(:,1:old_cap) = work(:,1:old_cap)
    end if
    call move_alloc(new_active, active)
    call move_alloc(new_work, work)
  end subroutine ensure_active_capacity

  subroutine oz_add_product_slice6_module(slices, dx, dy)
    real(dp), intent(inout) :: slices(6)
    real(dp), intent(in) :: dx, dy
    real(dp) :: prod, chk, cx, xh, xl, cy, yh, yl, err
    real(dp), parameter :: split = 134217729.0_dp

    prod = dx * dy
    cx = split * dx
    xh = cx - (cx - dx)
    xl = dx - xh
    cy = split * dy
    yh = cy - (cy - dy)
    yl = dy - yh
    chk = ((xh*yh - prod) + xh*yl + xl*yh) + xl*yl
    err = chk
    call oz_add_slice6_module(slices, prod)
    call oz_add_slice6_module(slices, err)
  end subroutine oz_add_product_slice6_module

  subroutine oz_add_slice6_module(slices, value)
    real(dp), intent(inout) :: slices(6)
    real(dp), intent(in) :: value
    real(dp) :: y, t, e
    integer :: k

    y = value
    do k = 1, 6
      t = slices(k) + y
      if (abs(slices(k)) >= abs(y)) then
        e = (slices(k) - t) + y
      else
        e = (y - t) + slices(k)
      end if
      slices(k) = t
      y = e
      if (y == 0.0_dp) exit
    end do
  end subroutine oz_add_slice6_module

end module sparse_triplet_builder
