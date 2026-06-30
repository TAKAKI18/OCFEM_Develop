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

  subroutine triplet_clear(mat)
    type(triplet_matrix), intent(inout) :: mat
    mat%nnz = 0
  end subroutine triplet_clear

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

    integer :: k, pos, i, start_pos, end_pos
    integer, allocatable :: order(:), next_pos(:), row_ptr(:)
    integer, allocatable :: sorted_row(:), sorted_col(:)
    real(dp), allocatable :: sorted_val(:)

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
    do k = 1, mat%nnz
      sorted_row(k) = mat%row(order(k))
      sorted_col(k) = mat%col(order(k))
      sorted_val(k) = mat%val(order(k))
    end do

    call sort_rows_by_column(mat%nrow, row_ptr, sorted_col, sorted_val)

    allocate(ic(mat%nnz), a(mat%nnz))
    allocate(ir(mat%nrow + 1))
    nnz_out = 0
    do i = 1, mat%nrow
      ir(i) = nnz_out + 1
      start_pos = row_ptr(i)
      end_pos = row_ptr(i + 1) - 1
      do k = start_pos, end_pos
        if (k == start_pos) then
          nnz_out = nnz_out + 1
          ic(nnz_out) = sorted_col(k)
          a(nnz_out) = sorted_val(k)
        else if (sorted_col(k) == ic(nnz_out)) then
          a(nnz_out) = a(nnz_out) + sorted_val(k)
        else
          nnz_out = nnz_out + 1
          ic(nnz_out) = sorted_col(k)
          a(nnz_out) = sorted_val(k)
        end if
      end do
    end do
    ir(mat%nrow + 1) = nnz_out + 1
  end subroutine triplet_to_csr

  subroutine sort_rows_by_column(nrow, ir, col, val)
    integer, intent(in) :: nrow
    integer, intent(in) :: ir(:)
    integer, intent(inout) :: col(:)
    real(dp), intent(inout) :: val(:)
    integer :: i, left, right

    do i = 1, nrow
      left = ir(i)
      right = ir(i + 1) - 1
      if (right > left) call insertion_sort_pair(col, val, left, right)
    end do
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

end module sparse_triplet_builder
