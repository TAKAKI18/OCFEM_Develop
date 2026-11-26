program test_sparse_triplet_builder
  use sparse_triplet_builder
  implicit none

  type(triplet_matrix) :: mat
  integer, allocatable :: ir(:), ic(:)
  real(dp), allocatable :: a(:)
  integer :: nnz

  call triplet_init(mat, 3, 3, 2)
  call triplet_add(mat, 2, 3, 4.0_dp)
  call triplet_add(mat, 1, 2, 1.0_dp)
  call triplet_add(mat, 1, 1, 2.0_dp)
  call triplet_add(mat, 1, 2, 3.0_dp)
  call triplet_add(mat, 3, 1, 5.0_dp)

  call triplet_to_csr(mat, ir, ic, a, nnz)

  if (nnz /= 4) stop 1
  if (any(ir /= [1, 3, 4, 5])) stop 2
  if (any(ic(1:nnz) /= [1, 2, 3, 1])) stop 3
  if (abs(a(1) - 2.0_dp) > 1.0e-12_dp) stop 4
  if (abs(a(2) - 4.0_dp) > 1.0e-12_dp) stop 5
  if (abs(a(3) - 4.0_dp) > 1.0e-12_dp) stop 6
  if (abs(a(4) - 5.0_dp) > 1.0e-12_dp) stop 7

  call triplet_init(mat, 4, 4, 1)
  call triplet_add(mat, 4, 4, 1.0_dp)
  call triplet_add(mat, 2, 1, -2.0_dp)
  call triplet_add(mat, 2, 1, 0.5_dp)
  call triplet_add(mat, 2, 3, 7.0_dp)
  call triplet_add(mat, 4, 2, 3.0_dp)
  call triplet_to_csr(mat, ir, ic, a, nnz)

  if (nnz /= 4) stop 8
  if (any(ir /= [1, 1, 3, 3, 5])) stop 9
  if (any(ic(1:nnz) /= [1, 3, 2, 4])) stop 10
  if (abs(a(1) + 1.5_dp) > 1.0e-12_dp) stop 11
  if (abs(a(2) - 7.0_dp) > 1.0e-12_dp) stop 12
  if (abs(a(3) - 3.0_dp) > 1.0e-12_dp) stop 13
  if (abs(a(4) - 1.0_dp) > 1.0e-12_dp) stop 14

  write(*,*) 'sparse_triplet_builder test passed'
end program test_sparse_triplet_builder
