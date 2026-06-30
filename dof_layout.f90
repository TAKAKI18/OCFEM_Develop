module dof_layout
  implicit none

  integer, parameter :: block_u = 1
  integer, parameter :: block_v = 2
  integer, parameter :: block_w = 3
  integer, parameter :: block_p = 4

contains

  integer function global_dof_index(itime, inode, iblock, nx_local, iat_local)
    integer, intent(in) :: itime, inode, iblock, nx_local, iat_local

    ! Unknown/correction block order: U, V, W, P.
    global_dof_index = itime + (inode - 1) * iat_local + (iblock - 1) * nx_local * iat_local
  end function global_dof_index

  integer function global_equation_block(local_block)
    integer, intent(in) :: local_block

    ! Local FFUVP/BBUVP order is continuity,U,V,W; global row order is continuity,V,W,U.
    select case (local_block)
    case (2)
      global_equation_block = block_p
    case (3)
      global_equation_block = block_v
    case (4)
      global_equation_block = block_w
    case default
      global_equation_block = local_block
    end select
  end function global_equation_block

end module dof_layout
