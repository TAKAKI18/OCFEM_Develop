!     CALCULATION OF OPERATOR MATRIX  Am (1-DIMENSIONAL)
!     Method  BY Takahiro Hasegawa
!            CODED BY T.OKUBO
!            HAKODATE NATIONAL COLLEGE OF TECHNOLOGY
!     PROGRAM [ COMV1T_TAKAKI_q.f90 ]
!
module collo_mod
      implicit none
      integer, parameter :: nta=50
      real(16) :: xx(0:nta),del(0:nta)
      integer :: nna,ian
      character fname1*100,fname2*100
end module collo_mod
!
!
program main
      use collo_mod
      implicit none
      real(16) :: temp,a(0:nta,0:nta),b(0:nta,0:nta)
      integer :: m,i,j
!     ------------ INPUT OF COLLOCATION POINT -------------------------
!
      write(*,*)'INPUT FILE NAME OF COLLOCATION POINT DATA : ?'
      read(*,*) fname1
      open(2,file=fname1)
        read(2,*) nna,nna
        ian=nna+1
        do i=0,ian
          read(2,*) xx(i)
            write(*,*) xx(i)
        end do
        do i=0,ian
          read(2,*) temp
        end do
        do i=0,ian-1
          read(2,*) del(i)
            write(*,*) del(i)
        end do
!
        close(2)
!
        m=1
        call collocation(m,a)
        m=2
        call collocation(m,b)
!
        write(*,*)'OUTPUT FILE NAME OF COLLOCATION MATRIX DATA : ?'
        read(*,*) fname2
        open(9,file=fname2)
           write(9,'(i5)') nna
           do i=0,ian
             write(9,'(d25.16)') xx(i)
           end do 
           do i=0,ian
             do j=0,ian
                write(9,'(2d25.16)') a(i,j),b(i,j)
             end do
           end do
        close(9)
!
      contains
!
!     ------------ CALCULATION OF MATRIX ------------------------------
!
!     ********************************************************
!
   real(16) function delta(i,j)
!
      use collo_mod
      implicit none
      integer :: i,j,k
      real(16) :: delt
!
      if(i > j) then
        delt=0.0q0
        do k=j,i-1
          delt=delt+del(k)
        end do
        delta=delt
      end if
!
      if(i < j) then
        delt=0.0q0
        do k=i,j-1
          delt=delt+del(k)
        end do
        delta=-delt
      end if
!
      if(i == j) then
        delta=0.0q0
      end if
!
   end function delta
!
!
! ***************************************************
!
   real(16) function T(i)
!
     use collo_mod 
     implicit none
     integer :: k,i
     real(16) :: Ti
   
     Ti=1.0q0
     do k=0,ian
       if(k /= i) then
          Ti=Ti*delta(i,k)
       end if
     end do
     T=Ti

  end function T
!
! ***************************************************
!
  subroutine collocation(m,a)
!
  use collo_mod
  implicit none
  integer :: i,j,l,N,m,imax,kk(0:50)
  real(16) :: Ri(0:nta),a(0:nta,0:nta)
  real(16) :: dl,Ti,Tj,aii,dsign
  
  N=nna
  dsign=(-1.0q0)**m
  imax=N/2+1
  do i=0,imax
     kk(0)=i
     Ti=T(i)
     do j=0,N+1
       if(i /= j) then
         kk(1)=j
         Ri(j)=recurs_R(kk,2,N,m)
         Tj=T(j)
         dl=delta(i,j)
         write(*,*) 'Ti=',Ti
         write(*,*) 'Tj=',Tj
         write(*,*) 'dl=',dl
         a(i,j)=m*Ti/(Tj*dl)*Ri(j)
       end if
     end do
     aii=0.0q0
     do l=0,N+1
       if(i /= l) then
         dl=delta(i,l)
         aii=aii+Ri(l)/dl
       end if
     end do
     a(i,i)=aii
   end do
   do i=0,imax
     do j=0,N+1
       a(N+1-i,N+1-j)=dsign*a(i,j)
     end do
   end do
!
   end subroutine collocation
!
! **********************************************************
!
   recursive function recurs_R(kk,r,N,m) result(v)
!
   implicit none
   integer :: kk(0:50),r,N,m,j,fg,k
   real(16) :: v,s
   
   if(r <= m) then
     s=0.0q0
     do k=0,N+1
       kk(r)=k
       fg=0
       do j=0,r-1
         if(kk(r) == kk(j)) then
           fg=fg+1
           exit
         end if
       end do
       if(fg == 0) then
         s=s+recurs_R(kk,r+1,N,m)/delta(kk(0),kk(r))
       end if
     end do
     v=s
   else
     v=1.0q0
   end if
   
   end function recurs_R
   
end program main



