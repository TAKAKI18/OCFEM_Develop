!   calculation of weight at integral point (calweight_TAKAKI_d.f90)
!                                 coded by T.Okubo
!                                 at 2015.7.2.(thu)
!
!     [ calweight_TAKAKI_d.f90]
!        FORTRAN90/95
!
!
!     *****************************************************************
!
program main
!
   implicit none
   real(16) :: zz1(0:210),zz2(0:210),delta(0:210)
   real(16) :: we11(0:210),we12(0:210),we21(0:210),we22(0:210)
   integer :: i,na1,na2,alp,beta
   character :: fname*100
!
   write(*,*)'input inner collocation number(na1,na2)= ?'
   read(*,*) na1,na2
   write(*,*)'Jacobi parameter (alp, beta)(0,0)= ?'
   read(*,*) alp,beta
!
   call collo_newton(na1,zz1,delta,alp,beta)
!
   call collo_weight(na1,zz1,alp,beta,we11,we12)
!
   call collo_newton(na2,zz2,delta,alp,beta)
!
   call collo_weight(na2,zz2,alp,beta,we21,we22)
!         
   write(*,*)'input file name of weight data at collocation: ?'
   read(*,*) fname
   open(2,file=fname)     
      write(2,'(2i5)') na1,na2
      do i=1,na1
        write(2,'(3d25.16)') zz1(i),we11(i),we12(i)
      end do
      do i=1,na2
        write(2,'(3d25.16)') zz2(i),we21(i),we22(i)
      end do
   close(2)
!
!
!          
      contains
!     *****************************************************************
!
subroutine collo_newton(na1,zz,delta,alp,beta)
!
      implicit none
      real(16),intent(out):: zz(0:210),delta(0:210)
	  real(16) :: z,E,del,dz
      integer,intent(in) :: na1,alp,beta
      integer :: i,k,ina
      zz(0)=0.0q0
      do i=0,na1/2
        z=zi0(i,na1,alp,beta)
        E=1.0q0
        do while(E>1.0d-30)
          del=z-zz(i)
          dz=-pp(z,na1,alp,beta)/dpp(z,na1,alp,beta)
          del=del+dz
          z=z+dz
          E=qabs(dz/del)
          write(*,*) 'del=',del
          write(*,*) 'z=',z
          write(*,*) 'E=',E
        end do
        zz(i+1)=z
        delta(i)=del
      end do
      if(mod(na1,2)/=0) then
        ina=na1/2+1
        zz(ina)=0.5q0
      end if
      do k=0,na1/2
        delta(na1-k)=delta(k)
      end do
      do k=1,na1
        zz(k)=zz(k-1)+delta(k-1)
      end do
      zz(na1+1)=1.0q0
!
end subroutine collo_newton
!
!     ***************************************************************
subroutine collo_weight(na1,zz,alp,beta,we1,we2)
!
      implicit none
	  real(16),intent(out):: we1(0:210),we2(0:210)
      real(16),intent(in):: zz(0:210)
      integer,intent(in) :: na1,alp,beta
      integer :: i,k,ina
!
      do i=1,na1
	  we1(i)=4.0q0/na1**2*zz(i)*(1.0q0-zz(i))/pp(zz(i),na1-1,alp,beta)**2
	  we2(i)=1.0q0/(zz(i)*(1.0q0-zz(i))*dpp(zz(i),na1,alp,beta)**2)
      end do
!
end subroutine collo_weight
!
!     ***************************************************************
!
      real(16) function zi0(i,na1,alp,beta)  
!
      implicit none
      real(16) :: pai
      integer :: i,na1,alp,beta
!
      pai=3.1415926535897932384626433832795028841971d0
!
      zi0=1.0q0/2.0d0*(1.0q0-qcos((i+1.0q0-(1.0q0-alp-beta)/4.0q0)&
          /(na1+(1.0q0+alp+beta)/2.0q0)*pai))
      write(*,*) zi0
!
      end function zi0
!
!     ****************************************************************
!
      real(16) function pp(z,na1,alp,beta)
!
      implicit none
      real(16),intent(in) :: z
      real(16) :: dr,cr,pr,qr,p(0:210)
      integer,intent(in) :: na1,alp,beta
      integer :: r
!
      p(0)=1.0q0
      p(1)=1.0q0-(alp+beta+2.0q0)/(alp+1.0q0)*z
        do r=1,na1-1
          dr=(2.0q0*r+alp+beta+2.0q0)/(r+alp+1.0q0)
          cr=1.0q0/(r+alp+beta+1.0q0)*dr
          pr=(r+beta)*r/(2.0q0*r+alp+beta)*cr
          qr=(2.0q0*r+alp+beta+1.0q0)*cr
          p(r+1)=(1.0q0+pr-qr*z)*p(r)-pr*p(r-1)
          write(*,*)'p(',r+1,')=',p(r+1)
        end do
      pp=p(na1)
      write(*,*)'pp=',pp
!
      end function pp
!
!     ****************************************************************
!
      real(16) function dpp(z,na1,alp,beta)
!
      implicit none
      real(16),intent(in) :: z
      real(16) :: dr,cr,pr,qr
      integer,intent(in) :: na1,alp,beta
!
        dr=(2.0q0*na1+alp+beta+2.0q0)/(na1+alp+1.0q0)
        cr=1.0q0/(na1+alp+beta+1.0q0)*dr
        pr=(na1+beta)*na1/(2.0q0*na1+alp+beta)*cr
        qr=(2.0q0*na1+alp+beta+1.0q0)*cr
        dpp=(pp(z,na1+1,alp,beta)-(1.0q0-dr*z)*pp(z,na1,alp,beta))/&
        (cr*z*(1.0q0-z))
      write(*,*)'dpp=',dpp
!
      end function dpp  
!
  end program main

