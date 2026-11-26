!   calculation of collocation by Newton Raphson(calcollQ_TKAKI_newton.f)
!                                 coded by T.Okubo
!                                 at 2015.7.2.(thu)
!
!     [ calcollQ_TAKAKI_newton_q.f]
!        FORTRAN90/95
!
!
!     *****************************************************************
!
program main
!
   implicit none
   real(16) :: zzx(0:210),zzy(0:210),deltax(0:210),deltay(0:210)
   integer :: i,nax,nay,alp,beta
   character :: fname*100
!
   write(*,*)'input inner collocation number(na1,na2)= ?'
   read(*,*) nax,nay
   write(*,*)'Jacobi parameter (alp, beta)(2)= ?'
   read(*,*) alp,beta
!     
   call collo_newton(nax,zzx,deltax,alp,beta)
!
   call collo_newton(nay,zzy,deltay,alp,beta)

   write(*,*)'input file name of collocation point data : ?'
   read(*,*) fname
   open(2,file=fname)     
      write(2,'(2i5)') nax,nay
      do i=0,nax+1
        write(2,'(d40.32)') zzx(i)
      end do
      do i=0,nay+1
        write(2,'(d40.32)') zzy(i)
      end do
      do i=0,nax
        write(2,'(d40.32)') deltax(i)
      end do 
      do i=0,nay
        write(2,'(d40.32)') deltay(i)
      end do 
   close(2)
!
!          
end program main
!
!
!     *****************************************************************
!
subroutine collo_newton(na1,zz,delta,alp,beta)
!
      implicit none
      real(16) :: zz(0:210),delta(0:210)
      real(16) :: z,del,E,dz
      integer :: na1,alp,beta
      integer :: i,k,ina
      zz(0)=0.0q0
      do i=0,na1/2
        z=zi0(i,na1,alp,beta)
        E=1.0q0
        do while(E>1.0q-30)
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
!
!
      contains
!     ***************************************************************
!
      real(16) function zi0(i,na1,alp,beta)  
!
      implicit none
      real(16) :: pai
      integer :: i,na1,alp,beta
!
      pai=3.1415926535897932384626433832795028841971q0
!
      zi0=1.0q0/2.0q0*(1.0q0-qcos((i+1.0q0-(1.0q0-alp-beta)/4.0q0)&
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
      real(16):: z
      real(16) :: dr,cr,pr,qr,p(0:210)
      integer :: na1,alp,beta
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
      real(16) :: z
      real(16) :: dr,cr,pr,qr
      integer :: na1,alp,beta
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
end subroutine collo_newton
