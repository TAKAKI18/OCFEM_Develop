!     CALCULATION OF OPERATOR MATRIX  Am( m=1,m=2), Cpq (p+q=2,p=1,q=1)
!     (3-DIMENSIONAL)
!     METHOD BY Takaki Ohkubo & Takahiro Hasegawa 
!       BY T.OKUBO (MODIFIED TO 2-DIMENSIONAL)
!            CODED BY T.OKUBO   2015.7.4(Sun)
!            HAKODATE NATIONAL COLLEGE OF TECHNOLOGY
!     PROGRAM [ COMV2d_TAKAKI_q.f90 ]
!
module collo_matrix_mod
      implicit none
      integer, parameter :: NTA=50,NT=50000
      real(16) ::  XX(1:NTA),YY(1:NTA),ZZ(1:NTA),DELX(0:NTA),DELY(0:NTA),DELZ(0:NTA)&
                  ,ax(0:NTA,0:NTA,1:3),ay(0:NTA,0:NTA,1:3),az(0:NTA,0:NTA,1:3),bx(0:NTA,0:NTA,1:3),by(0:NTA,0:NTA,1:3),bz(0:NTA,0:NTA,1:3)
      real(16), allocatable :: AAX(:,:),AAY(:,:),AAZ(:,:),BBX(:,:),BBY(:,:),BBZ(:,:),CCYZ(:,:),CCZX(:,:),CCXY(:,:)
      integer ::  NNX,NNY,NNZ,NNAX,NNAY,NNAZ,IAX,IAY,IAZ,IAN
      character :: FNAME1*100,FNAME2*100
!
end module collo_matrix_mod
!
module collo_mod
      implicit none
      integer :: nna1x,nna1y,nna1z,ian1x,ian1y,ian1z
end module collo_mod
!
program main
      use collo_matrix_mod
      use collo_mod 
      implicit none
      real(16) :: temp
      integer :: m,I,J,K,L,JK,IK,J1,J2,J3,I1,I2,I3,IA
!
!     ------------ INPUT OF COLLOCATION POINT -------------------------
!
      WRITE(*,*)'INPUT FILE NAME OF COLLOCATION POINT DATA : ?'
      READ(*,*) FNAME1
      OPEN(2,FILE=FNAME1)
        READ(2,*) NNAX,NNAY,NNAZ
        NNX=NNAX+1
        NNY=NNAY+1
        NNZ=NNAZ+1
        nna1x=NNAX
        nna1y=NNAY
        nna1z=NNAZ
        IAX=NNAX+2
        IAY=NNAY+2
        IAZ=NNAZ+2
        ian1x=NNAX+1
        ian1y=NNAY+1
        ian1z=NNAZ+1
        IAN=(NNAX+2)*(NNAY+2)*(NNAZ+2)
        ALLOCATE(AAX(1:IAN,1:IAN),AAY(1:IAN,1:IAN),AAZ(1:IAN,1:IAN),BBX(1:IAN,1:IAN),BBY(1:IAN,1:IAN),BBZ(1:IAN,1:IAN),&
                 CCYZ(1:IAN,1:IAN),CCZX(1:IAN,1:IAN),CCXY(1:IAN,1:IAN))
        DO I=1,IAX
          READ(2,*) XX(I)
        END DO
        DO I=1,IAY
          READ(2,*) YY(I)
        END DO
        DO I=1,IAZ
          READ(2,*) ZZ(I)
        END DO
        
        DO I=0,NNAX
          READ(2,*) DELX(I)
        END DO
        DO I=0,NNAY
          READ(2,*) DELY(I)
        END DO
        DO I=0,NNAZ
          READ(2,*) DELZ(I)
        END DO
        
      CLOSE(2)
!
!     ------------ CALCULATION OF MATRIX ------------------------------
!
      m=1
      call collocation(m,ax,nna1x,1)
      m=1
      call collocation(m,ay,nna1y,2)
      m=1
      call collocation(m,az,nna1z,3)      
      m=2
      call collocation(m,bx,nna1x,1)
      m=2
      call collocation(m,by,nna1y,2)
      m=2
      call collocation(m,bz,nna1z,3)
      
!      
      DO J1=1,IAX
      DO J2=1,IAY
      DO J3=1,IAZ
        DO I1=1,IAX
        DO I2=1,IAY
        DO I3=1,IAZ
!
          IF((I1.EQ.1).AND.(I2.EQ.1).AND.(I3.EQ.1)) THEN
            IK=1
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.1).AND.(I3.EQ.1)) THEN
            IK=2
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.IAY).AND.(I3.EQ.1)) THEN
            IK=3
          END IF
          IF((I1.EQ.1).AND.(I2.EQ.IAY).AND.(I3.EQ.1)) THEN
            IK=4
          ENDIF
          IF((I1.EQ.1).AND.(I2.EQ.1).AND.(I3.EQ.IAZ)) THEN
            IK=5
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.1).AND.(I3.EQ.IAZ)) THEN
            IK=6
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.IAY).AND.(I3.EQ.IAZ)) THEN
            IK=7
          END IF
          IF((I1.EQ.1).AND.(I2.EQ.IAY).AND.(I3.EQ.IAZ)) THEN
            IK=8
          ENDIF
!
          IF((J1.EQ.1).AND.(J2.EQ.1).AND.(J3.EQ.1)) THEN
            JK=1
          END IF
          IF((J1.EQ.IAX).AND.(J2.EQ.1).AND.(J3.EQ.1)) THEN
            JK=2
          END IF
          IF((J1.EQ.IAX).AND.(J2.EQ.IAY).AND.(J3.EQ.1)) THEN
            JK=3
          END IF
          IF((J1.EQ.1).AND.(J2.EQ.IAY).AND.(J3.EQ.1)) THEN
            JK=4
          ENDIF
          IF((J1.EQ.1).AND.(J2.EQ.1).AND.(J3.EQ.IAZ)) THEN
            JK=5
          END IF
          IF((J1.EQ.IAX).AND.(J2.EQ.1).AND.(J3.EQ.IAZ)) THEN
            JK=6
          END IF
          IF((J1.EQ.IAX).AND.(J2.EQ.IAY).AND.(J3.EQ.IAZ)) THEN
            JK=7
          END IF
          IF((J1.EQ.1).AND.(J2.EQ.IAY).AND.(J3.EQ.IAZ)) THEN
            JK=8
          ENDIF
!
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.(I2.EQ.1).AND.(I3.EQ.1)) THEN
            IK=8+(I1-1)
          END IF
          IF((I1.EQ.IAX).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.(I3.EQ.1)) THEN
            IK=8+NNAX+(I2-1)
          END IF
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.(I2.EQ.IAY).AND.(I3.EQ.1)) THEN
            IK=8+NNAX+NNAY+(I1-1)
          END IF
          IF((I1.EQ.1).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.(I3.EQ.1)) THEN
            IK=8+2*NNAX+NNAY+(I2-1)
          END IF
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.(I2.EQ.1).AND.(I3.EQ.IAZ)) THEN
            IK=8+2*(NNAX+NNAY)+(I1-1)
          END IF
          IF((I1.EQ.IAX).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.(I3.EQ.IAZ)) THEN
            IK=8+2*(NNAX+NNAY)+NNAX+(I2-1)
          END IF
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.(I2.EQ.IAY).AND.(I3.EQ.IAZ)) THEN
            IK=8+2*(NNAX+NNAY)+NNAX+NNAY+(I1-1)
          END IF
          IF((I1.EQ.1).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.(I3.EQ.IAZ)) THEN
            IK=8+2*(NNAX+NNAY)+2*NNAX+NNAY+(I2-1)
          END IF
          IF((I1.EQ.1).AND.(I2.EQ.1).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NNAX+NNAY)+(I3-1)
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.1).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NNAX+NNAY)+NNAZ+(I3-1)
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.IAY).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NNAX+NNAY)+2*NNAZ+(I3-1)
          END IF
          IF((I1.EQ.1).AND.(I2.EQ.IAY).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NNAX+NNAY)+3*NNAZ+(I3-1)
          END IF
!
          IF(((J1.GT.1).AND.(J1.LT.IAX)).AND.(J2.EQ.1).AND.(J3.EQ.1)) THEN
            JK=8+(J1-1)
          END IF
          IF((J1.EQ.IAX).AND.((J2.GT.1).AND.(J2.LT.IAY)).AND.(J3.EQ.1)) THEN
            JK=8+NNAX+(J2-1)
          END IF
          IF(((J1.GT.1).AND.(J1.LT.IAX)).AND.(J2.EQ.IAY).AND.(J3.EQ.1)) THEN
            JK=8+NNAX+NNAY+(J1-1)
          END IF
          IF((J1.EQ.1).AND.((J2.GT.1).AND.(J2.LT.IAY)).AND.(J3.EQ.1)) THEN
            JK=8+2*NNAX+NNAY+(J2-1)
          END IF
          IF(((J1.GT.1).AND.(J1.LT.IAX)).AND.(J2.EQ.1).AND.(J3.EQ.IAZ)) THEN
            JK=8+2*(NNAX+NNAY)+(J1-1)
          END IF
          IF((J1.EQ.IAX).AND.((J2.GT.1).AND.(J2.LT.IAY)).AND.(J3.EQ.IAZ)) THEN
            JK=8+2*(NNAX+NNAY)+NNAX+(J2-1)
          END IF
          IF(((J1.GT.1).AND.(J1.LT.IAX)).AND.(J2.EQ.IAY).AND.(J3.EQ.IAZ)) THEN
            JK=8+2*(NNAX+NNAY)+NNAX+NNAY+(J1-1)
          END IF
          IF((J1.EQ.1).AND.((J2.GT.1).AND.(J2.LT.IAY)).AND.(J3.EQ.IAZ)) THEN
            JK=8+2*(NNAX+NNAY)+2*NNAX+NNAY+(J2-1)
          END IF
          IF((J1.EQ.1).AND.(J2.EQ.1).AND.((J3.GT.1).AND.(J3.LT.IAZ))) THEN
            JK=8+4*(NNAX+NNAY)+(J3-1)
          END IF
          IF((J1.EQ.IAX).AND.(J2.EQ.1).AND.((J3.GT.1).AND.(J3.LT.IAZ))) THEN
            JK=8+4*(NNAX+NNAY)+NNAZ+(J3-1)
          END IF
          IF((J1.EQ.IAX).AND.(J2.EQ.IAY).AND.((J3.GT.1).AND.(J3.LT.IAZ))) THEN
            JK=8+4*(NNAX+NNAY)+2*NNAZ+(J3-1)
          END IF
          IF((J1.EQ.1).AND.(J2.EQ.IAY).AND.((J3.GT.1).AND.(J3.LT.IAZ))) THEN
            JK=8+4*(NNAX+NNAY)+3*NNAZ+(J3-1)
          END IF
!
          IF((((I1.GT.1).AND.(I1.LT.IAX)).AND.((I2.GT.1).AND.(I2.LT.IAY))).AND.(I3.EQ.1)) THEN
            IK=8+4*(NNAX+NNAY+NNAZ)+(I2-2)*NNAX+(I1-1)
          END IF
          IF((((I1.GT.1).AND.(I1.LT.IAX)).AND.((I2.GT.1).AND.(I2.LT.IAY))).AND.(I3.EQ.IAZ)) THEN
            IK=8+4*(NNAX+NNAY+NNAZ)+NNAX*NNAY+(I2-2)*NNAX+(I1-1)
          END IF
          IF((I2.EQ.1).AND.(((I1.GT.1).AND.(I1.LT.IAX)).AND.((I3.GT.1).AND.(I3.LT.IAZ)))) THEN
            IK=8+4*(NNAX+NNAY+NNAZ)+2*NNAX*NNAY+(I3-2)*NNAX+(I1-1)
          END IF
          IF((I1.EQ.IAX).AND.(((I2.GT.1).AND.(I2.LT.IAY)).AND.((I3.GT.1).AND.(I3.LT.IAZ)))) THEN
            IK=8+4*(NNAX+NNAY+NNAZ)+2*NNAX*NNAY+NNAX*NNAZ+(I3-2)*NNAY+(I2-1)
          END IF
          IF((I2.EQ.IAY).AND.(((I1.GT.1).AND.(I1.LT.IAX)).AND.((I3.GT.1).AND.(I3.LT.IAZ)))) THEN
            IK=8+4*(NNAX+NNAY+NNAZ)+2*NNAX*NNAY+NNAX*NNAZ+NNAY*NNAZ+(I3-2)*NNAX+(I1-1)
          END IF
          IF((I1.EQ.1).AND.(((I2.GT.1).AND.(I2.LT.IAY)).AND.((I3.GT.1).AND.(I3.LT.IAZ)))) THEN
            IK=8+4*(NNAX+NNAY+NNAZ)+2*NNAX*NNAY+2*NNAX*NNAZ+NNAY*NNAZ+(I3-2)*NNAY+(I2-1)
          END IF
!
          IF((((J1.GT.1).AND.(J1.LT.IAX)).AND.((J2.GT.1).AND.(J2.LT.IAY))).AND.(J3.EQ.1)) THEN
            JK=8+4*(NNAX+NNAY+NNAZ)+(J2-2)*NNAX+(J1-1)
          END IF
          IF((((J1.GT.1).AND.(J1.LT.IAX)).AND.((J2.GT.1).AND.(J2.LT.IAY))).AND.(J3.EQ.IAZ)) THEN
            JK=8+4*(NNAX+NNAY+NNAZ)+NNAX*NNAY+(J2-2)*NNAX+(J1-1)
          END IF
          IF((J2.EQ.1).AND.(((J1.GT.1).AND.(J1.LT.IAX)).AND.((J3.GT.1).AND.(J3.LT.IAZ)))) THEN
            JK=8+4*(NNAX+NNAY+NNAZ)+2*NNAX*NNAY+(J3-2)*NNAX+(J1-1)
          END IF
          IF((J1.EQ.IAX).AND.(((J2.GT.1).AND.(J2.LT.IAY)).AND.((J3.GT.1).AND.(J3.LT.IAZ)))) THEN
            JK=8+4*(NNAX+NNAY+NNAZ)+2*NNAX*NNAY+NNAX*NNAZ+(J3-2)*NNAY+(J2-1)
          END IF
          IF((J2.EQ.IAY).AND.(((J1.GT.1).AND.(J1.LT.IAX)).AND.((J3.GT.1).AND.(J3.LT.IAZ)))) THEN
            JK=8+4*(NNAX+NNAY+NNAZ)+2*NNAX*NNAY+NNAX*NNAZ+NNAY*NNAZ+(J3-2)*NNAX+(J1-1)
          END IF
          IF((J1.EQ.1).AND.(((J2.GT.1).AND.(J2.LT.IAY)).AND.((J3.GT.1).AND.(J3.LT.IAZ)))) THEN
            JK=8+4*(NNAX+NNAY+NNAZ)+2*NNAX*NNAY+2*NNAX*NNAZ+NNAY*NNAZ+(J3-2)*NNAY+(J2-1)
          END IF          
!
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NNAX+NNAY+NNAZ)+2*(NNAX*NNAY+NNAX*NNAZ+NNAY*NNAZ)+(I3-2)*NNAX*NNAY+(I2-2)*NNAX+(I1-1)
          END IF
!
          IF(((J1.GT.1).AND.(J1.LT.IAX)).AND.((J2.GT.1).AND.(J2.LT.IAY)).AND.((J3.GT.1).AND.(J3.LT.IAZ))) THEN
            JK=8+4*(NNAX+NNAY+NNAZ)+2*(NNAX*NNAY+NNAX*NNAZ+NNAY*NNAZ)+(J3-2)*NNAX*NNAY+(J2-2)*NNAX+(J1-1)
          END IF
!
!          
!
!
         AAX(JK,IK)=ALXYM(J1-1,I1-1,1,1)*ALXY(J2-1,I2-1,2)*ALXY(J3-1,I3-1,3)
         AAY(JK,IK)=ALXY(J1-1,I1-1,1)*ALXYM(J2-1,I2-1,1,2)*ALXY(J3-1,I3-1,3)
         AAZ(JK,IK)=ALXY(J1-1,I1-1,1)*ALXY(J2-1,I2-1,2)*ALXYM(J3-1,I3-1,1,3)
         BBX(JK,IK)=ALXYM(J1-1,I1-1,2,1)*ALXY(J2-1,I2-1,2)*ALXY(J3-1,I3-1,3)
         BBY(JK,IK)=ALXY(J1-1,I1-1,1)*ALXYM(J2-1,I2-1,2,2)*ALXY(J3-1,I3-1,3)
         BBZ(JK,IK)=ALXY(J1-1,I1-1,1)*ALXY(J2-1,I2-1,2)*ALXYM(J3-1,I3-1,2,3)
         CCXY(JK,IK)=ALXYM(J1-1,I1-1,1,1)*ALXYM(J2-1,I2-1,1,2)*ALXY(J3-1,I3-1,3)
         CCYZ(JK,IK)=ALXY(J1-1,I1-1,1)*ALXYM(J2-1,I2-1,1,2)*ALXYM(J3-1,I3-1,1,3)
         CCZX(JK,IK)=ALXYM(J1-1,I1-1,1,1)*ALXY(J2-1,I2-1,2)*ALXYM(J3-1,I3-1,1,3)
!         write(*,*)'jk,ik=', jk,ik
!       
!
!
        END DO
        END DO
        END DO
      END DO
      END DO
      END DO
!      write(*,*)'JK,IK,JK*IK,IAN,IAN*IAN',JK,IK,JK*IK,IAN,IAN*IAN       
!
!     -------- OUTPUT FILE OF A,B MATRIX ------------------------
! 
      WRITE(*,*)'OUTPUT FILE NAME OF A, B, C MATRIX:?'
	READ(*,*) FNAME2 
      OPEN(3,FILE=FNAME2)
        WRITE(3,'(3I10)') NNAX,NNAY,NNAZ
        DO I=1,IAX
          WRITE(3,'(D25.16)') XX(I)
        END DO
        DO I=1,IAY
          WRITE(3,'(D25.16)') YY(I)
        END DO
        DO I=1,IAZ
          WRITE(3,'(D25.16)') ZZ(I)
        END DO
        DO I=1,IAN
        DO J=1,IAN
          WRITE(3,'(9D25.16)') AAX(I,J),AAY(I,J),AAZ(I,J),BBX(I,J),BBY(I,J),BBZ(I,J),CCXY(I,J),CCYZ(I,J),CCZX(I,J)
        END DO
        END DO
      CLOSE(3)
!
      contains
!     ******************************************************************
!
   real(16) function alxy(i,j,mk)
   
      implicit none
      integer :: i,j,mk
      
      if(i == j) then
        alxy=1.0q0
      else
        alxy=0.0q0
      end if
      
   end function alxy
!
!     *******************************************************************
!
  real(16) function alxym(i,j,m,mk)

      use collo_matrix_mod
      implicit none
      integer :: i,j,m,mk
      
      if(mk.eq.1) then
      if(m == 1) then
        alxym=ax(i,j,mk)
      end if
      if(m == 2) then
        alxym=bx(i,j,mk)
      end if
      end if
      if(mk.eq.2) then
      if(m == 1) then
        alxym=ay(i,j,mk)
      end if
      if(m == 2) then
        alxym=by(i,j,mk)
      end if
      end if
      if(mk.eq.3) then
      if(m == 1) then
        alxym=az(i,j,mk)
      end if
      if(m == 2) then
        alxym=bz(i,j,mk)
      end if
      end if
      
  end function alxym
!
!     ********************************************************
!
   real(16) function delta(i,j,mk)
!
      use collo_matrix_mod, only : delx,dely,delz
      integer :: i,j,k,mk
      real(16) :: delt
!
      if(mk.eq.1) then
      if(i > j) then
        delt=0.0q0
        do k=j,i-1
          delt=delt+delx(k)
        end do
        delta=delt
      end if
!
      if(i < j) then
        delt=0.0q0
        do k=i,j-1
          delt=delt+delx(k)
        end do
        delta=-delt
      end if
!
      if(i == j) then
        delta=0.0q0
      end if
      end if
      
      if(mk.eq.2) then
      if(i > j) then
        delt=0.0q0
        do k=j,i-1
          delt=delt+dely(k)
        end do
        delta=delt
      end if
!
      if(i < j) then
        delt=0.0q0
        do k=i,j-1
          delt=delt+dely(k)
        end do
        delta=-delt
      end if
!
      if(i == j) then
        delta=0.0q0
      end if
      end if

      if(mk.eq.3) then
      if(i > j) then
        delt=0.0q0
        do k=j,i-1
          delt=delt+delz(k)
        end do
        delta=delt
      end if
!
      if(i < j) then
        delt=0.0q0
        do k=i,j-1
          delt=delt+delz(k)
        end do
        delta=-delt
      end if
!
      if(i == j) then
        delta=0.0q0
      end if
      end if
!
   end function delta
!
!
! ***************************************************
!
   real(16) function T(i,ian2,mk)
!
     use collo_mod 
     implicit none
     integer :: k,i,ian2,mk
     real(16) :: Ti
   
     
     Ti=1.0q0
     do k=0,ian2
       if(k /= i) then
          Ti=Ti*delta(i,k,mk)
       end if
     end do
     T=Ti

  end function T
!
! ***************************************************
!
  subroutine collocation(m,a,N,mk)

  use collo_mod
  implicit none
  integer :: i,j,l,N,m,imax,mk,kk(0:50)
  real(16) :: Ri(0:nt),a(0:nta,0:nta,1:3)
  real(16) :: dl,Ti,Tj,aii,dsign
  
  dsign=(-1.0q0)**m
  imax=N/2+1
  do i=0,imax
     kk(0)=i
     Ti=T(i,N+1,mk)
     do j=0,N+1
       if(i /= j) then
         kk(1)=j
         Ri(j)=recurs_R(kk,2,N,m,mk)
         Tj=T(j,N+1,mk)
         dl=delta(i,j,mk)
!         write(*,*) 'Ti=',Ti
!         write(*,*) 'Tj=',Tj
!         write(*,*) 'dl=',dl
         a(i,j,mk)=m*Ti/(Tj*dl)*Ri(j)
       end if
     end do
     aii=0.0q0
     do l=0,N+1
       if(i /= l) then
         dl=delta(i,l,mk)
         aii=aii+Ri(l)/dl
       end if
     end do
     a(i,i,mk)=aii
   end do
   do i=0,imax
     do j=0,N+1
       a(N+1-i,N+1-j,mk)=dsign*a(i,j,mk)
     end do
   end do
!
   end subroutine collocation
!
! **********************************************************
!
   recursive function recurs_R(kk,r,N,m,mk) result(v)
!
   implicit none
   integer :: kk(0:50),r,N,m,j,fg,k,mk
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
         s=s+recurs_R(kk,r+1,N,m,mk)/delta(kk(0),kk(r),mk)
       end if
     end do
     v=s
   else
     v=1.0q0
   end if
   
   end function recurs_R
   
end program main

