!     NUMERICAL ANALISYS OF UNSTEADY NAVIER_STOKES EQUATION
!         BY ORTHOGONAL COLLOCATION FINITE ELEMENT METHOD
!          
!
!         PROGRAM :  (UnSteady_Navire_Stokes_Straight_Duct_OCFEM_takaki_intel_LU_3D_MM_Bound2_OMP_lmn.f90)
!
!              CODED BY T.Okubo from 2011.3.31.(Thurs)
!              modified by T.Ohkubo  2022.7.29.(Fri)
!                                    2022.8.23.(Tues)
!                                    2025.1.9.(Thurs)
!                                    2026.
!     ************* Means of variable and parameter *******************
!
!     ------------  ALLAY PARAMETER OF DIMENSION  ----------------------
!     KNX           :MAXIMUM TOTAL COLLOCATION NUMBER
!     KMX           :MAXIMUM TOTAL ELEMENTS NUMBER
!     KNA           :MAXIMUM COLLOCATION NUMBER IN ONE DIMENSIONAL DIRECTION
!     KNT           :MAXIMUM COLLOCATION NUMBER IN TIME
!     KNAA          :MAXIMUM COLLOCATION NUMBER IN A ELEMENT
!     KNB           :MAXIMUM COLLOCATION NUMBER AT BOUNDARY OF SPACE
!     KNBT          :MAXIMUM COLLOCATION NUMBER IN INFLOW CONDITION
!     ------------------------------------------------------------------
!
!     ------------  PARAMETER ------------------------------------------
!     NA            :ONE DIMENSIONAL INNER COLLOCATION NUMBER
!     IAN           :THREE DIMENSIONAL TOTAL COLLOCATION NUMBER
!     IA            :ONE DIMENSIONAL TOTAL COLLOCATION NUMBER
!     NX            :TOTAL COLLOCATION NUMBER IN SPACE
!     MX            :ELEMENT NUMBER IN SPACE
!     ID            :TOTAL COLLOCATION NUMBER OF ELEMENT BOUNDARY CONNECTION
!     -----------------------------------------------------------------
!
!     ------------- PARAMETER OF FLUID MECHANICS ----------------------
!     RR            :REYNOLDS NUMBER
!     R1            :RADIAL LENGTH OF PIPE 
!     X1,X2,Y5,R1,Z1,D0,H0  :CHANNEL ELEMENT LENGTH 
!     U0            :MAXIMUM VEROCITY AT INFLOW
!     -----------------------------------------------------------------
!
!     -------------- VARIABLE -----------------------------------------
!     XXC(KNA)      :ORTHOGONAL COLLOCATION VALUE IN ONE DIMENSIONAL X_DIRECTION
!     YYC(KNA)      :ORTHOGONAL COLLOCATION VALUE IN ONE DIMENSIONAL Y_DIRECTION
!     ZZC(KNA)      :ORTHOGONAL COLLOCATION VALUE IN ONE DIMENSIONAL Z_DIRECTION
!     AX(KNAA,KNAA) :DIFFERENTIAL OPERATOR MATRIX OF ONE DIFFERENCIATE BY X
!     AY(KNAA,KNAA) :DIFFERENTIAL OPERATOR MATRIX OF ONE DIFFERENCIATE BY Y
!     AZ(KNAA,KNAA) :DIFFERENTIAL OPERATOR MATRIX OF ONE DIFFERENCIATE BY Z
!     BX(KNAA,KNAA) :DIFFERENTIAL OPERATOR MATRIX OF TWO DIFFERENCIATE BY X
!     BY(KNAA,KNAA) :DIFFERENTIAL OPERATOR MATRIX OF TWO DIFFERENCIATE BY Y
!     BZ(KNAA,KNAA) :DIFFERENTIAL OPERATOR MATRIX OF TWO DIFFERENCIATE BY Z
!     CXY(KNAA,KNAA):DIFFERENTIAL OPERATOR MATRIX OF TWO DIFFERENCIATE BY X,Y
!     CYZ(KNAA,KNAA):DIFFERENTIAL OPERATOR MATRIX OF TWO DIFFERENCIATE BY Y,Z
!     CZX(KNAA,KNAA):DIFFERENTIAL OPERATOR MATRIX OF TWO DIFFERENCIATE BY Z,X
!     IE(KNAA,KMX)  :RELATIONSHIP BETWEEN COLLOCAION AND ELEMENT
!     IDM(KMX)      :DIFFERENT DOMEAIN NUMBER
!     IT(KNX)       :SUPERPOSITION COUNT
!     XX(KNX)       :COORDINATE VALUE(X) AT COLLOCATION
!     YY(KNX)       :COORDINATE VALUE(Y) AT COLLOCATION
!     ZZ(KNX)       :COORDINATE VALUE(Z) AT COLLOCATION
!     U(KNT,KNX) :VALUE OF VEROCITY_U IN TIME & SPACE (YU1M:MODIFIED)
!     V(KNT,KNX) :VALUE OF VEROCITY_V IN TIME & SPACE (YV1M:MODIFIED)
!     W(KNT,KNX) :VALUE OF VEROCITY_V IN TIME & SPACE (YW1M:MODIFIED)
!     P(KNT,KNX) :VALUE OF PRESSURE_P IN TIME & SPACE (YP1M:MODIFIED)
!     YVW1(KNT,KNX) :VALUE OF VORTICITY_VW IN TIME & SPACE
!     YWU1(KNT,KNX) :VALUE OF VORTICITY_WU IN TIME & SPACE
!     YUV1(KNT,KNX) :VALUE OF VORTICITY_UV IN TIME & SPACE 
!     -----------------------------------------------------------------
!
!    -----------------------------------------------------------------
 
      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNT=10,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),TTC(1:KNT),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),U(KNT,KNX),V(KNT,KNX),W(KNT,KNX),P(KNT,KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA),UU0T(KNT,KNX),RR(KNT),UUT1(KNT),VVT1(KNT),WWT1(KNT)
      COMMON/PARAM1/ Re,U0,D0,H0,w8,X1,X2,ANYU,ROU,TT,DT,TC
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NAT,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX,NST,IAT
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),JD1(1:KNX),JD2(1:KNX),NSBU(1:KNB),NSBV(1:KNB),NSBW(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),NSBPO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB),IIT(1:KNX*4)
                     
      REAL(8),ALLOCATABLE :: AX(:,:),AY(:,:),AZ(:,:),BX(:,:),BY(:,:),BZ(:,:),CXY(:,:),CYZ(:,:),CZX(:,:),AT(:,:)
      ALLOCATE (AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA),AT(1:KNT,1:KNT))

      
      CALL INPUTDAT(AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX,AT)
      

      CALL SOLUTE(AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX,AT)

      STOP
      END


!     *****************************************************************


      SUBROUTINE INPUTDAT(AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX,AT)

      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNT=10,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),TTC(1:KNT),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),U(KNT,KNX),V(KNT,KNX),W(KNT,KNX),P(KNT,KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA),UU0T(KNT,KNX),RR(KNT),UUT1(KNT),VVT1(KNT),WWT1(KNT)
      COMMON/PARAM1/ Re,U0,D0,H0,w8,X1,X2,ANYU,ROU,TT,DT,TC
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NAT,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX,NST,IAT
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),JD1(1:KNX),JD2(1:KNX),NSBU(1:KNB),NSBV(1:KNB),NSBW(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),NSBPO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB),IIT(1:KNX*4)

      REAL(8)  AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA),AT(1:KNT,1:KNT)

      character FNAME*100


!     ------------- INPUT OF COLLOCATION VALUE OF SPACE-----------------

      WRITE(*,*)'FILE NAME OF DIFFERENTIAL OPERATOR MATRIX IN SPACE(COMV3D_  .DAT)?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)
        READ(2,*) NAX,NAY,NAZ
        NNAX=NAX
        NNAY=NAY
        NNAZ=NAZ
        IAX=NAX+2
        IAY=NAY+2
        IAZ=NAZ+2
        IAN=(NAX+2)*(NAY+2)*(NAZ+2)
        READ(2,2010)(XXC(I),I=1,IAX)
        READ(2,2010)(YYC(I),I=1,IAY)
        READ(2,2010)(ZZC(I),I=1,IAZ) 
        READ(2,2020)((AX(I,J),AY(I,J),AZ(I,J),BX(I,J),BY(I,J),BZ(I,J),CXY(I,J),CYZ(I,J),CZX(I,J),J=1,IAN),I=1,IAN)
      CLOSE(2)
 2010 FORMAT(F25.0)
 2020 FORMAT(9F25.0)

      WRITE(*,*)'FILE NAME OF DIFFERENTIAL OPERATOR MATRIX IN TIME(COMV1T_  .DAT) ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)
        READ(2,*) NAT
        IAT=NAT+2
        READ(2,2015)(TTC(I),I=1,IAT)
        READ(2,2025)((AT(I,J),temp,J=1,IAT),I=1,IAT)
        IAT=1
        AT(1,1)=0.0D0
      CLOSE(2)
 2015 FORMAT(F25.0)
 2025 FORMAT(2F25.0)

 
      WRITE(*,*)'FILE NAME OF COLLOCATE WEIGHT FOR FLOW RATE AT SECTION OF CHANNEL : ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)
        READ(2,*) NA1,NA2
        DO I=1,NA1
          READ(2,2100) TEMP,W1(I),TEMP
        END DO
        DO I=1,NA2
          READ(2,2100) TEMP,W2(I),TEMP
        END DO
      CLOSE(2)
 2100 FORMAT(3F25.0)

!     ------------ INPUT OF MESH ELEMET & COORDINATE DATA ------------

      WRITE(*,*)'INPUT FILE NAME OF NODE_ELEMENT DATA : ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)
        READ(2,1000) NX,MX
        DO I=1,MX
          READ(2,1010)(IE(J,I),J=1,IAN)
!             write(*,1010)(IE(J,I),J=1,IAN) 
        END DO
          
!          write(*,*) ie(14,1),ie(24,2)
        READ(2,1020)(IJ1,XX(I),YY(I),ZZ(I),I=1,NX)
      CLOSE(2)
 1000 FORMAT(2I10)
 1010 FORMAT(7I10)
 1020 FORMAT(I10,3F25.0)

!     ------------- INPUT OF BOUNDARY CONNECTION --------------------

      WRITE(*,*)'INPUT FILE NAME OF BOUNDARY CONNECTION : ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)
        READ(2,1030) ID
        READ(2,1040)(JD1(K),JD2(K),XL(JD1(K),JD2(K)),XM(JD1(K),JD2(K)),XN(JD1(K),JD2(K)),K=1,ID)
      CLOSE(2)
 1030 FORMAT(I10)
 1040 FORMAT(2I10,3F25.0)


!     ------------- INPUT FILE NAME OF BOUNDARY CONDITION -----------

      WRITE(*,*)'INPUT FILE NAME OF BOUNDARY CONDITION  U, V, P :?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)
     ! boudary condition u,v,w at wall
        READ(2,*) NSB1
        IF(NSB1.EQ.0) GO TO 10
          READ(2,3010)(NSBU(I),JSBU(I),ISBU(I),SBU(I),SBV(I),SBW(I),I=1,NSB1)
          do i=1,nsb1
 !           write(*,*)'nsbu(i)=',nsbu(i)
          end do

     ! boundary condition u,v,w at inflow
   10   READ(2,*) NSB2
        IF(NSB2.EQ.0) GO TO 20
          READ(2,3010)(NSBUI(I),JSBUI(I),ISBUI(I),SBUI(I),SBVI(I),SBWI(I),I=1,NSB2)

     ! boundary condition p at outflow
   20   READ(2,*) NSB3
        IF(NSB3.EQ.0) GO TO 25
          READ(2,3020)(NSBUO(I),JSBUO(I),ISBUO(I),SBPO(I),I=1,NSB3)

   25 CONTINUE

 3000 FORMAT(I10)
 3010 FORMAT(3I10,3F25.0)
 3020 FORMAT(3I10,F25.0)

      CLOSE(2)
   
 


!     ------------- INPUT OF PARAMETER -----------------------


      WRITE(*,*) 'ANYU,ROU= ?'
      READ(*,*) ANYU,ROU 
      
      WRITE(*,*)'Velocity U0,D0,H0,X1,X2 = ?'
      READ(*,*) U0,D0,H0,X1,X2
      Re=U0*D0/ANYU
      
      WRITE(*,*)'Time Trial number TT,DT,TC,=?'
      READ(*,*) TT,DT,TC
      NST=TT/DT
      
 


      RETURN
      END


!     *************************************************************

      SUBROUTINE TRANSF(JF,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
      
      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNT=10,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),TTC(1:KNT),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),U(KNT,KNX),V(KNT,KNX),W(KNT,KNX),P(KNT,KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA),UU0T(KNT,KNX),RR(KNT),UUT1(KNT),VVT1(KNT),WWT1(KNT)
      COMMON/PARAM1/ Re,U0,D0,H0,w8,X1,X2,ANYU,ROU,TT,DT,TC
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NAT,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX,NST,IAT
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),JD1(1:KNX),JD2(1:KNX),NSBU(1:KNB),NSBV(1:KNB),NSBW(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),NSBPO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB),IIT(1:KNX*4)
                     
      REAL(8)  AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA),AT(1:KNT,1:KNT),&
                 work(81)



!     --------------- COORDINATE TRANSFORMATION ----------------------

      PAI=3.14159265358979323846D0

     DO IIu=1,IAN
          DO iu1=1,9
          DO ju1=1,9
            DJJ(iu1,ju1)=0.0D0
            DJJA(iu1,ju1,IIu)=0.0D0
          END DO
          END DO
     END DO

      DO IIu=1,IAN
        DXDG=0.0d0
        DYDG=0.0d0
        DZDG=0.0D0
        DXDE=0.0d0
        DYDE=0.0d0
        DZDE=0.0D0
        DXDT=0.0d0
        DYDT=0.0d0
        DZDT=0.0D0
        D2XDG2=0.0d0
        D2YDG2=0.0d0
        D2ZDG2=0.0D0
        D2XDE2=0.0d0
        D2YDE2=0.0d0
        D2ZDE2=0.0D0
        D2XDT2=0.0d0
        D2YDT2=0.0d0
        D2ZDT2=0.0D0
        D2XDGDE=0.0d0
        D2YDGDE=0.0d0
        D2ZDGDE=0.0D0
        D2XDEDT=0.0d0
        D2YDEDT=0.0d0
        D2ZDEDT=0.0D0
        D2XDTDG=0.0d0
        D2YDTDG=0.0d0
        D2ZDTDG=0.0D0
        DO JJu=1,IAN
          Ju=IE(JJu,JF)
          DXDG=DXDG+AX(IIu,JJu)*XX(Ju)
          DYDG=DYDG+AX(IIu,JJu)*YY(Ju)
          DZDG=DZDG+AX(IIu,JJu)*ZZ(Ju)
          DXDE=DXDE+AY(IIu,JJu)*XX(Ju)
          DYDE=DYDE+AY(IIu,JJu)*YY(Ju)
          DZDE=DZDE+AY(IIu,JJu)*ZZ(Ju)
          DXDT=DXDT+AZ(IIu,JJu)*XX(Ju)
          DYDT=DYDT+AZ(IIu,JJu)*YY(Ju)
          DZDT=DZDT+AZ(IIu,JJu)*ZZ(Ju)
          D2XDG2=D2XDG2+BX(IIu,JJu)*XX(Ju)
          D2YDG2=D2YDG2+BX(IIu,JJu)*YY(Ju)
          D2ZDG2=D2ZDG2+BX(IIu,JJu)*ZZ(Ju)
          D2XDE2=D2XDE2+BY(IIu,JJu)*XX(Ju)
          D2YDE2=D2YDE2+BY(IIu,JJu)*YY(Ju)
          D2ZDE2=D2ZDE2+BY(IIu,JJu)*ZZ(Ju)
          D2XDT2=D2XDT2+BZ(IIu,JJu)*XX(Ju)
          D2YDT2=D2YDT2+BZ(IIu,JJu)*YY(Ju)
          D2ZDT2=D2ZDT2+BZ(IIu,JJu)*ZZ(Ju)
          D2XDGDE=D2XDGDE+CXY(IIu,JJu)*XX(Ju)
          D2YDGDE=D2YDGDE+CXY(IIu,JJu)*YY(Ju)
          D2ZDGDE=D2ZDGDE+CXY(IIu,JJu)*ZZ(Ju)
          D2XDEDT=D2XDEDT+CYZ(IIu,JJu)*XX(Ju)
          D2YDEDT=D2YDEDT+CYZ(IIu,JJu)*YY(Ju)
          D2ZDEDT=D2ZDEDT+CYZ(IIu,JJu)*ZZ(Ju)
          D2XDTDG=D2XDTDG+CZX(IIu,JJu)*XX(Ju)
          D2YDTDG=D2YDTDG+CZX(IIu,JJu)*YY(Ju)
          D2ZDTDG=D2ZDTDG+CZX(IIu,JJu)*ZZ(Ju)
        END DO
          
          
          DJJA(1,1,IIu)=DXDG
          DJJA(1,2,IIu)=DYDG
          DJJA(1,3,IIu)=DZDG
          DJJA(2,1,IIu)=DXDE
          DJJA(2,2,IIu)=DYDE
          DJJA(2,3,IIu)=DZDE
          DJJA(3,1,IIu)=DXDT
          DJJA(3,2,IIu)=DYDT
          DJJA(3,3,IIu)=DZDT
          DJJA(4,1,IIu)=D2XDG2
          DJJA(4,2,IIu)=D2YDG2
          DJJA(4,3,IIu)=D2ZDG2
          DJJA(4,4,IIu)=DXDG**2
          DJJA(4,5,IIu)=DYDG**2
          DJJA(4,6,IIu)=DZDG**2
          DJJA(4,7,IIu)=2*DXDG*DYDG
          DJJA(4,8,IIu)=2*DYDG*DZDG
          DJJA(4,9,IIu)=2*DZDG*DXDG
          DJJA(5,1,IIu)=D2XDE2
          DJJA(5,2,IIu)=D2YDE2
          DJJA(5,3,IIu)=D2ZDE2
          DJJA(5,4,IIu)=DXDE**2
          DJJA(5,5,IIu)=DYDE**2
          DJJA(5,6,IIu)=DZDE**2
          DJJA(5,7,IIu)=2*DXDE*DYDE
          DJJA(5,8,IIu)=2*DYDE*DZDE
          DJJA(5,9,IIu)=2*DZDE*DXDE
          DJJA(6,1,IIu)=D2XDT2
          DJJA(6,2,IIu)=D2YDT2
          DJJA(6,3,IIu)=D2ZDT2
          DJJA(6,4,IIu)=DXDT**2
          DJJA(6,5,IIu)=DYDT**2
          DJJA(6,6,IIu)=DZDT**2
          DJJA(6,7,IIu)=2*DXDT*DYDT
          DJJA(6,8,IIu)=2*DYDT*DZDT
          DJJA(6,9,IIu)=2*DZDT*DXDT
          DJJA(7,1,IIu)=D2XDGDE
          DJJA(7,2,IIu)=D2YDGDE
          DJJA(7,3,IIu)=D2ZDGDE          
          DJJA(7,4,IIu)=DXDG*DXDE
          DJJA(7,5,IIu)=DYDG*DYDE
          DJJA(7,6,IIu)=DZDG*DZDE
          DJJA(7,7,IIu)=DXDG*DYDE+DXDE*DYDG
          DJJA(7,8,IIu)=DYDG*DZDE+DYDE*DZDG
          DJJA(7,9,IIu)=DZDG*DXDE+DZDE*DXDG
          DJJA(8,1,IIu)=D2XDEDT
          DJJA(8,2,IIu)=D2YDEDT
          DJJA(8,3,IIu)=D2ZDEDT          
          DJJA(8,4,IIu)=DXDE*DXDT
          DJJA(8,5,IIu)=DYDE*DYDT
          DJJA(8,6,IIu)=DZDE*DZDT
          DJJA(8,7,IIu)=DXDE*DYDT+DXDT*DYDE
          DJJA(8,8,IIu)=DYDE*DZDT+DYDT*DZDE
          DJJA(8,9,IIu)=DZDE*DXDT+DZDT*DXDE
          DJJA(9,1,IIu)=D2XDTDG
          DJJA(9,2,IIu)=D2YDTDG
          DJJA(9,3,IIu)=D2ZDTDG          
          DJJA(9,4,IIu)=DXDT*DXDG
          DJJA(9,5,IIu)=DYDT*DYDG
          DJJA(9,6,IIu)=DZDT*DZDG
          DJJA(9,7,IIu)=DXDT*DYDG+DXDG*DYDT
          DJJA(9,8,IIu)=DYDT*DZDG+DYDG*DZDT
          DJJA(9,9,IIu)=DZDT*DXDG+DZDG*DXDT     
                    
        DO iu1=1,9
        DO ju1=1,9
          DJJ(iu1,ju1)=DJJA(iu1,ju1,IIu)
!          write(*,*)'clear ohkubo Djj',Jf,Iu1,Ju1,djja(iu1,ju1,iiu)
        END DO
        END DO
        
        call dgetrf( 9, 9, DJJ, 9, ipiv, info )
!        WRITE(*,*)'info=',info
        
        lwork=81
        call dgetri( 9, DJJ, 9, ipiv, work, lwork, info )
!        write(*,*)'info=',info

!        write(*,*)'jf=',jf
!
        DO iu1=1,9
        DO ju1=1,9
          DJJA(iu1,ju1,IIu)=DJJ(iu1,ju1)
          IF(JF==4) THEN
!            write(*,*)'clear ohkubo Djj',Jf,Iu1,Ju1,iiu,djja(iu1,ju1,iiu)
          END IF
        END DO
        END DO
!
!
!
      END DO
      
!      STOP
!
!
      RETURN
      END


!     *****************************************************************

      SUBROUTINE SOLUTE(AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX,AT)

      use sparse_triplet_builder, only: triplet_matrix, triplet_init, triplet_add, triplet_clear, triplet_to_csr
      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNT=10,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),TTC(1:KNT),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),U(KNT,KNX),V(KNT,KNX),W(KNT,KNX),P(KNT,KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA),UU0T(KNT,KNX),RR(KNT),UUT1(KNT),VVT1(KNT),WWT1(KNT)
      COMMON/PARAM1/ Re,U0,D0,H0,w8,X1,X2,ANYU,ROU,TT,DT,TC
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NAT,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX,NST,IAT
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),JD1(1:KNX),JD2(1:KNX),NSBU(1:KNB),NSBV(1:KNB),NSBW(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),NSBPO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB),IIT(1:KNX*4)
      
      REAL(8) AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA),AT(1:KNT,1:KNT)


      dimension abcxy(9),dj1(9),dj2(9),dj3(9),dj4(9),dj5(9),dj6(9),dj7(9),dj8(9),dj9(9),taxuvp1(knaa),tayuvp1(knaa),tazuvp1(knaa),&
                uu1(knaa),vv1(knaa),ww1(knaa),pp1(knaa),&
                tbxuvp1(knaa),tbyuvp1(knaa),tbzuvp1(knaa),tcxyuvp1(knaa),tcyzuvp1(knaa),tczxuvp1(knaa),tatuvp1(KNT),&
                divuvp1(knaa),ozwork(6,3)
      REAL(8) lapuvp1(knaa)



      character FNAME*100,FNAME_RESULT*100,trans*2

      INTEGER :: iparm(64)
      INTEGER(8) :: ipt(64)
      INTEGER GLOBAL_EQUATION_BLOCK,GLOBAL_DOF_INDEX,GLOBAL_EQUATION_INDEX
      INTEGER, PARAMETER :: BLOCK_U=1,BLOCK_V=2,BLOCK_W=3,BLOCK_P=4
      INTEGER,ALLOCATABLE :: IPIV(:),iperm(:),ir(:),ic(:),IS_CONSTRAINED(:),IS_CONSTRAINED_ROW(:),BC_ROW(:),IS_REPLACED_ROW(:),ISINGULAR_PIN(:)
      REAL(8),ALLOCATABLE :: FFS(:),FFS_PHYS(:),FFS_PHYS_BASE(:),FFS_LINEAR(:),FFS_REPLACED_ROW(:),YUVP1(:),x(:,:),a(:),FFUVP(:),FFUVP_LIN(:),FFUVP2(:),YUVWP(:),DDUVWP(:),DDUVWP_SD(:),FF1(:),BC_DELTA(:),ROW_SCALE(:),ROW_SCALE_BASE(:),COL_SCALE(:)
      REAL(8),ALLOCATABLE :: NORMAL_RHS(:),NORMAL_DIAG(:),SALS_RHS(:),SALS_DIAG(:)
      INTEGER, PARAMETER :: KOZ_WORK_THREADS=8
      REAL(8),ALLOCATABLE :: CG_X(:),CG_R(:),CG_Z(:),CG_P(:),CG_AP(:),CG_TMP(:),CG_DIAG(:),CG_WORK(:,:),CG_WORK_T(:,:,:)
      REAL(8),ALLOCATABLE :: FFS_FD_BASE(:),JCOL_FD_BASE(:),JCOL_FD_DIFF(:)
      REAL(8),ALLOCATABLE :: U_SAVE(:,:),V_SAVE(:,:),W_SAVE(:,:),P_SAVE(:,:)
      REAL(8),ALLOCATABLE :: U_LSCAN_BEST(:,:),V_LSCAN_BEST(:,:),W_LSCAN_BEST(:,:),P_LSCAN_BEST(:,:),DDUVWP_LSCAN_BEST(:)
      REAL(8),ALLOCATABLE :: BBUVP(:,:),BBUVP_LIN(:,:)
      REAL(8),ALLOCATABLE :: SSMMK(:,:),DENSE_A(:,:),DENSE_B(:,:)
      type(triplet_matrix) :: SPARSE_A, SPARSE_ELEMENT_A, SPARSE_LINEAR_A, SPARSE_PRESSURE_WALL_A, SPARSE_NORMAL_A
      N=4*NX*IAT
      IF(N.GT.20000) THEN
        WRITE(*,*)'DIRECT CSR NEWTON MATRIX ENABLED. DENSE DEBUG MATRIX SKIPPED.'
        WRITE(*,*)'N=',N,' DENSE MEMORY BYTES WOULD BE=',8.0D0*DBLE(N)*DBLE(N)
      END IF
      ALLOCATE (BBUVP(4*IAN*IAT,4*IAN*IAT),BBUVP_LIN(4*IAN*IAT,4*IAN*IAT))
      ALLOCATE (IPIV(N),iperm(N),FFS(N),FFS_PHYS(N),FFS_PHYS_BASE(N),FFS_LINEAR(N),YUVP1(N),x(N,1))
      ALLOCATE (FFUVP(4*IAN*IAT),FFUVP_LIN(4*IAN*IAT),FFUVP2(4*IAN*IAT),YUVWP(N),DDUVWP(N),DDUVWP_SD(N),FF1(N))
      ALLOCATE (U_SAVE(IAT,NX),V_SAVE(IAT,NX),W_SAVE(IAT,NX),P_SAVE(IAT,NX))
      ALLOCATE (U_LSCAN_BEST(IAT,NX),V_LSCAN_BEST(IAT,NX),W_LSCAN_BEST(IAT,NX),P_LSCAN_BEST(IAT,NX),DDUVWP_LSCAN_BEST(N))
      ALLOCATE (IS_CONSTRAINED(N),IS_CONSTRAINED_ROW(N),BC_ROW(N),IS_REPLACED_ROW(N),BC_DELTA(N),FFS_REPLACED_ROW(N))
      ALLOCATE (ROW_SCALE(N),ROW_SCALE_BASE(N),COL_SCALE(N))
      ALLOCATE (NORMAL_RHS(N),NORMAL_DIAG(N),SALS_RHS(N),SALS_DIAG(N))
      ALLOCATE (CG_X(N),CG_R(N),CG_Z(N),CG_P(N),CG_AP(N),CG_TMP(N),CG_DIAG(N),CG_WORK(6,N),CG_WORK_T(6,N,KOZ_WORK_THREADS))
      ALLOCATE (FFS_FD_BASE(N),JCOL_FD_BASE(N),JCOL_FD_DIFF(N))
      IF(N.LE.5000) ALLOCATE (ISINGULAR_PIN(N))
      IDENSE_DEBUG=0
      IF(IDENSE_DEBUG.NE.0.AND.N.LE.5000) ALLOCATE (SSMMK(N,N),DENSE_A(N,N),DENSE_B(N,1))
      PPE_CONV_SCALE=0.0D0
      PPE_JAC_SCALE=PPE_CONV_SCALE
      PPE_RHS_CONV_SCALE=PPE_CONV_SCALE
      IAPPLY_DGESV_PIN=0
      UPDATE_SIGN=1.0D0
      STEP_SIGN=UPDATE_SIGN
      DSD_PHYS_SIGN=UPDATE_SIGN
      IOPPOSITE_TRIED=0
      IPENDING_STEP=0
      ILINESEARCH_TRIES=0
      ILINESEARCH_ACCEPTED=0
      ILINESEARCH_EXTRA_MODE=0
      ILAMBDA_SCAN_ACTIVE=0
      ILAMBDA_SCAN_STAGE=0
      ILAMBDA_SCAN_COUNT=0
      IGAUSS_NEWTON_ONLY=1
      ILAMBDA_SCAN_MAX=3
      WW4_LAMBDA_BEST=1.0D300
      DLM_LAMBDA_SCAN_BEST=0.0D0
      DLM_LAMBDA_SCAN_CAND1=0.0D0
      DLM_LAMBDA_SCAN_CAND2=0.0D0
      DLM_LAMBDA_SCAN_CAND3=0.0D0
      DLM_LAMBDA_CURRENT=0.0D0
      DLM_LAMBDA_USED=0.0D0
      DLM_LAMBDA_MIN=0.0D0
      DLM_LAMBDA_MAX=1.0D8
      DLAMBDA_RATIO=0.0D0
      ISTAG_LAMBDA_COUNT=0
      ISTAG_LAMBDA_HOLD=0
      ISTAG_ACCEPT_STEP=-1
      ITRUST_REJECT_COUNT=0
!     Switch to the Gauss-Newton direction once Marquardt damping is small.
!     DGN_LAMBDA_SWITCH controls the threshold; the default is 1.0D-2.
      IGAUSS_NEWTON_TRIAL=1
      IGN_ACTIVE=0
      IGN_DECREASE_COUNT=0
      INEWTON_MAX=100
      DGN_LAMBDA_SWITCH=1.0D-2
      DGN_WW4_SWITCH=1.0D300
      DRELAX_START_SCALE=1.0D0
      IMATRIX_LINEAR_RESIDUAL=1
      IMATRIX_AXB_FINAL_RESIDUAL=0
      IDIAG_LOCAL_MOM=0
      IDIAG_TRIAL_RESIDUAL=0
      IPICARD_CONV_JAC=1
      DCONV_FACTOR=0.001D0
      WRITE(*,*)'CONVECTION CONTINUATION FACTOR=',DCONV_FACTOR
      WRITE(*,*)'PICARD CONVECTION MATRIX MODE=',IPICARD_CONV_JAC
      IPRESSURE_FD_CHECK=0
      IPRESSURE_FD_MODE=0
      IPRESSURE_FD_DONE=0
      IPRESSURE_FD_COL=0
      IPRESSURE_FD_NODE=0
      IPRESSURE_FD_TIME=0
      IPRESSURE_FD_TARGET_BLOCK=1
      IAPPLY_INLET_VELOCITY_BC=1
!     Primary solver path for the refactored steadyized problem.
!       0: disabled, 1: enabled
!     Keep Newton/PARDISO available as a diagnostic path, but use
!     Marquardt_CG as the standard nonlinear correction.
      INEWTON_PARDISO=0
      IMARQUARDT_CG=1
      IMARQUARDT_PCG=1
      IMARQUARDT_SD_FALLBACK=1
      IMARQUARDT_UNSCALED_CG=1
      IDIRECT_BICGSTAB=0
      IF(IMARQUARDT_CG.NE.0) THEN
        DLM_LAMBDA_MIN=0.0D0
        DLM_LAMBDA_MAX=1.0D4
        DLM_LAMBDA_CURRENT=1.0D-4
        IF(IMARQUARDT_UNSCALED_CG.NE.0) DLM_LAMBDA_CURRENT=1.0D0
        IF(Re.GE.1.0D5) THEN
          DLM_LAMBDA_MIN=1.0D-4
          DLM_LAMBDA_CURRENT=DMAX1(DLM_LAMBDA_CURRENT,DLM_LAMBDA_MIN)
          WRITE(*,*)'HIGH_RE MARQUARDT CONTROL GN_ON LAMBDA_MIN=',Re,DLM_LAMBDA_MIN
        END IF
      END IF
      IF(IGAUSS_NEWTON_ONLY.NE.0) THEN
        DLM_LAMBDA_MIN=0.0D0
        DLM_LAMBDA_MAX=0.0D0
        DLM_LAMBDA_CURRENT=0.0D0
        DLM_LAMBDA_USED=0.0D0
        IGAUSS_NEWTON_TRIAL=0
        ILAMBDA_SCAN_ACTIVE=0
        WRITE(*,*)'GAUSS_NEWTON_ONLY MODE: LAMBDA FIXED TO ZERO'
      END IF
      IF(IDIRECT_BICGSTAB.EQ.2.OR.IDIRECT_BICGSTAB.EQ.3) DLM_LAMBDA_CURRENT=1.0D0
      IPROJECT_CONTINUITY=0
      IPRESSURE_ONLY_CORRECTION=0
      DLM_SD_SWITCH=1.0D-2
      DPRESSURE_FD_EPS=1.0D-7
      DPRESSURE_FD_OLD=0.0D0
      FFS_FD_BASE=0.0D0
      JCOL_FD_BASE=0.0D0
      JCOL_FD_DIFF=0.0D0
      call triplet_init(SPARSE_ELEMENT_A, N, N, max(1, min(1000000, 20*N)))
      call triplet_init(SPARSE_LINEAR_A, N, N, max(1, min(1000000, 4*N)))
      call triplet_init(SPARSE_PRESSURE_WALL_A, N, N, max(1, min(1000000, 2*N)))
      
 
!     Do not force OpenMP to one thread; use OMP_NUM_THREADS or runtime default.

!     ---------------  JACOBIAN MATRIX FOR NEWTON RAPHSON METHOD ------
      PAI=3.1415926535897932384626433832795028841971D0

      WRITE(*,*)'OUTPUT FILE NAME OF TRANSIENT RESULT : ?'
      READ(*,*) FNAME
      FNAME_RESULT=FNAME
      OPEN(2,FILE=FNAME_RESULT,STATUS='REPLACE',ACTION='WRITE')

            T=0.0D0
 
    write(*,*) 'MX,NX',mx,nx
!       OMP disabled here: legacy assembly loops are not race-audited.
!       !$OMP parallel private(I,J)
!       !$OMP do
        DO J=1,IAT
        DO I=1,NX
          U(J,I)=0.0d0
          V(J,I)=0.0d0
          W(J,I)=0.0d0
          P(J,I)=0.0d0
    !      write(*,*)'clear 005! nx,iat',i,j
        END DO
        END DO
!       !$OMP end do
!       !$OMP end parallel
        write(*,*)'clear 006!'

        INX=IAN*IAT
        
      IPK1=0
      T=0.0D0
      DO ITT1=1,NST


      write(*,*)'clear 007!'

      WRITE(*,*)'INPUT FILE NAME OF STEADY STOKES INITIAL RESULT : ?'
      READ(*,*) FNAME
      OPEN(77,FILE=FNAME,STATUS='OLD')
      DO IB=1,NX
        READ(77,2100) XTMP,YTMP,ZTMP,UTMP,VTMP,WTMP,PTMP,RTMP
        U(1,IB)=UTMP
        V(1,IB)=VTMP
        W(1,IB)=WTMP
        P(1,IB)=PTMP
      END DO
      CLOSE(77)
      
      write(*,*)'clear 89 !'
      
      IPK=0

30    CONTINUE      
     
      DO I=1,NX*IAT*4
        IIT(I)=0
      END DO
      IF(ALLOCATED(SSMMK)) SSMMK=0.0D0
      call triplet_clear(SPARSE_ELEMENT_A)
      call triplet_clear(SPARSE_LINEAR_A)
      call triplet_clear(SPARSE_PRESSURE_WALL_A)
      IS_REPLACED_ROW=0
      IS_CONSTRAINED_ROW=0
      FFS_REPLACED_ROW=0.0D0
      FFS_LINEAR=0.0D0


!     OMP disabled here: keep residual initialization in the legacy order.
!     !$OMP parallel private(I)
!     !$OMP do
      DO I=1,NX*IAT*4
        FFS(I)=0.0D0
      END DO
!     !$OMP end do
!     !$OMP end parallel


      INX=IAN
      write(*,*)'clear 009!'
      


      DO IM=1,MX   
        DO I=1,4*IAN*IAT
        DO J=1,4*IAN*IAT
          BBUVP(I,J)=0.0D0
          BBUVP_LIN(I,J)=0.0D0
        END DO
        FFUVP(I)=0.0D0
        FFUVP_LIN(I)=0.0D0
        FFUVP2(I)=0.0D0
        END DO
     
        
        CALL TRANSF(IM,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)

        
!       OMP disabled here: element assembly writes shared work arrays.
!       !$OMP parallel private(I1,II,IQ1,IQ2,ITT,IBB,ISP,J1,JJ,JTT,JBB,OZWORK)
!       !$OMP do
        
        DO I1=1,IAN
        
          II=IE(I1,IM)
          DO IQ1=1,9
          DO IQ2=1,9
            DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,I1)
          END DO
          END DO
        
        DO ITT=1,IAT

          IIB=ITT+(I1-1)*IAT
          ISP=ITT+(II-1)*IAT

          IINTERIOR=0
          IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ)+NAX*NAY*NAZ)) THEN
            IINTERIOR=1
          END IF
          IF(IINTERIOR.EQ.1) THEN
            DO IJ=1,IAN
              TAXUVP1(IJ)=DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ)
              TAYUVP1(IJ)=DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ)
              TAZUVP1(IJ)=DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ)
              TBXUVP1(IJ)=DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)
              TBYUVP1(IJ)=DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)
              TBZUVP1(IJ)=DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ)
              TCXYUVP1(IJ)=DJJ(7,1)*AX(I1,IJ)+DJJ(7,2)*AY(I1,IJ)+DJJ(7,3)*AZ(I1,IJ)+DJJ(7,4)*BX(I1,IJ)+DJJ(7,5)*BY(I1,IJ)+DJJ(7,6)*BZ(I1,IJ)+DJJ(7,7)*CXY(I1,IJ)+DJJ(7,8)*CYZ(I1,IJ)+DJJ(7,9)*CZX(I1,IJ)
              TCYZUVP1(IJ)=DJJ(8,1)*AX(I1,IJ)+DJJ(8,2)*AY(I1,IJ)+DJJ(8,3)*AZ(I1,IJ)+DJJ(8,4)*BX(I1,IJ)+DJJ(8,5)*BY(I1,IJ)+DJJ(8,6)*BZ(I1,IJ)+DJJ(8,7)*CXY(I1,IJ)+DJJ(8,8)*CYZ(I1,IJ)+DJJ(8,9)*CZX(I1,IJ)
              TCZXUVP1(IJ)=DJJ(9,1)*AX(I1,IJ)+DJJ(9,2)*AY(I1,IJ)+DJJ(9,3)*AZ(I1,IJ)+DJJ(9,4)*BX(I1,IJ)+DJJ(9,5)*BY(I1,IJ)+DJJ(9,6)*BZ(I1,IJ)+DJJ(9,7)*CXY(I1,IJ)+DJJ(9,8)*CYZ(I1,IJ)+DJJ(9,9)*CZX(I1,IJ)
              UU1(IJ)=U(ITT,IE(IJ,IM))
              VV1(IJ)=V(ITT,IE(IJ,IM))
              WW1(IJ)=W(ITT,IE(IJ,IM))
            END DO
            TAXU=OZ_DDOT(IAN, TAXUVP1, 1, UU1, 1)
            TAYU=OZ_DDOT(IAN, TAYUVP1, 1, UU1, 1)
            TAZU=OZ_DDOT(IAN, TAZUVP1, 1, UU1, 1)
            TAXV=OZ_DDOT(IAN, TAXUVP1, 1, VV1, 1)
            TAYV=OZ_DDOT(IAN, TAYUVP1, 1, VV1, 1)
            TAZV=OZ_DDOT(IAN, TAZUVP1, 1, VV1, 1)
            TAXW=OZ_DDOT(IAN, TAXUVP1, 1, WW1, 1)
            TAYW=OZ_DDOT(IAN, TAYUVP1, 1, WW1, 1)
            TAZW=OZ_DDOT(IAN, TAZUVP1, 1, WW1, 1)
            TBXU=OZ_DDOT(IAN, TBXUVP1, 1, UU1, 1)
            TBYV=OZ_DDOT(IAN, TBYUVP1, 1, VV1, 1)
            TBZW=OZ_DDOT(IAN, TBZUVP1, 1, WW1, 1)
            TCXYU=OZ_DDOT(IAN, TCXYUVP1, 1, UU1, 1)
            TCZXU=OZ_DDOT(IAN, TCZXUVP1, 1, UU1, 1)
            TCYZV=OZ_DDOT(IAN, TCYZUVP1, 1, VV1, 1)
            TCXYV=OZ_DDOT(IAN, TCXYUVP1, 1, VV1, 1)
            TCYZW=OZ_DDOT(IAN, TCYZUVP1, 1, WW1, 1)
            TCZXW=OZ_DDOT(IAN, TCZXUVP1, 1, WW1, 1)
          END IF

          
          DO J1=1,IAN
          DO JTT=1,IAT
          
            JJ=IE(J1,IM)

            JJB=JTT+(J1-1)*IAT
          
            IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ).OR.(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.&
                I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
              BBUVP(IIB,JJB)=(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                  DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))*XL(I1,IM)+(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+&
                  DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))*XM(I1,IM)&
                  +(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                  DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))*XN(I1,IM)
              BBUVP(IIB+IAN*IAT,JJB+IAN*IAT)=(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                  DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))*XL(I1,IM)+(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+&
                  DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))*XM(I1,IM)&
                  +(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                  DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))*XN(I1,IM)
              BBUVP(IIB+2*IAN*IAT,JJB+2*IAN*IAT)=(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                  DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))*XL(I1,IM)+(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+&
                  DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))*XM(I1,IM)&
                  +(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                  DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))*XN(I1,IM)
              BBUVP(IIB+3*IAN*IAT,JJB+3*IAN*IAT)=(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                  DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))*XL(I1,IM)+(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+&
                  DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))*XM(I1,IM)&
                  +(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                  DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))*XN(I1,IM)
            END IF
            IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ)+NAX*NAY*NAZ)) THEN
                IF(ITT.EQ.JTT) THEN
                IF(.FALSE.) THEN
                DO IJ=1,IAN
                TAXUVP1(IJ)=DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ)
                TAYUVP1(IJ)=DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ)
                TAZUVP1(IJ)=DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ)
                TBXUVP1(IJ)=DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)
                TBYUVP1(IJ)=DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)
                TBZUVP1(IJ)=DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ)
                TCXYUVP1(IJ)=DJJ(7,1)*AX(I1,IJ)+DJJ(7,2)*AY(I1,IJ)+DJJ(7,3)*AZ(I1,IJ)+DJJ(7,4)*BX(I1,IJ)+DJJ(7,5)*BY(I1,IJ)+DJJ(7,6)*BZ(I1,IJ)+DJJ(7,7)*CXY(I1,IJ)+DJJ(7,8)*CYZ(I1,IJ)+DJJ(7,9)*CZX(I1,IJ)
                TCYZUVP1(IJ)=DJJ(8,1)*AX(I1,IJ)+DJJ(8,2)*AY(I1,IJ)+DJJ(8,3)*AZ(I1,IJ)+DJJ(8,4)*BX(I1,IJ)+DJJ(8,5)*BY(I1,IJ)+DJJ(8,6)*BZ(I1,IJ)+DJJ(8,7)*CXY(I1,IJ)+DJJ(8,8)*CYZ(I1,IJ)+DJJ(8,9)*CZX(I1,IJ)
                TCZXUVP1(IJ)=DJJ(9,1)*AX(I1,IJ)+DJJ(9,2)*AY(I1,IJ)+DJJ(9,3)*AZ(I1,IJ)+DJJ(9,4)*BX(I1,IJ)+DJJ(9,5)*BY(I1,IJ)+DJJ(9,6)*BZ(I1,IJ)+DJJ(9,7)*CXY(I1,IJ)+DJJ(9,8)*CYZ(I1,IJ)+DJJ(9,9)*CZX(I1,IJ)
                UU1(IJ)=U(ITT,IE(IJ,IM))
                VV1(IJ)=V(ITT,IE(IJ,IM))
                WW1(IJ)=W(ITT,IE(IJ,IM))
                END DO
                TAXU=OZ_DDOT(IAN, TAXUVP1, 1, UU1, 1)
                TAYU=OZ_DDOT(IAN, TAYUVP1, 1, UU1, 1)
                TAZU=OZ_DDOT(IAN, TAZUVP1, 1, UU1, 1)
                TAXV=OZ_DDOT(IAN, TAXUVP1, 1, VV1, 1)
                TAYV=OZ_DDOT(IAN, TAYUVP1, 1, VV1, 1)
                TAZV=OZ_DDOT(IAN, TAZUVP1, 1, VV1, 1)
                TAXW=OZ_DDOT(IAN, TAXUVP1, 1, WW1, 1)
                TAYW=OZ_DDOT(IAN, TAYUVP1, 1, WW1, 1)
                TAZW=OZ_DDOT(IAN, TAZUVP1, 1, WW1, 1)
                END IF
                IF(I1.EQ.J1) THEN
                BBUVP(IIB,JJB)=TAXUVP1(J1)
                BBUVP(IIB,JJB+IAN*IAT)=TAYUVP1(J1)
                BBUVP(IIB,JJB+2*IAN*IAT)=TAZUVP1(J1)
                BBUVP(IIB,JJB+3*IAN*IAT)=0.0D0
                
                IF(.FALSE.) THEN
                DO IJ=1,IAN
                TAXUVP1(IJ)=DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ)
                TAYUVP1(IJ)=DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ)
                TAZUVP1(IJ)=DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ)
                TBXUVP1(IJ)=DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)
                TBYUVP1(IJ)=DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)
                TBZUVP1(IJ)=DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ)
                TCXYUVP1(IJ)=DJJ(7,1)*AX(I1,IJ)+DJJ(7,2)*AY(I1,IJ)+DJJ(7,3)*AZ(I1,IJ)+DJJ(7,4)*BX(I1,IJ)+DJJ(7,5)*BY(I1,IJ)+DJJ(7,6)*BZ(I1,IJ)+DJJ(7,7)*CXY(I1,IJ)+DJJ(7,8)*CYZ(I1,IJ)+DJJ(7,9)*CZX(I1,IJ)
                TCYZUVP1(IJ)=DJJ(8,1)*AX(I1,IJ)+DJJ(8,2)*AY(I1,IJ)+DJJ(8,3)*AZ(I1,IJ)+DJJ(8,4)*BX(I1,IJ)+DJJ(8,5)*BY(I1,IJ)+DJJ(8,6)*BZ(I1,IJ)+DJJ(8,7)*CXY(I1,IJ)+DJJ(8,8)*CYZ(I1,IJ)+DJJ(8,9)*CZX(I1,IJ)
                TCZXUVP1(IJ)=DJJ(9,1)*AX(I1,IJ)+DJJ(9,2)*AY(I1,IJ)+DJJ(9,3)*AZ(I1,IJ)+DJJ(9,4)*BX(I1,IJ)+DJJ(9,5)*BY(I1,IJ)+DJJ(9,6)*BZ(I1,IJ)+DJJ(9,7)*CXY(I1,IJ)+DJJ(9,8)*CYZ(I1,IJ)+DJJ(9,9)*CZX(I1,IJ)
                UU1(IJ)=U(ITT,IE(IJ,IM))
                VV1(IJ)=V(ITT,IE(IJ,IM))
                WW1(IJ)=W(ITT,IE(IJ,IM))
                END DO
                TAXU=OZ_DDOT(IAN, TAXUVP1, 1, UU1, 1)
                TAYU=OZ_DDOT(IAN, TAYUVP1, 1, UU1, 1)
                TAZU=OZ_DDOT(IAN, TAZUVP1, 1, UU1, 1)
                TAXV=OZ_DDOT(IAN, TAXUVP1, 1, VV1, 1)
                TAYV=OZ_DDOT(IAN, TAYUVP1, 1, VV1, 1)
                TAZV=OZ_DDOT(IAN, TAZUVP1, 1, VV1, 1)
                TAXW=OZ_DDOT(IAN, TAXUVP1, 1, WW1, 1)
                TAYW=OZ_DDOT(IAN, TAYUVP1, 1, WW1, 1)
                TAZW=OZ_DDOT(IAN, TAZUVP1, 1, WW1, 1)
                TBXU=OZ_DDOT(IAN, TBXUVP1, 1, UU1, 1)
                TBYV=OZ_DDOT(IAN, TBYUVP1, 1, VV1, 1)
                TBZW=OZ_DDOT(IAN, TBZUVP1, 1, WW1, 1)
                TCXYU=OZ_DDOT(IAN, TCXYUVP1, 1, UU1, 1)
                TCZXU=OZ_DDOT(IAN, TCZXUVP1, 1, UU1, 1)
                TCYZV=OZ_DDOT(IAN, TCYZUVP1, 1, VV1, 1)
                TCXYV=OZ_DDOT(IAN, TCXYUVP1, 1, VV1, 1)
                TCYZW=OZ_DDOT(IAN, TCYZUVP1, 1, WW1, 1)
                TCZXW=OZ_DDOT(IAN, TCZXUVP1, 1, WW1, 1)
                END IF
                
                BBUVP(IIB,JJB)=TAXUVP1(J1)
                BBUVP(IIB,JJB+IAN*IAT)=TAYUVP1(J1)
                BBUVP(IIB,JJB+2*IAN*IAT)=TAZUVP1(J1)
                BBUVP(IIB+IAN*IAT,JJB)=0.0D0*AT(ITT,JTT)/DT+U(ITT,II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+TAXU+V(ITT,II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(ITT,II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))&
                     -ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+IAN*IAT,JJB+IAN*IAT)=TAYU
                BBUVP(IIB+IAN*IAT,JJB+2*IAN*IAT)=TAZU
                BBUVP(IIB+IAN*IAT,JJB+3*IAN*IAT)=1.0D0/ROU*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))
                BBUVP(IIB+2*IAN*IAT,JJB)=TAXV
                BBUVP(IIB+2*IAN*IAT,JJB+IAN*IAT)=0.0D0*AT(ITT,JTT)/DT+U(ITT,II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(ITT,II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+TAYV+W(ITT,II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))&
                     -ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+2*IAN*IAT,JJB+2*IAN*IAT)=TAZV
                BBUVP(IIB+2*IAN*IAT,JJB+3*IAN*IAT)=1.0D0/ROU*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+&
                     DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))
                BBUVP(IIB+3*IAN*IAT,JJB)=TAXW
                BBUVP(IIB+3*IAN*IAT,JJB+IAN*IAT)=TAYW
                BBUVP(IIB+3*IAN*IAT,JJB+2*IAN*IAT)=0.0D0*AT(ITT,JTT)/DT+U(ITT,II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(ITT,II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(ITT,II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))+TAZW-&
                     ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+3*IAN*IAT,JJB+3*IAN*IAT)=1.0D0/ROU*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                     DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))
                END IF
                IF(I1.NE.J1) THEN
                BBUVP(IIB,JJB)=TAXUVP1(J1)
                BBUVP(IIB,JJB+IAN*IAT)=TAYUVP1(J1)
                BBUVP(IIB,JJB+2*IAN*IAT)=TAZUVP1(J1)
                BBUVP(IIB,JJB+3*IAN*IAT)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB)=U(ITT,II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(ITT,II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(ITT,II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))&
                     -ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+IAN*IAT,JJB+IAN*IAT)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB+2*IAN*IAT)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB+3*IAN*IAT)=1.0D0/ROU*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))
                BBUVP(IIB+2*IAN*IAT,JJB)=0.0D0
                BBUVP(IIB+2*IAN*IAT,JJB+IAN*IAT)=U(ITT,II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(ITT,II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(ITT,II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))-&
                     ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+2*IAN*IAT,JJB+2*IAN*IAT)=0.0D0
                BBUVP(IIB+2*IAN*IAT,JJB+3*IAN*IAT)=1.0D0/ROU*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+&
                     DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))
                BBUVP(IIB+3*IAN*IAT,JJB)=0.0D0
                BBUVP(IIB+3*IAN*IAT,JJB+IAN*IAT)=0.0D0
                BBUVP(IIB+3*IAN*IAT,JJB+2*IAN*IAT)=U(ITT,II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(ITT,II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(ITT,II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))-&
                     ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+3*IAN*IAT,JJB+3*IAN*IAT)=1.0D0/ROU*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                     DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))
                END IF
                END IF
                              
                IF(ITT.NE.JTT) THEN
                IF(I1.EQ.J1) THEN
                BBUVP(IIB,JJB)=0.0D0
                BBUVP(IIB,JJB+IAN*IAT)=0.0D0
                BBUVP(IIB,JJB+2*IAN*IAT)=0.0D0
                BBUVP(IIB,JJB+3*IAN*IAT)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB)=0.0D0*AT(ITT,JTT)/DT
                BBUVP(IIB+IAN*IAT,JJB+IAN*IAT)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB+2*IAN*IAT)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB+3*IAN*IAT)=0.0D0
                BBUVP(IIB+2*IAN*IAT,JJB)=0.0D0
                BBUVP(IIB+2*IAN*IAT,JJB+IAN*IAT)=0.0D0*AT(ITT,JTT)/DT
                BBUVP(IIB+2*IAN*IAT,JJB+2*IAN*IAT)=0.0D0
                BBUVP(IIB+2*IAN*IAT,JJB+3*IAN*IAT)=0.0D0
                BBUVP(IIB+3*IAN*IAT,JJB)=0.0D0
                BBUVP(IIB+3*IAN*IAT,JJB+IAN*IAT)=0.0D0
                BBUVP(IIB+3*IAN*IAT,JJB+2*IAN*IAT)=0.0D0*AT(ITT,JTT)/DT
                BBUVP(IIB+3*IAN*IAT,JJB+3*IAN*IAT)=0.0D0
                END IF
                IF(I1.NE.J1) THEN
                BBUVP(IIB,JJB)=0.0D0
                BBUVP(IIB,JJB+IAN*IAT)=0.0D0
                BBUVP(IIB,JJB+2*IAN*IAT)=0.0D0
                BBUVP(IIB,JJB+3*IAN*IAT)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB+IAN*IAT)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB+2*IAN*IAT)=0.0D0
                BBUVP(IIB+IAN*IAT,JJB+3*IAN*IAT)=0.0D0
                BBUVP(IIB+2*IAN*IAT,JJB)=0.0D0
                BBUVP(IIB+2*IAN*IAT,JJB+IAN*IAT)=0.0D0
                BBUVP(IIB+2*IAN*IAT,JJB+2*IAN*IAT)=0.0D0
                BBUVP(IIB+2*IAN*IAT,JJB+3*IAN*IAT)=0.0D0
                BBUVP(IIB+3*IAN*IAT,JJB)=0.0D0
                BBUVP(IIB+3*IAN*IAT,JJB+IAN*IAT)=0.0D0
                BBUVP(IIB+3*IAN*IAT,JJB+2*IAN*IAT)=0.0D0
                BBUVP(IIB+3*IAN*IAT,JJB+3*IAN*IAT)=0.0D0
                END IF
                END IF
                IF(.FALSE.) THEN
                  BBUVP(IIB,JJB)=0.0D0
                  BBUVP(IIB,JJB+IAN*IAT)=0.0D0
                  BBUVP(IIB,JJB+2*IAN*IAT)=0.0D0
                  BBUVP(IIB,JJB+3*IAN*IAT)=0.0D0
                  IF(ITT.EQ.JTT) THEN
                    DDELTA_NODE=0.0D0
                    IF(I1.EQ.J1) DDELTA_NODE=1.0D0
                    DDX=TAXUVP1(J1)
                    DDY=TAYUVP1(J1)
                    DDZ=TAZUVP1(J1)
                    DDXX=TBXUVP1(J1)
                    DDYY=TBYUVP1(J1)
                    DDZZ=TBZUVP1(J1)
                    DDXY=TCXYUVP1(J1)
                    DDYZ=TCYZUVP1(J1)
                    DDZX=TCZXUVP1(J1)
                    BBUVP(IIB,JJB)=2.0D0*TAXU*DDX+U(ITT,II)*DDXX+2.0D0*TAXV*DDY+V(ITT,II)*DDXY+&
                         2.0D0*TAXW*DDZ+W(ITT,II)*DDZX+DDELTA_NODE*(TBXU+TCXYV+TCZXW)
                    BBUVP(IIB,JJB+IAN*IAT)=2.0D0*TAYU*DDX+U(ITT,II)*DDXY+2.0D0*TAYV*DDY+V(ITT,II)*DDYY+&
                         2.0D0*TAYW*DDZ+W(ITT,II)*DDYZ+DDELTA_NODE*(TCXYU+TBYV+TCYZW)
                    BBUVP(IIB,JJB+2*IAN*IAT)=2.0D0*TAZU*DDX+U(ITT,II)*DDZX+2.0D0*TAZV*DDY+V(ITT,II)*DDYZ+&
                         2.0D0*TAZW*DDZ+W(ITT,II)*DDZZ+DDELTA_NODE*(TCZXU+TCYZV+TBZW)
                    BBUVP(IIB,JJB+3*IAN*IAT)=(DDXX+DDYY+DDZZ)/ROU
                  END IF
                END IF
                IF(ITT.EQ.JTT) THEN
!                 Local equation row 1 is the continuity equation:
!                   du/dx + dv/dy + dw/dz = 0.
                  BBUVP(IIB,JJB)=TAXUVP1(J1)
                  BBUVP(IIB,JJB+IAN*IAT)=TAYUVP1(J1)
                  BBUVP(IIB,JJB+2*IAN*IAT)=TAZUVP1(J1)
                  BBUVP(IIB,JJB+3*IAN*IAT)=0.0D0
                END IF
                IF(IPICARD_CONV_JAC.NE.0.AND.ITT.EQ.JTT.AND.I1.EQ.J1) THEN
                  BBUVP(IIB+IAN*IAT,JJB)=BBUVP(IIB+IAN*IAT,JJB)-TAXU
                  BBUVP(IIB+IAN*IAT,JJB+IAN*IAT)=BBUVP(IIB+IAN*IAT,JJB+IAN*IAT)-TAYU
                  BBUVP(IIB+IAN*IAT,JJB+2*IAN*IAT)=BBUVP(IIB+IAN*IAT,JJB+2*IAN*IAT)-TAZU
                  BBUVP(IIB+2*IAN*IAT,JJB)=BBUVP(IIB+2*IAN*IAT,JJB)-TAXV
                  BBUVP(IIB+2*IAN*IAT,JJB+IAN*IAT)=BBUVP(IIB+2*IAN*IAT,JJB+IAN*IAT)-TAYV
                  BBUVP(IIB+2*IAN*IAT,JJB+2*IAN*IAT)=BBUVP(IIB+2*IAN*IAT,JJB+2*IAN*IAT)-TAZV
                  BBUVP(IIB+3*IAN*IAT,JJB)=BBUVP(IIB+3*IAN*IAT,JJB)-TAXW
                  BBUVP(IIB+3*IAN*IAT,JJB+IAN*IAT)=BBUVP(IIB+3*IAN*IAT,JJB+IAN*IAT)-TAYW
                  BBUVP(IIB+3*IAN*IAT,JJB+2*IAN*IAT)=BBUVP(IIB+3*IAN*IAT,JJB+2*IAN*IAT)-TAZW
                END IF
                IF(ITT.EQ.JTT) THEN
                  DLIN_LAP=DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+&
                       DJJ(4,5)*BY(I1,J1)+DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+&
                       DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+DJJ(5,4)*BX(I1,J1)+&
                       DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+&
                       DJJ(6,1)*AX(I1,J1)+DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+&
                       DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+DJJ(6,9)*CZX(I1,J1)
                  DLIN_GX=DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+&
                       DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1)
                  DLIN_GY=DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                       DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1)
                  DLIN_GZ=DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+&
                       DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1)
                  BBUVP_LIN(IIB,JJB)=DLIN_GX
                  BBUVP_LIN(IIB,JJB+IAN*IAT)=DLIN_GY
                  BBUVP_LIN(IIB,JJB+2*IAN*IAT)=DLIN_GZ
                  BBUVP_LIN(IIB,JJB+3*IAN*IAT)=0.0D0
                  BBUVP_LIN(IIB+IAN*IAT,JJB)=-ANYU*DLIN_LAP
                  BBUVP_LIN(IIB+IAN*IAT,JJB+3*IAN*IAT)=DLIN_GX/ROU
                  BBUVP_LIN(IIB+2*IAN*IAT,JJB+IAN*IAT)=-ANYU*DLIN_LAP
                  BBUVP_LIN(IIB+2*IAN*IAT,JJB+3*IAN*IAT)=DLIN_GY/ROU
                  BBUVP_LIN(IIB+3*IAN*IAT,JJB+2*IAN*IAT)=-ANYU*DLIN_LAP
                  BBUVP_LIN(IIB+3*IAN*IAT,JJB+3*IAN*IAT)=DLIN_GZ/ROU
                END IF
             END IF
             
         !   write(*,*)'IM,j1,jtt,jjb=',iM,j1,jtt,jjb
                
          END DO
          END DO
          
          
          IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ).OR.(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ&
              .AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
!            write(*,*)'clear 885!'

          DO IJ=1,IAN
            TAXUVP1(IJ)=DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ)
            TAYUVP1(IJ)=DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ)
            TAZUVP1(IJ)=DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ)
            UU1(IJ)=U(ITT,IE(IJ,IM))
            VV1(IJ)=V(ITT,IE(IJ,IM))
            WW1(IJ)=W(ITT,IE(IJ,IM))
            PP1(IJ)=P(ITT,IE(IJ,IM))
          END DO

          TAXU=OZ_DDOT(IAN, TAXUVP1, 1, UU1, 1)
          TAYU=OZ_DDOT(IAN, TAYUVP1, 1, UU1, 1)
          TAZU=OZ_DDOT(IAN, TAZUVP1, 1, UU1, 1)
          TAXV=OZ_DDOT(IAN, TAXUVP1, 1, VV1, 1)
          TAYV=OZ_DDOT(IAN, TAYUVP1, 1, VV1, 1)
          TAZV=OZ_DDOT(IAN, TAZUVP1, 1, VV1, 1)
          TAXW=OZ_DDOT(IAN, TAXUVP1, 1, WW1, 1)
          TAYW=OZ_DDOT(IAN, TAYUVP1, 1, WW1, 1)
          TAZW=OZ_DDOT(IAN, TAZUVP1, 1, WW1, 1)
          TAXP=OZ_DDOT(IAN, TAXUVP1, 1, PP1, 1)
          TAYP=OZ_DDOT(IAN, TAYUVP1, 1, PP1, 1)
          TAZP=OZ_DDOT(IAN, TAZUVP1, 1, PP1, 1)
          FFUVP(IIB)=TAXU*XL(I1,IM)+TAYU*XM(I1,IM)+TAZU*XN(I1,IM)
          FFUVP(IIB+IAN*IAT)=TAXV*XL(I1,IM)+TAYV*XM(I1,IM)+TAZV*XN(I1,IM)
          FFUVP(IIB+IAN*IAT*2)=TAXW*XL(I1,IM)+TAYW*XM(I1,IM)+TAZW*XN(I1,IM)
          FFUVP(IIB+IAN*IAT*3)=TAXP*XL(I1,IM)+TAYP*XM(I1,IM)+TAZP*XN(I1,IM)
        END IF
        
!        write(*,*)'clear 005!'

        IF(I1.GE.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ)+NAX*NAY*NAZ) THEN
          DO IJ=1,IAN
            TAXUVP1(IJ)=DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ)
            TAYUVP1(IJ)=DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ)
            TAZUVP1(IJ)=DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ)
            TBXUVP1(IJ)=DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)
            TBYUVP1(IJ)=DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)
            TBZUVP1(IJ)=DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ)
            TCXYUVP1(IJ)=DJJ(7,1)*AX(I1,IJ)+DJJ(7,2)*AY(I1,IJ)+DJJ(7,3)*AZ(I1,IJ)+DJJ(7,4)*BX(I1,IJ)+DJJ(7,5)*BY(I1,IJ)+DJJ(7,6)*BZ(I1,IJ)+DJJ(7,7)*CXY(I1,IJ)+DJJ(7,8)*CYZ(I1,IJ)+DJJ(7,9)*CZX(I1,IJ)
            TCYZUVP1(IJ)=DJJ(8,1)*AX(I1,IJ)+DJJ(8,2)*AY(I1,IJ)+DJJ(8,3)*AZ(I1,IJ)+DJJ(8,4)*BX(I1,IJ)+DJJ(8,5)*BY(I1,IJ)+DJJ(8,6)*BZ(I1,IJ)+DJJ(8,7)*CXY(I1,IJ)+DJJ(8,8)*CYZ(I1,IJ)+DJJ(8,9)*CZX(I1,IJ)
            TCZXUVP1(IJ)=DJJ(9,1)*AX(I1,IJ)+DJJ(9,2)*AY(I1,IJ)+DJJ(9,3)*AZ(I1,IJ)+DJJ(9,4)*BX(I1,IJ)+DJJ(9,5)*BY(I1,IJ)+DJJ(9,6)*BZ(I1,IJ)+DJJ(9,7)*CXY(I1,IJ)+DJJ(9,8)*CYZ(I1,IJ)+DJJ(9,9)*CZX(I1,IJ)
            UU1(IJ)=U(ITT,IE(IJ,IM))
            VV1(IJ)=V(ITT,IE(IJ,IM))
            WW1(IJ)=W(ITT,IE(IJ,IM))
            PP1(IJ)=P(ITT,IE(IJ,IM))
          END DO
          DO IJ=1,IAT
            TATUVP1(IJ)=0.0D0*AT(ITT,IJ)/DT
            UUT1(IJ)=U(IJ,IE(I1,IM))
            VVT1(IJ)=V(IJ,IE(I1,IM))
            WWT1(IJ)=W(IJ,IE(I1,IM))
          END DO
          TAXU=OZ_DDOT(IAN, TAXUVP1, 1, UU1, 1)
          TAYU=OZ_DDOT(IAN, TAYUVP1, 1, UU1, 1)
          TAZU=OZ_DDOT(IAN, TAZUVP1, 1, UU1, 1)
          TAXV=OZ_DDOT(IAN, TAXUVP1, 1, VV1, 1)
          TAYV=OZ_DDOT(IAN, TAYUVP1, 1, VV1, 1)
          TAZV=OZ_DDOT(IAN, TAZUVP1, 1, VV1, 1)
          TAXW=OZ_DDOT(IAN, TAXUVP1, 1, WW1, 1)
          TAYW=OZ_DDOT(IAN, TAYUVP1, 1, WW1, 1)
          TAZW=OZ_DDOT(IAN, TAZUVP1, 1, WW1, 1)
          TBXU=OZ_DDOT(IAN, TBXUVP1, 1, UU1, 1)
          TBYU=OZ_DDOT(IAN, TBYUVP1, 1, UU1, 1)
          TBZU=OZ_DDOT(IAN, TBZUVP1, 1, UU1, 1)
          TBXV=OZ_DDOT(IAN, TBXUVP1, 1, VV1, 1)
          TBYV=OZ_DDOT(IAN, TBYUVP1, 1, VV1, 1)
          TBZV=OZ_DDOT(IAN, TBZUVP1, 1, VV1, 1)
          TBXW=OZ_DDOT(IAN, TBXUVP1, 1, WW1, 1)
          TBYW=OZ_DDOT(IAN, TBYUVP1, 1, WW1, 1)
          TBZW=OZ_DDOT(IAN, TBZUVP1, 1, WW1, 1)
          TAXP=OZ_DDOT(IAN, TAXUVP1, 1, PP1, 1)
          TAYP=OZ_DDOT(IAN, TAYUVP1, 1, PP1, 1)
          TAZP=OZ_DDOT(IAN, TAZUVP1, 1, PP1, 1)
          TBXP=OZ_DDOT(IAN, TBXUVP1, 1, PP1, 1)
          TBYP=OZ_DDOT(IAN, TBYUVP1, 1, PP1, 1)
          TBZP=OZ_DDOT(IAN, TBZUVP1, 1, PP1, 1)
          TCXYU=OZ_DDOT(IAN, TCXYUVP1, 1, UU1, 1)
          TCYZU=OZ_DDOT(IAN, TCYZUVP1, 1, UU1, 1)
          TCZXU=OZ_DDOT(IAN, TCZXUVP1, 1, UU1, 1)
          TCXYV=OZ_DDOT(IAN, TCXYUVP1, 1, VV1, 1)
          TCYZV=OZ_DDOT(IAN, TCYZUVP1, 1, VV1, 1)
          TCZXV=OZ_DDOT(IAN, TCZXUVP1, 1, VV1, 1)
          TCXYW=OZ_DDOT(IAN, TCXYUVP1, 1, WW1, 1)
          TCYZW=OZ_DDOT(IAN, TCYZUVP1, 1, WW1, 1)
          TCZXW=OZ_DDOT(IAN, TCZXUVP1, 1, WW1, 1)
          TATU=OZ_DDOT(IAT, TATUVP1, 1, UUT1,1)
          TATV=OZ_DDOT(IAT, TATUVP1, 1, VVT1,1)
          TATW=OZ_DDOT(IAT, TATUVP1, 1, WWT1,1)

!           Continuity equation.
            FFUVP(IIB)=TAXU+TAYV+TAZW
            FFUVP(IIB+IAN*IAT)=TATU+U(ITT,II)*TAXU+V(ITT,II)*TAYU+W(ITT,II)*TAZU-ANYU*(TBYU+TBXU+TBZU)+1.0D0/ROU*TAXP
            FFUVP(IIB+IAN*IAT*2)=TATV+U(ITT,II)*TAXV+V(ITT,II)*TAYV+W(ITT,II)*TAZV-ANYU*(TBXV+TBYV+TBZV)+1.0D0/ROU*TAYP
            FFUVP(IIB+IAN*IAT*3)=TATW+U(ITT,II)*TAXW+V(ITT,II)*TAYW+W(ITT,II)*TAZW-ANYU*(TBXW+TBYW+TBZW)+1.0D0/ROU*TAZP
            FFUVP_LIN(IIB)=TAXU+TAYV+TAZW
            FFUVP_LIN(IIB+IAN*IAT)=-ANYU*(TBYU+TBXU+TBZU)+1.0D0/ROU*TAXP
            FFUVP_LIN(IIB+IAN*IAT*2)=-ANYU*(TBXV+TBYV+TBZV)+1.0D0/ROU*TAYP
            FFUVP_LIN(IIB+IAN*IAT*3)=-ANYU*(TBXW+TBYW+TBZW)+1.0D0/ROU*TAZP
            IF(IDIAG_LOCAL_MOM.NE.0.AND.IPK.GE.30.AND.ITT.EQ.1.AND.(II.EQ.344.OR.II.EQ.349.OR.II.EQ.379.OR.II.EQ.399.OR.II.EQ.404.OR.II.EQ.1759)) THEN
              WRITE(*,*)'LOCAL_MOM_SPLIT NODE/XYZ/U/V/W/P=',II,XX(II),YY(II),ZZ(II),U(ITT,II),V(ITT,II),W(ITT,II),P(ITT,II)
              WRITE(*,*)'LOCAL_CONT DIV/TAXU/TAYV/TAZW=',FFUVP(IIB),TAXU,TAYV,TAZW
              WRITE(*,*)'LOCAL_U TIME/CONV/VISC/PGRAD/LIN/TOTAL=',TATU,U(ITT,II)*TAXU+V(ITT,II)*TAYU+W(ITT,II)*TAZU,&
                   -ANYU*(TBYU+TBXU+TBZU),1.0D0/ROU*TAXP,FFUVP_LIN(IIB+IAN*IAT),FFUVP(IIB+IAN*IAT)
              WRITE(*,*)'LOCAL_V TIME/CONV/VISC/PGRAD/LIN/TOTAL=',TATV,U(ITT,II)*TAXV+V(ITT,II)*TAYV+W(ITT,II)*TAZV,&
                   -ANYU*(TBXV+TBYV+TBZV),1.0D0/ROU*TAYP,FFUVP_LIN(IIB+IAN*IAT*2),FFUVP(IIB+IAN*IAT*2)
              WRITE(*,*)'LOCAL_W TIME/CONV/VISC/PGRAD/LIN/TOTAL=',TATW,U(ITT,II)*TAXW+V(ITT,II)*TAYW+W(ITT,II)*TAZW,&
                   -ANYU*(TBXW+TBYW+TBZW),1.0D0/ROU*TAZP,FFUVP_LIN(IIB+IAN*IAT*3),FFUVP(IIB+IAN*IAT*3)
              WRITE(*,*)'LOCAL_EFFECTIVE_MOM DCONV/U/V/W=',DCONV_FACTOR,&
                   FFUVP_LIN(IIB+IAN*IAT)+DCONV_FACTOR*(FFUVP(IIB+IAN*IAT)-FFUVP_LIN(IIB+IAN*IAT)),&
                   FFUVP_LIN(IIB+IAN*IAT*2)+DCONV_FACTOR*(FFUVP(IIB+IAN*IAT*2)-FFUVP_LIN(IIB+IAN*IAT*2)),&
                   FFUVP_LIN(IIB+IAN*IAT*3)+DCONV_FACTOR*(FFUVP(IIB+IAN*IAT*3)-FFUVP_LIN(IIB+IAN*IAT*3))
              WRITE(*,*)'LOCAL_GRAD U TAX/TAY/TAZ V TAX/TAY/TAZ W TAX/TAY/TAZ=',TAXU,TAYU,TAZU,TAXV,TAYV,TAZV,TAXW,TAYW,TAZW
              WRITE(*,*)'LOCAL_LAP U/V/W=',TBYU+TBXU+TBZU,TBXV+TBYV+TBZV,TBXW+TBYW+TBZW
              WRITE(*,*)'LOCAL_PGRAD TAX/TAY/TAZ=',TAXP,TAYP,TAZP
!!$OMP CRITICAL(LOCAL_MOM_DEBUG)
              OPEN(96,FILE='results\local_mom_split_debug.dat',STATUS='UNKNOWN',POSITION='APPEND')
              WRITE(96,'(A,2I8,1P,28E24.15)')'LOCAL_SPLIT',II,ITT,XX(II),YY(II),ZZ(II),U(ITT,II),V(ITT,II),W(ITT,II),P(ITT,II),&
                   FFUVP(IIB),FFUVP(IIB+IAN*IAT),FFUVP(IIB+IAN*IAT*2),FFUVP(IIB+IAN*IAT*3),&
                   TATU,U(ITT,II)*TAXU+V(ITT,II)*TAYU+W(ITT,II)*TAZU,-ANYU*(TBYU+TBXU+TBZU),1.0D0/ROU*TAXP,&
                   TATV,U(ITT,II)*TAXV+V(ITT,II)*TAYV+W(ITT,II)*TAZV,-ANYU*(TBXV+TBYV+TBZV),1.0D0/ROU*TAYP,&
                   TATW,U(ITT,II)*TAXW+V(ITT,II)*TAYW+W(ITT,II)*TAZW,-ANYU*(TBXW+TBYW+TBZW),1.0D0/ROU*TAZP,&
                   TAXU,TAYV,TAZW,TAXP,TAYP,TAZP
              CLOSE(96)
!!$OMP END CRITICAL(LOCAL_MOM_DEBUG)
            END IF
 !           WRITE(*,*)'FFUVP=',FFUVP(IIB),FFUVP(IIB+IAN*IAT),FFUVP(IIB+2*IAN*IAT),FFUVP(IIB+3*IAN*IAT)
            END IF
            
!            write(*,*)'IM,i1,itt,iib=',IM,i1,itt,iib
              
        END DO
        END DO
!       !$OMP end do
!       !$OMP end parallel

        IF(DCONV_FACTOR.NE.1.0D0) THEN
          DO IPL1=1,4*IAN*IAT
            IBLOCK=(IPL1-1)/(IAN*IAT)+1
            IF(IBLOCK.GE.2.AND.IBLOCK.LE.4) THEN
              FFUVP(IPL1)=FFUVP_LIN(IPL1)+DCONV_FACTOR*(FFUVP(IPL1)-FFUVP_LIN(IPL1))
              DO IPL2=1,4*IAN*IAT
                BBUVP(IPL1,IPL2)=BBUVP_LIN(IPL1,IPL2)+DCONV_FACTOR*(BBUVP(IPL1,IPL2)-BBUVP_LIN(IPL1,IPL2))
              END DO
            END IF
          END DO
        END IF
        

        write(*,*)'clear ohkubo 66!'

        
!       ----------------- SUPERPOSITION OF BBUVP- AND TIME -------------

        DO IK1=1,4
!       Keep this loop serial while triplet assembly uses a shared append-only buffer.
!       Parallel calls to triplet_add race on SPARSE_ELEMENT_A%nnz and corrupt CSR input.
        DO JA1=1,IAN
            IF(JA1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ.AND.&
               JA1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ+NAX*NAY*NAZ) THEN
        DO ITT=1,IAT 
            II1=IE(JA1,IM)
!           Source equation order in FFUVP/BBUVP is CONT,U,V,W.
!           Global residual row order is CONT,V,W,U.
            IK1G=GLOBAL_EQUATION_BLOCK(IK1)
            ISP1=GLOBAL_EQUATION_INDEX(ITT,II1,IK1G,NX,IAT)
            IPL1=ITT+(JA1-1)*IAT+IAN*IAT*(IK1-1)
            IIT(ISP1)=IIT(ISP1)+1
!             Unknown and Newton correction column order is always U,V,W,P.
!             Unknown and Newton correction column order is U,V,W,P.
              DO IK2=1,4
              DO JA2=1,IAN
              DO JTT=1,IAT
                JJ2=IE(JA2,IM)
                ISP2=GLOBAL_DOF_INDEX(JTT,JJ2,IK2,NX,IAT)
                IPL2=JTT+(JA2-1)*IAT+IAN*IAT*(IK2-1)
                call triplet_add(SPARSE_ELEMENT_A, ISP1, ISP2, BBUVP(IPL1,IPL2))
                IF(BBUVP_LIN(IPL1,IPL2).NE.0.0D0) THEN
                  call triplet_add(SPARSE_LINEAR_A, ISP1, ISP2, BBUVP_LIN(IPL1,IPL2))
                END IF
              END DO             
              END DO
              END DO 
            FFS(ISP1)=FFS(ISP1)+FFUVP(IPL1)
            FFS_LINEAR(ISP1)=FFS_LINEAR(ISP1)+FFUVP_LIN(IPL1)
            IF(.FALSE..AND.IDIAG_LOCAL_MOM.NE.0.AND.ITT.EQ.1.AND.(II1.EQ.202.OR.II1.EQ.622).AND.IK1.GE.1.AND.IK1.LE.4) THEN
              WRITE(*,*)'SUPERPOSE_CONTRIB NODE/IM/JA1/EQ/ROW/FFUVP/FFLIN=',II1,IM,JA1,IK1,ISP1,FFUVP(IPL1),FFUVP_LIN(IPL1)
            END IF
        END DO
        END IF
        END DO
        END DO

      END DO
!
      MN_ELEMENT=SPARSE_ELEMENT_A%nnz
      WRITE(*,*)'ELEMENT TRIPLET CSR NNZ=',MN_ELEMENT
      IF(.FALSE..AND.IDIAG_LOCAL_MOM.NE.0) THEN
      DO ITARGET=1,4
        IF(ITARGET.EQ.1) ITNODE=202
        IF(ITARGET.EQ.2) ITNODE=344
        IF(ITARGET.EQ.3) ITNODE=374
        IF(ITARGET.EQ.4) ITNODE=622
        IF(IPK.LE.12) THEN
        DO IBLOCK=1,4
          IROWCHK=GLOBAL_EQUATION_INDEX(1,ITNODE,IBLOCK,NX,IAT)
          DROW_DOT=0.0D0
          DO I=1,SPARSE_LINEAR_A%nnz
            IF(SPARSE_LINEAR_A%row(I).EQ.IROWCHK) THEN
              ICOLCHK=SPARSE_LINEAR_A%col(I)
              ILOCAL=(ICOLCHK-1)/(NX*IAT)+1
              IPNODE=ICOLCHK-(ILOCAL-1)*NX*IAT
              IF(ILOCAL.EQ.1) DROW_DOT=DROW_DOT+SPARSE_LINEAR_A%val(I)*U(1,IPNODE)
              IF(ILOCAL.EQ.2) DROW_DOT=DROW_DOT+SPARSE_LINEAR_A%val(I)*V(1,IPNODE)
              IF(ILOCAL.EQ.3) DROW_DOT=DROW_DOT+SPARSE_LINEAR_A%val(I)*W(1,IPNODE)
              IF(ILOCAL.EQ.4) DROW_DOT=DROW_DOT+SPARSE_LINEAR_A%val(I)*P(1,IPNODE)
            END IF
          END DO
          WRITE(*,*)'POST_SUPERPOSE_ROW NODE/BLOCK/ROW/IIT/FFS/FFS_LINEAR/LIN_AX=',ITNODE,IBLOCK,IROWCHK,IIT(IROWCHK),FFS(IROWCHK),FFS_LINEAR(IROWCHK),DROW_DOT
        END DO
        END IF
      END DO
      END IF

              write(*,*)'clear ohkubo 1TC',TC

!     ---------------- FIRST X-DIFFERENTIAL .EQ. 0,  U,V,W  at B3 OUT FLOW --------------------
!     SLOPE U ---------
      
!     triplet_add updates a shared nnz counter; keep boundary triplet assembly serial.
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        DO ITT=1,IAT
        ISP1=GLOBAL_EQUATION_INDEX(ITT,II,4,NX,IAT)
        
        CALL TRANSF(5,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
          DO IJ=1,IAN
            ISP2=GLOBAL_DOF_INDEX(ITT,IE(IJ,5),BLOCK_U,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, (DJJ(1,1)*AX(J1,IJ)+DJJ(1,2)*AY(J1,IJ)+DJJ(1,3)*AZ(J1,IJ)+DJJ(1,4)*BX(J1,IJ)+DJJ(1,5)*BY(J1,IJ)+DJJ(1,6)*BZ(J1,IJ)+&
                 DJJ(1,7)*CXY(J1,IJ)+DJJ(1,8)*CYZ(J1,IJ)+DJJ(1,9)*CZX(J1,IJ)))
                 
          END DO
          
          TBYU=0.0D0
          DO IJ=1,IAN
          JJ=IE(IJ,5)
          TBYU=TBYU+(DJJ(1,1)*AX(J1,IJ)+DJJ(1,2)*AY(J1,IJ)+DJJ(1,3)*AZ(J1,IJ)+DJJ(1,4)*BX(J1,IJ)+DJJ(1,5)*BY(J1,IJ)+DJJ(1,6)*BZ(J1,IJ)+&
                 DJJ(1,7)*CXY(J1,IJ)+DJJ(1,8)*CYZ(J1,IJ)+DJJ(1,9)*CZX(J1,IJ))*U(ITT,JJ)
          END DO
          FFS(ISP1)=TBYU
          IS_REPLACED_ROW(ISP1)=1
          FFS_REPLACED_ROW(ISP1)=TBYU
        END DO
        END IF      
      END DO

!     SLOPE V -----------------
 
      
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        DO ITT=1,IAT
        ISP1=GLOBAL_EQUATION_INDEX(ITT,II,2,NX,IAT)
        
        CALL TRANSF(5,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
!         triplet_add updates a shared nnz counter; keep boundary triplet assembly serial.
          DO IJ=1,IAN
            ISP2=GLOBAL_DOF_INDEX(ITT,IE(IJ,5),BLOCK_V,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, (DJJ(1,1)*AX(J1,IJ)+DJJ(1,2)*AY(J1,IJ)+DJJ(1,3)*AZ(J1,IJ)+DJJ(1,4)*BX(J1,IJ)+DJJ(1,5)*BY(J1,IJ)+DJJ(1,6)*BZ(J1,IJ)+&
                 DJJ(1,7)*CXY(J1,IJ)+DJJ(1,8)*CYZ(J1,IJ)+DJJ(1,9)*CZX(J1,IJ)))
          END DO
          
          TBYV=0.0D0
!         !$OMP parallel private(IJ,JJ)
!         !$OMP do reduction(+ : TBYV)
          DO IJ=1,IAN
          JJ=IE(IJ,5)
          TBYV=TBYV+(DJJ(1,1)*AX(J1,IJ)+DJJ(1,2)*AY(J1,IJ)+DJJ(1,3)*AZ(J1,IJ)+DJJ(1,4)*BX(J1,IJ)+DJJ(1,5)*BY(J1,IJ)+DJJ(1,6)*BZ(J1,IJ)+&
                 DJJ(1,7)*CXY(J1,IJ)+DJJ(1,8)*CYZ(J1,IJ)+DJJ(1,9)*CZX(J1,IJ))*V(ITT,JJ)
          END DO
!         !$OMP end do
!         !$OMP end parallel
          FFS(ISP1)=TBYV
          IS_REPLACED_ROW(ISP1)=1
          FFS_REPLACED_ROW(ISP1)=TBYV
        END DO
        END IF      
      END DO
!     SLOPE W -----------------

      
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        DO ITT=1,IAT
        ISP1=GLOBAL_EQUATION_INDEX(ITT,II,3,NX,IAT)
        
        CALL TRANSF(5,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
!         triplet_add updates a shared nnz counter; keep boundary triplet assembly serial.
          DO IJ=1,IAN
            ISP2=GLOBAL_DOF_INDEX(ITT,IE(IJ,5),BLOCK_W,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, (DJJ(1,1)*AX(J1,IJ)+DJJ(1,2)*AY(J1,IJ)+DJJ(1,3)*AZ(J1,IJ)+DJJ(1,4)*BX(J1,IJ)+DJJ(1,5)*BY(J1,IJ)+DJJ(1,6)*BZ(J1,IJ)+&
                 DJJ(1,7)*CXY(J1,IJ)+DJJ(1,8)*CYZ(J1,IJ)+DJJ(1,9)*CZX(J1,IJ)))
          END DO
          
          TBYW=0.0D0
!         !$OMP parallel private(IJ,JJ)
!         !$OMP do reduction(+ : TBYW)
          DO IJ=1,IAN
          JJ=IE(IJ,5)
          TBYW=TBYW+(DJJ(1,1)*AX(J1,IJ)+DJJ(1,2)*AY(J1,IJ)+DJJ(1,3)*AZ(J1,IJ)+DJJ(1,4)*BX(J1,IJ)+DJJ(1,5)*BY(J1,IJ)+DJJ(1,6)*BZ(J1,IJ)+&
                 DJJ(1,7)*CXY(J1,IJ)+DJJ(1,8)*CYZ(J1,IJ)+DJJ(1,9)*CZX(J1,IJ))*W(ITT,JJ)
          END DO
!         !$OMP end do
!         !$OMP end parallel
          FFS(ISP1)=TBYW
          IS_REPLACED_ROW(ISP1)=1
          FFS_REPLACED_ROW(ISP1)=TBYW
        END DO
        END IF      
      END DO

!     ------------ SECOND X-DIFFERENTIAL OF PRESSURE AT INFLOW ----------------

      J=1
      CALL TRANSF(J,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
      
      DO I1=1,IAN
        IF((I1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ)) THEN
        DO ITT=1,IAT
                
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,I1)
        END DO
        END DO        
        
          JJ=IE(I1,J)
          ISP1=GLOBAL_EQUATION_INDEX(ITT,JJ,4,NX,IAT)
          I_SKIP_WALL=0
          DO I3=1,NSB1
            IF(NSBU(I3).EQ.JJ) I_SKIP_WALL=1
          END DO
          IF(I_SKIP_WALL.NE.0) CYCLE
!         triplet_add updates a shared nnz counter; keep boundary triplet assembly serial.
          DO IJ=1,IAN
            ISP2=GLOBAL_DOF_INDEX(ITT,IE(IJ,J),BLOCK_P,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, (DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)))
          END DO
          TBXP=0.0D0
!         !$OMP parallel private(IJ,JJ)
!         !$OMP do reduction(+ : TBXP)
          DO IJ=1,IAN
            JJ=IE(IJ,J)
            TBXP=TBXP+(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ))*P(ITT,JJ)
          END DO
!         !$OMP end do
!         !$OMP end parallel
          FFS(ISP1)=TBXP
          IS_REPLACED_ROW(ISP1)=1
          FFS_REPLACED_ROW(ISP1)=TBXP
      END DO
      END IF
    END DO


! -------- PRESSURE P BOUNDARY CONDITION AT WALL
      
      DO J1=1,5
        
        CALL TRANSF(J1,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        
        DO I1=1,IAN
        IF((I1.GE.9+4*(NAX+NAY+NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY).OR.&
            (I1.GE.1.AND.I1.LE.8).OR.(I1.GE.9.AND.I1.LE.8+4*(NAX+NAY+NAZ)).OR.&
            (I1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ).OR.&
            (I1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ)) THEN
          DO ITT=1,IAT
        
          DO IQ1=1,9
          DO IQ2=1,9
            DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,I1)
          END DO
          END DO        
        
          JJ=IE(I1,J1)
          ISP1=GLOBAL_EQUATION_INDEX(ITT,JJ,4,NX,IAT)
          I_SKIP_POUT=0
          DO I3=1,NSB3
            IF(NSBUO(I3).EQ.JJ) I_SKIP_POUT=1
          END DO
          IF(I_SKIP_POUT.NE.0) CYCLE
!         triplet_add updates a shared nnz counter; keep boundary triplet assembly serial.
          DO IJ=1,IAN
            DO JTT=1,IAT
            ISP2=GLOBAL_DOF_INDEX(JTT,IE(IJ,J1),BLOCK_U,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, -ANYU*XL(I1,J1)*(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+&
                 DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)+&
                 DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+&
                 DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)+DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+&
                 DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ)))
            ISP2=GLOBAL_DOF_INDEX(JTT,IE(IJ,J1),BLOCK_V,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, -ANYU*XM(I1,J1)*(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+&
                 DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)+&
                 DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+&
                 DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)+DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+&
                 DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ)))
            ISP2=GLOBAL_DOF_INDEX(JTT,IE(IJ,J1),BLOCK_W,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, -ANYU*XN(I1,J1)*(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+&
                 DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)+&
                 DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+&
                 DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)+DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+&
                 DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ)))
            ISP2=GLOBAL_DOF_INDEX(JTT,IE(IJ,J1),BLOCK_P,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, ((DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+&
                 DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ))*XL(I1,J1)+(DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+&
                 DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ))*XM(I1,J1)+&
                 (DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+&
                 DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ))*XN(I1,J1))/ROU)
            END DO
          END DO
          TAXYZP=0.0D0
          TNLAPU=0.0D0
          TNLAPV=0.0D0
          TNLAPW=0.0D0
!         !$OMP parallel private(IJ,JJ4)
!         !$OMP do reduction(+ : TAXYZP,TNLAPU,TNLAPV,TNLAPW)
          DO IJ=1,IAN
            DO JTT=1,IAT
            JJ4=IE(IJ,J1)
            TAXYZP=TAXYZP+((DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+&
                 DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ))*XL(I1,J1)+(DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+&
                 DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ))*XM(I1,J1)+&
                 (DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+&
                 DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ))*XN(I1,J1))*P(JTT,JJ4)
            TNLAPU=TNLAPU+(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)+DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+&
                 DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)+&
                 DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+&
                 DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ))*U(JTT,JJ4)
            TNLAPV=TNLAPV+(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)+DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+&
                 DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)+&
                 DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+&
                 DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ))*V(JTT,JJ4)
            TNLAPW=TNLAPW+(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)+DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+&
                 DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)+&
                 DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+&
                 DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ))*W(JTT,JJ4)
            END DO
          END DO
!         !$OMP end do
!         !$OMP end parallel
          TAXYZP=TAXYZP/ROU-ANYU*(XL(I1,J1)*TNLAPU+XM(I1,J1)*TNLAPV+XN(I1,J1)*TNLAPW)
          TDXP=0.0D0
          TDYP=0.0D0
          TDZP=0.0D0
!         !$OMP parallel private(IJ,JJ4)
!         !$OMP do reduction(+ : TDXP,TDYP,TDZP)
          DO IJ=1,IAN
            DO JTT=1,IAT
            JJ4=IE(IJ,J1)
            TDXP=TDXP+(DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+&
                 DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ))*P(JTT,JJ4)
            TDYP=TDYP+(DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+&
                 DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ))*P(JTT,JJ4)
            TDZP=TDZP+(DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+&
                 DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ))*P(JTT,JJ4)
            END DO
          END DO
!         !$OMP end do
!         !$OMP end parallel
          IF(JJ.EQ.1786.AND.ITT.EQ.1) THEN
            OPEN(97,FILE='results\pressure_wall_node1786_debug.dat',STATUS='UNKNOWN',POSITION='APPEND')
            WRITE(97,'(A)') '# node j1 i1 time x y z p nx ny nz dpdx dpdy dpdz wall_normal_momentum'
            WRITE(97,'(4I10,1P,11E24.14)') JJ,J1,I1,ITT,XX(JJ),YY(JJ),ZZ(JJ),P(ITT,JJ),XL(I1,J1),XM(I1,J1),XN(I1,J1),TDXP,TDYP,TDZP,TAXYZP
            CLOSE(97)
            WRITE(*,*)'DEBUG_P_WALL_1786 J1/I1/XYZ/P=',J1,I1,XX(JJ),YY(JJ),ZZ(JJ),P(ITT,JJ)
            WRITE(*,*)'DEBUG_P_WALL_1786 NORMAL/DP/NM=',XL(I1,J1),XM(I1,J1),XN(I1,J1),TDXP,TDYP,TDZP,TAXYZP
          END IF
          IS_REPLACED_ROW(ISP1)=1
          FFS_REPLACED_ROW(ISP1)=TAXYZP
          END DO
          END IF      
        END DO
      END DO

!     Recompute all replaced boundary residual rows from the assembled boundary
!     matrix. This keeps linear boundary FFS rows exactly consistent with A*x.
      DREPLACED_AXB_MAX=0.0D0
      DREPLACED_AXB_SUM=0.0D0
      NREPLACED_AXB=0
      IREPLACED_AXB_MAX=0
      DO I=1,N
        IF(IS_REPLACED_ROW(I).NE.0) FFS_REPLACED_ROW(I)=0.0D0
      END DO
      DO I=1,SPARSE_PRESSURE_WALL_A%nnz
        IROW=SPARSE_PRESSURE_WALL_A%row(I)
        IF(IS_REPLACED_ROW(IROW).NE.0) THEN
          ICOLCHK=SPARSE_PRESSURE_WALL_A%col(I)
          IBLOCK=(ICOLCHK-1)/(NX*IAT)+1
          ILOCAL=ICOLCHK-(IBLOCK-1)*NX*IAT
          IPNODE=(ILOCAL-1)/IAT+1
          IPTIME=ILOCAL-(IPNODE-1)*IAT
          DROW_VAR=0.0D0
          IF(IBLOCK.EQ.1) DROW_VAR=U(IPTIME,IPNODE)
          IF(IBLOCK.EQ.2) DROW_VAR=V(IPTIME,IPNODE)
          IF(IBLOCK.EQ.3) DROW_VAR=W(IPTIME,IPNODE)
          IF(IBLOCK.EQ.4) DROW_VAR=P(IPTIME,IPNODE)
          FFS_REPLACED_ROW(IROW)=FFS_REPLACED_ROW(IROW)+SPARSE_PRESSURE_WALL_A%val(I)*DROW_VAR
        END IF
      END DO
      DO I=1,N
        IF(IS_REPLACED_ROW(I).NE.0) THEN
          IF(DABS(FFS(I)-FFS_REPLACED_ROW(I)).GT.DREPLACED_AXB_MAX) THEN
            DREPLACED_AXB_MAX=DABS(FFS(I)-FFS_REPLACED_ROW(I))
            IREPLACED_AXB_MAX=I
          END IF
          DREPLACED_AXB_SUM=DREPLACED_AXB_SUM+DABS(FFS(I)-FFS_REPLACED_ROW(I))
          NREPLACED_AXB=NREPLACED_AXB+1
          FFS(I)=FFS_REPLACED_ROW(I)
        END IF
      END DO
      WRITE(*,*)'REPLACED_ROW_AX_MINUS_B_RECOMPUTE COUNT/MAXDIFF/SUMDIFF=',NREPLACED_AXB,DREPLACED_AXB_MAX,DREPLACED_AXB_SUM
      IF(IREPLACED_AXB_MAX.GT.0) THEN
        IBLOCK=(IREPLACED_AXB_MAX-1)/(NX*IAT)+1
        ILOCAL=IREPLACED_AXB_MAX-(IBLOCK-1)*NX*IAT
        IPNODE=(ILOCAL-1)/IAT+1
        IPTIME=ILOCAL-(IPNODE-1)*IAT
        WRITE(*,*)'REPLACED_ROW_AX_MINUS_B_MAX ROW/BLOCK/NODE/TIME/XYZ=',IREPLACED_AXB_MAX,IBLOCK,IPNODE,IPTIME,XX(IPNODE),YY(IPNODE),ZZ(IPNODE)
      END IF
      DINLET_P_RES_MAX=0.0D0
      DINLET_P_RES_SUM=0.0D0
      DINLET_P_RES_L2=0.0D0
      NINLET_P_RES=0
      IINLET_P_RES_MAX=0
      DO I3=1,NSB2
        II=NSBUI(I3)
        DO ITT=1,IAT
          IROW=GLOBAL_EQUATION_INDEX(ITT,II,4,NX,IAT)
          IF(IS_REPLACED_ROW(IROW).NE.0) THEN
            DROW_DOT=0.0D0
            DO I=1,SPARSE_PRESSURE_WALL_A%nnz
              IF(SPARSE_PRESSURE_WALL_A%row(I).EQ.IROW) THEN
                ICOLCHK=SPARSE_PRESSURE_WALL_A%col(I)
                IBLOCK=(ICOLCHK-1)/(NX*IAT)+1
                ILOCAL=ICOLCHK-(IBLOCK-1)*NX*IAT
                IPNODE=(ILOCAL-1)/IAT+1
                IPTIME=ILOCAL-(IPNODE-1)*IAT
                DROW_VAR=0.0D0
                IF(IBLOCK.EQ.1) DROW_VAR=U(IPTIME,IPNODE)
                IF(IBLOCK.EQ.2) DROW_VAR=V(IPTIME,IPNODE)
                IF(IBLOCK.EQ.3) DROW_VAR=W(IPTIME,IPNODE)
                IF(IBLOCK.EQ.4) DROW_VAR=P(IPTIME,IPNODE)
                DROW_DOT=DROW_DOT+SPARSE_PRESSURE_WALL_A%val(I)*DROW_VAR
              END IF
            END DO
            DAX_DIFF=DABS(DROW_DOT)
            DINLET_P_RES_SUM=DINLET_P_RES_SUM+DAX_DIFF
            DINLET_P_RES_L2=DINLET_P_RES_L2+DAX_DIFF*DAX_DIFF
            NINLET_P_RES=NINLET_P_RES+1
            IF(DAX_DIFF.GT.DINLET_P_RES_MAX) THEN
              DINLET_P_RES_MAX=DAX_DIFF
              IINLET_P_RES_MAX=IROW
              DINLET_P_RES_DOT=DROW_DOT
              IINLET_P_RES_NODE=II
              IINLET_P_RES_TIME=ITT
            END IF
          END IF
        END DO
      END DO
      DINLET_P_RES_L2=DSQRT(DINLET_P_RES_L2)
      WRITE(*,*)'INLET_PRESSURE_D2_ROW_FULL_RES COUNT/MAX/SUM/L2=',NINLET_P_RES,DINLET_P_RES_MAX,DINLET_P_RES_SUM,DINLET_P_RES_L2
      IF(IINLET_P_RES_MAX.GT.0) WRITE(*,*)'INLET_PRESSURE_D2_ROW_FULL_MAX ROW/NODE/TIME/XYZ/DOT=',IINLET_P_RES_MAX,IINLET_P_RES_NODE,IINLET_P_RES_TIME,XX(IINLET_P_RES_NODE),YY(IINLET_P_RES_NODE),ZZ(IINLET_P_RES_NODE),DINLET_P_RES_DOT

!     --------------  BUILD DIRICHLET CONSTRAINT LIST -----------------

      IS_CONSTRAINED=0
      BC_ROW=0
      BC_DELTA=0.0D0
      NCONSTRAINT=0
      NCONFLICT=0
      NPOUT_CONSTRAINT=0

      DO J=1,NSB1
        DO ITT=1,IAT
          JJ=NSBU(J)

          ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_U,NX,IAT)
          ISP1=GLOBAL_EQUATION_INDEX(ITT,JJ,4,NX,IAT)
          BCVAL=SBU(J)-U(ITT,JJ)
          IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
          IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
          IS_CONSTRAINED(ISP2)=1
          BC_ROW(ISP2)=ISP1
          IS_CONSTRAINED_ROW(ISP1)=1
          BC_DELTA(ISP2)=BCVAL

          ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_V,NX,IAT)
          ISP1=GLOBAL_EQUATION_INDEX(ITT,JJ,2,NX,IAT)
          BCVAL=SBV(J)-V(ITT,JJ)
          IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
          IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
          IS_CONSTRAINED(ISP2)=1
          BC_ROW(ISP2)=ISP1
          IS_CONSTRAINED_ROW(ISP1)=1
          BC_DELTA(ISP2)=BCVAL

          ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_W,NX,IAT)
          ISP1=GLOBAL_EQUATION_INDEX(ITT,JJ,3,NX,IAT)
          BCVAL=SBW(J)-W(ITT,JJ)
          IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
          IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
          IS_CONSTRAINED(ISP2)=1
          BC_ROW(ISP2)=ISP1
          IS_CONSTRAINED_ROW(ISP1)=1
          BC_DELTA(ISP2)=BCVAL
        END DO
      END DO

      DO I3=1,NSB3
        DO ITT=1,IAT
          II=NSBUO(I3)
          ISP2=GLOBAL_DOF_INDEX(ITT,II,BLOCK_P,NX,IAT)
          ISP1=GLOBAL_EQUATION_INDEX(ITT,II,4,NX,IAT)
         BCVAL=SBPO(I3)-P(ITT,II)
          IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
          IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
          IS_CONSTRAINED(ISP2)=1
          BC_ROW(ISP2)=ISP1
          IS_CONSTRAINED_ROW(ISP1)=1
          BC_DELTA(ISP2)=BCVAL
          NPOUT_CONSTRAINT=NPOUT_CONSTRAINT+1
          IF(ITT.EQ.1.AND.I3.LE.8) WRITE(*,*)'P_OUT_DIRICHLET NODE/VALUE=',II,SBPO(I3)
        END DO
      END DO

      IF(IAPPLY_INLET_VELOCITY_BC.NE.0) THEN
        DO J=1,NSB2
          DO ITT=1,IAT
            JJ=NSBUI(J)

            ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_U,NX,IAT)
            ISP1=GLOBAL_EQUATION_INDEX(ITT,JJ,4,NX,IAT)
            BCVAL=SBUI(J)-U(ITT,JJ)
            IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
            IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
            IS_CONSTRAINED(ISP2)=1
            BC_ROW(ISP2)=ISP1
            IS_CONSTRAINED_ROW(ISP1)=1
            BC_DELTA(ISP2)=BCVAL

            ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_V,NX,IAT)
            ISP1=GLOBAL_EQUATION_INDEX(ITT,JJ,2,NX,IAT)
            BCVAL=SBVI(J)-V(ITT,JJ)
            IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
            IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
            IS_CONSTRAINED(ISP2)=1
            BC_ROW(ISP2)=ISP1
            IS_CONSTRAINED_ROW(ISP1)=1
            BC_DELTA(ISP2)=BCVAL

            ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_W,NX,IAT)
            ISP1=GLOBAL_EQUATION_INDEX(ITT,JJ,3,NX,IAT)
            BCVAL=SBWI(J)-W(ITT,JJ)
            IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
            IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
            IS_CONSTRAINED(ISP2)=1
            BC_ROW(ISP2)=ISP1
            IS_CONSTRAINED_ROW(ISP1)=1
            BC_DELTA(ISP2)=BCVAL
          END DO
        END DO
      ELSE
        WRITE(*,*)'INLET VELOCITY DIRICHLET BC SKIPPED; USING STEADY STOKES INITIAL FIELD'
      END IF

      WRITE(*,*)'DIRICHLET DOF COUNT, CONFLICTS=',NCONSTRAINT,NCONFLICT
      WRITE(*,*)'P_OUT_DIRICHLET DOF COUNT=',NPOUT_CONSTRAINT
      NREPLACED_CONSTRAINED=0
      NRC1=0
      NRC2=0
      NRC3=0
      NRC4=0
      NRC_PRINT=0
      DO I=1,N
        IF(IS_REPLACED_ROW(I).NE.0.AND.IS_CONSTRAINED_ROW(I).NE.0) THEN
          NREPLACED_CONSTRAINED=NREPLACED_CONSTRAINED+1
          IF(I.LE.NX*IAT) THEN
            NRC1=NRC1+1
          ELSE IF(I.LE.2*NX*IAT) THEN
            NRC2=NRC2+1
          ELSE IF(I.LE.3*NX*IAT) THEN
            NRC3=NRC3+1
          ELSE
            NRC4=NRC4+1
          END IF
          IF(NRC_PRINT.LT.40) THEN
            IBLOCK=(I-1)/(NX*IAT)+1
            ILOCAL=I-(IBLOCK-1)*NX*IAT
            IPNODE=(ILOCAL-1)/IAT+1
            IPTIME=ILOCAL-(IPNODE-1)*IAT
            WRITE(*,*)'REPLACED_CONSTRAINED_DETAIL ROW/BLOCK/NODE/TIME/XYZ=',I,IBLOCK,IPNODE,IPTIME,XX(IPNODE),YY(IPNODE),ZZ(IPNODE)
            NRC_PRINT=NRC_PRINT+1
          END IF
        END IF
      END DO
      WRITE(*,*)'REPLACED_ROWS_CONSTRAINED TOTAL/CONT/V/W/U=',NREPLACED_CONSTRAINED,NRC1,NRC2,NRC3,NRC4
      NINLET_CONT_REPLACED=0
      NINLET_CONT_NORMAL=0
      NINLET_CONT_CONSTRAINED=0
      DO I3=1,NSB2
        II=NSBUI(I3)
        DO ITT=1,IAT
          IROW=GLOBAL_EQUATION_INDEX(ITT,II,1,NX,IAT)
          IF(IS_CONSTRAINED_ROW(IROW).NE.0) THEN
            NINLET_CONT_CONSTRAINED=NINLET_CONT_CONSTRAINED+1
          ELSE IF(IS_REPLACED_ROW(IROW).NE.0) THEN
            NINLET_CONT_REPLACED=NINLET_CONT_REPLACED+1
          ELSE
            NINLET_CONT_NORMAL=NINLET_CONT_NORMAL+1
            IF(NINLET_CONT_NORMAL.LE.8) WRITE(*,*)'INLET_CONT_NORMAL_ROW NODE/TIME/ROW/XYZ=',II,ITT,IROW,XX(II),YY(II),ZZ(II)
          END IF
        END DO
      END DO
      WRITE(*,*)'INLET_CONT_ROW_STATUS REPLACED/NORMAL/CONSTRAINED=',NINLET_CONT_REPLACED,NINLET_CONT_NORMAL,NINLET_CONT_CONSTRAINED
      FFS_PHYS=FFS
      DO ISP2=1,N
        IF(IS_CONSTRAINED(ISP2).NE.0) THEN
          IROW=BC_ROW(ISP2)
          IF(IROW.GT.0) FFS_PHYS(IROW)=BC_DELTA(ISP2)
        END IF
      END DO

      IF(IMATRIX_LINEAR_RESIDUAL.NE.0) THEN
        DMATRIX_RES_MAX=0.0D0
        DMATRIX_RES_SUM=0.0D0
        NMATRIX_RES=0
        IMATRIX_RES_MAX=0
        DO I=1,NX*IAT
          IF(IS_REPLACED_ROW(I).EQ.0.AND.IS_CONSTRAINED_ROW(I).EQ.0) THEN
            FFS_REPLACED_ROW(I)=0.0D0
          END IF
        END DO
        DO I=1,SPARSE_ELEMENT_A%nnz
          IROW=SPARSE_ELEMENT_A%row(I)
          IF(IROW.GE.1.AND.IROW.LE.NX*IAT) THEN
            IF(IS_REPLACED_ROW(IROW).EQ.0.AND.IS_CONSTRAINED_ROW(IROW).EQ.0) THEN
              ICOLCHK=SPARSE_ELEMENT_A%col(I)
              IBLOCK=(ICOLCHK-1)/(NX*IAT)+1
              ILOCAL=ICOLCHK-(IBLOCK-1)*NX*IAT
              IPNODE=(ILOCAL-1)/IAT+1
              IPTIME=ILOCAL-(IPNODE-1)*IAT
              DROW_VAR=0.0D0
              IF(IBLOCK.EQ.1) DROW_VAR=U(IPTIME,IPNODE)
              IF(IBLOCK.EQ.2) DROW_VAR=V(IPTIME,IPNODE)
              IF(IBLOCK.EQ.3) DROW_VAR=W(IPTIME,IPNODE)
              IF(IBLOCK.EQ.4) DROW_VAR=P(IPTIME,IPNODE)
              FFS_REPLACED_ROW(IROW)=FFS_REPLACED_ROW(IROW)+SPARSE_ELEMENT_A%val(I)*DROW_VAR
            END IF
          END IF
        END DO
        DO I=1,NX*IAT
          IF(IS_REPLACED_ROW(I).EQ.0.AND.IS_CONSTRAINED_ROW(I).EQ.0) THEN
            IF(DABS(FFS(I)-FFS_REPLACED_ROW(I)).GT.DMATRIX_RES_MAX) THEN
              DMATRIX_RES_MAX=DABS(FFS(I)-FFS_REPLACED_ROW(I))
              IMATRIX_RES_MAX=I
            END IF
            DMATRIX_RES_SUM=DMATRIX_RES_SUM+DABS(FFS(I)-FFS_REPLACED_ROW(I))
            NMATRIX_RES=NMATRIX_RES+1
            FFS(I)=FFS_REPLACED_ROW(I)
          END IF
        END DO
        WRITE(*,*)'MATRIX_LINEAR_RESIDUAL CONT COUNT/MAXDIFF/SUMDIFF=',NMATRIX_RES,DMATRIX_RES_MAX,DMATRIX_RES_SUM
        IF(IMATRIX_RES_MAX.GT.0) THEN
          IPNODE=(IMATRIX_RES_MAX-1)/IAT+1
          IPTIME=IMATRIX_RES_MAX-(IPNODE-1)*IAT
          WRITE(*,*)'MATRIX_LINEAR_RESIDUAL CONT MAX ROW/NODE/TIME/XYZ=',IMATRIX_RES_MAX,IPNODE,IPTIME,XX(IPNODE),YY(IPNODE),ZZ(IPNODE)
        END IF
        DO I=NX*IAT+1,4*NX*IAT
          IF(IS_REPLACED_ROW(I).EQ.0.AND.IS_CONSTRAINED_ROW(I).EQ.0) THEN
            FFS_REPLACED_ROW(I)=0.0D0
          END IF
        END DO
        DO I=1,SPARSE_LINEAR_A%nnz
          IROW=SPARSE_LINEAR_A%row(I)
          IF(IROW.GT.NX*IAT.AND.IROW.LE.4*NX*IAT) THEN
            IF(IS_REPLACED_ROW(IROW).EQ.0.AND.IS_CONSTRAINED_ROW(IROW).EQ.0) THEN
              ICOLCHK=SPARSE_LINEAR_A%col(I)
              IBLOCK=(ICOLCHK-1)/(NX*IAT)+1
              ILOCAL=ICOLCHK-(IBLOCK-1)*NX*IAT
              IPNODE=(ILOCAL-1)/IAT+1
              IPTIME=ILOCAL-(IPNODE-1)*IAT
              DROW_VAR=0.0D0
              IF(IBLOCK.EQ.1) DROW_VAR=U(IPTIME,IPNODE)
              IF(IBLOCK.EQ.2) DROW_VAR=V(IPTIME,IPNODE)
              IF(IBLOCK.EQ.3) DROW_VAR=W(IPTIME,IPNODE)
              IF(IBLOCK.EQ.4) DROW_VAR=P(IPTIME,IPNODE)
              FFS_REPLACED_ROW(IROW)=FFS_REPLACED_ROW(IROW)+SPARSE_LINEAR_A%val(I)*DROW_VAR
            END IF
          END IF
        END DO
        DMOM_AX_MAX_U=0.0D0
        DMOM_AX_MAX_V=0.0D0
        DMOM_AX_MAX_W=0.0D0
        DMOM_AX_SUM_U=0.0D0
        DMOM_AX_SUM_V=0.0D0
        DMOM_AX_SUM_W=0.0D0
        NMOM_AX_U=0
        NMOM_AX_V=0
        NMOM_AX_W=0
        IMOM_AX_U=0
        IMOM_AX_V=0
        IMOM_AX_W=0
        DO I=NX*IAT+1,4*NX*IAT
          IF(IS_REPLACED_ROW(I).EQ.0.AND.IS_CONSTRAINED_ROW(I).EQ.0) THEN
            DAX_DIFF=DABS(FFS_LINEAR(I)-FFS_REPLACED_ROW(I))
            IF(I.LE.2*NX*IAT) THEN
              NMOM_AX_U=NMOM_AX_U+1
              DMOM_AX_SUM_U=DMOM_AX_SUM_U+DAX_DIFF
              IF(DAX_DIFF.GT.DMOM_AX_MAX_U) THEN
                DMOM_AX_MAX_U=DAX_DIFF
                IMOM_AX_U=I
              END IF
            ELSE IF(I.LE.3*NX*IAT) THEN
              NMOM_AX_V=NMOM_AX_V+1
              DMOM_AX_SUM_V=DMOM_AX_SUM_V+DAX_DIFF
              IF(DAX_DIFF.GT.DMOM_AX_MAX_V) THEN
                DMOM_AX_MAX_V=DAX_DIFF
                IMOM_AX_V=I
              END IF
            ELSE
              NMOM_AX_W=NMOM_AX_W+1
              DMOM_AX_SUM_W=DMOM_AX_SUM_W+DAX_DIFF
              IF(DAX_DIFF.GT.DMOM_AX_MAX_W) THEN
                DMOM_AX_MAX_W=DAX_DIFF
                IMOM_AX_W=I
              END IF
            END IF
            FFS(I)=FFS(I)-FFS_LINEAR(I)+FFS_REPLACED_ROW(I)
          END IF
        END DO
        WRITE(*,*)'MATRIX_LINEAR_RESIDUAL MOM U COUNT/MAXDIFF/SUMDIFF=',NMOM_AX_U,DMOM_AX_MAX_U,DMOM_AX_SUM_U
        WRITE(*,*)'MATRIX_LINEAR_RESIDUAL MOM V COUNT/MAXDIFF/SUMDIFF=',NMOM_AX_V,DMOM_AX_MAX_V,DMOM_AX_SUM_V
        WRITE(*,*)'MATRIX_LINEAR_RESIDUAL MOM W COUNT/MAXDIFF/SUMDIFF=',NMOM_AX_W,DMOM_AX_MAX_W,DMOM_AX_SUM_W
        DO ITARGET=1,3
          IF(ITARGET.EQ.1) IMATRIX_RES_MAX=IMOM_AX_U
          IF(ITARGET.EQ.2) IMATRIX_RES_MAX=IMOM_AX_V
          IF(ITARGET.EQ.3) IMATRIX_RES_MAX=IMOM_AX_W
          IF(IMATRIX_RES_MAX.GT.0) THEN
            IBLOCK=(IMATRIX_RES_MAX-1)/(NX*IAT)+1
            ILOCAL=IMATRIX_RES_MAX-(IBLOCK-1)*NX*IAT
            IPNODE=(ILOCAL-1)/IAT+1
            IPTIME=ILOCAL-(IPNODE-1)*IAT
            WRITE(*,*)'MATRIX_LINEAR_RESIDUAL MOM MAX EQBLOCK/ROW/NODE/TIME/XYZ=',IBLOCK,IMATRIX_RES_MAX,IPNODE,IPTIME,XX(IPNODE),YY(IPNODE),ZZ(IPNODE)
          END IF
        END DO
        DO I=NX*IAT+1,4*NX*IAT
          IF(IS_REPLACED_ROW(I).EQ.0.AND.IS_CONSTRAINED_ROW(I).EQ.0) THEN
            FFS_REPLACED_ROW(I)=0.0D0
          END IF
        END DO
        DO I=1,SPARSE_ELEMENT_A%nnz
          IROW=SPARSE_ELEMENT_A%row(I)
          IF(IROW.GT.NX*IAT.AND.IROW.LE.4*NX*IAT) THEN
            IF(IS_REPLACED_ROW(IROW).EQ.0.AND.IS_CONSTRAINED_ROW(IROW).EQ.0) THEN
              ICOLCHK=SPARSE_ELEMENT_A%col(I)
              IBLOCK=(ICOLCHK-1)/(NX*IAT)+1
              ILOCAL=ICOLCHK-(IBLOCK-1)*NX*IAT
              IPNODE=(ILOCAL-1)/IAT+1
              IPTIME=ILOCAL-(IPNODE-1)*IAT
              DROW_VAR=0.0D0
              IF(IBLOCK.EQ.1) DROW_VAR=U(IPTIME,IPNODE)
              IF(IBLOCK.EQ.2) DROW_VAR=V(IPTIME,IPNODE)
              IF(IBLOCK.EQ.3) DROW_VAR=W(IPTIME,IPNODE)
              IF(IBLOCK.EQ.4) DROW_VAR=P(IPTIME,IPNODE)
              FFS_REPLACED_ROW(IROW)=FFS_REPLACED_ROW(IROW)+SPARSE_ELEMENT_A%val(I)*DROW_VAR
            END IF
          END IF
        END DO
        DMOM_AX_MAX_U=0.0D0
        DMOM_AX_MAX_V=0.0D0
        DMOM_AX_MAX_W=0.0D0
        DMOM_AX_SUM_U=0.0D0
        DMOM_AX_SUM_V=0.0D0
        DMOM_AX_SUM_W=0.0D0
        NMOM_AX_U=0
        NMOM_AX_V=0
        NMOM_AX_W=0
        IMOM_AX_U=0
        IMOM_AX_V=0
        IMOM_AX_W=0
        DO I=NX*IAT+1,4*NX*IAT
          IF(IS_REPLACED_ROW(I).EQ.0.AND.IS_CONSTRAINED_ROW(I).EQ.0) THEN
            DAX_DIFF=DABS(FFS(I)-FFS_REPLACED_ROW(I))
            IF(I.LE.2*NX*IAT) THEN
              NMOM_AX_U=NMOM_AX_U+1
              DMOM_AX_SUM_U=DMOM_AX_SUM_U+DAX_DIFF
              IF(DAX_DIFF.GT.DMOM_AX_MAX_U) THEN
                DMOM_AX_MAX_U=DAX_DIFF
                IMOM_AX_U=I
              END IF
            ELSE IF(I.LE.3*NX*IAT) THEN
              NMOM_AX_V=NMOM_AX_V+1
              DMOM_AX_SUM_V=DMOM_AX_SUM_V+DAX_DIFF
              IF(DAX_DIFF.GT.DMOM_AX_MAX_V) THEN
                DMOM_AX_MAX_V=DAX_DIFF
                IMOM_AX_V=I
              END IF
            ELSE
              NMOM_AX_W=NMOM_AX_W+1
              DMOM_AX_SUM_W=DMOM_AX_SUM_W+DAX_DIFF
              IF(DAX_DIFF.GT.DMOM_AX_MAX_W) THEN
                DMOM_AX_MAX_W=DAX_DIFF
                IMOM_AX_W=I
              END IF
            END IF
          END IF
        END DO
        WRITE(*,*)'MOMENTUM_FULL_JAC_AX_DIFF U COUNT/MAX/SUM=',NMOM_AX_U,DMOM_AX_MAX_U,DMOM_AX_SUM_U
        WRITE(*,*)'MOMENTUM_FULL_JAC_AX_DIFF V COUNT/MAX/SUM=',NMOM_AX_V,DMOM_AX_MAX_V,DMOM_AX_SUM_V
        WRITE(*,*)'MOMENTUM_FULL_JAC_AX_DIFF W COUNT/MAX/SUM=',NMOM_AX_W,DMOM_AX_MAX_W,DMOM_AX_SUM_W
        DO ITARGET=1,3
          IF(ITARGET.EQ.1) IMATRIX_RES_MAX=IMOM_AX_U
          IF(ITARGET.EQ.2) IMATRIX_RES_MAX=IMOM_AX_V
          IF(ITARGET.EQ.3) IMATRIX_RES_MAX=IMOM_AX_W
          IF(IMATRIX_RES_MAX.GT.0) THEN
            IBLOCK=(IMATRIX_RES_MAX-1)/(NX*IAT)+1
            ILOCAL=IMATRIX_RES_MAX-(IBLOCK-1)*NX*IAT
            IPNODE=(ILOCAL-1)/IAT+1
            IPTIME=ILOCAL-(IPNODE-1)*IAT
            WRITE(*,*)'MOMENTUM_FULL_JAC_AX_DIFF MAX EQBLOCK/ROW/NODE/TIME/XYZ=',IBLOCK,IMATRIX_RES_MAX,IPNODE,IPTIME,XX(IPNODE),YY(IPNODE),ZZ(IPNODE)
          END IF
        END DO
      END IF

      DO ITARGET=1,4
        IF(ITARGET.EQ.1) ITNODE=386
        IF(ITARGET.EQ.2) ITNODE=548
        IF(ITARGET.EQ.3) ITNODE=550
        IF(ITARGET.EQ.4) ITNODE=1786
        IF(ITNODE.GE.1.AND.ITNODE.LE.NX) THEN
          IROWCHK=GLOBAL_EQUATION_INDEX(1,ITNODE,1,NX,IAT)
          IF(IS_REPLACED_ROW(IROWCHK).NE.0) THEN
            DROW_DOT=0.0D0
            DROW_ABS=0.0D0
            NROW_NNZ=0
            DO I=1,SPARSE_PRESSURE_WALL_A%nnz
              IF(SPARSE_PRESSURE_WALL_A%row(I).EQ.IROWCHK) THEN
                ICOLCHK=SPARSE_PRESSURE_WALL_A%col(I)
                IBLOCK=(ICOLCHK-1)/(NX*IAT)+1
                ILOCAL=ICOLCHK-(IBLOCK-1)*NX*IAT
                IPNODE=(ILOCAL-1)/IAT+1
                IPTIME=ILOCAL-(IPNODE-1)*IAT
                DROW_VAR=0.0D0
                IF(IBLOCK.EQ.1) DROW_VAR=U(IPTIME,IPNODE)
                IF(IBLOCK.EQ.2) DROW_VAR=V(IPTIME,IPNODE)
                IF(IBLOCK.EQ.3) DROW_VAR=W(IPTIME,IPNODE)
                IF(IBLOCK.EQ.4) DROW_VAR=P(IPTIME,IPNODE)
                DROW_DOT=DROW_DOT+SPARSE_PRESSURE_WALL_A%val(I)*DROW_VAR
                DROW_ABS=DROW_ABS+DABS(SPARSE_PRESSURE_WALL_A%val(I))
                NROW_NNZ=NROW_NNZ+1
              END IF
            END DO
            WRITE(*,*)'WALL_PRESSURE_ROW_CHECK NODE/ROW/NNZ=',ITNODE,IROWCHK,NROW_NNZ
            WRITE(*,*)'WALL_PRESSURE_ROW_CHECK RHS/DOT/DIFF/ABSSUM=',FFS_REPLACED_ROW(IROWCHK),DROW_DOT,FFS_REPLACED_ROW(IROWCHK)-DROW_DOT,DROW_ABS
          END IF
        END IF
      END DO

!     --------------  APPLY DIRICHLET CONDITIONS (DIRECT CSR PATH) ---

!     Raw equation consistency before Dirichlet column elimination.
!     This is the correct place to check r = A*x - b for all assembled rows.
      IF(IPK.LE.1) THEN
      DRAW_AX_DIFF_MAX=0.0D0
      DRAW_AX_DIFF_SUM=0.0D0
      DRAW_AX_DIFF_L2=0.0D0
      IRAW_AX_MAX_ROW=0
      DO IBLOCK=1,4
        IF(IBLOCK.EQ.1) THEN
          DRAW_BLOCK_MAX1=0.0D0
          DRAW_BLOCK_SUM1=0.0D0
          DRAW_BLOCK_L21=0.0D0
        ELSE IF(IBLOCK.EQ.2) THEN
          DRAW_BLOCK_MAX2=0.0D0
          DRAW_BLOCK_SUM2=0.0D0
          DRAW_BLOCK_L22=0.0D0
        ELSE IF(IBLOCK.EQ.3) THEN
          DRAW_BLOCK_MAX3=0.0D0
          DRAW_BLOCK_SUM3=0.0D0
          DRAW_BLOCK_L23=0.0D0
        ELSE
          DRAW_BLOCK_MAX4=0.0D0
          DRAW_BLOCK_SUM4=0.0D0
          DRAW_BLOCK_L24=0.0D0
        END IF
      END DO
      DRAW_KIND_MAX0=0.0D0
      DRAW_KIND_MAX1=0.0D0
      DRAW_KIND_MAX2=0.0D0
      DRAW_KIND_SUM0=0.0D0
      DRAW_KIND_SUM1=0.0D0
      DRAW_KIND_SUM2=0.0D0
      NRAW_KIND0=0
      NRAW_KIND1=0
      NRAW_KIND2=0
      DO IROW=1,N
        DROW_DOT=0.0D0
        IF(IS_REPLACED_ROW(IROW).NE.0) THEN
          DO I=1,SPARSE_PRESSURE_WALL_A%nnz
            IF(SPARSE_PRESSURE_WALL_A%row(I).EQ.IROW) THEN
              ICOLCHK=SPARSE_PRESSURE_WALL_A%col(I)
              ILOCAL=(ICOLCHK-1)-(ICOLCHK-1)/(NX*IAT)*(NX*IAT)+1
              IPNODE=(ILOCAL-1)/IAT+1
              IPTIME=ILOCAL-(IPNODE-1)*IAT
              JBLOCK=(ICOLCHK-1)/(NX*IAT)+1
              IF(JBLOCK.EQ.1) DROW_DOT=DROW_DOT+SPARSE_PRESSURE_WALL_A%val(I)*U(IPTIME,IPNODE)
              IF(JBLOCK.EQ.2) DROW_DOT=DROW_DOT+SPARSE_PRESSURE_WALL_A%val(I)*V(IPTIME,IPNODE)
              IF(JBLOCK.EQ.3) DROW_DOT=DROW_DOT+SPARSE_PRESSURE_WALL_A%val(I)*W(IPTIME,IPNODE)
              IF(JBLOCK.EQ.4) DROW_DOT=DROW_DOT+SPARSE_PRESSURE_WALL_A%val(I)*P(IPTIME,IPNODE)
            END IF
          END DO
          DRAW_TARGET=FFS_REPLACED_ROW(IROW)
          IROW_KIND=1
        ELSE IF(IS_CONSTRAINED_ROW(IROW).NE.0) THEN
          DRAW_TARGET=FFS(IROW)
          DROW_DOT=FFS(IROW)
          IROW_KIND=2
        ELSE
          DO I=1,SPARSE_ELEMENT_A%nnz
            IF(SPARSE_ELEMENT_A%row(I).EQ.IROW) THEN
              ICOLCHK=SPARSE_ELEMENT_A%col(I)
              ILOCAL=(ICOLCHK-1)-(ICOLCHK-1)/(NX*IAT)*(NX*IAT)+1
              IPNODE=(ILOCAL-1)/IAT+1
              IPTIME=ILOCAL-(IPNODE-1)*IAT
              JBLOCK=(ICOLCHK-1)/(NX*IAT)+1
              IF(JBLOCK.EQ.1) DROW_DOT=DROW_DOT+SPARSE_ELEMENT_A%val(I)*U(IPTIME,IPNODE)
              IF(JBLOCK.EQ.2) DROW_DOT=DROW_DOT+SPARSE_ELEMENT_A%val(I)*V(IPTIME,IPNODE)
              IF(JBLOCK.EQ.3) DROW_DOT=DROW_DOT+SPARSE_ELEMENT_A%val(I)*W(IPTIME,IPNODE)
              IF(JBLOCK.EQ.4) DROW_DOT=DROW_DOT+SPARSE_ELEMENT_A%val(I)*P(IPTIME,IPNODE)
            END IF
          END DO
          DRAW_TARGET=FFS(IROW)
          IROW_KIND=0
        END IF
        DAX_DIFF=DABS(DROW_DOT-DRAW_TARGET)
        DRAW_AX_DIFF_SUM=DRAW_AX_DIFF_SUM+DAX_DIFF
        DRAW_AX_DIFF_L2=DRAW_AX_DIFF_L2+DAX_DIFF*DAX_DIFF
        IF(DAX_DIFF.GT.DRAW_AX_DIFF_MAX) THEN
          DRAW_AX_DIFF_MAX=DAX_DIFF
          IRAW_AX_MAX_ROW=IROW
          DRAW_AX_MAX_DOT=DROW_DOT
          DRAW_AX_MAX_FFS=DRAW_TARGET
        END IF
        IBLOCK=(IROW-1)/(NX*IAT)+1
        IF(IBLOCK.EQ.1) THEN
          DRAW_BLOCK_MAX1=DMAX1(DRAW_BLOCK_MAX1,DAX_DIFF)
          DRAW_BLOCK_SUM1=DRAW_BLOCK_SUM1+DAX_DIFF
          DRAW_BLOCK_L21=DRAW_BLOCK_L21+DAX_DIFF*DAX_DIFF
        ELSE IF(IBLOCK.EQ.2) THEN
          DRAW_BLOCK_MAX2=DMAX1(DRAW_BLOCK_MAX2,DAX_DIFF)
          DRAW_BLOCK_SUM2=DRAW_BLOCK_SUM2+DAX_DIFF
          DRAW_BLOCK_L22=DRAW_BLOCK_L22+DAX_DIFF*DAX_DIFF
        ELSE IF(IBLOCK.EQ.3) THEN
          DRAW_BLOCK_MAX3=DMAX1(DRAW_BLOCK_MAX3,DAX_DIFF)
          DRAW_BLOCK_SUM3=DRAW_BLOCK_SUM3+DAX_DIFF
          DRAW_BLOCK_L23=DRAW_BLOCK_L23+DAX_DIFF*DAX_DIFF
        ELSE
          DRAW_BLOCK_MAX4=DMAX1(DRAW_BLOCK_MAX4,DAX_DIFF)
          DRAW_BLOCK_SUM4=DRAW_BLOCK_SUM4+DAX_DIFF
          DRAW_BLOCK_L24=DRAW_BLOCK_L24+DAX_DIFF*DAX_DIFF
        END IF
        IF(IROW_KIND.EQ.0) THEN
          NRAW_KIND0=NRAW_KIND0+1
          DRAW_KIND_MAX0=DMAX1(DRAW_KIND_MAX0,DAX_DIFF)
          DRAW_KIND_SUM0=DRAW_KIND_SUM0+DAX_DIFF
        ELSE IF(IROW_KIND.EQ.1) THEN
          NRAW_KIND1=NRAW_KIND1+1
          DRAW_KIND_MAX1=DMAX1(DRAW_KIND_MAX1,DAX_DIFF)
          DRAW_KIND_SUM1=DRAW_KIND_SUM1+DAX_DIFF
        ELSE
          NRAW_KIND2=NRAW_KIND2+1
          DRAW_KIND_MAX2=DMAX1(DRAW_KIND_MAX2,DAX_DIFF)
          DRAW_KIND_SUM2=DRAW_KIND_SUM2+DAX_DIFF
        END IF
      END DO
      DRAW_AX_DIFF_L2=DSQRT(DRAW_AX_DIFF_L2)
      WRITE(*,*)'RAW_ALL_ROW_AX_CHECK_BEFORE_BC COUNT/MAX/SUM/L2=',N,DRAW_AX_DIFF_MAX,DRAW_AX_DIFF_SUM,DRAW_AX_DIFF_L2
      IF(IRAW_AX_MAX_ROW.NE.0) THEN
        ILOCAL=(IRAW_AX_MAX_ROW-1)-(IRAW_AX_MAX_ROW-1)/(NX*IAT)*(NX*IAT)+1
        IPNODE=(ILOCAL-1)/IAT+1
        IPTIME=ILOCAL-(IPNODE-1)*IAT
        IBLOCK=(IRAW_AX_MAX_ROW-1)/(NX*IAT)+1
        WRITE(*,*)'RAW_ALL_ROW_AX_MAX BLOCK/NODE/TIME/ROW/AX/FFS/DIFF=',IBLOCK,IPNODE,IPTIME,IRAW_AX_MAX_ROW,DRAW_AX_MAX_DOT,DRAW_AX_MAX_FFS,DRAW_AX_DIFF_MAX
      END IF
      WRITE(*,*)'RAW_ALL_ROW_AX_BLOCK_MAX CONT/V/W/U=',DRAW_BLOCK_MAX1,DRAW_BLOCK_MAX2,DRAW_BLOCK_MAX3,DRAW_BLOCK_MAX4
      WRITE(*,*)'RAW_ALL_ROW_AX_BLOCK_SUM CONT/V/W/U=',DRAW_BLOCK_SUM1,DRAW_BLOCK_SUM2,DRAW_BLOCK_SUM3,DRAW_BLOCK_SUM4
      WRITE(*,*)'RAW_ALL_ROW_AX_BLOCK_L2 CONT/V/W/U=',DSQRT(DRAW_BLOCK_L21),DSQRT(DRAW_BLOCK_L22),DSQRT(DRAW_BLOCK_L23),DSQRT(DRAW_BLOCK_L24)
      WRITE(*,*)'RAW_ALL_ROW_AX_KIND COUNT NORMAL/REPLACED/CONSTRAINED=',NRAW_KIND0,NRAW_KIND1,NRAW_KIND2
      WRITE(*,*)'RAW_ALL_ROW_AX_KIND_MAX NORMAL/REPLACED/CONSTRAINED=',DRAW_KIND_MAX0,DRAW_KIND_MAX1,DRAW_KIND_MAX2
      WRITE(*,*)'RAW_ALL_ROW_AX_KIND_SUM NORMAL/REPLACED/CONSTRAINED=',DRAW_KIND_SUM0,DRAW_KIND_SUM1,DRAW_KIND_SUM2
      END IF

      IF(ALLOCATED(ir)) DEALLOCATE(ir)
      IF(ALLOCATED(ic)) DEALLOCATE(ic)
      IF(ALLOCATED(a)) DEALLOCATE(a)
      call triplet_init(SPARSE_A, n, n, max(1, min(4000000, SPARSE_ELEMENT_A%nnz+SPARSE_PRESSURE_WALL_A%nnz+NCONSTRAINT)))
      DO ISP2=1,N
        IF(IS_CONSTRAINED(ISP2).NE.0) THEN
          BCVAL=BC_DELTA(ISP2)
          DO I=1,SPARSE_ELEMENT_A%nnz
            IF(SPARSE_ELEMENT_A%col(I).EQ.ISP2.AND.IS_REPLACED_ROW(SPARSE_ELEMENT_A%row(I)).EQ.0.AND.IS_CONSTRAINED_ROW(SPARSE_ELEMENT_A%row(I)).EQ.0) THEN
              FFS(SPARSE_ELEMENT_A%row(I))=FFS(SPARSE_ELEMENT_A%row(I))-SPARSE_ELEMENT_A%val(I)*BCVAL
            END IF
          END DO
          DO I=1,SPARSE_PRESSURE_WALL_A%nnz
            IF(SPARSE_PRESSURE_WALL_A%col(I).EQ.ISP2) THEN
              FFS_REPLACED_ROW(SPARSE_PRESSURE_WALL_A%row(I))=FFS_REPLACED_ROW(SPARSE_PRESSURE_WALL_A%row(I))-SPARSE_PRESSURE_WALL_A%val(I)*BCVAL
            END IF
          END DO
        END IF
      END DO
      DO I=1,SPARSE_ELEMENT_A%nnz
        IF(IS_REPLACED_ROW(SPARSE_ELEMENT_A%row(I)).EQ.0.AND.IS_CONSTRAINED_ROW(SPARSE_ELEMENT_A%row(I)).EQ.0.AND.IS_CONSTRAINED(SPARSE_ELEMENT_A%col(I)).EQ.0) THEN
          call triplet_add(SPARSE_A, SPARSE_ELEMENT_A%row(I), SPARSE_ELEMENT_A%col(I), SPARSE_ELEMENT_A%val(I))
        END IF
      END DO
      DO I=1,SPARSE_PRESSURE_WALL_A%nnz
        IF(IS_CONSTRAINED_ROW(SPARSE_PRESSURE_WALL_A%row(I)).EQ.0.AND.IS_CONSTRAINED(SPARSE_PRESSURE_WALL_A%col(I)).EQ.0) THEN
          call triplet_add(SPARSE_A, SPARSE_PRESSURE_WALL_A%row(I), SPARSE_PRESSURE_WALL_A%col(I), SPARSE_PRESSURE_WALL_A%val(I))
        END IF
      END DO
      DO I=1,N
        IF(IS_REPLACED_ROW(I).NE.0.AND.IS_CONSTRAINED_ROW(I).EQ.0) FFS(I)=FFS_REPLACED_ROW(I)
      END DO
      DO ISP2=1,N
        IF(IS_CONSTRAINED(ISP2).NE.0) THEN
          IROW=BC_ROW(ISP2)
          call triplet_add(SPARSE_A, IROW, ISP2, 1.0D0)
          FFS(IROW)=BC_DELTA(ISP2)
        END IF
      END DO
      call triplet_to_csr(SPARSE_A, ir, ic, a, mn)
      mnr=n+1
      write(*,*)'direct csr !mnr,mn',mnr,mn

!     Effective right hand side for final CSR diagnostics.  Most rows are
!     residual rows with b=0 after Dirichlet column elimination; constrained
!     rows keep their prescribed value.
      NORMAL_RHS=0.0D0
      DO IROW=1,N
        IF(IS_CONSTRAINED_ROW(IROW).NE.0) NORMAL_RHS(IROW)=FFS(IROW)
      END DO

!     Matrix-based residual path:
!     use the final CSR matrix itself to build the residual r = A*x - b
!     for every row.  This makes the stored FFS vector exactly consistent
!     with the matrix rows used by the linear solver.
      IF(IMATRIX_AXB_FINAL_RESIDUAL.NE.0) THEN
        DAX_DIFF=0.0D0
        DO IROW=1,N
          DROW_DOT=0.0D0
          DO I=ir(IROW),ir(IROW+1)-1
            ICOLCHK=ic(I)
            ILOCAL=(ICOLCHK-1)-(ICOLCHK-1)/(NX*IAT)*(NX*IAT)+1
            IPNODE=(ILOCAL-1)/IAT+1
            IPTIME=ILOCAL-(IPNODE-1)*IAT
            JBLOCK=(ICOLCHK-1)/(NX*IAT)+1
            IF(JBLOCK.EQ.1) DROW_DOT=DROW_DOT+a(I)*U(IPTIME,IPNODE)
            IF(JBLOCK.EQ.2) DROW_DOT=DROW_DOT+a(I)*V(IPTIME,IPNODE)
            IF(JBLOCK.EQ.3) DROW_DOT=DROW_DOT+a(I)*W(IPTIME,IPNODE)
            IF(JBLOCK.EQ.4) DROW_DOT=DROW_DOT+a(I)*P(IPTIME,IPNODE)
          END DO
          FFS(IROW)=DROW_DOT-NORMAL_RHS(IROW)
          DAX_DIFF=DMAX1(DAX_DIFF,DABS(FFS(IROW)))
        END DO
        FFS_PHYS=FFS
        WRITE(*,*)'FINAL CSR MATRIX-BASED RESIDUAL ENABLED MAXABS=',DAX_DIFF
      END IF

!     Final CSR diagnostics:
!     normal nonlinear rows report |f(x)|, while replaced linear rows report
!     the Ax-b consistency error.  Constrained rows are delta-equation rows
!     and are not checked against the current state vector here.
!     The Newton linear solve itself is checked later by NEWTON_PARDISO_REL_RESIDUAL.
      DALL_AX_DIFF_MAX=0.0D0
      DALL_AX_DIFF_SUM=0.0D0
      DALL_AX_DIFF_L2=0.0D0
      IALL_AX_MAX_ROW=0
      DO IBLOCK=1,4
        IF(IBLOCK.EQ.1) THEN
          DALL_BLOCK_MAX1=0.0D0
          DALL_BLOCK_SUM1=0.0D0
          DALL_BLOCK_L21=0.0D0
        ELSE IF(IBLOCK.EQ.2) THEN
          DALL_BLOCK_MAX2=0.0D0
          DALL_BLOCK_SUM2=0.0D0
          DALL_BLOCK_L22=0.0D0
        ELSE IF(IBLOCK.EQ.3) THEN
          DALL_BLOCK_MAX3=0.0D0
          DALL_BLOCK_SUM3=0.0D0
          DALL_BLOCK_L23=0.0D0
        ELSE
          DALL_BLOCK_MAX4=0.0D0
          DALL_BLOCK_SUM4=0.0D0
          DALL_BLOCK_L24=0.0D0
        END IF
      END DO
      DALL_KIND_MAX0=0.0D0
      DALL_KIND_MAX1=0.0D0
      DALL_KIND_MAX2=0.0D0
      DALL_KIND_SUM0=0.0D0
      DALL_KIND_SUM1=0.0D0
      DALL_KIND_SUM2=0.0D0
      NALL_KIND0=0
      NALL_KIND1=0
      NALL_KIND2=0
      DO IROW=1,N
        DROW_DOT=0.0D0
        DO I=ir(IROW),ir(IROW+1)-1
          ICOLCHK=ic(I)
          ILOCAL=(ICOLCHK-1)-(ICOLCHK-1)/(NX*IAT)*(NX*IAT)+1
          IPNODE=(ILOCAL-1)/IAT+1
          IPTIME=ILOCAL-(IPNODE-1)*IAT
          JBLOCK=(ICOLCHK-1)/(NX*IAT)+1
          IF(JBLOCK.EQ.1) DROW_DOT=DROW_DOT+a(I)*U(IPTIME,IPNODE)
          IF(JBLOCK.EQ.2) DROW_DOT=DROW_DOT+a(I)*V(IPTIME,IPNODE)
          IF(JBLOCK.EQ.3) DROW_DOT=DROW_DOT+a(I)*W(IPTIME,IPNODE)
          IF(JBLOCK.EQ.4) DROW_DOT=DROW_DOT+a(I)*P(IPTIME,IPNODE)
        END DO
        IF(IS_CONSTRAINED_ROW(IROW).NE.0) THEN
          DRHS_EFFECTIVE=NORMAL_RHS(IROW)
          IROW_KIND=2
        ELSE
          DRHS_EFFECTIVE=0.0D0
          IF(IS_REPLACED_ROW(IROW).NE.0) THEN
            IROW_KIND=1
          ELSE
            IROW_KIND=0
          END IF
        END IF
        IF(IMATRIX_AXB_FINAL_RESIDUAL.NE.0) THEN
          DAX_DIFF=DABS((DROW_DOT-DRHS_EFFECTIVE)-FFS(IROW))
        ELSE IF(IROW_KIND.EQ.0) THEN
          DAX_DIFF=DABS(FFS(IROW))
          DROW_DOT=0.0D0
          DRHS_EFFECTIVE=0.0D0
        ELSE IF(IROW_KIND.EQ.1) THEN
          DAX_DIFF=DABS(DROW_DOT-FFS(IROW))
        ELSE
          DAX_DIFF=0.0D0
        END IF
        DALL_AX_DIFF_SUM=DALL_AX_DIFF_SUM+DAX_DIFF
        DALL_AX_DIFF_L2=DALL_AX_DIFF_L2+DAX_DIFF*DAX_DIFF
        IF(DAX_DIFF.GT.DALL_AX_DIFF_MAX) THEN
          DALL_AX_DIFF_MAX=DAX_DIFF
          IALL_AX_MAX_ROW=IROW
          DALL_AX_MAX_DOT=DROW_DOT
          DALL_AX_MAX_RHS=DRHS_EFFECTIVE
          DALL_AX_MAX_FFS=FFS(IROW)
        END IF
        IBLOCK=(IROW-1)/(NX*IAT)+1
        IF(IBLOCK.EQ.1) THEN
          DALL_BLOCK_MAX1=DMAX1(DALL_BLOCK_MAX1,DAX_DIFF)
          DALL_BLOCK_SUM1=DALL_BLOCK_SUM1+DAX_DIFF
          DALL_BLOCK_L21=DALL_BLOCK_L21+DAX_DIFF*DAX_DIFF
        ELSE IF(IBLOCK.EQ.2) THEN
          DALL_BLOCK_MAX2=DMAX1(DALL_BLOCK_MAX2,DAX_DIFF)
          DALL_BLOCK_SUM2=DALL_BLOCK_SUM2+DAX_DIFF
          DALL_BLOCK_L22=DALL_BLOCK_L22+DAX_DIFF*DAX_DIFF
        ELSE IF(IBLOCK.EQ.3) THEN
          DALL_BLOCK_MAX3=DMAX1(DALL_BLOCK_MAX3,DAX_DIFF)
          DALL_BLOCK_SUM3=DALL_BLOCK_SUM3+DAX_DIFF
          DALL_BLOCK_L23=DALL_BLOCK_L23+DAX_DIFF*DAX_DIFF
        ELSE
          DALL_BLOCK_MAX4=DMAX1(DALL_BLOCK_MAX4,DAX_DIFF)
          DALL_BLOCK_SUM4=DALL_BLOCK_SUM4+DAX_DIFF
          DALL_BLOCK_L24=DALL_BLOCK_L24+DAX_DIFF*DAX_DIFF
        END IF
        IF(IROW_KIND.EQ.0) THEN
          NALL_KIND0=NALL_KIND0+1
          DALL_KIND_MAX0=DMAX1(DALL_KIND_MAX0,DAX_DIFF)
          DALL_KIND_SUM0=DALL_KIND_SUM0+DAX_DIFF
        ELSE IF(IROW_KIND.EQ.1) THEN
          NALL_KIND1=NALL_KIND1+1
          DALL_KIND_MAX1=DMAX1(DALL_KIND_MAX1,DAX_DIFF)
          DALL_KIND_SUM1=DALL_KIND_SUM1+DAX_DIFF
        ELSE
          NALL_KIND2=NALL_KIND2+1
          DALL_KIND_MAX2=DMAX1(DALL_KIND_MAX2,DAX_DIFF)
          DALL_KIND_SUM2=DALL_KIND_SUM2+DAX_DIFF
        END IF
      END DO
      DALL_AX_DIFF_L2=DSQRT(DALL_AX_DIFF_L2)
      WRITE(*,*)'FINAL_CSR_MIXED_DIAG COUNT/MAX/SUM/L2=',N,DALL_AX_DIFF_MAX,DALL_AX_DIFF_SUM,DALL_AX_DIFF_L2
      IF(IALL_AX_MAX_ROW.NE.0) THEN
        ILOCAL=(IALL_AX_MAX_ROW-1)-(IALL_AX_MAX_ROW-1)/(NX*IAT)*(NX*IAT)+1
        IPNODE=(ILOCAL-1)/IAT+1
        IPTIME=ILOCAL-(IPNODE-1)*IAT
        IBLOCK=(IALL_AX_MAX_ROW-1)/(NX*IAT)+1
        WRITE(*,*)'FINAL_CSR_MIXED_MAX BLOCK/NODE/TIME/ROW/AX/B/FFS/DIFF=',IBLOCK,IPNODE,IPTIME,IALL_AX_MAX_ROW,DALL_AX_MAX_DOT,DALL_AX_MAX_RHS,DALL_AX_MAX_FFS,DALL_AX_DIFF_MAX
      END IF
      WRITE(*,*)'FINAL_CSR_MIXED_BLOCK_MAX CONT/V/W/U=',DALL_BLOCK_MAX1,DALL_BLOCK_MAX2,DALL_BLOCK_MAX3,DALL_BLOCK_MAX4
      WRITE(*,*)'FINAL_CSR_MIXED_BLOCK_SUM CONT/V/W/U=',DALL_BLOCK_SUM1,DALL_BLOCK_SUM2,DALL_BLOCK_SUM3,DALL_BLOCK_SUM4
      WRITE(*,*)'FINAL_CSR_MIXED_BLOCK_L2 CONT/V/W/U=',DSQRT(DALL_BLOCK_L21),DSQRT(DALL_BLOCK_L22),DSQRT(DALL_BLOCK_L23),DSQRT(DALL_BLOCK_L24)
      WRITE(*,*)'FINAL_CSR_ROW_KIND COUNT NORMAL/REPLACED/CONSTRAINED=',NALL_KIND0,NALL_KIND1,NALL_KIND2
      WRITE(*,*)'FINAL_CSR_NORMAL_PHYS_RESIDUAL MAX/SUM=',DALL_KIND_MAX0,DALL_KIND_SUM0
      WRITE(*,*)'FINAL_CSR_AXB_CONSISTENCY MAX REPLACED=',DALL_KIND_MAX1
      WRITE(*,*)'FINAL_CSR_AXB_CONSISTENCY SUM REPLACED=',DALL_KIND_SUM1
      WRITE(*,*)'FINAL_CSR_CONSTRAINED_ROWS SKIPPED_FOR_CURRENT_STATE=',NALL_KIND2

!     Check consistency of outflow velocity rows after CSR assembly.
      DOUT_DIFF_MAX=0.0D0
      DOUT_DIFF_SUM=0.0D0
      NOUT_CHECK=0
      IOUT_MAX_ROW=0
      IOUT_MAX_BLOCK=0
      IOUT_MAX_NODE=0
      DOUT_MAX_DOT=0.0D0
      DOUT_MAX_FFS=0.0D0
      DO I3=1,NSB3
        II=NSBUO(I3)
        DO ITT=1,IAT
          DO IBLOCK=2,4
            IROWCHK=GLOBAL_EQUATION_INDEX(ITT,II,IBLOCK,NX,IAT)
            IF(IS_CONSTRAINED_ROW(IROWCHK).EQ.0) THEN
              DROW_DOT=0.0D0
              DO I=ir(IROWCHK),ir(IROWCHK+1)-1
                ICOLCHK=ic(I)
                ILOCAL=(ICOLCHK-1)-(ICOLCHK-1)/(NX*IAT)*(NX*IAT)+1
                IPNODE=(ILOCAL-1)/IAT+1
                IPTIME=ILOCAL-(IPNODE-1)*IAT
                JBLOCK=(ICOLCHK-1)/(NX*IAT)+1
                IF(JBLOCK.EQ.1) DROW_DOT=DROW_DOT+a(I)*U(IPTIME,IPNODE)
                IF(JBLOCK.EQ.2) DROW_DOT=DROW_DOT+a(I)*V(IPTIME,IPNODE)
                IF(JBLOCK.EQ.3) DROW_DOT=DROW_DOT+a(I)*W(IPTIME,IPNODE)
                IF(JBLOCK.EQ.4) DROW_DOT=DROW_DOT+a(I)*P(IPTIME,IPNODE)
              END DO
              DOUT_DIFF=DABS(FFS(IROWCHK)-DROW_DOT)
              DOUT_DIFF_SUM=DOUT_DIFF_SUM+DOUT_DIFF
              NOUT_CHECK=NOUT_CHECK+1
              IF(DOUT_DIFF.GT.DOUT_DIFF_MAX) THEN
                DOUT_DIFF_MAX=DOUT_DIFF
                IOUT_MAX_ROW=IROWCHK
                IOUT_MAX_BLOCK=IBLOCK
                IOUT_MAX_NODE=II
                DOUT_MAX_DOT=DROW_DOT
                DOUT_MAX_FFS=FFS(IROWCHK)
              END IF
              IF(I3.LE.2.AND.ITT.EQ.1) WRITE(*,*)'OUTFLOW_ROW_AX_SAMPLE BLOCK/NODE/ROW/FFS/DOT/DIFF=',IBLOCK,II,IROWCHK,FFS(IROWCHK),DROW_DOT,FFS(IROWCHK)-DROW_DOT
            END IF
          END DO
        END DO
      END DO
      WRITE(*,*)'OUTFLOW_ROW_AX_CHECK COUNT/MAXDIFF/SUMDIFF=',NOUT_CHECK,DOUT_DIFF_MAX,DOUT_DIFF_SUM
      IF(IOUT_MAX_ROW.NE.0) WRITE(*,*)'OUTFLOW_ROW_AX_MAX BLOCK/NODE/ROW/FFS/DOT/DIFF=',IOUT_MAX_BLOCK,IOUT_MAX_NODE,IOUT_MAX_ROW,DOUT_MAX_FFS,DOUT_MAX_DOT,DOUT_MAX_FFS-DOUT_MAX_DOT
      DAXB_TOL=1.0D-6
      IF(DALL_KIND_MAX1.LE.DAXB_TOL.AND.DOUT_DIFF_MAX.LE.DAXB_TOL) THEN
        WRITE(*,*)'LINEAR_AXB_SUMMARY PASS TOL/MAX_REPLACED/MAX_OUTFLOW=',DAXB_TOL,DALL_KIND_MAX1,DOUT_DIFF_MAX
      ELSE
        WRITE(*,*)'LINEAR_AXB_SUMMARY FAIL TOL/MAX_REPLACED/MAX_OUTFLOW=',DAXB_TOL,DALL_KIND_MAX1,DOUT_DIFF_MAX
      END IF

!     --------------  FINITE-DIFFERENCE CHECK OF PRESSURE JACOBIAN COLUMN
      IF(IPRESSURE_FD_CHECK.NE.0.AND.IPRESSURE_FD_DONE.EQ.0) THEN
        IF(IPRESSURE_FD_MODE.EQ.0) THEN
          IPRESSURE_FD_COL=3*NX*IAT+1
          DMAX_FD=0.0D0
          DO I=3*NX*IAT+1,4*NX*IAT
            IF(DABS(FFS(I)).GT.DMAX_FD) THEN
              DMAX_FD=DABS(FFS(I))
              IPRESSURE_FD_COL=I
            END IF
          END DO
          ILOCAL=IPRESSURE_FD_COL-3*NX*IAT
          IPRESSURE_FD_NODE=(ILOCAL-1)/IAT+1
          IPRESSURE_FD_TIME=ILOCAL-(IPRESSURE_FD_NODE-1)*IAT
          IPRESSURE_FD_COL=GLOBAL_DOF_INDEX(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE,IPRESSURE_FD_TARGET_BLOCK,NX,IAT)
          FFS_FD_BASE=FFS
          JCOL_FD_BASE=0.0D0
          DO I=1,N
            DO J=ir(I),ir(I+1)-1
              IF(ic(J).EQ.IPRESSURE_FD_COL) JCOL_FD_BASE(I)=JCOL_FD_BASE(I)+a(J)
            END DO
          END DO
          IF(IPRESSURE_FD_TARGET_BLOCK.EQ.1) THEN
            DPRESSURE_FD_OLD=U(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)
            U(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD+DPRESSURE_FD_EPS
          ELSE IF(IPRESSURE_FD_TARGET_BLOCK.EQ.2) THEN
            DPRESSURE_FD_OLD=V(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)
            V(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD+DPRESSURE_FD_EPS
          ELSE IF(IPRESSURE_FD_TARGET_BLOCK.EQ.3) THEN
            DPRESSURE_FD_OLD=W(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)
            W(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD+DPRESSURE_FD_EPS
          ELSE
            DPRESSURE_FD_OLD=P(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)
            P(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD+DPRESSURE_FD_EPS
          END IF
          IPRESSURE_FD_MODE=1
          WRITE(*,*)'PRESSURE FD JACOBIAN CHECK START BLOCK/COL/NODE/TIME/EPS=',IPRESSURE_FD_TARGET_BLOCK,IPRESSURE_FD_COL,IPRESSURE_FD_NODE,IPRESSURE_FD_TIME,DPRESSURE_FD_EPS
          GO TO 30
        ELSE IF(IPRESSURE_FD_MODE.EQ.1) THEN
          OPEN(98,FILE='results\pressure_jacobian_fd_check.dat',STATUS='REPLACE')
          WRITE(98,'(A)') '# row block node time analytic fd diff absdiff rhs_base rhs_pert'
          DMAX_FD=0.0D0
          DSUM_FD=0.0D0
          DMAX_FD_REL=0.0D0
          IMAX_FD=0
          DO I=1,N
            DFD=(FFS(I)-FFS_FD_BASE(I))/DPRESSURE_FD_EPS
            JCOL_FD_DIFF(I)=DFD-JCOL_FD_BASE(I)
            DABSFD=DABS(JCOL_FD_DIFF(I))
            DSUM_FD=DSUM_FD+DABSFD
            DDEN=DMAX1(1.0D-30,DABS(DFD),DABS(JCOL_FD_BASE(I)))
            DMAX_FD_REL=DMAX1(DMAX_FD_REL,DABSFD/DDEN)
            IF(DABSFD.GT.DMAX_FD) THEN
              DMAX_FD=DABSFD
              IMAX_FD=I
            END IF
          END DO
          DO IPRINT=1,40
            DMAX_PRINT=-1.0D0
            IMAX_PRINT=0
            DO I=1,N
              IF(DABS(JCOL_FD_DIFF(I)).GT.DMAX_PRINT) THEN
                ISEEN=0
                IF(ISEEN.EQ.0) THEN
                  DMAX_PRINT=DABS(JCOL_FD_DIFF(I))
                  IMAX_PRINT=I
                END IF
              END IF
            END DO
            IF(IMAX_PRINT.LE.0) EXIT
            IBLOCK=(IMAX_PRINT-1)/(NX*IAT)+1
            ILOCAL=IMAX_PRINT-(IBLOCK-1)*NX*IAT
            IPNODE=(ILOCAL-1)/IAT+1
            IPTIME=ILOCAL-(IPNODE-1)*IAT
            DFD=(FFS(IMAX_PRINT)-FFS_FD_BASE(IMAX_PRINT))/DPRESSURE_FD_EPS
            WRITE(98,'(4I10,1P,6E24.14)') IMAX_PRINT,IBLOCK,IPNODE,IPTIME,JCOL_FD_BASE(IMAX_PRINT),DFD,JCOL_FD_DIFF(IMAX_PRINT),&
                 DABS(JCOL_FD_DIFF(IMAX_PRINT)),FFS_FD_BASE(IMAX_PRINT),FFS(IMAX_PRINT)
            JCOL_FD_DIFF(IMAX_PRINT)=0.0D0
          END DO
          CLOSE(98)
          IBLOCK=(IMAX_FD-1)/(NX*IAT)+1
          ILOCAL=IMAX_FD-(IBLOCK-1)*NX*IAT
          IPNODE=(ILOCAL-1)/IAT+1
          IPTIME=ILOCAL-(IPNODE-1)*IAT
          WRITE(*,*)'PRESSURE FD JACOBIAN CHECK DONE BLOCK/COL/NODE/TIME=',IPRESSURE_FD_TARGET_BLOCK,IPRESSURE_FD_COL,IPRESSURE_FD_NODE,IPRESSURE_FD_TIME
          WRITE(*,*)'PRESSURE FD MAX ABS DIFF ROW/BLOCK/NODE/TIME=',IMAX_FD,IBLOCK,IPNODE,IPTIME
          WRITE(*,*)'PRESSURE FD MAX_ABS_DIFF SUM_ABS_DIFF MAX_REL_DIFF=',DMAX_FD,DSUM_FD,DMAX_FD_REL
          WRITE(*,*)'PRESSURE FD DETAILS FILE=results\pressure_jacobian_fd_check.dat'
          IF(IPRESSURE_FD_TARGET_BLOCK.EQ.1) THEN
            U(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD
          ELSE IF(IPRESSURE_FD_TARGET_BLOCK.EQ.2) THEN
            V(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD
          ELSE IF(IPRESSURE_FD_TARGET_BLOCK.EQ.3) THEN
            W(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD
          ELSE
            P(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD
          END IF
          IPRESSURE_FD_DONE=1
          STOP
        END IF
      END IF

      GO TO 125

           
!     --------------  BOUNDARY CONDITION AT WALL --------------

        DO J=1,NSB1
        DO ITT=1,IAT
          JJ=NSBU(J)
          ISP2=ITT+(JJ-1)*IAT
          BCVAL=U(ITT,JJ)-SBU(J)
          DO IS1=1,NX*4
          DO JTT=1,IAT
            ISP1=JTT+(IS1-1)*IAT
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*BCVAL
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=BCVAL
        END DO
      END DO
      
      DO J=1,NSB1
        DO ITT=1,IAT
          JJ=NSBU(J)
          ISP2=ITT+(JJ-1)*IAT+NX*IAT
          BCVAL=V(ITT,JJ)-SBV(J)
          DO IS1=1,NX*4
          DO JTT=1,IAT
            ISP1=JTT+(IS1-1)*IAT
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*BCVAL
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=BCVAL
        END DO
      END DO

      DO J=1,NSB1
        DO ITT=1,IAT
          JJ=NSBU(J)
          ISP2=ITT+(JJ-1)*IAT+2*NX*IAT
          BCVAL=W(ITT,JJ)-SBW(J)
          DO IS1=1,NX*4
          DO JTT=1,IAT
            ISP1=JTT+(IS1-1)*IAT
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*BCVAL
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=BCVAL
        END DO
      END DO
      
      
! -------- PRESSURE AT OUTOFLOW CONDITION 
    J1=5
     DO I3=1,NSB3
       DO ITT=1,IAT      
         II=NSBUO(I3)
         ISP2=ITT+(II-1)*IAT+3*NX*IAT
         BCVAL=P(ITT,II)-SBPO(I3)
         DO IS1=1,NX*4
         DO JTT=1,IAT
           ISP1=JTT+(IS1-1)*IAT
           FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*BCVAL
           SSMMK(ISP1,ISP2)=0.0D0
           SSMMK(ISP2,ISP1)=0.0D0
         END DO
         END DO
         SSMMK(ISP2,ISP2)=1.0D0
         FFS(ISP2)=BCVAL
       END DO
     END DO
       

!     ------------ INFLOW BOUNDARY CONDITION ------------------   
      DO J=1,NSB2
        DO ITT=1,IAT
          JJ=NSBUI(J)
          ISP2=ITT+(JJ-1)*IAT
          BCVAL=U(ITT,JJ)-SBUI(J)
          DO IS1=1,NX*4
          DO JTT=1,IAT
            ISP1=JTT+(IS1-1)*IAT
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*BCVAL
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=BCVAL
        END DO
      END DO

      DO J=1,NSB2
        DO ITT=1,IAT
          JJ=NSBUI(J)
          ISP2=ITT+(JJ-1)*IAT+NX*IAT
          BCVAL=V(ITT,JJ)-SBVI(J)
          DO IS1=1,NX*4
          DO JTT=1,IAT
            ISP1=JTT+(IS1-1)*IAT
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*BCVAL
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=BCVAL
        END DO
      END DO
      
      DO J=1,NSB2
        DO ITT=1,IAT
          JJ=NSBUI(J)
          ISP2=ITT+(JJ-1)*IAT+2*NX*IAT
          BCVAL=W(ITT,JJ)-SBWI(J)
          DO IS1=1,NX*4
          DO JTT=1,IAT
            ISP1=JTT+(IS1-1)*IAT
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*BCVAL
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=BCVAL
        END DO
      END DO

  125 CONTINUE
      WRITE(*,*)'CLEAR 103 !'

      WW4_SUM=0.0D0
      WW4_RHS_SUM=0.0D0
!     !$OMP parallel private(I)
!     !$OMP do reduction(+ : WW4_SUM,WW4_RHS_SUM)
      DO I=1,4*NX*IAT
        WW4_SUM=WW4_SUM+FFS_PHYS(I)**2
        WW4_RHS_SUM=WW4_RHS_SUM+FFS(I)**2
      END DO
!     !$OMP end do
!     !$OMP end parallel
      WW4=DSQRT(WW4_SUM)
      WW4_LINEAR_RHS=DSQRT(WW4_RHS_SUM)
      
      WRITE(*,*)'WW4-',WW4
      WRITE(*,*)'WW4_SUM=',WW4_SUM
      WRITE(*,*)'WW4_NORM=',WW4
      WRITE(*,*)'WW4_PHYS_NORM=',WW4
      WRITE(*,*)'WW4_LINEAR_RHS_NORM=',WW4_LINEAR_RHS
      DRU=0.0D0
      DRV=0.0D0
      DRW=0.0D0
      DRP=0.0D0
      DUMAX=0.0D0
      DVMAX=0.0D0
      DWMAX=0.0D0
      DPMAX=0.0D0
      IUMAX=0
      IVMAX=0
      IWMAX=0
      IPMAX=0
      DO I=1,NX*IAT
        DRU=DRU+FFS_PHYS(I)**2
        DRV=DRV+FFS_PHYS(I+NX*IAT)**2
        DRW=DRW+FFS_PHYS(I+2*NX*IAT)**2
        DRP=DRP+FFS_PHYS(I+3*NX*IAT)**2
        IF(DABS(FFS_PHYS(I)).GT.DUMAX) THEN
          DUMAX=DABS(FFS_PHYS(I))
          IUMAX=I
        END IF
        IF(DABS(FFS_PHYS(I+NX*IAT)).GT.DVMAX) THEN
          DVMAX=DABS(FFS_PHYS(I+NX*IAT))
          IVMAX=I+NX*IAT
        END IF
        IF(DABS(FFS_PHYS(I+2*NX*IAT)).GT.DWMAX) THEN
          DWMAX=DABS(FFS_PHYS(I+2*NX*IAT))
          IWMAX=I+2*NX*IAT
        END IF
        IF(DABS(FFS_PHYS(I+3*NX*IAT)).GT.DPMAX) THEN
          DPMAX=DABS(FFS_PHYS(I+3*NX*IAT))
          IPMAX=I+3*NX*IAT
        END IF
      END DO
      DRESMAX_CONT=DUMAX
      DRESMAX_V=DVMAX
      DRESMAX_W=DWMAX
      DRESMAX_U=DPMAX
      WRITE(*,*)'RESIDUAL_ROW_NORM CONT/V-MOM/W-MOM/U-MOM=',DSQRT(DRU),DSQRT(DRV),DSQRT(DRW),DSQRT(DRP)
      DO IBLOCK=1,4
        IF(IBLOCK.EQ.1) ILOCAL=IUMAX
        IF(IBLOCK.EQ.2) ILOCAL=IVMAX-NX*IAT
        IF(IBLOCK.EQ.3) ILOCAL=IWMAX-2*NX*IAT
        IF(IBLOCK.EQ.4) ILOCAL=IPMAX-3*NX*IAT
        IPNODE=(ILOCAL-1)/IAT+1
        IPTIME=ILOCAL-(IPNODE-1)*IAT
        IF(IBLOCK.EQ.1) WRITE(*,*)'MAX_RESIDUAL CONT ROW/NODE/TIME/VAL=',IUMAX,IPNODE,IPTIME,FFS_PHYS(IUMAX)
        IF(IBLOCK.EQ.2) WRITE(*,*)'MAX_RESIDUAL V-MOM ROW/NODE/TIME/VAL=',IVMAX,IPNODE,IPTIME,FFS_PHYS(IVMAX)
        IF(IBLOCK.EQ.3) WRITE(*,*)'MAX_RESIDUAL W-MOM ROW/NODE/TIME/VAL=',IWMAX,IPNODE,IPTIME,FFS_PHYS(IWMAX)
        IF(IBLOCK.EQ.4) WRITE(*,*)'MAX_RESIDUAL U-MOM ROW/NODE/TIME/VAL=',IPMAX,IPNODE,IPTIME,FFS_PHYS(IPMAX)
        ION_WALL=0
        ION_INLET=0
        ION_POUT=0
        DO I3=1,NSB1
          IF(NSBU(I3).EQ.IPNODE) ION_WALL=1
        END DO
        DO I3=1,NSB2
          IF(NSBUI(I3).EQ.IPNODE) ION_INLET=1
        END DO
        DO I3=1,NSB3
          IF(NSBUO(I3).EQ.IPNODE) ION_POUT=1
        END DO
        WRITE(*,*)'MAX_RESIDUAL_NODE_INFO EQBLOCK/NODE/XYZ/ROWCONSTR/WALL/INLET/POUT=',IBLOCK,IPNODE,XX(IPNODE),YY(IPNODE),ZZ(IPNODE),IS_CONSTRAINED_ROW(ILOCAL+NX*IAT*(IBLOCK-1)),ION_WALL,ION_INLET,ION_POUT
      END DO

      IF(IPENDING_STEP.NE.0.AND.IDIAG_TRIAL_RESIDUAL.NE.0) THEN
        DO IBLOCK=1,4
          DROWMAX=-1.0D300
          IOUT_MAX_ROW=0
          DROWSUM=0.0D0
          DO I=1,NX*IAT
            IROW=I+(IBLOCK-1)*NX*IAT
            DINCR=DABS(FFS_PHYS(IROW))-DABS(FFS_PHYS_BASE(IROW))
            IF(DINCR.GT.0.0D0) DROWSUM=DROWSUM+DINCR*DINCR
            IF(DINCR.GT.DROWMAX) THEN
              DROWMAX=DINCR
              IOUT_MAX_ROW=IROW
            END IF
          END DO
          ILOCAL=IOUT_MAX_ROW-(IBLOCK-1)*NX*IAT
          IPNODE=(ILOCAL-1)/IAT+1
          IPTIME=ILOCAL-(IPNODE-1)*IAT
          ION_WALL=0
          ION_INLET=0
          ION_POUT=0
          DO I3=1,NSB1
            IF(NSBU(I3).EQ.IPNODE) ION_WALL=1
          END DO
          DO I3=1,NSB2
            IF(NSBUI(I3).EQ.IPNODE) ION_INLET=1
          END DO
          DO I3=1,NSB3
            IF(NSBUO(I3).EQ.IPNODE) ION_POUT=1
          END DO
          WRITE(*,*)'TRIAL_RESIDUAL_INCREASE BLOCK/MAXINC/SQRTSUMINC/ROW/NODE/TIME=',IBLOCK,DROWMAX,DSQRT(DROWSUM),IOUT_MAX_ROW,IPNODE,IPTIME
          WRITE(*,*)'TRIAL_RESIDUAL_INCREASE OLD/NEW/XYZ/WALL/INLET/POUT=',FFS_PHYS_BASE(IOUT_MAX_ROW),FFS_PHYS(IOUT_MAX_ROW),XX(IPNODE),YY(IPNODE),ZZ(IPNODE),ION_WALL,ION_INLET,ION_POUT
          IUROW=IPTIME+(IPNODE-1)*IAT
          IVROW=IUROW+NX*IAT
          IWROW=IUROW+2*NX*IAT
          IPROW=IUROW+3*NX*IAT
          WRITE(*,*)'TRIAL_RESIDUAL_INCREASE UVWP_SAVE=',U_SAVE(IPTIME,IPNODE),V_SAVE(IPTIME,IPNODE),W_SAVE(IPTIME,IPNODE),P_SAVE(IPTIME,IPNODE)
          WRITE(*,*)'TRIAL_RESIDUAL_INCREASE UVWP_TRIAL=',U(IPTIME,IPNODE),V(IPTIME,IPNODE),W(IPTIME,IPNODE),P(IPTIME,IPNODE)
          WRITE(*,*)'TRIAL_RESIDUAL_INCREASE RAW_DUVWP=',DDUVWP(IUROW),DDUVWP(IVROW),DDUVWP(IWROW),DDUVWP(IPROW)
        END DO
      END IF

      IF(IPENDING_STEP.NE.0) THEN
        IF(ILINESEARCH_EXTRA_MODE.EQ.11) THEN
          DO IBLOCK=1,4
            DS_OLD=0.0D0
            DS_NEW=0.0D0
            DROWMAX=0.0D0
            IOUT_MAX_ROW=0
            DO I=1,NX*IAT
              IROW=I+(IBLOCK-1)*NX*IAT
              DS_OLD=DS_OLD+FFS_PHYS_BASE(IROW)*FFS_PHYS_BASE(IROW)
              DS_NEW=DS_NEW+FFS_PHYS(IROW)*FFS_PHYS(IROW)
              IF(DABS(FFS_PHYS(IROW))-DABS(FFS_PHYS_BASE(IROW)).GT.DROWMAX) THEN
                DROWMAX=DABS(FFS_PHYS(IROW))-DABS(FFS_PHYS_BASE(IROW))
                IOUT_MAX_ROW=IROW
              END IF
            END DO
            IF(IOUT_MAX_ROW.GT.0) THEN
              ILOCAL=IOUT_MAX_ROW-(IBLOCK-1)*NX*IAT
              IPNODE=(ILOCAL-1)/IAT+1
              IPTIME=ILOCAL-(IPNODE-1)*IAT
            ELSE
              IPNODE=0
              IPTIME=0
            END IF
            WRITE(*,*)'STAG_SD_BLOCK_TRIAL OLD/NEW/DIFF BLOCK/MAXINC/ROW/NODE/TIME=',DSQRT(DS_OLD),DSQRT(DS_NEW),DSQRT(DS_NEW)-DSQRT(DS_OLD),IBLOCK,DROWMAX,IOUT_MAX_ROW,IPNODE,IPTIME
          END DO
        END IF
        IF(WW4.GE.WW4_PREV) THEN
          DO IBLOCK=1,4
            DS_OLD=0.0D0
            DS_NEW=0.0D0
            DROWMAX=0.0D0
            IOUT_MAX_ROW=0
            DO I=1,NX*IAT
              IROW=I+(IBLOCK-1)*NX*IAT
              DS_OLD=DS_OLD+FFS_PHYS_BASE(IROW)*FFS_PHYS_BASE(IROW)
              DS_NEW=DS_NEW+FFS_PHYS(IROW)*FFS_PHYS(IROW)
              IF(DABS(FFS_PHYS(IROW))-DABS(FFS_PHYS_BASE(IROW)).GT.DROWMAX) THEN
                DROWMAX=DABS(FFS_PHYS(IROW))-DABS(FFS_PHYS_BASE(IROW))
                IOUT_MAX_ROW=IROW
              END IF
            END DO
            IF(IOUT_MAX_ROW.GT.0) THEN
              ILOCAL=IOUT_MAX_ROW-(IBLOCK-1)*NX*IAT
              IPNODE=(ILOCAL-1)/IAT+1
              IPTIME=ILOCAL-(IPNODE-1)*IAT
            ELSE
              IPNODE=0
              IPTIME=0
            END IF
            WRITE(*,*)'TRIAL_BLOCK_PHYS_NORM OLD/NEW/DIFF BLOCK/MAXINC/ROW/NODE/TIME=',DSQRT(DS_OLD),DSQRT(DS_NEW),DSQRT(DS_NEW)-DSQRT(DS_OLD),IBLOCK,DROWMAX,IOUT_MAX_ROW,IPNODE,IPTIME
          END DO
        END IF
        IF(IPK.LE.12) THEN
        DO IBLOCK=1,4
          DRESNORM=0.0D0
          DAXNORM=0.0D0
          DROWSUM=0.0D0
          DROWMAX=0.0D0
          IOUT_MAX_ROW=0
          DO I=1,NX*IAT
            IROW=I+(IBLOCK-1)*NX*IAT
            DAXI=STEP_SIGN*(RELAX_RETRY_U*CG_AP(IROW)+RELAX_RETRY_VW*CG_R(IROW)+RELAX_RETRY_P*CG_Z(IROW))
            DFSCALED=FFS(IROW)-CG_P(IROW)
            DRESNORM=DRESNORM+DFSCALED*DFSCALED
            DAXNORM=DAXNORM+DAXI*DAXI
            DROWSUM=DROWSUM+(DFSCALED-DAXI)*(DFSCALED-DAXI)
            IF(DABS(DFSCALED-DAXI).GT.DROWMAX) THEN
              DROWMAX=DABS(DFSCALED-DAXI)
              IOUT_MAX_ROW=IROW
            END IF
          END DO
          DRESNORM=DSQRT(DRESNORM)
          DAXNORM=DSQRT(DAXNORM)
          DROWSUM=DSQRT(DROWSUM)
          IF(DRESNORM.GT.0.0D0) THEN
            WRITE(*,*)'DIRECTIONAL_RAW_CHECK BLOCK/FD/JD/DIFF/REL/MAXROW/MAXDIFF=',IBLOCK,DRESNORM,DAXNORM,DROWSUM,DROWSUM/DRESNORM,IOUT_MAX_ROW,DROWMAX
          ELSE
            WRITE(*,*)'DIRECTIONAL_RAW_CHECK BLOCK/FD/JD/DIFF/REL/MAXROW/MAXDIFF=',IBLOCK,DRESNORM,DAXNORM,DROWSUM,0.0D0,IOUT_MAX_ROW,DROWMAX
          END IF
        END DO
        DO IBLOCK=1,4
          DRESNORM=0.0D0
          DAXNORM=0.0D0
          DROWSUM=0.0D0
          DROWMAX=0.0D0
          IOUT_MAX_ROW=0
          DO I=1,NX*IAT
            IROW=I+(IBLOCK-1)*NX*IAT
            DAXI=STEP_SIGN*(RELAX_RETRY_U*CG_AP(IROW)+RELAX_RETRY_VW*CG_R(IROW)+RELAX_RETRY_P*CG_Z(IROW))
            DFSCALED=FFS_PHYS(IROW)-FFS_PHYS_BASE(IROW)
            DRESNORM=DRESNORM+DFSCALED*DFSCALED
            DAXNORM=DAXNORM+DAXI*DAXI
            DROWSUM=DROWSUM+(DFSCALED-DAXI)*(DFSCALED-DAXI)
            IF(DABS(DFSCALED-DAXI).GT.DROWMAX) THEN
              DROWMAX=DABS(DFSCALED-DAXI)
              IOUT_MAX_ROW=IROW
            END IF
          END DO
          DRESNORM=DSQRT(DRESNORM)
          DAXNORM=DSQRT(DAXNORM)
          DROWSUM=DSQRT(DROWSUM)
          IF(DRESNORM.GT.0.0D0) THEN
            WRITE(*,*)'DIRECTIONAL_PHYS_CHECK BLOCK/FD/JD/DIFF/REL/MAXROW/MAXDIFF=',IBLOCK,DRESNORM,DAXNORM,DROWSUM,DROWSUM/DRESNORM,IOUT_MAX_ROW,DROWMAX
          ELSE
            WRITE(*,*)'DIRECTIONAL_PHYS_CHECK BLOCK/FD/JD/DIFF/REL/MAXROW/MAXDIFF=',IBLOCK,DRESNORM,DAXNORM,DROWSUM,0.0D0,IOUT_MAX_ROW,DROWMAX
          END IF
        END DO
        DO IBLOCK=1,4
          DRESNORM=0.0D0
          DAXNORM=0.0D0
          DROWSUM=0.0D0
          DROWMAX=0.0D0
          IOUT_MAX_ROW=0
          DO I=1,NX*IAT
            IROW=I+(IBLOCK-1)*NX*IAT
            IF(ROW_SCALE(IROW).GT.0.0D0) THEN
              DAXI=STEP_SIGN*(RELAX_RETRY_U*CG_AP(IROW)+RELAX_RETRY_VW*CG_R(IROW)+RELAX_RETRY_P*CG_Z(IROW))*ROW_SCALE(IROW)
              DFSCALED=(FFS(IROW)-CG_P(IROW))*ROW_SCALE(IROW)
            ELSE
              DAXI=0.0D0
              DFSCALED=0.0D0
            END IF
            DRESNORM=DRESNORM+DFSCALED*DFSCALED
            DAXNORM=DAXNORM+DAXI*DAXI
            DROWSUM=DROWSUM+(DFSCALED-DAXI)*(DFSCALED-DAXI)
            IF(DABS(DFSCALED-DAXI).GT.DROWMAX) THEN
              DROWMAX=DABS(DFSCALED-DAXI)
              IOUT_MAX_ROW=IROW
            END IF
          END DO
          DRESNORM=DSQRT(DRESNORM)
          DAXNORM=DSQRT(DAXNORM)
          DROWSUM=DSQRT(DROWSUM)
          IF(DRESNORM.GT.0.0D0) THEN
            WRITE(*,*)'DIRECTIONAL_SCALED_CHECK BLOCK/FD/JD/DIFF/REL/MAXROW/MAXDIFF=',IBLOCK,DRESNORM,DAXNORM,DROWSUM,DROWSUM/DRESNORM,IOUT_MAX_ROW,DROWMAX
          ELSE
            WRITE(*,*)'DIRECTIONAL_SCALED_CHECK BLOCK/FD/JD/DIFF/REL/MAXROW/MAXDIFF=',IBLOCK,DRESNORM,DAXNORM,DROWSUM,0.0D0,IOUT_MAX_ROW,DROWMAX
          END IF
        END DO
        END IF
        IF(INEWTON_PARDISO.NE.0) THEN
          IF(WW4.GE.WW4_PREV) THEN
            IF(0.NE.0.AND.ILINESEARCH_TRIES.GE.8.AND.ILINESEARCH_EXTRA_MODE.EQ.0) THEN
              DO I=1,NX
              DO ITT=1,IAT
                U(ITT,I)=U_SAVE(ITT,I)
                V(ITT,I)=V_SAVE(ITT,I)
                W(ITT,I)=W_SAVE(ITT,I)
                P(ITT,I)=P_SAVE(ITT,I)
              END DO
              END DO
              DO I=1,N
                DDUVWP(I)=DDUVWP_SD(I)
              END DO
              RELAX_RETRY_U=1.0D0
              RELAX_RETRY_VW=1.0D0
              RELAX_RETRY_P=1.0D0
              ILINESEARCH_TRIES=0
              ILINESEARCH_EXTRA_MODE=1
              WRITE(*,*)'NEWTON_PARDISO SWITCH_TO_STEEPEST_DESCENT WW4_PREV/WW4=',WW4_PREV,WW4
              DO IK=1,4
              DO I=1,NX
              DO ITT=1,IAT
                ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
                IF(IS_CONSTRAINED(ISP).NE.0) THEN
                  IF(ISP.LE.NX*IAT) THEN
                    U(ITT,I)=U(ITT,I)+RELAX_RETRY_U*DDUVWP(ISP)
                  ELSE IF(ISP.LE.2*NX*IAT) THEN
                    V(ITT,I)=V(ITT,I)+RELAX_RETRY_VW*DDUVWP(ISP)
                  ELSE IF(ISP.LE.3*NX*IAT) THEN
                    W(ITT,I)=W(ITT,I)+RELAX_RETRY_VW*DDUVWP(ISP)
                  ELSE
                    P(ITT,I)=P(ITT,I)+RELAX_RETRY_P*DDUVWP(ISP)
                  END IF
                ELSE IF(ISP.LE.NX*IAT) THEN
                  U(ITT,I)=U(ITT,I)+RELAX_RETRY_U*DDUVWP(ISP)
                ELSE IF(ISP.LE.2*NX*IAT) THEN
                  V(ITT,I)=V(ITT,I)+RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE IF(ISP.LE.3*NX*IAT) THEN
                  W(ITT,I)=W(ITT,I)+RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE
                  P(ITT,I)=P(ITT,I)+RELAX_RETRY_P*DDUVWP(ISP)
                END IF
              END DO
              END DO
              END DO
              IPENDING_STEP=1
              GOTO 30
            END IF
            IF(ILINESEARCH_TRIES.LT.24.AND.DMAX1(RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P).GT.1.0D-10) THEN
              DO I=1,NX
              DO ITT=1,IAT
                U(ITT,I)=U_SAVE(ITT,I)
                V(ITT,I)=V_SAVE(ITT,I)
                W(ITT,I)=W_SAVE(ITT,I)
                P(ITT,I)=P_SAVE(ITT,I)
              END DO
              END DO
              RELAX_RETRY_U=0.5D0*RELAX_RETRY_U
              RELAX_RETRY_VW=0.5D0*RELAX_RETRY_VW
              RELAX_RETRY_P=0.5D0*RELAX_RETRY_P
              ILINESEARCH_TRIES=ILINESEARCH_TRIES+1
              WRITE(*,*)'NEWTON_PARDISO BACKTRACK WW4_PREV/WW4/TRY/RELAX_U/VW/P=',WW4_PREV,WW4,ILINESEARCH_TRIES,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
              DO IK=1,4
              DO I=1,NX
              DO ITT=1,IAT
                ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
                IF(IS_CONSTRAINED(ISP).NE.0) THEN
                  IF(ISP.LE.NX*IAT) THEN
                    U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                  ELSE IF(ISP.LE.2*NX*IAT) THEN
                    V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                  ELSE IF(ISP.LE.3*NX*IAT) THEN
                    W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                  ELSE
                    P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                  END IF
                ELSE IF(ISP.LE.NX*IAT) THEN
                  U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                ELSE IF(ISP.LE.2*NX*IAT) THEN
                  V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE IF(ISP.LE.3*NX*IAT) THEN
                  W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE
                  P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                END IF
              END DO
              END DO
              END DO
              IPENDING_STEP=1
              GOTO 30
            END IF
            WRITE(*,*)'NEWTON_PARDISO LINESEARCH REJECT LIMIT WW4_PREV/WW4=',WW4_PREV,WW4
            STOP
          END IF
          DLAMBDA_RATIO=0.0D0
          DLM_LAMBDA_CURRENT=0.0D0
          WRITE(*,*)'NEWTON_PARDISO LINESEARCH ACCEPT WW4_PREV/WW4/TRY/RELAX_U/VW/P=',WW4_PREV,WW4,ILINESEARCH_TRIES,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
          GOTO 7788
        END IF
        DWW4_DROP=WW4_PREV-WW4
        DWW4_REL_DROP=DWW4_DROP/DMAX1(1.0D-30,DABS(WW4_PREV))
        IF(Re.GE.1.0D5.AND.IMARQUARDT_CG.NE.0.AND.ILINESEARCH_TRIES.EQ.0.AND.ILINESEARCH_EXTRA_MODE.EQ.0) THEN
          IF(ILAMBDA_SCAN_ACTIVE.NE.0) THEN
            IF(WW4.LT.WW4_LAMBDA_BEST) THEN
              WW4_LAMBDA_BEST=WW4
              DLM_LAMBDA_SCAN_BEST=DLM_LAMBDA_USED
              DO I=1,NX
              DO ITT=1,IAT
                U_LSCAN_BEST(ITT,I)=U(ITT,I)
                V_LSCAN_BEST(ITT,I)=V(ITT,I)
                W_LSCAN_BEST(ITT,I)=W(ITT,I)
                P_LSCAN_BEST(ITT,I)=P(ITT,I)
              END DO
              END DO
              DO I=1,N
                DDUVWP_LSCAN_BEST(I)=DDUVWP(I)
              END DO
            END IF
            WRITE(*,*)'LAMBDA_SCAN CANDIDATE STAGE/LAMBDA/WW4/BEST=',ILAMBDA_SCAN_STAGE,DLM_LAMBDA_USED,WW4,WW4_LAMBDA_BEST
            IF(ILAMBDA_SCAN_STAGE.LT.ILAMBDA_SCAN_COUNT) THEN
              ILAMBDA_SCAN_STAGE=ILAMBDA_SCAN_STAGE+1
              DO I=1,NX
              DO ITT=1,IAT
                U(ITT,I)=U_SAVE(ITT,I)
                V(ITT,I)=V_SAVE(ITT,I)
                W(ITT,I)=W_SAVE(ITT,I)
                P(ITT,I)=P_SAVE(ITT,I)
              END DO
              END DO
              IF(ILAMBDA_SCAN_STAGE.EQ.2) THEN
                DLM_LAMBDA_CURRENT=DLM_LAMBDA_SCAN_CAND2
              ELSE
                DLM_LAMBDA_CURRENT=DLM_LAMBDA_SCAN_CAND3
              END IF
              IPENDING_STEP=0
              ILINESEARCH_TRIES=0
              ILINESEARCH_ACCEPTED=0
              ILINESEARCH_EXTRA_MODE=0
              WRITE(*,*)'LAMBDA_SCAN NEXT STAGE/LAMBDA=',ILAMBDA_SCAN_STAGE,DLM_LAMBDA_CURRENT
              GOTO 30
            ELSE
              DO I=1,NX
              DO ITT=1,IAT
                U(ITT,I)=U_LSCAN_BEST(ITT,I)
                V(ITT,I)=V_LSCAN_BEST(ITT,I)
                W(ITT,I)=W_LSCAN_BEST(ITT,I)
                P(ITT,I)=P_LSCAN_BEST(ITT,I)
              END DO
              END DO
              DO I=1,N
                DDUVWP(I)=DDUVWP_LSCAN_BEST(I)
              END DO
              WW4=WW4_LAMBDA_BEST
              DLM_LAMBDA_CURRENT=DLM_LAMBDA_SCAN_BEST
              DLM_LAMBDA_USED=DLM_LAMBDA_SCAN_BEST
              ILAMBDA_SCAN_ACTIVE=0
              ILAMBDA_SCAN_STAGE=0
              DWW4_DROP=WW4_PREV-WW4
              DWW4_REL_DROP=DWW4_DROP/DMAX1(1.0D-30,DABS(WW4_PREV))
              WRITE(*,*)'LAMBDA_SCAN ACCEPT BEST LAMBDA/WW4/REL_DROP=',DLM_LAMBDA_CURRENT,WW4,DWW4_REL_DROP
            END IF
          ELSE IF(DWW4_DROP.GT.0.0D0.AND.DWW4_REL_DROP.LT.3.0D-3) THEN
            ILAMBDA_SCAN_ACTIVE=1
            ILAMBDA_SCAN_STAGE=1
            ILAMBDA_SCAN_COUNT=3
            WW4_LAMBDA_BEST=WW4
            DLM_LAMBDA_SCAN_CAND1=DMAX1(DLM_LAMBDA_MIN,DLM_LAMBDA_USED)
            DLM_LAMBDA_SCAN_CAND2=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,10.0D0*DLM_LAMBDA_SCAN_CAND1))
            DLM_LAMBDA_SCAN_CAND3=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,0.1D0*DLM_LAMBDA_SCAN_CAND1))
            DLM_LAMBDA_SCAN_BEST=DLM_LAMBDA_SCAN_CAND1
            DO I=1,NX
            DO ITT=1,IAT
              U_LSCAN_BEST(ITT,I)=U(ITT,I)
              V_LSCAN_BEST(ITT,I)=V(ITT,I)
              W_LSCAN_BEST(ITT,I)=W(ITT,I)
              P_LSCAN_BEST(ITT,I)=P(ITT,I)
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            DO I=1,N
              DDUVWP_LSCAN_BEST(I)=DDUVWP(I)
            END DO
            DLM_LAMBDA_CURRENT=DLM_LAMBDA_SCAN_CAND2
            ILAMBDA_SCAN_STAGE=2
            IPENDING_STEP=0
            ILINESEARCH_TRIES=0
            ILINESEARCH_ACCEPTED=0
            ILINESEARCH_EXTRA_MODE=0
            WRITE(*,*)'LAMBDA_SCAN START REL_DROP/CAND1/CAND2/CAND3=',DWW4_REL_DROP,DLM_LAMBDA_SCAN_CAND1,DLM_LAMBDA_SCAN_CAND2,DLM_LAMBDA_SCAN_CAND3
            WRITE(*,*)'LAMBDA_SCAN NEXT STAGE/LAMBDA=',2,DLM_LAMBDA_CURRENT
            GOTO 30
          END IF
        END IF
        IF(DWW4_DROP.GT.0.0D0) THEN
          IGN_DECREASE_COUNT=IGN_DECREASE_COUNT+1
        ELSE
          IGN_DECREASE_COUNT=0
          IF(IGN_ACTIVE.NE.0) THEN
            DGN_LAMBDA_SWITCH=DMAX1(1.0D-12,DGN_LAMBDA_SWITCH/10.0D0)
            WRITE(*,*)'GAUSS_NEWTON_TRIAL FAILED LOWER_SWITCH=',DGN_LAMBDA_SWITCH
          END IF
        END IF
        DLM_LAMBDA_OLD=DLM_LAMBDA_CURRENT
        ISALS_EVAL=0
!       Lambda is controlled by the measured nonlinear residual drop only:
!       DWW4_DROP = WW4_PREV - WW4.  The SALS ratio below is kept as a
!       diagnostic of the linearized model, but it does not choose lambda.
!       Backtracking and block trials must not change lambda.
        IF(ILINESEARCH_TRIES.EQ.0.AND.ILINESEARCH_EXTRA_MODE.EQ.0) THEN
        ISALS_EVAL=1
        DS_OLD=0.0D0
        DS_NEW=0.0D0
        DS_PRED=0.0D0
        DPHYS_OLD=0.0D0
        DPHYS_NEW=0.0D0
        DPHYS_PRED=0.0D0
        DLM_FDOTJD=0.0D0
        DLM_JD_NORM2=0.0D0
        DLM_GRAD_DROP=0.0D0
        DLM_DAMP_QUAD=0.0D0
        DLM_TRACE=0.0D0
        DLM_TRACE_INV=0.0D0
        DO I=1,N
!         DAXI gives R*J*d for diagnostics.  The SALS ratio itself uses
!         d^T*(J^T*R^T*R*f)-lambda*d^T*D*d below.
          DAXI=STEP_SIGN*(CG_AP(I)+CG_R(I)+CG_Z(I))
          DZ_SALS=0.0D0
          IF(COL_SCALE(I).GT.0.0D0) DZ_SALS=STEP_SIGN*DDUVWP(I)/COL_SCALE(I)
          DAXI_SCALED=DAXI*ROW_SCALE_BASE(I)
          DFSCALED=FFS_PHYS_BASE(I)*ROW_SCALE_BASE(I)
          DS_OLD=DS_OLD+DFSCALED*DFSCALED
          DS_NEW=DS_NEW+(FFS_PHYS(I)*ROW_SCALE_BASE(I))*(FFS_PHYS(I)*ROW_SCALE_BASE(I))
          DS_PRED=DS_PRED+(DFSCALED+DAXI_SCALED)*(DFSCALED+DAXI_SCALED)
          DLM_TRACE=DLM_TRACE+SALS_DIAG(I)
          IF(SALS_DIAG(I).GT.1.0D-300) DLM_TRACE_INV=DLM_TRACE_INV+1.0D0/SALS_DIAG(I)
          DPHYS_JD=DAXI
          DPHYS_OLD=DPHYS_OLD+FFS_PHYS_BASE(I)*FFS_PHYS_BASE(I)
          DPHYS_NEW=DPHYS_NEW+FFS_PHYS(I)*FFS_PHYS(I)
          DPHYS_PRED=DPHYS_PRED+(FFS_PHYS_BASE(I)+DPHYS_JD)*(FFS_PHYS_BASE(I)+DPHYS_JD)
          DLM_FDOTJD=DLM_FDOTJD+DFSCALED*DAXI_SCALED
          DLM_JD_NORM2=DLM_JD_NORM2+DAXI_SCALED*DAXI_SCALED
          DLM_GRAD_DROP=DLM_GRAD_DROP+DZ_SALS*SALS_RHS(I)
          DLM_DAMP_QUAD=DLM_DAMP_QUAD+SALS_DIAG(I)*DZ_SALS*DZ_SALS
        END DO
        DLM_PRED_DROP_NORM=DPHYS_OLD-DPHYS_PRED
        DLM_PRED_DROP=DLM_GRAD_DROP-DLM_LAMBDA_OLD*DLM_DAMP_QUAD
        DLM_ACT_DROP=DS_OLD-DS_NEW
        IF(DABS(DLM_PRED_DROP).GT.1.0D-300) THEN
          DLAMBDA_RATIO=DLM_ACT_DROP/DLM_PRED_DROP
        ELSE
          DLAMBDA_RATIO=-50.0D0
          WRITE(*,*)'MARQUARDT RATIO INVALID PRED_DROP_ZERO ACT_DROP=',DLM_PRED_DROP,DLM_ACT_DROP
        END IF
        DTRUST_TARGET=1.0D-4
        DTRUST_LOW=1.0D-5
        DTRUST_HIGH=1.0D-3
        DLAMBDA_TRUST_FACTOR=1.0D0
        ILAMBDA_REASON=12
        DLM_LAMBDA_C=1.0D0/DMAX1(DLM_TRACE_INV,1.0D-300)
        IF(DWW4_DROP.GT.0.0D0) THEN
          IF(Re.GE.1.0D5.AND.DWW4_REL_DROP.LT.1.0D-4) THEN
            DLAMBDA_TRUST_FACTOR=1.0D0
            ILAMBDA_REASON=34
          ELSE
            DLAMBDA_TRUST_FACTOR=1.0D0/1.2D0
            DLM_LAMBDA_CURRENT=DLM_LAMBDA_CURRENT*DLAMBDA_TRUST_FACTOR
            ILAMBDA_REASON=31
          END IF
        ELSE
          DLAMBDA_TRUST_FACTOR=2.0D0
          IF(DLM_LAMBDA_CURRENT.LE.0.0D0) DLM_LAMBDA_CURRENT=DMAX1(DLM_LAMBDA_C,1.0D0)
          DLM_LAMBDA_CURRENT=DLM_LAMBDA_CURRENT*DLAMBDA_TRUST_FACTOR
          ILAMBDA_REASON=32
          DRELAX_START_SCALE=1.0D0
        END IF
        DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,DLM_LAMBDA_CURRENT))
        ELSE
          DLAMBDA_TRUST_FACTOR=1.0D0
          ILAMBDA_REASON=0
          WRITE(*,*)'MARQUARDT SALS SKIP LINESEARCH TRY/BLOCKMODE=',ILINESEARCH_TRIES,ILINESEARCH_EXTRA_MODE
        END IF
        IF(.FALSE..AND.DWW4_REL_DROP.LT.1.0D-6.AND.DMAX1(RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P).LT.1.0D-5.AND.ISTAG_ACCEPT_STEP.NE.IPK) THEN
          ISTAG_LAMBDA_COUNT=ISTAG_LAMBDA_COUNT+1
          ISTAG_ACCEPT_STEP=IPK
          IF(ISTAG_LAMBDA_COUNT.LE.1) THEN
            DLAMBDA_FACTOR=4.0D0
          ELSE IF(ISTAG_LAMBDA_COUNT.EQ.2) THEN
            DLAMBDA_FACTOR=10.0D0
          ELSE
            DLAMBDA_FACTOR=100.0D0
          END IF
          DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,DLAMBDA_FACTOR*DMAX1(DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT)))
          DRELAX_START_SCALE=1.0D0
          ILAMBDA_REASON=25
          ISTAG_LAMBDA_HOLD=3
          WRITE(*,*)'MARQUARDT STAGNATION_ACCEPT INCREASE_LAMBDA COUNT/FACTOR/REL_DROP/RELAX=',ISTAG_LAMBDA_COUNT,DLAMBDA_FACTOR,DWW4_REL_DROP,DMAX1(RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P)
        END IF
        IF(ISALS_EVAL.NE.0) THEN
          WRITE(*,*)'MARQUARDT SALS THEORY ACT/PRED/R/S_OLD/S_NEW/S_PRED=',DLM_ACT_DROP,DLM_PRED_DROP,DLAMBDA_RATIO,DS_OLD,DS_NEW,DS_PRED
          WRITE(*,*)'MARQUARDT SALS THEORY PHYS ACT/PRED/R/S_OLD/S_NEW/S_PRED=',DLM_ACT_DROP,DLM_PRED_DROP,DLAMBDA_RATIO,DPHYS_OLD,DPHYS_NEW,DPHYS_PRED
          WRITE(*,*)'MARQUARDT SALS PRED_FORMULA JTRF-LAMBDA_DTD/JTRF/DTD=',DLM_PRED_DROP,DLM_GRAD_DROP,DLM_DAMP_QUAD
          WRITE(*,*)'MARQUARDT LINEAR_MODEL_DIAG -2FTJD-JD2/NORMDIFF/FTJD/JD2=',-2.0D0*DLM_FDOTJD-DLM_JD_NORM2,DLM_PRED_DROP_NORM,DLM_FDOTJD,DLM_JD_NORM2
          WRITE(*,*)'MARQUARDT TRACE DIAG/INV_APPROX=',DLM_TRACE,DLM_TRACE_INV
          WRITE(*,*)'MARQUARDT LAMBDA_C_INV_TRACE_APPROX=',1.0D0/DMAX1(DLM_TRACE_INV,1.0D-300)
          WRITE(*,*)'MARQUARDT SALS_RATIO_DIAGNOSTIC/MEASURED_FACTOR=',DLAMBDA_RATIO,DLAMBDA_TRUST_FACTOR
          WRITE(*,*)'MARQUARDT MEASURED_LAMBDA_CONTROL DROP/REL_DROP=',DWW4_DROP,DWW4_REL_DROP
          WRITE(*,*)'GAUSS_NEWTON_TRIAL READINESS COUNT/LAMBDA/SWITCH=',IGN_DECREASE_COUNT,DLM_LAMBDA_CURRENT,DGN_LAMBDA_SWITCH
          WRITE(*,*)'MARQUARDT LAMBDA UPDATE OLD/NEW/REASON/REL_DROP=',DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT,ILAMBDA_REASON,DWW4_REL_DROP
        ELSE
          WRITE(*,*)'MARQUARDT LAMBDA KEEP DURING LINESEARCH OLD/NEW=',DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT
        END IF
        DAXI=DMAX1(RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P,1.0D-30)
        WRITE(*,*)'MARQUARDT_ACTUAL_DPHI STEP_SIGN/RELAX/DPHI_PER_REL/OLD/NEW=',STEP_SIGN,DAXI,(WW4*WW4-WW4_PREV*WW4_PREV)/DAXI,WW4_PREV,WW4
        IFAST_RESCUE=0
        IF(.FALSE..AND.ISALS_EVAL.NE.0.AND.DLAMBDA_RATIO.LT.0.0D0) IFAST_RESCUE=1
        IF(DWW4_DROP.GT.0.0D0.AND.DWW4_REL_DROP.LE.1.0D-6) IFAST_RESCUE=1
        IBLOCK_TRY_LIMIT=24
        IF(IFAST_RESCUE.NE.0) IBLOCK_TRY_LIMIT=4
        WRITE(*,*)'MARQUARDT EARLY_BLOCK_CANDIDATE FLAG/LIMIT/R/REL_RISE=',IFAST_RESCUE,IBLOCK_TRY_LIMIT,DLAMBDA_RATIO,(WW4-WW4_PREV)/DMAX1(1.0D-30,DABS(WW4_PREV))
        IF(WW4.LT.WW4_LINE_BEST) THEN
          WW4_LINE_BEST=WW4
          RELAX_BEST_U=RELAX_RETRY_U
          RELAX_BEST_VW=RELAX_RETRY_VW
          RELAX_BEST_P=RELAX_RETRY_P
          DBEST_STEP_SIGN=STEP_SIGN
          IBEST_BACKTRACKS=ILINESEARCH_TRIES
          WRITE(*,*)'MARQUARDT CANDIDATE BEST WW4_PREV/BEST/SIGN/TRY/RELAX_U/VW/P=',WW4_PREV,WW4_LINE_BEST,DBEST_STEP_SIGN,IBEST_BACKTRACKS,RELAX_BEST_U,RELAX_BEST_VW,RELAX_BEST_P
        END IF
        IF(.FALSE..AND.IPK.GE.30.AND.DWW4_REL_DROP.GT.0.0D0.AND.DWW4_REL_DROP.LT.3.0D-3.AND.ISTAG_SD_TRIED.EQ.0) THEN
          DO I=1,NX
          DO ITT=1,IAT
            U(ITT,I)=U_SAVE(ITT,I)
            V(ITT,I)=V_SAVE(ITT,I)
            W(ITT,I)=W_SAVE(ITT,I)
            P(ITT,I)=P_SAVE(ITT,I)
          END DO
          END DO
          DO I=1,N
            DDUVWP(I)=DDUVWP_SD(I)
          END DO
          DU_SCALE=1.0D-12
          DV_SCALE=1.0D-12
          DW_SCALE=1.0D-12
          DP_SCALE=1.0D-12
          DO I=1,NX
          DO ITT=1,IAT
            DU_SCALE=DMAX1(DU_SCALE,DABS(U_SAVE(ITT,I)))
            DV_SCALE=DMAX1(DV_SCALE,DABS(V_SAVE(ITT,I)))
            DW_SCALE=DMAX1(DW_SCALE,DABS(W_SAVE(ITT,I)))
            DP_SCALE=DMAX1(DP_SCALE,DABS(P_SAVE(ITT,I)))
          END DO
          END DO
          DSD_MAX_REL=0.0D0
          DO I=1,N
            IF(I.LE.NX*IAT) THEN
              DSD_MAX_REL=DMAX1(DSD_MAX_REL,DABS(DDUVWP(I))/DU_SCALE)
            ELSE IF(I.LE.2*NX*IAT) THEN
              DSD_MAX_REL=DMAX1(DSD_MAX_REL,DABS(DDUVWP(I))/DV_SCALE)
            ELSE IF(I.LE.3*NX*IAT) THEN
              DSD_MAX_REL=DMAX1(DSD_MAX_REL,DABS(DDUVWP(I))/DW_SCALE)
            ELSE
              DSD_MAX_REL=DMAX1(DSD_MAX_REL,DABS(DDUVWP(I))/DP_SCALE)
            END IF
          END DO
          DSD_SCALE=DMAX1(1.0D0,DSD_MAX_REL/1.0D-2)
          DO I=1,N
            DDUVWP(I)=DDUVWP(I)/DSD_SCALE
          END DO
          IF(DDIR_STEEPEST.LT.0.0D0) THEN
            STEP_SIGN=1.0D0
          ELSE
            STEP_SIGN=-1.0D0
          END IF
          ISTAG_SD_TRIED=1
          ILINESEARCH_TRIES=0
          ILINESEARCH_EXTRA_MODE=11
          IOPPOSITE_TRIED=1
          RELAX_RETRY_U=1.0D0
          RELAX_RETRY_VW=1.0D0
          RELAX_RETRY_P=1.0D0
          WRITE(*,*)'MARQUARDT STAGNATION SMALL_DROP TRY SCALED_STEEPEST SIGN/MAXREL/SCALE/TARGET_REL/REL_DROP=',STEP_SIGN,DSD_MAX_REL,DSD_SCALE,1.0D-2,DWW4_REL_DROP
          DO IK=1,4
          DO I=1,NX
          DO ITT=1,IAT
            ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
            IF(ISP.LE.NX*IAT) THEN
              U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
            ELSE IF(ISP.LE.2*NX*IAT) THEN
              V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
            ELSE IF(ISP.LE.3*NX*IAT) THEN
              W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
            ELSE
              P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
            END IF
          END DO
          END DO
          END DO
          IPENDING_STEP=1
          GOTO 30
        END IF
        IF(WW4.GE.WW4_PREV) THEN
          IF(ILAMBDA_REASON.EQ.25) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            IPENDING_STEP=0
            ILINESEARCH_TRIES=0
            IOPPOSITE_TRIED=0
            RELAX_RETRY_U=RELAX_BASE_U
            RELAX_RETRY_VW=RELAX_BASE_VW
            RELAX_RETRY_P=RELAX_BASE_P
            WRITE(*,*)'MARQUARDT STAGNATION RECOMPUTE_DIRECTION LAMBDA=',DLM_LAMBDA_CURRENT
            GOTO 30
          END IF
          IF(.FALSE..AND.DLAMBDA_RATIO.LT.0.0D0.AND.ILINESEARCH_TRIES.GE.2.AND.ISTAG_SD_TRIED.EQ.0) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            DO I=1,N
              DDUVWP(I)=DDUVWP_SD(I)
            END DO
            DU_SCALE=1.0D-12
            DV_SCALE=1.0D-12
            DW_SCALE=1.0D-12
            DP_SCALE=1.0D-12
            DO I=1,NX
            DO ITT=1,IAT
              DU_SCALE=DMAX1(DU_SCALE,DABS(U_SAVE(ITT,I)))
              DV_SCALE=DMAX1(DV_SCALE,DABS(V_SAVE(ITT,I)))
              DW_SCALE=DMAX1(DW_SCALE,DABS(W_SAVE(ITT,I)))
              DP_SCALE=DMAX1(DP_SCALE,DABS(P_SAVE(ITT,I)))
            END DO
            END DO
            DSD_MAX_REL=0.0D0
            DO I=1,N
              IF(I.LE.NX*IAT) THEN
                DSD_MAX_REL=DMAX1(DSD_MAX_REL,DABS(DDUVWP(I))/DU_SCALE)
              ELSE IF(I.LE.2*NX*IAT) THEN
                DSD_MAX_REL=DMAX1(DSD_MAX_REL,DABS(DDUVWP(I))/DV_SCALE)
              ELSE IF(I.LE.3*NX*IAT) THEN
                DSD_MAX_REL=DMAX1(DSD_MAX_REL,DABS(DDUVWP(I))/DW_SCALE)
              ELSE
                DSD_MAX_REL=DMAX1(DSD_MAX_REL,DABS(DDUVWP(I))/DP_SCALE)
              END IF
            END DO
            DSD_SCALE=DMAX1(1.0D0,DSD_MAX_REL/1.0D-2)
            DO I=1,N
              DDUVWP(I)=DDUVWP(I)/DSD_SCALE
            END DO
            STEP_SIGN=DSD_PHYS_SIGN
            ISTAG_SD_TRIED=1
            ILINESEARCH_TRIES=0
            ILINESEARCH_EXTRA_MODE=11
            IOPPOSITE_TRIED=1
            RELAX_RETRY_U=1.0D0
            RELAX_RETRY_VW=1.0D0
            RELAX_RETRY_P=1.0D0
            WRITE(*,*)'MARQUARDT STAGNATION TRY SCALED_STEEPEST SIGN/MAXREL/SCALE/TARGET_REL=',STEP_SIGN,DSD_MAX_REL,DSD_SCALE,1.0D-2
            DO IK=1,4
            DO I=1,NX
            DO ITT=1,IAT
              ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
              IF(ISP.LE.NX*IAT) THEN
                U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
              ELSE IF(ISP.LE.2*NX*IAT) THEN
                V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE IF(ISP.LE.3*NX*IAT) THEN
                W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE
                P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
              END IF
            END DO
            END DO
            END DO
            IPENDING_STEP=1
            GOTO 30
          END IF
          IF(.FALSE..AND.IOPPOSITE_TRIED.NE.0.AND.ILINESEARCH_EXTRA_MODE.GE.10.AND.ILINESEARCH_TRIES.GE.IBLOCK_TRY_LIMIT.AND.DLAMBDA_RATIO.LT.0.0D0.AND.DLM_LAMBDA_OLD.GT.1.01D0*DLM_LAMBDA_MIN) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,10.0D0*DLM_LAMBDA_OLD))
            IPENDING_STEP=0
            ILINESEARCH_TRIES=0
            IOPPOSITE_TRIED=0
            RELAX_RETRY_U=RELAX_BASE_U
            RELAX_RETRY_VW=RELAX_BASE_VW
            RELAX_RETRY_P=RELAX_BASE_P
            WRITE(*,*)'MARQUARDT NEGATIVE_R INCREASE_LAMBDA OLD/NEW/R/WW4_PREV/WW4=',DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT,DLAMBDA_RATIO,WW4_PREV,WW4
            GOTO 30
          END IF
          IF(.FALSE..AND.IOPPOSITE_TRIED.NE.0.AND.ILINESEARCH_EXTRA_MODE.GE.10.AND.ILINESEARCH_TRIES.GE.IBLOCK_TRY_LIMIT.AND.ERR1_REL_MAX.LT.1.0D-8.AND.WW4.GT.1.0D-8.AND.DLM_LAMBDA_OLD.GT.1.01D0*DLM_LAMBDA_MIN) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,10.0D0*DLM_LAMBDA_OLD))
            IPENDING_STEP=0
            ILINESEARCH_TRIES=0
            IOPPOSITE_TRIED=0
            RELAX_RETRY_U=RELAX_BASE_U
            RELAX_RETRY_VW=RELAX_BASE_VW
            RELAX_RETRY_P=RELAX_BASE_P
            WRITE(*,*)'MARQUARDT SMALL_UPDATE INCREASE_LAMBDA OLD/NEW/ERR1/WW4_PREV/WW4=',DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT,ERR1_REL_MAX,WW4_PREV,WW4
            GOTO 30
          END IF
          IF(.FALSE..AND.DLAMBDA_RATIO.LT.0.0D0.AND.ILINESEARCH_TRIES.GE.2.AND.IOPPOSITE_TRIED.NE.0.AND.ISTAG_SD_TRIED.EQ.0.AND.DLM_LAMBDA_CURRENT.GT.DLM_LAMBDA_OLD) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            IPENDING_STEP=0
            ILINESEARCH_TRIES=0
            IOPPOSITE_TRIED=0
            RELAX_RETRY_U=RELAX_BASE_U
            RELAX_RETRY_VW=RELAX_BASE_VW
            RELAX_RETRY_P=RELAX_BASE_P
            WRITE(*,*)'MARQUARDT NEGATIVE_R RECOMPUTE_DIRECTION LAMBDA_OLD/NEW/R=',DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT,DLAMBDA_RATIO
            GOTO 30
          END IF
          IF(ILINESEARCH_TRIES.LT.IBLOCK_TRY_LIMIT.AND.DMAX1(RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P).GT.1.0D-10) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            RELAX_RETRY_U=0.5D0*RELAX_RETRY_U
            RELAX_RETRY_VW=0.5D0*RELAX_RETRY_VW
            RELAX_RETRY_P=0.5D0*RELAX_RETRY_P
            ILINESEARCH_TRIES=ILINESEARCH_TRIES+1
            WRITE(*,*)'MARQUARDT BACKTRACK TRY WW4_PREV/WW4/TRY/RELAX_U/VW/P=',WW4_PREV,WW4,ILINESEARCH_TRIES,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
            DO IK=1,4
            DO I=1,NX
            DO ITT=1,IAT
              ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
              IF(IS_CONSTRAINED(ISP).NE.0) THEN
                IF(ISP.LE.NX*IAT) THEN
                  U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                ELSE IF(ISP.LE.2*NX*IAT) THEN
                  V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE IF(ISP.LE.3*NX*IAT) THEN
                  W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE
                  P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                END IF
              ELSE IF(ISP.LE.NX*IAT) THEN
                U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
              ELSE IF(ISP.LE.2*NX*IAT) THEN
                V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE IF(ISP.LE.3*NX*IAT) THEN
                W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE
                P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
              END IF
            END DO
            END DO
            END DO
            IPENDING_STEP=1
            GOTO 30
          END IF
          IF(.FALSE..AND.ILINESEARCH_TRIES.GE.IBLOCK_TRY_LIMIT.AND.ILINESEARCH_EXTRA_MODE.LT.10) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            ILINESEARCH_EXTRA_MODE=ILINESEARCH_EXTRA_MODE+1
            ILINESEARCH_TRIES=0
            RELAX_RETRY_U=RELAX_BASE_U
            RELAX_RETRY_VW=RELAX_BASE_VW
            RELAX_RETRY_P=RELAX_BASE_P
            IF(ILINESEARCH_EXTRA_MODE.EQ.1) THEN
              RELAX_RETRY_U=0.0D0
              RELAX_RETRY_VW=0.0D0
              RELAX_RETRY_P=RELAX_BASE_P
            ELSE IF(ILINESEARCH_EXTRA_MODE.EQ.2) THEN
              RELAX_RETRY_U=0.0D0
              RELAX_RETRY_VW=0.0D0
              RELAX_RETRY_P=0.3D0*RELAX_BASE_P
            ELSE IF(ILINESEARCH_EXTRA_MODE.EQ.3) THEN
              RELAX_RETRY_P=0.3D0*RELAX_BASE_P
            ELSE IF(ILINESEARCH_EXTRA_MODE.EQ.4) THEN
              RELAX_RETRY_P=0.1D0*RELAX_BASE_P
            ELSE IF(ILINESEARCH_EXTRA_MODE.EQ.5) THEN
              RELAX_RETRY_U=0.3D0*RELAX_BASE_U
              RELAX_RETRY_VW=0.3D0*RELAX_BASE_VW
              RELAX_RETRY_P=0.0D0
            ELSE IF(ILINESEARCH_EXTRA_MODE.EQ.6) THEN
              RELAX_RETRY_U=0.0D0
              RELAX_RETRY_VW=0.0D0
              RELAX_RETRY_P=1.0D-2*RELAX_BASE_P
            ELSE IF(ILINESEARCH_EXTRA_MODE.EQ.7) THEN
              RELAX_RETRY_U=0.0D0
              RELAX_RETRY_VW=0.0D0
              RELAX_RETRY_P=1.0D-3*RELAX_BASE_P
            ELSE IF(ILINESEARCH_EXTRA_MODE.EQ.8) THEN
              RELAX_RETRY_U=1.0D-2*RELAX_BASE_U
              RELAX_RETRY_VW=1.0D-2*RELAX_BASE_VW
              RELAX_RETRY_P=0.0D0
            ELSE IF(ILINESEARCH_EXTRA_MODE.EQ.9) THEN
              RELAX_RETRY_U=1.0D-3*RELAX_BASE_U
              RELAX_RETRY_VW=1.0D-3*RELAX_BASE_VW
              RELAX_RETRY_P=0.0D0
            ELSE
              RELAX_RETRY_U=1.0D-3*RELAX_BASE_U
              RELAX_RETRY_VW=1.0D-3*RELAX_BASE_VW
              RELAX_RETRY_P=1.0D-3*RELAX_BASE_P
            END IF
            WRITE(*,*)'MARQUARDT TRY BLOCK_DAMPED MODE/RELAX_U/VW/P=',ILINESEARCH_EXTRA_MODE,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
            DO IK=1,4
            DO I=1,NX
            DO ITT=1,IAT
              ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
              IF(IS_CONSTRAINED(ISP).NE.0) THEN
                IF(ISP.LE.NX*IAT) THEN
                  U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                ELSE IF(ISP.LE.2*NX*IAT) THEN
                  V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE IF(ISP.LE.3*NX*IAT) THEN
                  W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE
                  P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                END IF
              ELSE IF(ISP.LE.NX*IAT) THEN
                U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
              ELSE IF(ISP.LE.2*NX*IAT) THEN
                V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE IF(ISP.LE.3*NX*IAT) THEN
                W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE
                P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
              END IF
            END DO
            END DO
            END DO
            IPENDING_STEP=1
            GOTO 30
          END IF
          IF(.FALSE..AND.ILINESEARCH_TRIES.GE.IBLOCK_TRY_LIMIT.AND.DLAMBDA_RATIO.GE.0.0D0.AND.(WW4-WW4_PREV).LE.1.0D-6*DMAX1(1.0D-30,DABS(WW4_PREV))) THEN
            DLM_LAMBDA_CURRENT=DLM_LAMBDA_OLD
            WRITE(*,*)'MARQUARDT NONMONOTONE ACCEPT SMALL_RISE WW4_PREV/WW4/REL_RISE/LAMBDA=',WW4_PREV,WW4,(WW4-WW4_PREV)/DMAX1(1.0D-30,DABS(WW4_PREV)),DLM_LAMBDA_CURRENT
            GOTO 7788
          END IF
          RELAX_RETRY_U=RELAX_BASE_U
          RELAX_RETRY_VW=RELAX_BASE_VW
          RELAX_RETRY_P=RELAX_BASE_P
          DO I=1,NX
          DO ITT=1,IAT
            U(ITT,I)=U_SAVE(ITT,I)
            V(ITT,I)=V_SAVE(ITT,I)
            W(ITT,I)=W_SAVE(ITT,I)
            P(ITT,I)=P_SAVE(ITT,I)
          END DO
          END DO
          IF(IOPPOSITE_TRIED.EQ.0) THEN
            STEP_SIGN=-STEP_SIGN
            IOPPOSITE_TRIED=1
            WRITE(*,*)'MARQUARDT TRUST TRY OPPOSITE SIGN R/LAMBDA/WW4_PREV/WW4=',DLAMBDA_RATIO,DLM_LAMBDA_OLD,WW4_PREV,WW4
            DO IK=1,4
            DO I=1,NX
            DO ITT=1,IAT
              ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
              IF(IS_CONSTRAINED(ISP).NE.0) THEN
                IF(ISP.LE.NX*IAT) THEN
                  U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                ELSE IF(ISP.LE.2*NX*IAT) THEN
                  V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE IF(ISP.LE.3*NX*IAT) THEN
                  W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE
                  P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                END IF
              ELSE IF(ISP.LE.NX*IAT) THEN
                U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
              ELSE IF(ISP.LE.2*NX*IAT) THEN
                V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE IF(ISP.LE.3*NX*IAT) THEN
                W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE
                P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
              END IF
            END DO
            END DO
            END DO
            IPENDING_STEP=1
            GOTO 30
          END IF
          IF(IPK.GE.30.AND.WW4_LINE_BEST.LT.1.0D299.AND.WW4_LINE_BEST.LE.WW4_PREV*(1.0D0+1.0D-6)) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            STEP_SIGN=DBEST_STEP_SIGN
            RELAX_RETRY_U=RELAX_BEST_U
            RELAX_RETRY_VW=RELAX_BEST_VW
            RELAX_RETRY_P=RELAX_BEST_P
            DO IK=1,4
            DO I=1,NX
            DO ITT=1,IAT
              ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
              IF(IS_CONSTRAINED(ISP).NE.0) THEN
                IF(ISP.LE.NX*IAT) THEN
                  U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                ELSE IF(ISP.LE.2*NX*IAT) THEN
                  V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE IF(ISP.LE.3*NX*IAT) THEN
                  W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE
                  P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                END IF
              ELSE IF(ISP.LE.NX*IAT) THEN
                U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
              ELSE IF(ISP.LE.2*NX*IAT) THEN
                V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE IF(ISP.LE.3*NX*IAT) THEN
                W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE
                P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
              END IF
            END DO
            END DO
            END DO
            WW4=WW4_LINE_BEST
            ILINESEARCH_TRIES=IBEST_BACKTRACKS
            WRITE(*,*)'MARQUARDT ACCEPT BEST CANDIDATE WW4_PREV/BEST/SIGN/TRY/RELAX_U/VW/P=',WW4_PREV,WW4,STEP_SIGN,ILINESEARCH_TRIES,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
            GOTO 7788
          END IF
          ITRUST_REJECT_COUNT=ITRUST_REJECT_COUNT+1
          DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,10.0D0*DMAX1(DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT,1.0D0/DMAX1(DLM_TRACE_INV,1.0D-300))))
          WRITE(*,*)'MARQUARDT TRUST REJECT STEP R/LAMBDA_OLD/LAMBDA_NEW=',DLAMBDA_RATIO,DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT
          WRITE(*,*)'MARQUARDT TRUST REJECT COUNT WW4_PREV/WW4=',ITRUST_REJECT_COUNT,WW4_PREV,WW4
          IF(ITRUST_REJECT_COUNT.GE.30) THEN
            WRITE(*,*)'MARQUARDT TRUST REJECT LIMIT REACHED. CHECK JTF SIGN OR J/F CONSISTENCY.'
            STOP
          END IF
          IPENDING_STEP=0
          GOTO 30
        END IF
7788   CONTINUE
        ITRUST_REJECT_COUNT=0
        ERR1_ABS=0.0D0
        ERR1_REL_L2=0.0D0
        ERR1_REL_MAX=0.0D0
        DU_SCALE=0.0D0
        DV_SCALE=0.0D0
        DW_SCALE=0.0D0
        DP_SCALE=0.0D0
        DO I=1,NX
        DO ITT=1,IAT
          DU_SCALE=DMAX1(DU_SCALE,DABS(U(ITT,I)))
          DV_SCALE=DMAX1(DV_SCALE,DABS(V(ITT,I)))
          DW_SCALE=DMAX1(DW_SCALE,DABS(W(ITT,I)))
          DP_SCALE=DMAX1(DP_SCALE,DABS(P(ITT,I)))
        END DO
        END DO
        DU_SCALE=DMAX1(DU_SCALE,1.0D-12)
        DV_SCALE=DMAX1(DV_SCALE,1.0D-12)
        DW_SCALE=DMAX1(DW_SCALE,1.0D-12)
        DP_SCALE=DMAX1(DP_SCALE,1.0D-12)
        DO I=1,N
          IF(I.LE.IAT*NX) THEN
            DDU_USED=RELAX_RETRY_U*DDUVWP(I)
            DDU_REL=DABS(DDU_USED)/DU_SCALE
          ELSE IF(I.LE.2*IAT*NX) THEN
            DDU_USED=RELAX_RETRY_VW*DDUVWP(I)
            DDU_REL=DABS(DDU_USED)/DV_SCALE
          ELSE IF(I.LE.3*IAT*NX) THEN
            DDU_USED=RELAX_RETRY_VW*DDUVWP(I)
            DDU_REL=DABS(DDU_USED)/DW_SCALE
          ELSE
            DDU_USED=RELAX_RETRY_P*DDUVWP(I)
            DDU_REL=DABS(DDU_USED)/DP_SCALE
          END IF
          ERR1_ABS=ERR1_ABS+DDU_USED*DDU_USED
          ERR1_REL_L2=ERR1_REL_L2+DDU_REL*DDU_REL
          ERR1_REL_MAX=DMAX1(ERR1_REL_MAX,DDU_REL)
        END DO
        ERR1_ABS=DSQRT(ERR1_ABS)
        ERR1_REL_L2=DSQRT(ERR1_REL_L2)
        ERR1=ERR1_REL_MAX
        IPK=IPK+1
        DUMIN=1.0D300
        DUMAXV=-1.0D300
        DVMIN=1.0D300
        DVMAXV=-1.0D300
        DWMIN=1.0D300
        DWMAXV=-1.0D300
        DPMIN=1.0D300
        DPMAXV=-1.0D300
        IUMIN=0
        ITUMIN=0
        DO I=1,NX
        DO ITT=1,IAT
          IF(U(ITT,I).LT.DUMIN) THEN
            DUMIN=U(ITT,I)
            IUMIN=I
            ITUMIN=ITT
          END IF
          DUMAXV=DMAX1(DUMAXV,U(ITT,I))
          DVMIN=DMIN1(DVMIN,V(ITT,I))
          DVMAXV=DMAX1(DVMAXV,V(ITT,I))
          DWMIN=DMIN1(DWMIN,W(ITT,I))
          DWMAXV=DMAX1(DWMAXV,W(ITT,I))
          DPMIN=DMIN1(DPMIN,P(ITT,I))
          DPMAXV=DMAX1(DPMAXV,P(ITT,I))
        END DO
        END DO
        WRITE(*,*)'STEP_U_MIN_MAX=',IPK,DUMIN,DUMAXV
        WRITE(*,*)'STEP_V_MIN_MAX=',IPK,DVMIN,DVMAXV
        WRITE(*,*)'STEP_W_MIN_MAX=',IPK,DWMIN,DWMAXV
        WRITE(*,*)'STEP_P_MIN_MAX=',IPK,DPMIN,DPMAXV
        IF(IGAUSS_NEWTON_ONLY.NE.0) THEN
          DLM_LAMBDA_CURRENT=0.0D0
          DLM_LAMBDA_USED=0.0D0
        END IF
        WRITE(*,*)'STEP_LAMBDA_NEXT=',IPK,DLM_LAMBDA_CURRENT
        WRITE(*,*)'STEP_LAMBDA_USED_NEXT=',IPK,DLM_LAMBDA_USED,DLM_LAMBDA_CURRENT
        WRITE(*,*)'STEP_RESIDUAL_SUM=',IPK,WW4**2
        WRITE(*,*)'STEP_RESIDUAL_NORM=',IPK,WW4
        WRITE(*,*)'STEP_BLOCK_RESIDUAL_NORM CONT/V/W/U=',IPK,DSQRT(DRU),DSQRT(DRV),DSQRT(DRW),DSQRT(DRP)
        WRITE(*,*)'STEP_BLOCK_RESIDUAL_MAX CONT/V/W/U=',IPK,DRESMAX_CONT,DRESMAX_V,DRESMAX_W,DRESMAX_U
        WRITE(*,*)'STEP_ERR1_WW4=',IPK,ERR1,WW4
        WRITE(*,*)'STEP_ERR1_ABS_REL_L2_REL_MAX=',IPK,ERR1_ABS,ERR1_REL_L2,ERR1_REL_MAX
        WRITE(*,*)'STEP_SUMMARY=',IPK,' LAMBDA_USED=',DLM_LAMBDA_USED,' LAMBDA_NEXT=',DLM_LAMBDA_CURRENT,' ERR1=',ERR1,' WW4=',WW4
        WRITE(*,*)'STEP_SUMMARY_BLOCK_NORM CONT/V/W/U=',IPK,DSQRT(DRU),DSQRT(DRV),DSQRT(DRW),DSQRT(DRP)
        WRITE(*,*)'STEP_SUMMARY_UVWP_MINMAX=',IPK,DUMIN,DUMAXV,DVMIN,DVMAXV,DWMIN,DWMAXV,DPMIN,DPMAXV
        IF(WW4.LE.1.0D-8) THEN
          WRITE(*,*)'NEWTON_CONVERGED_BY_WW4 IPK/WW4=',IPK,WW4
          STOP
        END IF
        IF(IPK.GE.INEWTON_MAX) THEN
          WRITE(*,*)'NEWTON_MAX_REACHED IPK/MAX/WW4=',IPK,INEWTON_MAX,WW4
          STOP
        END IF
        IPENDING_STEP=0
      END IF

      IF(0.NE.0.AND.IPENDING_STEP.NE.0) THEN
        DUMIN=1.0D300
        IUMIN=0
        ITUMIN=0
        DO I=1,NX
        DO ITT=1,IAT
          IF(U(ITT,I).LT.DUMIN) THEN
            DUMIN=U(ITT,I)
            IUMIN=I
            ITUMIN=ITT
          END IF
        END DO
        END DO
        ION_WALL=0
        ION_INLET=0
        ION_POUT=0
        DO I3=1,NSB1
          IF(NSBU(I3).EQ.IUMIN) ION_WALL=1
        END DO
        DO I3=1,NSB2
          IF(NSBUI(I3).EQ.IUMIN) ION_INLET=1
        END DO
        DO I3=1,NSB3
          IF(NSBUO(I3).EQ.IUMIN) ION_POUT=1
        END DO
        IF(DUMIN.LT.0.0D0) THEN
          WRITE(*,*)'LINESEARCH_TRIAL_U_MIN NODE/TIME/XYZ/WALL/INLET/POUT=',IUMIN,ITUMIN,XX(IUMIN),YY(IUMIN),ZZ(IUMIN),ION_WALL,ION_INLET,ION_POUT
          WRITE(*,*)'LINESEARCH_TRIAL_U_MIN_VALUE=',DUMIN
        END IF
        DWW4_TRIAL_REL_DROP=(WW4_PREV-WW4)/DMAX1(1.0D0,DABS(WW4_PREV))
        IF(WW4.LT.WW4_PREV) THEN
          ILINESEARCH_ACCEPTED=ILINESEARCH_ACCEPTED+1
          IF(ILINESEARCH_ACCEPTED.EQ.1.OR.WW4.LT.WW4_LINE_BEST) THEN
            WW4_LINE_BEST=WW4
            RELAX_BEST_U=RELAX_RETRY_U
            RELAX_BEST_VW=RELAX_RETRY_VW
            RELAX_BEST_P=RELAX_RETRY_P
            IBEST_BACKTRACKS=ILINESEARCH_TRIES
          END IF
          WRITE(*,*)'LINESEARCH ACCEPT TRY WW4_PREV/WW4/BEST/RELAX_U/VW/P=',WW4_PREV,WW4,WW4_LINE_BEST,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
          IF(ILINESEARCH_ACCEPTED.LT.3.AND.ILINESEARCH_TRIES.LT.12) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            IF(ILINESEARCH_ACCEPTED.EQ.1.AND.ILINESEARCH_TRIES.GT.0.AND.ILINESEARCH_EXTRA_MODE.EQ.0) THEN
              RELAX_RETRY_U=1.5D0*RELAX_BEST_U
              RELAX_RETRY_VW=1.5D0*RELAX_BEST_VW
              RELAX_RETRY_P=1.5D0*RELAX_BEST_P
              ILINESEARCH_EXTRA_MODE=1
            ELSE
              RELAX_RETRY_U=0.5D0*RELAX_RETRY_U
              RELAX_RETRY_VW=0.5D0*RELAX_RETRY_VW
              RELAX_RETRY_P=0.5D0*RELAX_RETRY_P
              ILINESEARCH_EXTRA_MODE=2
            END IF
            ILINESEARCH_TRIES=ILINESEARCH_TRIES+1
            WRITE(*,*)'LINESEARCH EXTRA TRY MODE/BEST/RELAX_U/VW/P=',ILINESEARCH_EXTRA_MODE,WW4_LINE_BEST,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
            DO IK=1,4
            DO I=1,NX
            DO ITT=1,IAT
              ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
              IF(IS_CONSTRAINED(ISP).NE.0) THEN
                IF(ISP.LE.NX*IAT) THEN
                  U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                ELSE IF(ISP.LE.2*NX*IAT) THEN
                  V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE IF(ISP.LE.3*NX*IAT) THEN
                  W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE
                  P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                END IF
              ELSE IF(ISP.LE.NX*IAT) THEN
                U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
              ELSE IF(ISP.LE.2*NX*IAT) THEN
                V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE IF(ISP.LE.3*NX*IAT) THEN
                W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE
                P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
              END IF
            END DO
            END DO
            END DO
            GO TO 30
          END IF
          DO I=1,NX
          DO ITT=1,IAT
            U(ITT,I)=U_SAVE(ITT,I)
            V(ITT,I)=V_SAVE(ITT,I)
            W(ITT,I)=W_SAVE(ITT,I)
            P(ITT,I)=P_SAVE(ITT,I)
          END DO
          END DO
          RELAX_RETRY_U=RELAX_BEST_U
          RELAX_RETRY_VW=RELAX_BEST_VW
          RELAX_RETRY_P=RELAX_BEST_P
          DO IK=1,4
          DO I=1,NX
          DO ITT=1,IAT
            ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
            IF(IS_CONSTRAINED(ISP).NE.0) THEN
              IF(ISP.LE.NX*IAT) THEN
                U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
              ELSE IF(ISP.LE.2*NX*IAT) THEN
                V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE IF(ISP.LE.3*NX*IAT) THEN
                W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE
                P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
              END IF
            ELSE IF(ISP.LE.NX*IAT) THEN
              U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
            ELSE IF(ISP.LE.2*NX*IAT) THEN
              V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
            ELSE IF(ISP.LE.3*NX*IAT) THEN
              W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
            ELSE
              P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
            END IF
          END DO
          END DO
          END DO
          WW4=WW4_LINE_BEST
          ILINESEARCH_TRIES=IBEST_BACKTRACKS
          WRITE(*,*)'LINESEARCH ACCEPT BEST WW4_PREV/WW4/RELAX_U/VW/P=',WW4_PREV,WW4,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
          IF(RELAX_BASE_U.GT.0.0D0) THEN
            DRELAX_START_SCALE=DMIN1(1.0D0,DMAX1(1.0D-5,2.0D0*RELAX_RETRY_U/RELAX_BASE_U))
          END IF
          DWW4_DROP=WW4_PREV-WW4
          DWW4_REL_DROP=DWW4_DROP/DMAX1(1.0D-30,DABS(WW4_PREV))
          IF(DWW4_DROP.GT.0.0D0) THEN
            IGN_DECREASE_COUNT=IGN_DECREASE_COUNT+1
          ELSE
            IGN_DECREASE_COUNT=0
            IF(IGN_ACTIVE.NE.0) THEN
              DGN_LAMBDA_SWITCH=DMAX1(1.0D-12,DGN_LAMBDA_SWITCH/10.0D0)
              WRITE(*,*)'GAUSS_NEWTON_TRIAL FAILED LOWER_SWITCH=',DGN_LAMBDA_SWITCH
            END IF
          END IF
          DLM_LAMBDA_OLD=DLM_LAMBDA_CURRENT
          DS_OLD=0.0D0
          DS_NEW=0.0D0
          DS_PRED=0.0D0
          DPHYS_OLD=0.0D0
          DPHYS_NEW=0.0D0
          DPHYS_PRED=0.0D0
          DLM_FDOTJD=0.0D0
          DLM_JD_NORM2=0.0D0
          DLM_TRACE=0.0D0
          DLM_TRACE_INV=0.0D0
          DO I=1,N
            DAXI=STEP_SIGN*(RELAX_RETRY_U*CG_AP(I)+RELAX_RETRY_VW*CG_R(I)+RELAX_RETRY_P*CG_Z(I))
            DAXI_SCALED=DAXI*ROW_SCALE_BASE(I)
            DFSCALED=FFS_PHYS_BASE(I)*ROW_SCALE_BASE(I)
            DS_OLD=DS_OLD+DFSCALED*DFSCALED
            DS_NEW=DS_NEW+(FFS_PHYS(I)*ROW_SCALE_BASE(I))*(FFS_PHYS(I)*ROW_SCALE_BASE(I))
            DS_PRED=DS_PRED+(DFSCALED+DAXI_SCALED)*(DFSCALED+DAXI_SCALED)
            DLM_TRACE=DLM_TRACE+NORMAL_DIAG(I)
            IF(NORMAL_DIAG(I).GT.1.0D-300) DLM_TRACE_INV=DLM_TRACE_INV+1.0D0/NORMAL_DIAG(I)
            DPHYS_JD=DAXI
            DPHYS_OLD=DPHYS_OLD+FFS_PHYS_BASE(I)*FFS_PHYS_BASE(I)
            DPHYS_NEW=DPHYS_NEW+FFS_PHYS(I)*FFS_PHYS(I)
            DPHYS_PRED=DPHYS_PRED+(FFS_PHYS_BASE(I)+DPHYS_JD)*(FFS_PHYS_BASE(I)+DPHYS_JD)
            DLM_FDOTJD=DLM_FDOTJD+DFSCALED*DAXI_SCALED
            DLM_JD_NORM2=DLM_JD_NORM2+DAXI_SCALED*DAXI_SCALED
          END DO
          DLM_PRED_DROP_NORM=DPHYS_OLD-DPHYS_PRED
          DLM_PRED_DROP=-2.0D0*DLM_FDOTJD-DLM_JD_NORM2
          DLM_ACT_DROP=DS_OLD-DS_NEW
          IF(DABS(DLM_PRED_DROP).GT.1.0D-300) THEN
            DLAMBDA_RATIO=DLM_ACT_DROP/DLM_PRED_DROP
          ELSE
            DLAMBDA_RATIO=-50.0D0
            WRITE(*,*)'MARQUARDT SALS_RATIO_DIAG INVALID PRED_DROP_ZERO ACT_DROP=',DLM_PRED_DROP,DLM_ACT_DROP
          END IF
          ILAMBDA_REASON=12
          DLM_LAMBDA_C=1.0D0/DMAX1(DLM_TRACE_INV,1.0D-300)
          IF(DWW4_DROP.GT.0.0D0) THEN
            IF(Re.GE.1.0D5.AND.DWW4_REL_DROP.LT.1.0D-4) THEN
              DLAMBDA_FACTOR=1.0D0
              ILAMBDA_REASON=34
            ELSE
              DLAMBDA_FACTOR=1.0D0/1.2D0
              DLM_LAMBDA_CURRENT=DLM_LAMBDA_CURRENT*DLAMBDA_FACTOR
              DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,DLM_LAMBDA_CURRENT))
              ISTAG_LAMBDA_COUNT=0
              ILAMBDA_REASON=31
            END IF
          ELSE
            DLAMBDA_FACTOR=2.0D0
            IF(DLM_LAMBDA_CURRENT.LE.0.0D0) DLM_LAMBDA_CURRENT=DMAX1(DLM_LAMBDA_C,1.0D0)
            DLM_LAMBDA_CURRENT=DLM_LAMBDA_CURRENT*DLAMBDA_FACTOR
            DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,DLM_LAMBDA_CURRENT))
            DRELAX_START_SCALE=1.0D0
            ILAMBDA_REASON=32
          END IF
          IF(.FALSE..AND.DWW4_REL_DROP.LT.1.0D-6.AND.ILINESEARCH_TRIES.GE.12.AND.ISTAG_ACCEPT_STEP.NE.IPK) THEN
            ISTAG_LAMBDA_COUNT=ISTAG_LAMBDA_COUNT+1
            ISTAG_ACCEPT_STEP=IPK
            IF(ISTAG_LAMBDA_COUNT.LE.1) THEN
              DLAMBDA_FACTOR=4.0D0
            ELSE IF(ISTAG_LAMBDA_COUNT.EQ.2) THEN
              DLAMBDA_FACTOR=10.0D0
            ELSE
              DLAMBDA_FACTOR=100.0D0
            END IF
            DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,DLAMBDA_FACTOR*DMAX1(DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT)))
            DRELAX_START_SCALE=1.0D0
            ILAMBDA_REASON=25
            ISTAG_LAMBDA_HOLD=3
            WRITE(*,*)'MARQUARDT STAGNATION_ACCEPT INCREASE_LAMBDA COUNT/FACTOR/REL_DROP/BACKTRACKS=',ISTAG_LAMBDA_COUNT,DLAMBDA_FACTOR,DWW4_REL_DROP,ILINESEARCH_TRIES
          END IF
          WRITE(*,*)'MARQUARDT SALS_RATIO_DIAG ACT/PRED/R/S_OLD/S_NEW/S_PRED=',DLM_ACT_DROP,DLM_PRED_DROP,DLAMBDA_RATIO,DS_OLD,DS_NEW,DS_PRED
          WRITE(*,*)'MARQUARDT SALS_RATIO_DIAG_PHYS ACT/PRED/R/S_OLD/S_NEW/S_PRED=',DLM_ACT_DROP,DLM_PRED_DROP,DLAMBDA_RATIO,DPHYS_OLD,DPHYS_NEW,DPHYS_PRED
          WRITE(*,*)'MARQUARDT LINEAR_MODEL_DIAG PRED/NORMDIFF/FTJD/JDNORM2=',DLM_PRED_DROP,DLM_PRED_DROP_NORM,DLM_FDOTJD,DLM_JD_NORM2
          WRITE(*,*)'MARQUARDT MEASURED_LAMBDA_CONTROL DROP/REL_DROP=',DWW4_DROP,DWW4_REL_DROP
          WRITE(*,*)'GAUSS_NEWTON_TRIAL READINESS COUNT/LAMBDA/SWITCH=',IGN_DECREASE_COUNT,DLM_LAMBDA_CURRENT,DGN_LAMBDA_SWITCH
          WRITE(*,*)'MARQUARDT LAMBDA UPDATE OLD/NEW/REASON/REL_DROP/BACKTRACKS=',DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT,ILAMBDA_REASON,DWW4_REL_DROP,ILINESEARCH_TRIES
          WRITE(*,*)'MARQUARDT NEXT RELAX START SCALE=',DRELAX_START_SCALE
          IPENDING_STEP=0
          ILINESEARCH_TRIES=0
          ILINESEARCH_ACCEPTED=0
          ILINESEARCH_EXTRA_MODE=0
        ELSE
          IF(ILINESEARCH_ACCEPTED.GT.0) THEN
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            RELAX_RETRY_U=RELAX_BEST_U
            RELAX_RETRY_VW=RELAX_BEST_VW
            RELAX_RETRY_P=RELAX_BEST_P
            DO IK=1,4
            DO I=1,NX
            DO ITT=1,IAT
              ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
              IF(IS_CONSTRAINED(ISP).NE.0) THEN
                IF(ISP.LE.NX*IAT) THEN
                  U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                ELSE IF(ISP.LE.2*NX*IAT) THEN
                  V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE IF(ISP.LE.3*NX*IAT) THEN
                  W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE
                  P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                END IF
              ELSE IF(ISP.LE.NX*IAT) THEN
                U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
              ELSE IF(ISP.LE.2*NX*IAT) THEN
                V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE IF(ISP.LE.3*NX*IAT) THEN
                W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE
                P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
              END IF
            END DO
            END DO
            END DO
            WW4=WW4_LINE_BEST
            ILINESEARCH_TRIES=IBEST_BACKTRACKS
            WRITE(*,*)'LINESEARCH ACCEPT BEST AFTER RISE WW4_PREV/WW4/RELAX_U/VW/P=',WW4_PREV,WW4,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
            IF(RELAX_BASE_U.GT.0.0D0) THEN
              DRELAX_START_SCALE=DMIN1(1.0D0,DMAX1(1.0D-5,2.0D0*RELAX_RETRY_U/RELAX_BASE_U))
            END IF
            DWW4_DROP=WW4_PREV-WW4
            DWW4_REL_DROP=DWW4_DROP/DMAX1(1.0D-30,DABS(WW4_PREV))
            IF(DWW4_DROP.GT.0.0D0) THEN
              IGN_DECREASE_COUNT=IGN_DECREASE_COUNT+1
            ELSE
              IGN_DECREASE_COUNT=0
              IF(IGN_ACTIVE.NE.0) THEN
                DGN_LAMBDA_SWITCH=DMAX1(1.0D-12,DGN_LAMBDA_SWITCH/10.0D0)
                WRITE(*,*)'GAUSS_NEWTON_TRIAL FAILED LOWER_SWITCH=',DGN_LAMBDA_SWITCH
              END IF
            END IF
            DLM_LAMBDA_OLD=DLM_LAMBDA_CURRENT
            DS_OLD=0.0D0
            DS_NEW=0.0D0
            DS_PRED=0.0D0
            DPHYS_OLD=0.0D0
            DPHYS_NEW=0.0D0
            DPHYS_PRED=0.0D0
            DLM_FDOTJD=0.0D0
            DLM_JD_NORM2=0.0D0
            DLM_TRACE=0.0D0
            DLM_TRACE_INV=0.0D0
            DO I=1,N
              DAXI=STEP_SIGN*(RELAX_RETRY_U*CG_AP(I)+RELAX_RETRY_VW*CG_R(I)+RELAX_RETRY_P*CG_Z(I))
              DAXI_SCALED=DAXI*ROW_SCALE_BASE(I)
              DFSCALED=FFS_PHYS_BASE(I)*ROW_SCALE_BASE(I)
              DS_OLD=DS_OLD+DFSCALED*DFSCALED
              DS_NEW=DS_NEW+(FFS_PHYS(I)*ROW_SCALE_BASE(I))*(FFS_PHYS(I)*ROW_SCALE_BASE(I))
              DS_PRED=DS_PRED+(DFSCALED+DAXI_SCALED)*(DFSCALED+DAXI_SCALED)
              DLM_TRACE=DLM_TRACE+NORMAL_DIAG(I)
              IF(NORMAL_DIAG(I).GT.1.0D-300) DLM_TRACE_INV=DLM_TRACE_INV+1.0D0/NORMAL_DIAG(I)
              DPHYS_JD=DAXI
              DPHYS_OLD=DPHYS_OLD+FFS_PHYS_BASE(I)*FFS_PHYS_BASE(I)
              DPHYS_NEW=DPHYS_NEW+FFS_PHYS(I)*FFS_PHYS(I)
              DPHYS_PRED=DPHYS_PRED+(FFS_PHYS_BASE(I)+DPHYS_JD)*(FFS_PHYS_BASE(I)+DPHYS_JD)
              DLM_FDOTJD=DLM_FDOTJD+DFSCALED*DAXI_SCALED
              DLM_JD_NORM2=DLM_JD_NORM2+DAXI_SCALED*DAXI_SCALED
            END DO
            DLM_PRED_DROP_NORM=DPHYS_OLD-DPHYS_PRED
            DLM_PRED_DROP=-2.0D0*DLM_FDOTJD-DLM_JD_NORM2
            DLM_ACT_DROP=DS_OLD-DS_NEW
            IF(DABS(DLM_PRED_DROP).GT.1.0D-300) THEN
              DLAMBDA_RATIO=DLM_ACT_DROP/DLM_PRED_DROP
            ELSE
              DLAMBDA_RATIO=-50.0D0
              WRITE(*,*)'MARQUARDT SALS RISE RATIO INVALID PRED_DROP_ZERO ACT_DROP=',DLM_PRED_DROP,DLM_ACT_DROP
            END IF
            ILAMBDA_REASON=12
            DLM_LAMBDA_C=1.0D0/DMAX1(DLM_TRACE_INV,1.0D-300)
            IF(DWW4_DROP.GT.0.0D0) THEN
              IF(Re.GE.1.0D5.AND.DWW4_REL_DROP.LT.1.0D-4) THEN
                DLAMBDA_FACTOR=1.0D0
                ILAMBDA_REASON=34
              ELSE
                DLAMBDA_FACTOR=1.0D0/1.2D0
                DLM_LAMBDA_CURRENT=DLM_LAMBDA_CURRENT*DLAMBDA_FACTOR
                DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,DLM_LAMBDA_CURRENT))
                ISTAG_LAMBDA_COUNT=0
                ILAMBDA_REASON=31
              END IF
            ELSE
              DLAMBDA_FACTOR=2.0D0
              IF(DLM_LAMBDA_CURRENT.LE.0.0D0) DLM_LAMBDA_CURRENT=DMAX1(DLM_LAMBDA_C,1.0D0)
              DLM_LAMBDA_CURRENT=DLM_LAMBDA_CURRENT*DLAMBDA_FACTOR
              DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,DLM_LAMBDA_CURRENT))
              DRELAX_START_SCALE=1.0D0
              ILAMBDA_REASON=32
            END IF
            IF(.FALSE..AND.DWW4_REL_DROP.LT.1.0D-6.AND.ILINESEARCH_TRIES.GE.12.AND.ISTAG_ACCEPT_STEP.NE.IPK) THEN
              ISTAG_LAMBDA_COUNT=ISTAG_LAMBDA_COUNT+1
              ISTAG_ACCEPT_STEP=IPK
              IF(ISTAG_LAMBDA_COUNT.LE.1) THEN
                DLAMBDA_FACTOR=4.0D0
              ELSE IF(ISTAG_LAMBDA_COUNT.EQ.2) THEN
                DLAMBDA_FACTOR=10.0D0
              ELSE
                DLAMBDA_FACTOR=100.0D0
              END IF
              DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,DLAMBDA_FACTOR*DMAX1(DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT)))
              DRELAX_START_SCALE=1.0D0
              ILAMBDA_REASON=25
              ISTAG_LAMBDA_HOLD=3
              WRITE(*,*)'MARQUARDT STAGNATION_ACCEPT INCREASE_LAMBDA COUNT/FACTOR/REL_DROP/BACKTRACKS=',ISTAG_LAMBDA_COUNT,DLAMBDA_FACTOR,DWW4_REL_DROP,ILINESEARCH_TRIES
            END IF
            WRITE(*,*)'MARQUARDT SALS_RATIO_DIAG ACT/PRED/R/S_OLD/S_NEW/S_PRED=',DLM_ACT_DROP,DLM_PRED_DROP,DLAMBDA_RATIO,DS_OLD,DS_NEW,DS_PRED
            WRITE(*,*)'MARQUARDT SALS_RATIO_DIAG_PHYS ACT/PRED/R/S_OLD/S_NEW/S_PRED=',DLM_ACT_DROP,DLM_PRED_DROP,DLAMBDA_RATIO,DPHYS_OLD,DPHYS_NEW,DPHYS_PRED
            WRITE(*,*)'MARQUARDT LINEAR_MODEL_DIAG PRED/NORMDIFF/FTJD/JDNORM2=',DLM_PRED_DROP,DLM_PRED_DROP_NORM,DLM_FDOTJD,DLM_JD_NORM2
            WRITE(*,*)'MARQUARDT MEASURED_LAMBDA_CONTROL DROP/REL_DROP=',DWW4_DROP,DWW4_REL_DROP
            WRITE(*,*)'GAUSS_NEWTON_TRIAL READINESS COUNT/LAMBDA/SWITCH=',IGN_DECREASE_COUNT,DLM_LAMBDA_CURRENT,DGN_LAMBDA_SWITCH
            WRITE(*,*)'MARQUARDT LAMBDA UPDATE OLD/NEW/REASON/REL_DROP/BACKTRACKS=',DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT,ILAMBDA_REASON,DWW4_REL_DROP,ILINESEARCH_TRIES
            WRITE(*,*)'MARQUARDT NEXT RELAX START SCALE=',DRELAX_START_SCALE
            IPENDING_STEP=0
            ILINESEARCH_TRIES=0
            ILINESEARCH_ACCEPTED=0
            ILINESEARCH_EXTRA_MODE=0
          ELSE
          IF(ILINESEARCH_TRIES.GE.12) THEN
            IF(IOPPOSITE_TRIED.EQ.0) THEN
              DO I=1,NX
              DO ITT=1,IAT
                U(ITT,I)=U_SAVE(ITT,I)
                V(ITT,I)=V_SAVE(ITT,I)
                W(ITT,I)=W_SAVE(ITT,I)
                P(ITT,I)=P_SAVE(ITT,I)
              END DO
              END DO
              STEP_SIGN=-STEP_SIGN
              IOPPOSITE_TRIED=1
              ILINESEARCH_TRIES=0
              RELAX_RETRY_U=RELAX_BASE_U*DRELAX_START_SCALE
              RELAX_RETRY_VW=RELAX_BASE_VW*DRELAX_START_SCALE
              RELAX_RETRY_P=RELAX_BASE_P*DRELAX_START_SCALE
              WRITE(*,*)'LINESEARCH TRY OPPOSITE SIGN RELAX_U/VW/P=',RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
              DO IK=1,4
              DO I=1,NX
              DO ITT=1,IAT
                ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
                IF(IS_CONSTRAINED(ISP).NE.0) THEN
                  IF(ISP.LE.NX*IAT) THEN
                    U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                  ELSE IF(ISP.LE.2*NX*IAT) THEN
                    V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                  ELSE IF(ISP.LE.3*NX*IAT) THEN
                    W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                  ELSE
                    P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                  END IF
                ELSE IF(ISP.LE.NX*IAT) THEN
                  U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
                ELSE IF(ISP.LE.2*NX*IAT) THEN
                  V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE IF(ISP.LE.3*NX*IAT) THEN
                  W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
                ELSE
                  P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
                END IF
              END DO
              END DO
              END DO
              GO TO 30
            END IF
            DO I=1,NX
            DO ITT=1,IAT
              U(ITT,I)=U_SAVE(ITT,I)
              V(ITT,I)=V_SAVE(ITT,I)
              W(ITT,I)=W_SAVE(ITT,I)
              P(ITT,I)=P_SAVE(ITT,I)
            END DO
            END DO
            DLM_LAMBDA_OLD=DLM_LAMBDA_CURRENT
            DLM_TRACE=0.0D0
            DO I=1,N
              DLM_TRACE=DLM_TRACE+NORMAL_DIAG(I)
            END DO
            IF(DLM_LAMBDA_CURRENT.LE.0.0D0) THEN
              DLM_LAMBDA_CURRENT=1.0D0/DMAX1(DLM_TRACE,1.0D-300)
            ELSE
              DLM_LAMBDA_CURRENT=10.0D0*DLM_LAMBDA_CURRENT
            END IF
            DLM_LAMBDA_CURRENT=DMIN1(DLM_LAMBDA_MAX,DMAX1(DLM_LAMBDA_MIN,DLM_LAMBDA_CURRENT))
            DRELAX_START_SCALE=1.0D0
            IPENDING_STEP=0
            ILINESEARCH_TRIES=0
            WRITE(*,*)'MARQUARDT LAMBDA UPDATE OLD/NEW/REASON/REL_DROP/BACKTRACKS=',DLM_LAMBDA_OLD,DLM_LAMBDA_CURRENT,9,0.0D0,ILINESEARCH_TRIES
            WRITE(*,*)'MARQUARDT NEXT RELAX START SCALE=',DRELAX_START_SCALE
            GOTO 30
          END IF
          DO I=1,NX
          DO ITT=1,IAT
            U(ITT,I)=U_SAVE(ITT,I)
            V(ITT,I)=V_SAVE(ITT,I)
            W(ITT,I)=W_SAVE(ITT,I)
            P(ITT,I)=P_SAVE(ITT,I)
          END DO
          END DO
          RELAX_RETRY_U=0.5D0*RELAX_RETRY_U
          RELAX_RETRY_VW=0.5D0*RELAX_RETRY_VW
          RELAX_RETRY_P=0.5D0*RELAX_RETRY_P
          ILINESEARCH_TRIES=ILINESEARCH_TRIES+1
          WRITE(*,*)'LINESEARCH REJECT WW4_PREV/WW4 TRY RELAX_U/VW/P=',WW4_PREV,WW4,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
          DO IK=1,4
          DO I=1,NX
          DO ITT=1,IAT
            ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
            IF(IS_CONSTRAINED(ISP).NE.0) THEN
              IF(ISP.LE.NX*IAT) THEN
                U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
              ELSE IF(ISP.LE.2*NX*IAT) THEN
                V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE IF(ISP.LE.3*NX*IAT) THEN
                W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
              ELSE
                P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
              END IF
            ELSE IF(ISP.LE.NX*IAT) THEN
              U(ITT,I)=U(ITT,I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(ISP)
            ELSE IF(ISP.LE.2*NX*IAT) THEN
              V(ITT,I)=V(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
            ELSE IF(ISP.LE.3*NX*IAT) THEN
              W(ITT,I)=W(ITT,I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(ISP)
            ELSE
              P(ITT,I)=P(ITT,I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(ISP)
            END IF
          END DO
          END DO
          END DO
          GO TO 30
          END IF
        END IF
      END IF

      GO TO 126

      IZERO_ROW=0
      IZERO_COL=0
      ISMALL_DIAG=0
      IP_ZERO_ROW=0
      IP_ZERO_COL=0
      IP_SMALL_DIAG=0
      IU_ZERO_ROW=0
      IV_ZERO_ROW=0
      IW_ZERO_ROW=0
      IU_SMALL_DIAG=0
      IV_SMALL_DIAG=0
      IW_SMALL_DIAG=0
      DMINDIAG=1.0D300
      DMAXDIAG=0.0D0
      DMINROWSUM=1.0D300
      DMINCOLSUM=1.0D300
      DO I=1,4*NX*IAT
        DROWSUM=0.0D0
        DCOLSUM=0.0D0
        DO J=1,4*NX*IAT
          DROWSUM=DROWSUM+DABS(SSMMK(I,J))
          DCOLSUM=DCOLSUM+DABS(SSMMK(J,I))
        END DO
        DMINDIAG=DMIN1(DMINDIAG,DABS(SSMMK(I,I)))
        DMAXDIAG=DMAX1(DMAXDIAG,DABS(SSMMK(I,I)))
        DMINROWSUM=DMIN1(DMINROWSUM,DROWSUM)
        DMINCOLSUM=DMIN1(DMINCOLSUM,DCOLSUM)
        IF(DROWSUM.LT.1.0D-12) IZERO_ROW=IZERO_ROW+1
        IF(DCOLSUM.LT.1.0D-12) IZERO_COL=IZERO_COL+1
        IF(DABS(SSMMK(I,I)).LT.1.0D-12) ISMALL_DIAG=ISMALL_DIAG+1
        IF(I.LE.NX*IAT) THEN
          IF(DROWSUM.LT.1.0D-12) IU_ZERO_ROW=IU_ZERO_ROW+1
          IF(DABS(SSMMK(I,I)).LT.1.0D-12) IU_SMALL_DIAG=IU_SMALL_DIAG+1
        ELSE IF(I.LE.2*NX*IAT) THEN
          IF(DROWSUM.LT.1.0D-12) IV_ZERO_ROW=IV_ZERO_ROW+1
          IF(DABS(SSMMK(I,I)).LT.1.0D-12) IV_SMALL_DIAG=IV_SMALL_DIAG+1
        ELSE IF(I.LE.3*NX*IAT) THEN
          IF(DROWSUM.LT.1.0D-12) IW_ZERO_ROW=IW_ZERO_ROW+1
          IF(DABS(SSMMK(I,I)).LT.1.0D-12) IW_SMALL_DIAG=IW_SMALL_DIAG+1
        ELSE
          IF(DROWSUM.LT.1.0D-12) IP_ZERO_ROW=IP_ZERO_ROW+1
          IF(DCOLSUM.LT.1.0D-12) IP_ZERO_COL=IP_ZERO_COL+1
          IF(DABS(SSMMK(I,I)).LT.1.0D-12) IP_SMALL_DIAG=IP_SMALL_DIAG+1
          IF(DROWSUM.LT.1.0D-12) THEN
            ILOCAL=I-3*NX*IAT
            IPNODE=(ILOCAL-1)/IAT+1
            IPTIME=ILOCAL-(IPNODE-1)*IAT
            WRITE(*,*)'ZERO PRESSURE ROW INDEX/NODE/TIME/RHS=',I,IPNODE,IPTIME,FFS(I)
          END IF
        END IF
      END DO
      WRITE(*,*)'MATRIX ZERO ROW/COL/SMALL DIAG=',IZERO_ROW,IZERO_COL,ISMALL_DIAG
      WRITE(*,*)'U/V/W/P ZERO ROWS=',IU_ZERO_ROW,IV_ZERO_ROW,IW_ZERO_ROW,IP_ZERO_ROW
      WRITE(*,*)'U/V/W/P SMALL DIAG=',IU_SMALL_DIAG,IV_SMALL_DIAG,IW_SMALL_DIAG,IP_SMALL_DIAG
      WRITE(*,*)'PRESSURE ZERO ROW/COL/SMALL DIAG=',IP_ZERO_ROW,IP_ZERO_COL,IP_SMALL_DIAG
      WRITE(*,*)'MATRIX MIN/MAX DIAG MIN ROW/COL SUM=',DMINDIAG,DMAXDIAG,DMINROWSUM,DMINCOLSUM

      IPINNED_ROW=0
      DO I=1,4*NX*IAT
        DROWSUM=0.0D0
        DO J=1,4*NX*IAT
          DROWSUM=DROWSUM+DABS(SSMMK(I,J))
        END DO
        IF(DROWSUM.LT.1.0D-12) THEN
          DO J=1,4*NX*IAT
            SSMMK(I,J)=0.0D0
            SSMMK(J,I)=0.0D0
          END DO
          SSMMK(I,I)=1.0D0
          FFS(I)=0.0D0
          IPINNED_ROW=IPINNED_ROW+1
        END IF
      END DO
      IF(IPINNED_ROW.GT.0) WRITE(*,*)'PINNED ZERO MATRIX ROW/COL=',IPINNED_ROW
       
      n=nx*iat*4
  
      IF(ALLOCATED(ir)) DEALLOCATE(ir)
      IF(ALLOCATED(ic)) DEALLOCATE(ic)
      IF(ALLOCATED(a)) DEALLOCATE(a)
      call triplet_init(SPARSE_A, n, n, max(1, min(1000000, 20*n)))
      DO I=1,N
        DO J=1,N
          IF(SSMMK(I,J).NE.0.0D0) THEN
            call triplet_add(SPARSE_A, I, J, SSMMK(I,J))
          END IF
        END DO
      END DO
      call triplet_to_csr(SPARSE_A, ir, ic, a, mn)
      mnr=n+1
      
      write(*,*)'clear 99 !mnr,mn',mnr,mn    

      IF(ALLOCATED(DENSE_A)) THEN
        DENSE_A=0.0D0
        DO I=1,N
          DO J=ir(I),ir(I+1)-1
            DENSE_A(I,ic(J))=DENSE_A(I,ic(J))+a(J)
          END DO
        END DO
        DCSRMAX=0.0D0
        DCSRSUM=0.0D0
        DO I=1,N
          DO J=1,N
            DCSRMAX=DMAX1(DCSRMAX,DABS(DENSE_A(I,J)-SSMMK(I,J)))
            DCSRSUM=DCSRSUM+DABS(DENSE_A(I,J)-SSMMK(I,J))
          END DO
        END DO
        WRITE(*,*)'CSR_RECONSTRUCT_MAX_SUM_DIFF=',DCSRMAX,DCSRSUM
      END IF
      
  126 CONTINUE
      do i=1,64
        ipt(i)=0
        iparm(i)=0
      end do
      iparm(1)=1
      iparm(2)=2
      iparm(8)=50
      iparm(10)=13
      iparm(11)=1
      iparm(13)=1
      maxfct=1
      mnum=1
      mtype=11
      iphase=13 
      nrhs=1     
      msglvl=0
      
      WRITE(*,*)'clear 777 !'
      
      DO I=1,N
        FF1(I)=FFS(I)
      END DO

!     Marquardt_CG is evaluated on the physical, unscaled equations.
!     Ozaki dot products and CSR matvecs are used for accuracy, but no row,
!     column, block, or boundary weighting is applied.
      ROW_SCALE=1.0D0
      COL_SCALE=1.0D0
      DROW_SCALE_CAP=1.0D4
      DRS_MIN1=1.0D300
      DRS_MIN2=1.0D300
      DRS_MIN3=1.0D300
      DRS_MIN4=1.0D300
      DRS_MAX1=0.0D0
      DRS_MAX2=0.0D0
      DRS_MAX3=0.0D0
      DRS_MAX4=0.0D0
      ROW_SCALE=1.0D0
      DPRE_SCALE_MAX=0.0D0
      IPRE_SCALE_ROW=0
      DO I=1,N
        ROW_SCALE(I)=1.0D0
        IF(I.LE.NX*IAT.AND.DABS(FFS(I)).GT.DPRE_SCALE_MAX) THEN
          DPRE_SCALE_MAX=DABS(FFS(I))
          IPRE_SCALE_ROW=I
        END IF
        IBLOCK=(I-1)/(NX*IAT)+1
        IF(IBLOCK.EQ.1) THEN
          DRS_MIN1=DMIN1(DRS_MIN1,ROW_SCALE(I))
          DRS_MAX1=DMAX1(DRS_MAX1,ROW_SCALE(I))
        ELSE IF(IBLOCK.EQ.2) THEN
          DRS_MIN2=DMIN1(DRS_MIN2,ROW_SCALE(I))
          DRS_MAX2=DMAX1(DRS_MAX2,ROW_SCALE(I))
        ELSE IF(IBLOCK.EQ.3) THEN
          DRS_MIN3=DMIN1(DRS_MIN3,ROW_SCALE(I))
          DRS_MAX3=DMAX1(DRS_MAX3,ROW_SCALE(I))
        ELSE
          DRS_MIN4=DMIN1(DRS_MIN4,ROW_SCALE(I))
          DRS_MAX4=DMAX1(DRS_MAX4,ROW_SCALE(I))
        END IF
      END DO
      WRITE(*,*)'ROW_SCALE_MINMAX CONT/V/W/U CAP=',DRS_MIN1,DRS_MAX1,DRS_MIN2,DRS_MAX2,DRS_MIN3,DRS_MAX3,DRS_MIN4,DRS_MAX4,DROW_SCALE_CAP
      IF(IPRE_SCALE_ROW.GT.0) THEN
        IPNODE=(IPRE_SCALE_ROW-1)/IAT+1
        IPTIME=IPRE_SCALE_ROW-(IPNODE-1)*IAT
        WRITE(*,*)'MARQUARDT_PRE_SCALE_CONT_MAX ROW/NODE/TIME/FFS/ROWSCALE/CONSTR/REPLACED=',&
             IPRE_SCALE_ROW,IPNODE,IPTIME,FFS(IPRE_SCALE_ROW),ROW_SCALE(IPRE_SCALE_ROW),&
             IS_CONSTRAINED_ROW(IPRE_SCALE_ROW),IS_REPLACED_ROW(IPRE_SCALE_ROW)
      END IF
      DBLOCK_RMS1=0.0D0
      DBLOCK_RMS2=0.0D0
      DBLOCK_RMS3=0.0D0
      DBLOCK_RMS4=0.0D0
      DBLOCK_TARGET=0.0D0
      DBLOCK_W1=1.0D0
      DBLOCK_W2=1.0D0
      DBLOCK_W3=1.0D0
      DBLOCK_W4=1.0D0
      IF(.FALSE..AND.IMARQUARDT_CG.NE.0.AND.DBLOCK_TARGET.GT.0.0D0) THEN
        DBLOCK_W1=DMIN1(5.0D0,DMAX1(0.2D0,DBLOCK_TARGET/DMAX1(DBLOCK_RMS1,1.0D-30)))
        DBLOCK_W2=DMIN1(5.0D0,DMAX1(0.2D0,DBLOCK_TARGET/DMAX1(DBLOCK_RMS2,1.0D-30)))
        DBLOCK_W3=DMIN1(5.0D0,DMAX1(0.2D0,DBLOCK_TARGET/DMAX1(DBLOCK_RMS3,1.0D-30)))
        DBLOCK_W4=DMIN1(5.0D0,DMAX1(0.2D0,DBLOCK_TARGET/DMAX1(DBLOCK_RMS4,1.0D-30)))
      END IF
      WRITE(*,*)'MARQUARDT PHYSICAL_SCALE BLOCK WEIGHTS CONT/V/W/U=',DBLOCK_W1,DBLOCK_W2,DBLOCK_W3,DBLOCK_W4
      DOUTFLOW_VEL_ROW_WEIGHT=1.0D0
      IF(.FALSE..AND.IMARQUARDT_CG.NE.0) THEN
        DO I3=1,NSB3
          II=NSBUO(I3)
          DO ITT=1,IAT
            IROW=GLOBAL_EQUATION_INDEX(ITT,II,2,NX,IAT)
            ROW_SCALE(IROW)=ROW_SCALE(IROW)*DOUTFLOW_VEL_ROW_WEIGHT
            IROW=GLOBAL_EQUATION_INDEX(ITT,II,3,NX,IAT)
            ROW_SCALE(IROW)=ROW_SCALE(IROW)*DOUTFLOW_VEL_ROW_WEIGHT
            IROW=GLOBAL_EQUATION_INDEX(ITT,II,4,NX,IAT)
            ROW_SCALE(IROW)=ROW_SCALE(IROW)*DOUTFLOW_VEL_ROW_WEIGHT
          END DO
        END DO
      END IF
      IF(IMARQUARDT_CG.NE.0) WRITE(*,*)'MARQUARDT OUTFLOW VELOCITY ROW EXTRA WEIGHT=',DOUTFLOW_VEL_ROW_WEIGHT
      DINLET_P_D2_ROW_WEIGHT=1.0D0
      IF(.FALSE.) THEN
        DO I3=1,NSB2
          II=NSBUI(I3)
          DO ITT=1,IAT
            IROW=GLOBAL_EQUATION_INDEX(ITT,II,4,NX,IAT)
            IF(IS_REPLACED_ROW(IROW).NE.0.AND.IS_CONSTRAINED_ROW(IROW).EQ.0) THEN
              ROW_SCALE(IROW)=ROW_SCALE(IROW)*DINLET_P_D2_ROW_WEIGHT
            END IF
          END DO
        END DO
      END IF
      WRITE(*,*)'INLET PRESSURE D2 ROW EXTRA WEIGHT=',DINLET_P_D2_ROW_WEIGHT
      DO I=1,N
        ROW_SCALE(I)=1.0D0
      END DO
      DPOISSON_ROW_EXTRA_WEIGHT=1.0D0
      IF(IMARQUARDT_CG.NE.0) WRITE(*,*)'MARQUARDT UNSCALED_CG FLAG=',IMARQUARDT_UNSCALED_CG
      IF(IMARQUARDT_UNSCALED_CG.NE.0) WRITE(*,*)'MARQUARDT ROW/COL SCALING DISABLED FOR CG AND SALS_R'
      DSU=0.0D0
      DSV=0.0D0
      DSW=0.0D0
      DSP=0.0D0
      DSMAXU=0.0D0
      DSMAXV=0.0D0
      DSMAXW=0.0D0
      DSMAXP=0.0D0
      ISCALED_MAX_CONT_ROW=0
      DO I=1,NX*IAT
        DSU=DSU+FFS(I)**2
        DSV=DSV+FFS(I+NX*IAT)**2
        DSW=DSW+FFS(I+2*NX*IAT)**2
        DSP=DSP+FFS(I+3*NX*IAT)**2
        IF(DABS(FFS(I)).GT.DSMAXU) THEN
          DSMAXU=DABS(FFS(I))
          ISCALED_MAX_CONT_ROW=I
        END IF
        DSMAXV=DMAX1(DSMAXV,DABS(FFS(I+NX*IAT)))
        DSMAXW=DMAX1(DSMAXW,DABS(FFS(I+2*NX*IAT)))
        DSMAXP=DMAX1(DSMAXP,DABS(FFS(I+3*NX*IAT)))
      END DO
      WRITE(*,*)'MARQUARDT_SCALED_ROW_NORM CONT/V/W/U=',DSQRT(DSU),DSQRT(DSV),DSQRT(DSW),DSQRT(DSP)
      WRITE(*,*)'MARQUARDT_SCALED_ROW_MAX CONT/V/W/U=',DSMAXU,DSMAXV,DSMAXW,DSMAXP
      IF(ISCALED_MAX_CONT_ROW.GT.0) THEN
        IPNODE=(ISCALED_MAX_CONT_ROW-1)/IAT+1
        IPTIME=ISCALED_MAX_CONT_ROW-(IPNODE-1)*IAT
        WRITE(*,*)'MARQUARDT_SCALED_CONT_MAX ROW/NODE/TIME/FFS/ROWSCALE/CONSTR/REPLACED/XYZ=',&
             ISCALED_MAX_CONT_ROW,IPNODE,IPTIME,FFS(ISCALED_MAX_CONT_ROW),ROW_SCALE(ISCALED_MAX_CONT_ROW),&
             IS_CONSTRAINED_ROW(ISCALED_MAX_CONT_ROW),IS_REPLACED_ROW(ISCALED_MAX_CONT_ROW),XX(IPNODE),YY(IPNODE),ZZ(IPNODE)
      END IF
      DO I=1,N
        COL_SCALE(I)=1.0D0
      END DO
!     COL_SCALE is identically one in the physical-scale path; skip the
!     matrix multiply loop.

      IF(INEWTON_PARDISO.NE.0) THEN
        NORMAL_RHS=0.0D0
        DO I=1,N
          NORMAL_RHS(I)=-FFS(I)
          FFS(I)=NORMAL_RHS(I)
          x(I,1)=0.0D0
        END DO
        WRITE(*,*)'NEWTON_PARDISO RHS_SIGN_MINUS'
        call pardiso (ipt, maxfct, mnum, mtype, iphase, n, a, ir, ic, iperm, nrhs, iparm, msglvl, FFS, x, ierror)
        WRITE(*,*)'NEWTON_PARDISO IERROR=',ierror
        WRITE(*,*)'PARDISO IPARM 14 18 19 20 22=',iparm(14),iparm(18),iparm(19),iparm(20),iparm(22)
        IF(ierror.NE.0) THEN
          WRITE(*,*)'NEWTON_PARDISO FAILED. STOP NEWTON ITERATION.'
          STOP
        END IF
        DRESNORM=0.0D0
        DRHSNORM=0.0D0
        DO I=1,N
          DAXI=0.0D0
          DO J=ir(I),ir(I+1)-1
            DAXI=DAXI+a(J)*x(ic(J),1)
          END DO
          DRESNORM=DRESNORM+(DAXI-NORMAL_RHS(I))**2
          DRHSNORM=DRHSNORM+NORMAL_RHS(I)**2
        END DO
        DRESNORM=DSQRT(DRESNORM)
        DRHSNORM=DSQRT(DRHSNORM)
        IF(DRHSNORM.GT.0.0D0) THEN
          DLINEAR_REL=DRESNORM/DRHSNORM
          WRITE(*,*)'NEWTON_PARDISO_REL_RESIDUAL CSR=',DLINEAR_REL
        ELSE
          DLINEAR_REL=0.0D0
          WRITE(*,*)'NEWTON_PARDISO_ABS_RESIDUAL CSR=',DRESNORM
        END IF
        DO I=1,N
          x(I,1)=COL_SCALE(I)*x(I,1)
        END DO
        ICG_DONE=0
        ICG_MAX=0
        GO TO 127
      END IF

      IF(IDIRECT_BICGSTAB.EQ.1.OR.IDIRECT_BICGSTAB.EQ.2) THEN
        DLM_LAMBDA=DLM_LAMBDA_CURRENT
        CG_X=0.0D0
        NORMAL_DIAG=0.0D0
        DO I=1,N
          DO J=ir(I),ir(I+1)-1
            IF(ic(J).EQ.I) NORMAL_DIAG(I)=NORMAL_DIAG(I)+a(J)
          END DO
          IF(IDIRECT_BICGSTAB.EQ.2) NORMAL_DIAG(I)=NORMAL_DIAG(I)+DLM_LAMBDA
        END DO
        DO I=1,N
          IF(DABS(NORMAL_DIAG(I)).GT.1.0D-30) THEN
            NORMAL_DIAG(I)=1.0D0/NORMAL_DIAG(I)
          ELSE
            NORMAL_DIAG(I)=1.0D0
          END IF
        END DO
        DO I=1,N
          CG_R(I)=-FFS(I)
          CG_Z(I)=CG_R(I)
          CG_P(I)=0.0D0
          CG_AP(I)=0.0D0
        END DO
        DRHSNORM=0.0D0
        DO I=1,N
          DRHSNORM=DRHSNORM+FFS(I)*FFS(I)
        END DO
        DRHSNORM=DSQRT(DRHSNORM)
        DBICG_RHO_OLD=1.0D0
        DBICG_ALPHA=1.0D0
        DBICG_OMEGA=1.0D0
        IBICG_MAX=2000
        ICG_MAX=IBICG_MAX
        DCG_TOL=1.0D-6
        DRESNORM=DRHSNORM
        DO IBICG_DONE=1,IBICG_MAX
          DBICG_RHO_NEW=0.0D0
          DO I=1,N
            DBICG_RHO_NEW=DBICG_RHO_NEW+CG_Z(I)*CG_R(I)
          END DO
          IF(DABS(DBICG_RHO_NEW).LT.1.0D-300) EXIT
          IF(IBICG_DONE.EQ.1) THEN
            DO I=1,N
              CG_P(I)=CG_R(I)
            END DO
          ELSE
            DBICG_BETA=(DBICG_RHO_NEW/DBICG_RHO_OLD)*(DBICG_ALPHA/DBICG_OMEGA)
            DO I=1,N
              CG_P(I)=CG_R(I)+DBICG_BETA*(CG_P(I)-DBICG_OMEGA*CG_AP(I))
            END DO
          END IF
          DO I=1,N
            NORMAL_RHS(I)=NORMAL_DIAG(I)*CG_P(I)
          END DO
          CG_AP=0.0D0
          DO I=1,N
            DO J=ir(I),ir(I+1)-1
              CG_AP(I)=CG_AP(I)+a(J)*NORMAL_RHS(ic(J))
            END DO
            IF(IDIRECT_BICGSTAB.EQ.2) CG_AP(I)=CG_AP(I)+DLM_LAMBDA*NORMAL_RHS(I)
          END DO
          DBICG_DEN=0.0D0
          DO I=1,N
            DBICG_DEN=DBICG_DEN+CG_Z(I)*CG_AP(I)
          END DO
          IF(DABS(DBICG_DEN).LT.1.0D-300) EXIT
          DBICG_ALPHA=DBICG_RHO_NEW/DBICG_DEN
          DO I=1,N
            CG_TMP(I)=CG_R(I)-DBICG_ALPHA*CG_AP(I)
          END DO
          DRESNORM=0.0D0
          DO I=1,N
            DRESNORM=DRESNORM+CG_TMP(I)*CG_TMP(I)
          END DO
          DRESNORM=DSQRT(DRESNORM)
          IF(DRHSNORM.GT.0.0D0.AND.DRESNORM/DRHSNORM.LT.DCG_TOL) THEN
            DO I=1,N
              CG_X(I)=CG_X(I)+DBICG_ALPHA*NORMAL_RHS(I)
            END DO
            EXIT
          END IF
          DO I=1,N
            NORMAL_RHS(I)=NORMAL_DIAG(I)*CG_TMP(I)
          END DO
          CG_DIAG=0.0D0
          DO I=1,N
            DO J=ir(I),ir(I+1)-1
              CG_DIAG(I)=CG_DIAG(I)+a(J)*NORMAL_RHS(ic(J))
            END DO
            IF(IDIRECT_BICGSTAB.EQ.2) CG_DIAG(I)=CG_DIAG(I)+DLM_LAMBDA*NORMAL_RHS(I)
          END DO
          DBICG_NUM=0.0D0
          DBICG_DEN=0.0D0
          DO I=1,N
            DBICG_NUM=DBICG_NUM+CG_DIAG(I)*CG_TMP(I)
            DBICG_DEN=DBICG_DEN+CG_DIAG(I)*CG_DIAG(I)
          END DO
          IF(DABS(DBICG_DEN).LT.1.0D-300) EXIT
          DBICG_OMEGA=DBICG_NUM/DBICG_DEN
          DO I=1,N
            CG_X(I)=CG_X(I)+DBICG_ALPHA*NORMAL_DIAG(I)*CG_P(I)+DBICG_OMEGA*NORMAL_RHS(I)
            CG_R(I)=CG_TMP(I)-DBICG_OMEGA*CG_DIAG(I)
          END DO
          DRESNORM=0.0D0
          DO I=1,N
            DRESNORM=DRESNORM+CG_R(I)*CG_R(I)
          END DO
          DRESNORM=DSQRT(DRESNORM)
          IF(DRHSNORM.GT.0.0D0.AND.DRESNORM/DRHSNORM.LT.DCG_TOL) EXIT
          IF(DABS(DBICG_OMEGA).LT.1.0D-300) EXIT
          DBICG_RHO_OLD=DBICG_RHO_NEW
        END DO
        ICG_DONE=IBICG_DONE
        IF(DRHSNORM.GT.0.0D0) THEN
          DLINEAR_REL=DRESNORM/DRHSNORM
        ELSE
          DLINEAR_REL=DRESNORM
        END IF
        IF(IDIRECT_BICGSTAB.EQ.2) THEN
          WRITE(*,*)'SHIFTED_JACOBIAN_BICGSTAB ITER/MAX/RELRES/TOL/LAMBDA=',ICG_DONE,IBICG_MAX,DLINEAR_REL,DCG_TOL,DLM_LAMBDA
        ELSE
          WRITE(*,*)'DIRECT_BICGSTAB ITER/MAX/RELRES/TOL=',ICG_DONE,IBICG_MAX,DLINEAR_REL,DCG_TOL
        END IF
        DO I=1,N
          x(I,1)=COL_SCALE(I)*CG_X(I)
        END DO
        GO TO 127
      END IF

      IF(IDIRECT_BICGSTAB.EQ.3) THEN
        DLM_LAMBDA=DLM_LAMBDA_CURRENT
        NORMAL_RHS=0.0D0
        DO I=1,N
          NORMAL_RHS(I)=FFS(I)
        END DO
        NDIAG_MISSING=0
        DO I=1,N
          IDIAG=0
          DO J=ir(I),ir(I+1)-1
            IF(ic(J).EQ.I) THEN
              a(J)=a(J)+DLM_LAMBDA
              IDIAG=1
              EXIT
            END IF
          END DO
          IF(IDIAG.EQ.0) NDIAG_MISSING=NDIAG_MISSING+1
        END DO
        DO I=1,N
          FFS(I)=NORMAL_RHS(I)
          x(I,1)=0.0D0
        END DO
        WRITE(*,*)'SHIFTED_JACOBIAN_PARDISO RHS_SIGN_PLUS LAMBDA/DIAG_MISSING=',DLM_LAMBDA,NDIAG_MISSING
        call pardiso (ipt, maxfct, mnum, mtype, iphase, n, a, ir, ic, iperm, nrhs, iparm, msglvl, FFS, x, ierror)
        WRITE(*,*)'SHIFTED_JACOBIAN_PARDISO IERROR=',ierror
        IF(ierror.NE.0) THEN
          WRITE(*,*)'SHIFTED_JACOBIAN_PARDISO FAILED. STOP NEWTON ITERATION.'
          STOP
        END IF
        DRESNORM=0.0D0
        DRHSNORM=0.0D0
        DO I=1,N
          DAXI=0.0D0
          DO J=ir(I),ir(I+1)-1
            DAXI=DAXI+a(J)*x(ic(J),1)
          END DO
          DRESNORM=DRESNORM+(DAXI-NORMAL_RHS(I))**2
          DRHSNORM=DRHSNORM+NORMAL_RHS(I)**2
        END DO
        DRESNORM=DSQRT(DRESNORM)
        DRHSNORM=DSQRT(DRHSNORM)
        IF(DRHSNORM.GT.0.0D0) THEN
          DLINEAR_REL=DRESNORM/DRHSNORM
          WRITE(*,*)'SHIFTED_JACOBIAN_PARDISO_REL_RESIDUAL=',DLINEAR_REL
        ELSE
          DLINEAR_REL=0.0D0
          WRITE(*,*)'SHIFTED_JACOBIAN_PARDISO_ABS_RESIDUAL=',DRESNORM
        END IF
        DO I=1,N
          x(I,1)=COL_SCALE(I)*x(I,1)
        END DO
        ICG_DONE=0
        ICG_MAX=0
        GO TO 127
      END IF

!     -------- MARQUARDT DIAGONAL NORMAL EQUATION -------------------
!     Solve (diag(J^T J) + lambda*diag(J^T J)) dy = J^T F
!     for the scaled variable dy, then recover dx = COL_SCALE*dy.
!     This is the lightweight Marquardt step; the full J^T J path was too
!     heavy for the 5,6,6 case with the current triplet CSR builder.
      DLM_LAMBDA=DLM_LAMBDA_CURRENT
      IF(IGAUSS_NEWTON_ONLY.NE.0) THEN
        DLM_LAMBDA=0.0D0
        DLM_LAMBDA_CURRENT=0.0D0
        DLM_LAMBDA_USED=0.0D0
      END IF
      IGN_ACTIVE=0
      IF(IGAUSS_NEWTON_TRIAL.NE.0.AND.IGN_DECREASE_COUNT.GE.3.AND.DLM_LAMBDA_CURRENT.LE.DGN_LAMBDA_SWITCH.AND.WW4.LE.DGN_WW4_SWITCH) THEN
        DLM_LAMBDA=0.0D0
        IGN_ACTIVE=1
        WRITE(*,*)'GAUSS_NEWTON_TRIAL ACTIVE WW4/LAMBDA_CURRENT/COUNT=',WW4,DLM_LAMBDA_CURRENT,IGN_DECREASE_COUNT
      END IF
      IF(IMARQUARDT_CG.NE.0) THEN
        NORMAL_RHS=0.0D0
        CG_DIAG=0.0D0
        NORMAL_DIAG=0.0D0
        CALL OZ_CSR_TMV_OMP(N,ir,ic,a,FFS,NORMAL_RHS,CG_WORK_T)
        DO I=1,N
          NORMAL_RHS(I)=-NORMAL_RHS(I)
        END DO
        DO I=1,N
          DO K=1,6
            CG_WORK(K,I)=0.0D0
          END DO
        END DO
        DO I=1,N
          DO J=ir(I),ir(I+1)-1
            CALL OZ_ADD_PRODUCT_SLICE6(CG_WORK(1,ic(J)),a(J),a(J))
          END DO
        END DO
        DO I=1,N
          DO K=6,1,-1
            CG_DIAG(I)=CG_DIAG(I)+CG_WORK(K,I)
          END DO
        END DO
        DADD_MAX=0.0D0
        DADD_MIN=1.0D300
        DADD_SUM=0.0D0
        IADD_COUNT=0
        DO I=1,N
          IF(CG_DIAG(I).GT.0.0D0) THEN
            DADD_MAX=DMAX1(DADD_MAX,CG_DIAG(I))
            DADD_MIN=DMIN1(DADD_MIN,CG_DIAG(I))
            DADD_SUM=DADD_SUM+CG_DIAG(I)
            IADD_COUNT=IADD_COUNT+1
          END IF
        END DO
        IF(IADD_COUNT.GT.0) THEN
          DADD_MEAN=DADD_SUM/DBLE(IADD_COUNT)
          DADD_FLOOR=DMAX1(1.0D-12*DADD_MAX,1.0D-8*DADD_MEAN)
        ELSE
          DADD_MIN=1.0D0
          DADD_MAX=1.0D0
          DADD_MEAN=1.0D0
          DADD_FLOOR=1.0D0
        END IF
        !$OMP PARALLEL DO PRIVATE(I) IF(N.GT.1000)
        DO I=1,N
!         Marquardt damping matrix: D = I + diag(A^T A).
          NORMAL_DIAG(I)=1.0D0+CG_DIAG(I)
        END DO
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO PRIVATE(I) IF(N.GT.1000)
        DO I=1,N
!         NORMAL_RHS is -J^T*R^T*R*f for the CG equation.
!         Store the SALS b-vector J^T*R^T*R*f separately.
          SALS_RHS(I)=-NORMAL_RHS(I)
          SALS_DIAG(I)=NORMAL_DIAG(I)
        END DO
        !$OMP END PARALLEL DO
        IF(IPK.EQ.0.AND.IGAUSS_NEWTON_ONLY.EQ.0) THEN
          DLM_TRACE_INV_INIT=0.0D0
          DO I=1,N
            IF(SALS_DIAG(I).GT.1.0D-300) DLM_TRACE_INV_INIT=DLM_TRACE_INV_INIT+1.0D0/SALS_DIAG(I)
          END DO
          DLM_LAMBDA_C_INIT=1.0D0/DMAX1(DLM_TRACE_INV_INIT,1.0D-300)
          IF(DLM_LAMBDA_CURRENT.LT.DLM_LAMBDA_C_INIT) THEN
            WRITE(*,*)'MARQUARDT INITIAL LAMBDA RAISE OLD/LC=',DLM_LAMBDA_CURRENT,DLM_LAMBDA_C_INIT
            DLM_LAMBDA_CURRENT=DLM_LAMBDA_C_INIT
            DLM_LAMBDA=DLM_LAMBDA_CURRENT
          END IF
        END IF
        IF(IGAUSS_NEWTON_ONLY.NE.0) THEN
          DLM_LAMBDA=0.0D0
          DLM_LAMBDA_CURRENT=0.0D0
          DLM_LAMBDA_USED=0.0D0
          WRITE(*,*)'GAUSS_NEWTON_ONLY CG NORMAL EQUATION LAMBDA=',DLM_LAMBDA
        END IF
        WRITE(*,*)'MARQUARDT_ADD_DIAG MODE/MIN/MAX/MEAN/FLOOR=',2,DADD_MIN,DADD_MAX,DADD_MEAN,DADD_FLOOR
        WRITE(*,*)'MARQUARDT_DAMPING_DIAG FORM=I_PLUS_DIAG_ATA'
        WRITE(*,*)'MARQUARDT_ADD_DIAG_MIN=',DADD_MIN
        WRITE(*,*)'MARQUARDT_ADD_DIAG_MAX=',DADD_MAX
        WRITE(*,*)'MARQUARDT_ADD_DIAG_MEAN=',DADD_MEAN
        WRITE(*,*)'MARQUARDT_ADD_DIAG_FLOOR=',DADD_FLOOR
        CG_TMP=0.0D0
        DOUT_ROW_RNORM=0.0D0
        DOUT_ROW_JNORM=0.0D0
        DO I3=1,NSB3
          II=NSBUO(I3)
          DO ITT=1,IAT
            IROW=GLOBAL_EQUATION_INDEX(ITT,II,4,NX,IAT)
            DOUT_ROW_RNORM=DOUT_ROW_RNORM+FFS(IROW)*FFS(IROW)
            DO J=ir(IROW),ir(IROW+1)-1
              CG_TMP(ic(J))=CG_TMP(ic(J))-a(J)*FFS(IROW)
              DOUT_ROW_JNORM=DOUT_ROW_JNORM+a(J)*a(J)
            END DO
          END DO
        END DO
        DOUT_JTR_NORM=0.0D0
        DOUT_TOTAL_JTR_NORM=0.0D0
        DOUT_JTR_DOT=0.0D0
        DOUT_COL_DIAG_MEAN=0.0D0
        DOUT_COL_SCALE_MIN=1.0D300
        DOUT_COL_SCALE_MAX=0.0D0
        IOUT_COL_COUNT=0
        IOUT_MAX_ROW=0
        DOUT_DIFF_MAX=0.0D0
        DO I=1,NX*IAT
          DOUT_JTR_NORM=DOUT_JTR_NORM+CG_TMP(I)*CG_TMP(I)
          DOUT_TOTAL_JTR_NORM=DOUT_TOTAL_JTR_NORM+NORMAL_RHS(I)*NORMAL_RHS(I)
          DOUT_JTR_DOT=DOUT_JTR_DOT+CG_TMP(I)*NORMAL_RHS(I)
          IF(DABS(CG_TMP(I)).GT.DOUT_DIFF_MAX) THEN
            DOUT_DIFF_MAX=DABS(CG_TMP(I))
            IOUT_MAX_ROW=I
          END IF
          IF(CG_TMP(I).NE.0.0D0) THEN
            DOUT_COL_DIAG_MEAN=DOUT_COL_DIAG_MEAN+CG_DIAG(I)
            DOUT_COL_SCALE_MIN=DMIN1(DOUT_COL_SCALE_MIN,COL_SCALE(I))
            DOUT_COL_SCALE_MAX=DMAX1(DOUT_COL_SCALE_MAX,COL_SCALE(I))
            IOUT_COL_COUNT=IOUT_COL_COUNT+1
          END IF
        END DO
        IF(IOUT_COL_COUNT.GT.0) THEN
          DOUT_COL_DIAG_MEAN=DOUT_COL_DIAG_MEAN/DBLE(IOUT_COL_COUNT)
        ELSE
          DOUT_COL_SCALE_MIN=0.0D0
        END IF
        DOUT_JTR_NORM=DSQRT(DOUT_JTR_NORM)
        DOUT_TOTAL_JTR_NORM=DSQRT(DOUT_TOTAL_JTR_NORM)
        DOUT_ROW_RNORM=DSQRT(DOUT_ROW_RNORM)
        DOUT_ROW_JNORM=DSQRT(DOUT_ROW_JNORM)
        IF(DOUT_JTR_NORM.GT.0.0D0.AND.DOUT_TOTAL_JTR_NORM.GT.0.0D0) THEN
          DOUT_COS=DOUT_JTR_DOT/(DOUT_JTR_NORM*DOUT_TOTAL_JTR_NORM)
        ELSE
          DOUT_COS=0.0D0
        END IF
        WRITE(*,*)'OUTFLOW_U_NORMAL_DIAG ROW_RNORM/ROW_JNORM/JTR_NORM/TOTAL_U_JTR/COS=',DOUT_ROW_RNORM,DOUT_ROW_JNORM,DOUT_JTR_NORM,DOUT_TOTAL_JTR_NORM,DOUT_COS
        WRITE(*,*)'OUTFLOW_U_NORMAL_DIAG COL_COUNT/DIAG_MEAN/COL_SCALE_MINMAX=',IOUT_COL_COUNT,DOUT_COL_DIAG_MEAN,DOUT_COL_SCALE_MIN,DOUT_COL_SCALE_MAX
        IF(IOUT_MAX_ROW.NE.0) WRITE(*,*)'OUTFLOW_U_NORMAL_MAX COL/NODE/TIME/OUTJTR/TOTALJTR/CGDIAG/NORMALDIAG=',IOUT_MAX_ROW,(IOUT_MAX_ROW-1)/IAT+1,IOUT_MAX_ROW-((IOUT_MAX_ROW-1)/IAT)*IAT,CG_TMP(IOUT_MAX_ROW),NORMAL_RHS(IOUT_MAX_ROW),CG_DIAG(IOUT_MAX_ROW),NORMAL_DIAG(IOUT_MAX_ROW)
        CG_X=0.0D0
        CG_R=NORMAL_RHS
        CG_Z=0.0D0
        CG_P=0.0D0
        !$OMP PARALLEL DO PRIVATE(I) IF(N.GT.1000)
        DO I=1,N
          IF(IMARQUARDT_PCG.NE.0) THEN
            IF(CG_DIAG(I)+DLM_LAMBDA*NORMAL_DIAG(I).GT.0.0D0) THEN
              CG_Z(I)=CG_R(I)/(CG_DIAG(I)+DLM_LAMBDA*NORMAL_DIAG(I))
            ELSE
              CG_Z(I)=0.0D0
            END IF
          ELSE
            CG_Z(I)=CG_R(I)
          END IF
          CG_P(I)=CG_Z(I)
        END DO
        !$OMP END PARALLEL DO
        CG_RHSN=OZ_DDOT(N,NORMAL_RHS,1,NORMAL_RHS,1)
        CG_RHO=OZ_DDOT(N,CG_R,1,CG_Z,1)
        CG_RHSN=DSQRT(DMAX1(CG_RHSN,0.0D0))
        IF(CG_RHSN.LE.0.0D0) CG_RHSN=1.0D0
        ICG_MAX=3000
        DCG_TOL=1.0D-6
        CG_RHO_NEW=CG_RHO
        ICG_DONE=0
        DO ICG=1,ICG_MAX
          ICG_DONE=ICG
          CALL OZ_CSR_MV(N,ir,ic,a,CG_P,CG_TMP)
          CALL OZ_CSR_TMV_OMP(N,ir,ic,a,CG_TMP,CG_AP,CG_WORK_T)
          !$OMP PARALLEL DO PRIVATE(I,K) IF(N.GT.1000)
          DO I=1,N
            DO K=1,6
              CG_WORK(K,I)=0.0D0
            END DO
            CALL OZ_ADD_SLICE6(CG_WORK(1,I),CG_AP(I))
            CALL OZ_ADD_PRODUCT_SLICE6(CG_WORK(1,I),DLM_LAMBDA*NORMAL_DIAG(I),CG_P(I))
            CG_AP(I)=0.0D0
            DO K=6,1,-1
              CG_AP(I)=CG_AP(I)+CG_WORK(K,I)
            END DO
          END DO
          !$OMP END PARALLEL DO
          CG_PAP=OZ_DDOT(N,CG_P,1,CG_AP,1)
          IF(CG_PAP.LE.0.0D0) EXIT
          CG_ALPHA=CG_RHO/CG_PAP
          !$OMP PARALLEL DO PRIVATE(I) IF(N.GT.1000)
          DO I=1,N
            CG_X(I)=CG_X(I)+CG_ALPHA*CG_P(I)
            CG_R(I)=CG_R(I)-CG_ALPHA*CG_AP(I)
          END DO
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(I) IF(N.GT.1000)
          DO I=1,N
            IF(IMARQUARDT_PCG.NE.0) THEN
              IF(CG_DIAG(I)+DLM_LAMBDA*NORMAL_DIAG(I).GT.0.0D0) THEN
                CG_Z(I)=CG_R(I)/(CG_DIAG(I)+DLM_LAMBDA*NORMAL_DIAG(I))
              ELSE
                CG_Z(I)=0.0D0
              END IF
            ELSE
              CG_Z(I)=CG_R(I)
            END IF
          END DO
          !$OMP END PARALLEL DO
          CG_RHO_NEW=OZ_DDOT(N,CG_R,1,CG_Z,1)
          DRESNORM=OZ_DDOT(N,CG_R,1,CG_R,1)
          IF(DSQRT(DMAX1(DRESNORM,0.0D0))/CG_RHSN.LT.DCG_TOL) EXIT
          CG_BETA=CG_RHO_NEW/CG_RHO
          !$OMP PARALLEL DO PRIVATE(I) IF(N.GT.1000)
          DO I=1,N
            CG_P(I)=CG_Z(I)+CG_BETA*CG_P(I)
          END DO
          !$OMP END PARALLEL DO
          CG_RHO=CG_RHO_NEW
        END DO
        !$OMP PARALLEL DO PRIVATE(I) IF(N.GT.1000)
        DO I=1,N
          x(I,1)=COL_SCALE(I)*CG_X(I)
        END DO
        !$OMP END PARALLEL DO
        DRESNORM=OZ_DDOT(N,CG_R,1,CG_R,1)
        DLINEAR_REL=DSQRT(DMAX1(DRESNORM,0.0D0))/CG_RHSN
        DLM_LAMBDA_USED=DLM_LAMBDA
        WRITE(*,*)'MARQUARDT_CG ITER/MAX/RELRES/TOL/LAMBDA/PCG=',ICG_DONE,ICG_MAX,DLINEAR_REL,DCG_TOL,DLM_LAMBDA,IMARQUARDT_PCG
        IF(IGN_ACTIVE.NE.0) WRITE(*,*)'GAUSS_NEWTON_TRIAL RESULT ITER/RELRES=',ICG_DONE,DLINEAR_REL
        GOTO 127
      END IF

      NORMAL_RHS=0.0D0
      NORMAL_DIAG=0.0D0
      CALL OZ_CSR_TMV_OMP(N,ir,ic,a,FFS,NORMAL_RHS,CG_WORK_T)
      DO I=1,N
        NORMAL_RHS(I)=-NORMAL_RHS(I)
        DO K=1,6
          CG_WORK(K,I)=0.0D0
        END DO
      END DO
      DO I=1,N
        DO J=ir(I),ir(I+1)-1
          CALL OZ_ADD_PRODUCT_SLICE6(CG_WORK(1,ic(J)),a(J),a(J))
        END DO
      END DO
      DO I=1,N
        DO K=6,1,-1
          NORMAL_DIAG(I)=NORMAL_DIAG(I)+CG_WORK(K,I)
        END DO
      END DO
      call triplet_init(SPARSE_NORMAL_A, n, n, n)
      DO I=1,N
        IF(NORMAL_DIAG(I).GT.0.0D0) THEN
          call triplet_add(SPARSE_NORMAL_A, I, I, (1.0D0+DLM_LAMBDA)*NORMAL_DIAG(I))
        ELSE
          call triplet_add(SPARSE_NORMAL_A, I, I, DLM_LAMBDA)
        END IF
      END DO
      IF(ALLOCATED(ir)) DEALLOCATE(ir)
      IF(ALLOCATED(ic)) DEALLOCATE(ic)
      IF(ALLOCATED(a)) DEALLOCATE(a)
      call triplet_to_csr(SPARSE_NORMAL_A, ir, ic, a, mn)
      mnr=n+1
      DO I=1,N
        FFS(I)=NORMAL_RHS(I)
      END DO
      WRITE(*,*)'MARQUARDT DIAG NORMAL CSR NNZ/LAMBDA=',mn,DLM_LAMBDA

      IF(ALLOCATED(DENSE_A)) THEN
        IF(ALLOCATED(ISINGULAR_PIN)) ISINGULAR_PIN=0
        DENSE_A=SSMMK
        DO I=1,N
          DENSE_B(I,1)=FF1(I)
        END DO
        CALL DGESV(N,NRHS,DENSE_A,N,IPIV,DENSE_B,N,INFO)
        WRITE(*,*)'DGESV INFO=',INFO
        IF(INFO.GT.0) THEN
          ISINGULAR_PIN=0
          INFO2=INFO
          DO IDIAG=1,8
            IBLOCK=(INFO2-1)/(NX*IAT)+1
            ILOCAL=INFO2-(IBLOCK-1)*NX*IAT
            IPNODE=(ILOCAL-1)/IAT+1
            IPTIME=ILOCAL-(IPNODE-1)*IAT
            DROWSUM=0.0D0
            DCOLSUM=0.0D0
            DO J=1,N
              DROWSUM=DROWSUM+DABS(SSMMK(INFO2,J))
              DCOLSUM=DCOLSUM+DABS(SSMMK(J,INFO2))
            END DO
            WRITE(*,*)'DGESV SINGULAR DOF BLOCK/NODE/TIME=',IBLOCK,IPNODE,IPTIME
            WRITE(*,*)'DGESV SINGULAR DIAG/ROWSUM/COLSUM/RHS=',SSMMK(INFO2,INFO2),DROWSUM,DCOLSUM,FF1(INFO2)
            ISINGULAR_PIN(INFO2)=1
            DENSE_A=SSMMK
            DO I=1,N
              DENSE_B(I,1)=FF1(I)
            END DO
            DO I=1,N
              IF(ISINGULAR_PIN(I).NE.0) THEN
                DO J=1,N
                  DENSE_A(I,J)=0.0D0
                  DENSE_A(J,I)=0.0D0
                END DO
                DENSE_A(I,I)=1.0D0
                DENSE_B(I,1)=0.0D0
              END IF
            END DO
            CALL DGESV(N,NRHS,DENSE_A,N,IPIV,DENSE_B,N,INFO2)
            WRITE(*,*)'DGESV INFO AFTER PIN COUNT=',IDIAG,INFO2
            IF(INFO2.EQ.0) EXIT
          END DO
        END IF
        IF(INFO.EQ.0) THEN
          DRESNORM_DENSE=0.0D0
          DRHSNORM=0.0D0
          DO I=1,N
            DAXI=0.0D0
            DO J=1,N
              DAXI=DAXI+SSMMK(I,J)*DENSE_B(J,1)
            END DO
            DRESNORM_DENSE=DRESNORM_DENSE+(DAXI-FF1(I))**2
            DRHSNORM=DRHSNORM+FF1(I)**2
          END DO
          DRESNORM_DENSE=DSQRT(DRESNORM_DENSE)
          DRHSNORM=DSQRT(DRHSNORM)
          IF(DRHSNORM.GT.0.0D0) THEN
            WRITE(*,*)'DGESV_REL_RESIDUAL=',DRESNORM_DENSE/DRHSNORM
          ELSE
            WRITE(*,*)'DGESV_ABS_RESIDUAL=',DRESNORM_DENSE
          END IF
        END IF
      END IF

      IF(ALLOCATED(ISINGULAR_PIN).AND.IAPPLY_DGESV_PIN.NE.0) THEN
        IPINNED_ROW=0
        DO I=1,N
          IF(ISINGULAR_PIN(I).NE.0) THEN
            DO J=1,N
              SSMMK(I,J)=0.0D0
              SSMMK(J,I)=0.0D0
            END DO
            SSMMK(I,I)=1.0D0
            FFS(I)=0.0D0
            IPINNED_ROW=IPINNED_ROW+1
          END IF
        END DO
        IF(IPINNED_ROW.GT.0) THEN
          WRITE(*,*)'APPLIED DGESV SINGULAR PIN TO CSR MATRIX=',IPINNED_ROW
          IF(ALLOCATED(ir)) DEALLOCATE(ir)
          IF(ALLOCATED(ic)) DEALLOCATE(ic)
          IF(ALLOCATED(a)) DEALLOCATE(a)
          call triplet_init(SPARSE_A, n, n, max(1, min(1000000, 20*n)))
          DO I=1,N
            DO J=1,N
              IF(SSMMK(I,J).NE.0.0D0) THEN
                call triplet_add(SPARSE_A, I, J, SSMMK(I,J))
              END IF
            END DO
          END DO
          call triplet_to_csr(SPARSE_A, ir, ic, a, mn)
          mnr=n+1
          write(*,*)'rebuilt pinned csr !mnr,mn',mnr,mn
        END IF
      END IF

      call pardiso (ipt, maxfct, mnum, mtype, iphase, n, a, ir, ic, iperm, nrhs, iparm, msglvl, FFS, x, ierror)

      write(*,*)'ierror=',ierror
      WRITE(*,*)'PARDISO IPARM 14 18 19 20 22=',iparm(14),iparm(18),iparm(19),iparm(20),iparm(22)
      IF(ierror.NE.0) THEN
        WRITE(*,*)'PARDISO FAILED. STOP NEWTON ITERATION.'
        STOP
      END IF

      DO I=1,N
        x(I,1)=COL_SCALE(I)*x(I,1)
      END DO

      DRESNORM=0.0D0
      DRHSNORM=0.0D0
      DO I=1,N
        DAXI=0.0D0
        DO J=ir(I),ir(I+1)-1
          DAXI=DAXI+a(J)*(x(ic(J),1)/COL_SCALE(ic(J)))
        END DO
        DRESNORM=DRESNORM+(DAXI-NORMAL_RHS(I))**2
        DRHSNORM=DRHSNORM+NORMAL_RHS(I)**2
      END DO
      DRESNORM=DSQRT(DRESNORM)
      DRHSNORM=DSQRT(DRHSNORM)
      IF(DRHSNORM.GT.0.0D0) THEN
        DLINEAR_REL=DRESNORM/DRHSNORM
        WRITE(*,*)'MARQUARDT_LINEAR_REL_RESIDUAL CSR=',DLINEAR_REL
      ELSE
        DLINEAR_REL=0.0D0
        WRITE(*,*)'MARQUARDT_LINEAR_ABS_RESIDUAL CSR=',DRESNORM
      END IF
127   CONTINUE
      DDUVWP_SD=0.0D0
      DO I=1,N
        DFSCALED=FF1(I)*ROW_SCALE(I)
        DO J=ir(I),ir(I+1)-1
          DDUVWP_SD(ic(J))=DDUVWP_SD(ic(J))-COL_SCALE(ic(J))*a(J)*DFSCALED
        END DO
      END DO
      IF(ALLOCATED(DENSE_A).AND.INFO.EQ.0) THEN
        DO I=1,N
          X(I,1)=DENSE_B(I,1)
        END DO
        WRITE(*,*)'USING_DGESV_SOLUTION_FOR_NEWTON_UPDATE'
      END IF

!     -------- SOLUTION OF LINEAR SYSTEM ------------
!      trans='N'
!      n=nx*4
!      nrhs=1
!      lda=200000
!      ldb=200000
!        write(*,*)' Clear ohkubo 110!'   
!      CALL dgetrf(n,n,ssmmk,lda,ipiv,info)
!
!      CALL dgetrs(trans,n,nrhs,ssmmk,lda,ipiv,ffs,ldb,info)
!      write(*,*) '----------info=',info
!
!
!     The linear system row order is the equation order, but the solution
!     vector X and the Newton correction DDUVWP are ordered as U,V,W,P.
     DO I=1,4*IAT*NX
        DDUVWP(I)=X(I,1)
      END DO
      DS_OLD=0.0D0
      DS_PRED=0.0D0
      DROWSUM=0.0D0
      DAXNORM=0.0D0
      DRESNORM=0.0D0
      DO I=1,N
        DFSCALED=FF1(I)*ROW_SCALE(I)
        DJD_NEWTON=0.0D0
        DJD_STEEPEST=0.0D0
        DO J=ir(I),ir(I+1)-1
          DJD_NEWTON=DJD_NEWTON+a(J)*(DDUVWP(ic(J))/COL_SCALE(ic(J)))
          DJD_STEEPEST=DJD_STEEPEST+a(J)*(DDUVWP_SD(ic(J))/COL_SCALE(ic(J)))
        END DO
        DS_OLD=DS_OLD+DFSCALED*DFSCALED
        DS_PRED=DS_PRED+(DFSCALED+DJD_NEWTON)*(DFSCALED+DJD_NEWTON)
        DROWSUM=DROWSUM+(DFSCALED-DJD_NEWTON)*(DFSCALED-DJD_NEWTON)
        DAXNORM=DAXNORM+(DFSCALED+DJD_STEEPEST)*(DFSCALED+DJD_STEEPEST)
        DRESNORM=DRESNORM+(DFSCALED-DJD_STEEPEST)*(DFSCALED-DJD_STEEPEST)
      END DO
      WRITE(*,*)'MARQUARDT_VS_SD_PRED OLD/MARQ_PLUS/MARQ_MINUS/SD_PLUS/SD_MINUS=',DSQRT(DS_OLD),DSQRT(DS_PRED),DSQRT(DROWSUM),DSQRT(DAXNORM),DSQRT(DRESNORM)
      DS_OLD=0.0D0
      DS_PRED=0.0D0
      DROWSUM=0.0D0
      DAXNORM=0.0D0
      DRESNORM=0.0D0
      DO I=1,N
        DFSCALED=FFS_PHYS(I)
        DJD_NEWTON=0.0D0
        DJD_STEEPEST=0.0D0
        DO J=ir(I),ir(I+1)-1
          DJD_NEWTON=DJD_NEWTON+a(J)*(DDUVWP(ic(J))/COL_SCALE(ic(J)))
          DJD_STEEPEST=DJD_STEEPEST+a(J)*(DDUVWP_SD(ic(J))/COL_SCALE(ic(J)))
        END DO
        IF(ROW_SCALE(I).GT.0.0D0) THEN
          DJD_NEWTON=DJD_NEWTON/ROW_SCALE(I)
          DJD_STEEPEST=DJD_STEEPEST/ROW_SCALE(I)
        ELSE
          DJD_NEWTON=0.0D0
          DJD_STEEPEST=0.0D0
        END IF
        DS_OLD=DS_OLD+DFSCALED*DFSCALED
        DS_PRED=DS_PRED+(DFSCALED+DJD_NEWTON)*(DFSCALED+DJD_NEWTON)
        DROWSUM=DROWSUM+(DFSCALED-DJD_NEWTON)*(DFSCALED-DJD_NEWTON)
        DAXNORM=DAXNORM+(DFSCALED+DJD_STEEPEST)*(DFSCALED+DJD_STEEPEST)
        DRESNORM=DRESNORM+(DFSCALED-DJD_STEEPEST)*(DFSCALED-DJD_STEEPEST)
      END DO
      WRITE(*,*)'MARQUARDT_VS_SD_PHYS_PRED OLD/MARQ_PLUS/MARQ_MINUS/SD_PLUS/SD_MINUS=',DSQRT(DS_OLD),DSQRT(DS_PRED),DSQRT(DROWSUM),DSQRT(DAXNORM),DSQRT(DRESNORM)
      IF(DRESNORM.LT.DAXNORM) THEN
        DSD_PHYS_SIGN=-UPDATE_SIGN
      ELSE
        DSD_PHYS_SIGN=UPDATE_SIGN
      END IF
      WRITE(*,*)'MARQUARDT_SD_PHYS_SIGN_SELECT SIGN/SD_PLUS/SD_MINUS=',DSD_PHYS_SIGN,DSQRT(DAXNORM),DSQRT(DRESNORM)
      IF(.FALSE..AND.DMIN1(DAXNORM,DRESNORM).LT.DMIN1(DS_PRED,DROWSUM)) THEN
        DO I=1,N
          DDUVWP(I)=DDUVWP_SD(I)
        END DO
        IUSE_STEEPEST=1
        WRITE(*,*)'MARQUARDT SWITCH_TO_STEEPEST_BY_PRED'
      END IF
      DDIR_NEWTON=0.0D0
      DDIR_STEEPEST=0.0D0
      DSD_JDNORM2=0.0D0
      DO I=1,N
        DFSCALED=FF1(I)*ROW_SCALE(I)
        DJD_NEWTON=0.0D0
        DJD_STEEPEST=0.0D0
        DO J=ir(I),ir(I+1)-1
          DJD_NEWTON=DJD_NEWTON+a(J)*(DDUVWP(ic(J))/COL_SCALE(ic(J)))
          DJD_STEEPEST=DJD_STEEPEST+a(J)*(DDUVWP_SD(ic(J))/COL_SCALE(ic(J)))
        END DO
        DDIR_NEWTON=DDIR_NEWTON+DFSCALED*DJD_NEWTON
        DDIR_STEEPEST=DDIR_STEEPEST+DFSCALED*DJD_STEEPEST
        DSD_JDNORM2=DSD_JDNORM2+DJD_STEEPEST*DJD_STEEPEST
      END DO
      IF(DSD_JDNORM2.GT.0.0D0) THEN
        DSD_ALPHA=-DDIR_STEEPEST/DSD_JDNORM2
      ELSE
        DSD_ALPHA=0.0D0
      END IF
      WRITE(*,*)'MARQUARDT_DIRECTION_CHECK STEP/GRAD/LINEAR_REL=',DDIR_NEWTON,DDIR_STEEPEST,DLINEAR_REL
      WRITE(*,*)'MARQUARDT_DIR_DPHI SELECTED_PLUS/MINUS=',2.0D0*DDIR_NEWTON,-2.0D0*DDIR_NEWTON
      WRITE(*,*)'MARQUARDT_SD_DPHI PLUS/MINUS=',2.0D0*DDIR_STEEPEST,-2.0D0*DDIR_STEEPEST
      WRITE(*,*)'MARQUARDT_PURE_SD_ALPHA DIR/JDNORM2/ALPHA=',DDIR_STEEPEST,DSD_JDNORM2,DSD_ALPHA
      IF(DDIR_NEWTON.LT.0.0D0) THEN
        WRITE(*,*)'MARQUARDT_DIR_DESCENT_SIGN SELECTED=PLUS'
      ELSE IF(DDIR_NEWTON.GT.0.0D0) THEN
        WRITE(*,*)'MARQUARDT_DIR_DESCENT_SIGN SELECTED=MINUS'
      ELSE
        WRITE(*,*)'MARQUARDT_DIR_DESCENT_SIGN SELECTED=ZERO'
      END IF
      IF(DDIR_STEEPEST.LT.0.0D0) THEN
        WRITE(*,*)'MARQUARDT_SD_DESCENT_SIGN=PLUS'
      ELSE IF(DDIR_STEEPEST.GT.0.0D0) THEN
        WRITE(*,*)'MARQUARDT_SD_DESCENT_SIGN=MINUS'
      ELSE
        WRITE(*,*)'MARQUARDT_SD_DESCENT_SIGN=ZERO'
      END IF
      DO IBLOCK=1,4
        DRESNORM=0.0D0
        DAXNORM=0.0D0
        DROWSUM=0.0D0
        DROWMAX=0.0D0
        DO I=1,NX*IAT
          IROW=I+(IBLOCK-1)*NX*IAT
          DAXI=0.0D0
          DO J=ir(IROW),ir(IROW+1)-1
            DAXI=DAXI+a(J)*(DDUVWP(ic(J))/COL_SCALE(ic(J)))
          END DO
          DFSCALED=FF1(IROW)*ROW_SCALE(IROW)
          DRESNORM=DRESNORM+DFSCALED*DFSCALED
          DAXNORM=DAXNORM+DAXI*DAXI
          DROWSUM=DROWSUM+(DFSCALED+DAXI)*(DFSCALED+DAXI)
          DROWMAX=DMAX1(DROWMAX,DABS(DFSCALED+DAXI))
        END DO
        WRITE(*,*)'PRED_LINEAR_BLOCK IBLOCK RNORM JDNORM POSTNORM POSTMAX=',IBLOCK,DSQRT(DRESNORM),DSQRT(DAXNORM),DSQRT(DROWSUM),DROWMAX
      END DO
      DO IBLOCK=2,4
        DRESNORM=0.0D0
        DAXNORM=0.0D0
        DROWSUM=0.0D0
        DROWMAX=0.0D0
        IOUT_MAX_NODE=0
        DO I3=1,NSB3
          II=NSBUO(I3)
          DO ITT=1,IAT
            IROW=GLOBAL_EQUATION_INDEX(ITT,II,IBLOCK,NX,IAT)
            DAXI=0.0D0
            DO J=ir(IROW),ir(IROW+1)-1
              DAXI=DAXI+a(J)*(DDUVWP(ic(J))/COL_SCALE(ic(J)))
            END DO
            DFSCALED=FF1(IROW)*ROW_SCALE(IROW)
            DRESNORM=DRESNORM+DFSCALED*DFSCALED
            DAXNORM=DAXNORM+DAXI*DAXI
            DROWSUM=DROWSUM+(DFSCALED+DAXI)*(DFSCALED+DAXI)
            IF(DABS(DFSCALED+DAXI).GT.DROWMAX) THEN
              DROWMAX=DABS(DFSCALED+DAXI)
              IOUT_MAX_NODE=II
              IOUT_MAX_ROW=IROW
              DOUT_MAX_FFS=DFSCALED
              DOUT_MAX_DOT=DAXI
            END IF
          END DO
        END DO
        WRITE(*,*)'PRED_OUTFLOW_BLOCK IBLOCK RNORM JDNORM POSTNORM POSTMAX=',IBLOCK,DSQRT(DRESNORM),DSQRT(DAXNORM),DSQRT(DROWSUM),DROWMAX
        WRITE(*,*)'PRED_OUTFLOW_MAX IBLOCK/NODE/ROW/R/JD/POST=',IBLOCK,IOUT_MAX_NODE,IOUT_MAX_ROW,DOUT_MAX_FFS,DOUT_MAX_DOT,DOUT_MAX_FFS+DOUT_MAX_DOT
      END DO
      IUSE_STEEPEST=0
      IF(.FALSE..AND.IMARQUARDT_SD_FALLBACK.NE.0.AND.DLM_LAMBDA_CURRENT.GE.1.0D2.AND.WW4.GT.1.0D-8.AND.DSD_ALPHA.GT.0.0D0) THEN
        DO I=1,N
          DDUVWP(I)=DSD_ALPHA*DDUVWP_SD(I)
        END DO
        IUSE_STEEPEST=1
        RELAX_RETRY_U=1.0D0
        RELAX_RETRY_VW=1.0D0
        RELAX_RETRY_P=1.0D0
        WRITE(*,*)'MARQUARDT USING PURE STEEPEST DESCENT LAMBDA/ALPHA=',DLM_LAMBDA_CURRENT,DSD_ALPHA
      END IF
      DO I=1,N
        IF(IS_CONSTRAINED(I).NE.0) DDUVWP(I)=0.0D0
      END DO
      IF(IPRESSURE_ONLY_CORRECTION.NE.0) THEN
        DO I=1,3*NX*IAT
          DDUVWP(I)=0.0D0
        END DO
        WRITE(*,*)'PRESSURE_ONLY_CORRECTION ENABLED'
      END IF
      IF(IPROJECT_CONTINUITY.NE.0) THEN
        DO I=1,N
          CG_TMP(I)=DDUVWP(I)/COL_SCALE(I)
        END DO
        DPROJ_BEFORE=0.0D0
        DPROJ_AFTER=0.0D0
        DPROJ_MAX_BEFORE=0.0D0
        DPROJ_MAX_AFTER=0.0D0
        DO I=1,NX*IAT
          DAXI=FFS(I)
          DO J=ir(I),ir(I+1)-1
            IF(ic(J).LE.3*NX*IAT) DAXI=DAXI+a(J)*CG_TMP(ic(J))
          END DO
          DPROJ_BEFORE=DPROJ_BEFORE+DAXI*DAXI
          DPROJ_MAX_BEFORE=DMAX1(DPROJ_MAX_BEFORE,DABS(DAXI))
        END DO
        DPROJ_OMEGA=1.0D-3
        DPROJ_CAP=1.0D-8
        DO IPROJ_SWEEP=1,1
          DO I=1,NX*IAT
            DAXI=FFS(I)
            DROWSUM=0.0D0
            DO J=ir(I),ir(I+1)-1
              IF(ic(J).LE.3*NX*IAT.AND.IS_CONSTRAINED(ic(J)).EQ.0) THEN
                DAXI=DAXI+a(J)*CG_TMP(ic(J))
                DROWSUM=DROWSUM+a(J)*a(J)
              END IF
            END DO
            IF(DROWSUM.GT.1.0D-8) THEN
              DPROJ_ALPHA=-DPROJ_OMEGA*DAXI/DROWSUM
              DPROJ_STEP_NORM=DABS(DPROJ_ALPHA)*DSQRT(DROWSUM)
              IF(DPROJ_STEP_NORM.GT.DPROJ_CAP) DPROJ_ALPHA=DPROJ_ALPHA*(DPROJ_CAP/DPROJ_STEP_NORM)
              DO J=ir(I),ir(I+1)-1
                IF(ic(J).LE.3*NX*IAT.AND.IS_CONSTRAINED(ic(J)).EQ.0) THEN
                  CG_TMP(ic(J))=CG_TMP(ic(J))+DPROJ_ALPHA*a(J)
                END IF
              END DO
            END IF
          END DO
        END DO
        DO I=1,N
          DDUVWP(I)=COL_SCALE(I)*CG_TMP(I)
          IF(IS_CONSTRAINED(I).NE.0) DDUVWP(I)=0.0D0
        END DO
        DO I=1,NX*IAT
          DAXI=FFS(I)
          DO J=ir(I),ir(I+1)-1
            IF(ic(J).LE.3*NX*IAT) DAXI=DAXI+a(J)*DDUVWP(ic(J))
          END DO
          DPROJ_AFTER=DPROJ_AFTER+DAXI*DAXI
          DPROJ_MAX_AFTER=DMAX1(DPROJ_MAX_AFTER,DABS(DAXI))
        END DO
        WRITE(*,*)'CONTINUITY_PROJECTION NORM BEFORE/AFTER MAX BEFORE/AFTER=',DSQRT(DPROJ_BEFORE),DSQRT(DPROJ_AFTER),DPROJ_MAX_BEFORE,DPROJ_MAX_AFTER
      END IF

!     Save the raw linear model f + J*delta used for WW4-based
!     Marquardt ratio control before the next trial rebuilds CSR.
      DO I=1,N
        CG_AP(I)=0.0D0
        CG_R(I)=0.0D0
        CG_Z(I)=0.0D0
        CG_P(I)=FF1(I)
        ROW_SCALE_BASE(I)=ROW_SCALE(I)
      END DO
      CG_P=0.0D0
      DO I=1,NX*IAT
        CG_P(I)=DDUVWP(I)/COL_SCALE(I)
      END DO
      CALL OZ_CSR_MV(N,ir,ic,a,CG_P,CG_TMP)
      DO I=1,N
        IF(ROW_SCALE(I).GT.0.0D0) CG_AP(I)=CG_TMP(I)/ROW_SCALE(I)
      END DO
      CG_P=0.0D0
      DO I=NX*IAT+1,3*NX*IAT
        CG_P(I)=DDUVWP(I)/COL_SCALE(I)
      END DO
      CALL OZ_CSR_MV(N,ir,ic,a,CG_P,CG_TMP)
      DO I=1,N
        IF(ROW_SCALE(I).GT.0.0D0) CG_R(I)=CG_TMP(I)/ROW_SCALE(I)
      END DO
      CG_P=0.0D0
      DO I=3*NX*IAT+1,N
        CG_P(I)=DDUVWP(I)/COL_SCALE(I)
      END DO
      CALL OZ_CSR_MV(N,ir,ic,a,CG_P,CG_TMP)
      DO I=1,N
        IF(ROW_SCALE(I).GT.0.0D0) CG_Z(I)=CG_TMP(I)/ROW_SCALE(I)
      END DO
      
      DO IK=1,4
      DO I=1,NX
      DO ITT=1,IAT
        ISP=ITT+(I-1)*IAT+(IK-1)*NX*IAT
        IF(ISP.GE.1.AND.ISP.LE.NX*IAT) THEN
          YUVWP(ISP)=U(ITT,I)
        END IF
        IF(ISP.GT.NX*IAT.AND.ISP.LE.2*NX*IAT) THEN
          YUVWP(ISP)=V(ITT,I)
        END IF
        IF(ISP.GT.2*NX*IAT.AND.ISP.LE.3*NX*IAT) THEN
          YUVWP(ISP)=W(ITT,I)
        END IF
        IF(ISP.GT.3*NX*IAT.AND.ISP.LE.4*NX*IAT) THEN
          YUVWP(ISP)=P(ITT,I)
        END IF
      END DO
      END DO
      END DO
      
!     Start from full relaxation, then apply a generic trust-region based
!     on the current variable scales.  This avoids problem-specific
!     constraints such as Umin while preventing an oversized Newton step.
      RELAX_U=1.0D0
      RELAX_VW=1.0D0
      RELAX_P=1.0D0
      DUPDATE_CAP_U=1.0D300
      DUPDATE_CAP_VW=1.0D300
      DUPDATE_CAP_P=1.0D300
      DRAW_MAX_U=0.0D0
      DRAW_MAX_V=0.0D0
      DRAW_MAX_W=0.0D0
      DRAW_MAX_P=0.0D0
      DSTATE_MAX_U=0.0D0
      DSTATE_MAX_V=0.0D0
      DSTATE_MAX_W=0.0D0
      DSTATE_MAX_P=0.0D0
      DO I=1,4*IAT*NX
        IF(I.LE.IAT*NX) THEN
          DRAW_MAX_U=DMAX1(DRAW_MAX_U,DABS(DDUVWP(I)))
          DSTATE_MAX_U=DMAX1(DSTATE_MAX_U,DABS(YUVWP(I)))
        ELSE IF(I.LE.2*IAT*NX) THEN
          DRAW_MAX_V=DMAX1(DRAW_MAX_V,DABS(DDUVWP(I)))
          DSTATE_MAX_V=DMAX1(DSTATE_MAX_V,DABS(YUVWP(I)))
        ELSE IF(I.LE.3*IAT*NX) THEN
          DRAW_MAX_W=DMAX1(DRAW_MAX_W,DABS(DDUVWP(I)))
          DSTATE_MAX_W=DMAX1(DSTATE_MAX_W,DABS(YUVWP(I)))
        ELSE
          DRAW_MAX_P=DMAX1(DRAW_MAX_P,DABS(DDUVWP(I)))
          DSTATE_MAX_P=DMAX1(DSTATE_MAX_P,DABS(YUVWP(I)))
        END IF
      END DO
      ITRUST_UPDATE_CAP=0
      IF(ITRUST_UPDATE_CAP.NE.0) THEN
        DTRUST_FRAC=1.0D-1
        DUPDATE_CAP_U=DMAX1(5.0D-3,DTRUST_FRAC*DMAX1(DSTATE_MAX_U,DABS(U0),1.0D-6))
        DUPDATE_CAP_VW=DMAX1(5.0D-3,DTRUST_FRAC*DMAX1(DSTATE_MAX_V,DSTATE_MAX_W,DABS(U0),1.0D-6))
        DUPDATE_CAP_P=DMAX1(1.0D-1,DTRUST_FRAC*DMAX1(DSTATE_MAX_P,1.0D0))
      END IF
      IF(DRAW_MAX_U*RELAX_U.GT.DUPDATE_CAP_U.AND.DRAW_MAX_U.GT.0.0D0) RELAX_U=DUPDATE_CAP_U/DRAW_MAX_U
      IF(DMAX1(DRAW_MAX_V,DRAW_MAX_W)*RELAX_VW.GT.DUPDATE_CAP_VW.AND.DMAX1(DRAW_MAX_V,DRAW_MAX_W).GT.0.0D0) THEN
        RELAX_VW=DUPDATE_CAP_VW/DMAX1(DRAW_MAX_V,DRAW_MAX_W)
      END IF
      IF(DRAW_MAX_P*RELAX_P.GT.DUPDATE_CAP_P.AND.DRAW_MAX_P.GT.0.0D0) RELAX_P=DUPDATE_CAP_P/DRAW_MAX_P
!     ERR1 is diagnostic only.  Accept/reject is controlled by WW4, so a
!     larger correction is allowed when it lowers the nonlinear residual.
      DERR1_CAP=1.0D300
      DERR1_TRIAL=0.0D0
      DO I=1,4*IAT*NX
        IF(I.LE.IAT*NX) THEN
          DERR1_TRIAL=DERR1_TRIAL+(RELAX_U*DDUVWP(I))**2
        ELSE IF(I.LE.3*IAT*NX) THEN
          DERR1_TRIAL=DERR1_TRIAL+(RELAX_VW*DDUVWP(I))**2
        ELSE
          DERR1_TRIAL=DERR1_TRIAL+(RELAX_P*DDUVWP(I))**2
        END IF
      END DO
      DERR1_TRIAL=DSQRT(DERR1_TRIAL)
      RELAX_BASE_U=RELAX_U
      RELAX_BASE_VW=RELAX_VW
      RELAX_BASE_P=RELAX_P
      DRELAX_START_SCALE=1.0D0
      WRITE(*,*)'MARQUARDT_RELAX_START U/VW/P=',RELAX_U,RELAX_VW,RELAX_P
      WRITE(*,*)'RAW_MAX_UPDATE U/V/W/P=',DRAW_MAX_U,DRAW_MAX_V,DRAW_MAX_W,DRAW_MAX_P
      WRITE(*,*)'TRUST_REGION_SCALE U/V/W/P=',DSTATE_MAX_U,DSTATE_MAX_V,DSTATE_MAX_W,DSTATE_MAX_P
      IF(ITRUST_UPDATE_CAP.NE.0) THEN
        WRITE(*,*)'TRUST_REGION_CAP U/VW/P=',DUPDATE_CAP_U,DUPDATE_CAP_VW,DUPDATE_CAP_P
      ELSE
        WRITE(*,*)'TRUST_REGION_CAP DISABLED'
      END IF
      WRITE(*,*)'TRIAL_ERR1_DIAGNOSTIC=',DERR1_TRIAL
      RELAX_RETRY_U=RELAX_U
      RELAX_RETRY_VW=RELAX_VW
      RELAX_RETRY_P=RELAX_P
      WW4_PREV=WW4
      DO I=1,4*IAT*NX
        FFS_PHYS_BASE(I)=FFS_PHYS(I)
      END DO
      WW4_LINE_BEST=1.0D300
      RELAX_BEST_U=RELAX_RETRY_U
      RELAX_BEST_VW=RELAX_RETRY_VW
      RELAX_BEST_P=RELAX_RETRY_P
      IBEST_BACKTRACKS=0
      IF(IPENDING_STEP.EQ.0) THEN
        ILINESEARCH_TRIES=0
        ILINESEARCH_ACCEPTED=0
        ILINESEARCH_EXTRA_MODE=0
        ISTAG_SD_TRIED=0
      END IF
      IPENDING_STEP=1
      STEP_SIGN=UPDATE_SIGN
      DS_OLD=0.0D0
      DS_PRED=0.0D0
      DROWSUM=0.0D0
      DPHYS_OLD=0.0D0
      DPHYS_PLUS=0.0D0
      DPHYS_MINUS=0.0D0
      DO I=1,N
        DFSCALED=FF1(I)*ROW_SCALE(I)
        DAXI=(RELAX_RETRY_U*CG_AP(I)+RELAX_RETRY_VW*CG_R(I)+RELAX_RETRY_P*CG_Z(I))*ROW_SCALE(I)
        DS_OLD=DS_OLD+DFSCALED*DFSCALED
        DS_PRED=DS_PRED+(DFSCALED+DAXI)*(DFSCALED+DAXI)
        DROWSUM=DROWSUM+(DFSCALED-DAXI)*(DFSCALED-DAXI)
        DPHYS_JD=DAXI
        DPHYS_OLD=DPHYS_OLD+FFS_PHYS(I)*FFS_PHYS(I)
        DPHYS_PLUS=DPHYS_PLUS+(FFS_PHYS(I)+DPHYS_JD)*(FFS_PHYS(I)+DPHYS_JD)
        DPHYS_MINUS=DPHYS_MINUS+(FFS_PHYS(I)-DPHYS_JD)*(FFS_PHYS(I)-DPHYS_JD)
      END DO
      WRITE(*,*)'MARQUARDT_SIGN_PREDICT NORM_OLD/PLUS/MINUS=',DSQRT(DS_OLD),DSQRT(DS_PRED),DSQRT(DROWSUM)
      WRITE(*,*)'MARQUARDT_SIGN_PHYS_PREDICT NORM_OLD/PLUS/MINUS=',DSQRT(DPHYS_OLD),DSQRT(DPHYS_PLUS),DSQRT(DPHYS_MINUS)
      IF(DPHYS_MINUS.LT.DPHYS_PLUS) THEN
        STEP_SIGN=-UPDATE_SIGN
        WRITE(*,*)'MARQUARDT_SIGN_SELECT MINUS PHYS'
      ELSE
        WRITE(*,*)'MARQUARDT_SIGN_SELECT PLUS PHYS'
      END IF
      IF(IUSE_STEEPEST.NE.0) THEN
        STEP_SIGN=DSD_PHYS_SIGN
        WRITE(*,*)'MARQUARDT_STEEPEST_SIGN_FORCE STEP_SIGN/DDIR=',STEP_SIGN,DDIR_STEEPEST
      END IF
      DBEST_STEP_SIGN=STEP_SIGN
      IOPPOSITE_TRIED=0
      DO I=1,NX
      DO ITT=1,IAT
        U_SAVE(ITT,I)=U(ITT,I)
        V_SAVE(ITT,I)=V(ITT,I)
        W_SAVE(ITT,I)=W(ITT,I)
        P_SAVE(ITT,I)=P(ITT,I)
      END DO
      END DO

      DO I=1,4*IAT*NX
        IF(IS_CONSTRAINED(I).NE.0) THEN
          IF(I.LE.IAT*NX) THEN
             YUVWP(I)=YUVWP(I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(I)
          ELSE IF(I.LE.3*IAT*NX) THEN
             YUVWP(I)=YUVWP(I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(I)
          ELSE
             YUVWP(I)=YUVWP(I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(I)
          END IF
        ELSE IF(I.GE.1.AND.I.LE.IAT*NX) THEN
           YUVWP(I)=YUVWP(I)+STEP_SIGN*RELAX_RETRY_U*DDUVWP(I)
        ELSE IF(I.GT.IAT*NX.AND.I.LE.3*IAT*NX) THEN
           YUVWP(I)=YUVWP(I)+STEP_SIGN*RELAX_RETRY_VW*DDUVWP(I)
        ELSE IF(I.GT.3*IAT*NX) THEN
           YUVWP(I)=YUVWP(I)+STEP_SIGN*RELAX_RETRY_P*DDUVWP(I)
        END IF
      END DO
           
      
      ERR1_ABS=0.0D0
      ERR1_REL_L2=0.0D0
      ERR1_REL_MAX=0.0D0
      DU_SCALE=0.0D0
      DV_SCALE=0.0D0
      DW_SCALE=0.0D0
      DP_SCALE=0.0D0
      DUMAX=0.0D0
      DVMAX=0.0D0
      DWMAX=0.0D0
      DPMAX=0.0D0
      IUMAX=0
      IVMAX=0
      IWMAX=0
      IPMAX=0
      DO I=1,4*IAT*NX
        IF(I.LE.IAT*NX) THEN
          DU_SCALE=DMAX1(DU_SCALE,DABS(YUVWP(I)))
        ELSE IF(I.LE.2*IAT*NX) THEN
          DV_SCALE=DMAX1(DV_SCALE,DABS(YUVWP(I)))
        ELSE IF(I.LE.3*IAT*NX) THEN
          DW_SCALE=DMAX1(DW_SCALE,DABS(YUVWP(I)))
        ELSE
          DP_SCALE=DMAX1(DP_SCALE,DABS(YUVWP(I)))
        END IF
      END DO
      DU_SCALE=DMAX1(DU_SCALE,1.0D-12)
      DV_SCALE=DMAX1(DV_SCALE,1.0D-12)
      DW_SCALE=DMAX1(DW_SCALE,1.0D-12)
      DP_SCALE=DMAX1(DP_SCALE,1.0D-12)
      DO I=1,4*IAT*NX
        IF(IS_CONSTRAINED(I).NE.0) THEN
          IF(I.LE.IAT*NX) THEN
            DDU_USED=RELAX_RETRY_U*DDUVWP(I)
          ELSE IF(I.LE.3*IAT*NX) THEN
            DDU_USED=RELAX_RETRY_VW*DDUVWP(I)
          ELSE
            DDU_USED=RELAX_RETRY_P*DDUVWP(I)
          END IF
        ELSE IF(I.GE.1.AND.I.LE.IAT*NX) THEN
          DDU_USED=RELAX_RETRY_U*DDUVWP(I)
        ELSE IF(I.GT.IAT*NX.AND.I.LE.3*IAT*NX) THEN
          DDU_USED=RELAX_RETRY_VW*DDUVWP(I)
        ELSE IF(I.GT.3*IAT*NX) THEN
          DDU_USED=RELAX_RETRY_P*DDUVWP(I)
        END IF
        ERR1_ABS=ERR1_ABS+DDU_USED**2
        IF(I.LE.IAT*NX) THEN
          DDU_REL=DABS(DDU_USED)/DU_SCALE
        ELSE IF(I.LE.2*IAT*NX) THEN
          DDU_REL=DABS(DDU_USED)/DV_SCALE
        ELSE IF(I.LE.3*IAT*NX) THEN
          DDU_REL=DABS(DDU_USED)/DW_SCALE
        ELSE
          DDU_REL=DABS(DDU_USED)/DP_SCALE
        END IF
        ERR1_REL_L2=ERR1_REL_L2+DDU_REL**2
        ERR1_REL_MAX=DMAX1(ERR1_REL_MAX,DDU_REL)
        IF(I.LE.IAT*NX) THEN
          IF(DABS(DDUVWP(I)).GT.DUMAX) THEN
            DUMAX=DABS(DDUVWP(I))
            IUMAX=I
          END IF
        ELSE IF(I.LE.2*IAT*NX) THEN
          IF(DABS(DDUVWP(I)).GT.DVMAX) THEN
            DVMAX=DABS(DDUVWP(I))
            IVMAX=I
          END IF
        ELSE IF(I.LE.3*IAT*NX) THEN
          IF(DABS(DDUVWP(I)).GT.DWMAX) THEN
            DWMAX=DABS(DDUVWP(I))
            IWMAX=I
          END IF
        ELSE
          IF(DABS(DDUVWP(I)).GT.DPMAX) THEN
            DPMAX=DABS(DDUVWP(I))
            IPMAX=I
          END IF
        END IF
      END DO
      ERR1_ABS=DSQRT(ERR1_ABS)
      ERR1_REL_L2=DSQRT(ERR1_REL_L2)
      ERR1=ERR1_REL_MAX
      WRITE(*,*)'ERR1_CURRENT=',ERR1
      WRITE(*,*)'ERR1_ABS_L2=',ERR1_ABS
      WRITE(*,*)'ERR1_REL_L2=',ERR1_REL_L2
      WRITE(*,*)'ERR1_REL_MAX=',ERR1_REL_MAX
      WRITE(*,*)'ERR1_SCALE U,V,W,P=',DU_SCALE,DV_SCALE,DW_SCALE,DP_SCALE
      WRITE(*,*)'MAX_UPDATE U,V,W,P=',DUMAX,DVMAX,DWMAX,DPMAX
      IF(IPMAX.GT.0) THEN
        ILOCAL=IPMAX-3*NX*IAT
        IPNODE=(ILOCAL-1)/IAT+1
        IPTIME=ILOCAL-(IPNODE-1)*IAT
        WRITE(*,*)'MAX_PRESSURE_UPDATE_INDEX=',IPMAX,'NODE=',IPNODE,'TIME=',IPTIME
        IF(ALLOCATED(SSMMK)) THEN
          DROWSUM=0.0D0
          DCOLSUM=0.0D0
          DO J=1,4*NX*IAT
            DROWSUM=DROWSUM+DABS(SSMMK(IPMAX,J))
            DCOLSUM=DCOLSUM+DABS(SSMMK(J,IPMAX))
          END DO
          WRITE(*,*)'MAX_PRESSURE_ROW_DIAG_ROWSUM_COLSUM_RHS=',SSMMK(IPMAX,IPMAX),DROWSUM,DCOLSUM,FFS(IPMAX)
        END IF
      END IF
      IF(DUMAX.GT.1.0D6.OR.DVMAX.GT.1.0D6.OR.DWMAX.GT.1.0D6.OR.DPMAX.GT.1.0D6) THEN
        WRITE(*,*)'NEWTON UPDATE IS TOO LARGE. STOP NEWTON ITERATION.'
        STOP
      END IF
      IF(ERR1.NE.ERR1.OR.WW4.NE.WW4) THEN
        WRITE(*,*)'NEWTON RESIDUAL BECAME NaN. STOP NEWTON ITERATION.'
        STOP
      END IF
      
      DO IK=1,4  
      DO I=1,NX
      DO ITT=1,IAT
        ISP=ITT+(I-1)*IAT+NX*IAT*(IK-1)
        IF(ISP.LE.NX*IAT) THEN
          U(ITT,I)=YUVWP(ISP)
        END IF
        IF((ISP.GT.NX*IAT).AND.(ISP.LE.2*NX*IAT)) THEN
          V(ITT,I)=YUVWP(ISP)
        END IF
        IF((ISP.GT.2*NX*IAT).AND.(ISP.LE.3*NX*IAT)) THEN
          W(ITT,I)=YUVWP(ISP)
        END IF
        IF(ISP.GT.3*NX*IAT) THEN
          P(ITT,I)=YUVWP(ISP)
        END IF
      END DO
      END DO
      END DO

  128 CONTINUE
      IF(IPENDING_STEP.NE.0) THEN
        WRITE(*,*)'MARQUARDT TRUST TRIAL PENDING RELAX_U/VW/P=',RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
        GO TO 30
      END IF

      IF(.TRUE.) THEN
        DUMIN=1.0D300
        DUMAXV=-1.0D300
        DVMIN=1.0D300
        DVMAXV=-1.0D300
        DWMIN=1.0D300
        DWMAXV=-1.0D300
        DPMIN=1.0D300
        DPMAXV=-1.0D300
        IUMIN=0
        ITUMIN=0
        IF(IPK.EQ.0) THEN
          OPEN(89,FILE='results\newton_uvwp_by_step.dat',STATUS='REPLACE')
          WRITE(89,*)'STEP ITT NODE X Y Z U V W P'
        ELSE
          OPEN(89,FILE='results\newton_uvwp_by_step.dat',STATUS='UNKNOWN',POSITION='APPEND')
        END IF
        DO I=1,NX
        DO ITT=1,IAT
          IF(U(ITT,I).LT.DUMIN) THEN
            DUMIN=U(ITT,I)
            IUMIN=I
            ITUMIN=ITT
          END IF
          DUMAXV=DMAX1(DUMAXV,U(ITT,I))
          DVMIN=DMIN1(DVMIN,V(ITT,I))
          DVMAXV=DMAX1(DVMAXV,V(ITT,I))
          DWMIN=DMIN1(DWMIN,W(ITT,I))
          DWMAXV=DMAX1(DWMAXV,W(ITT,I))
          DPMIN=DMIN1(DPMIN,P(ITT,I))
          DPMAXV=DMAX1(DPMAXV,P(ITT,I))
          WRITE(89,'(I8,1X,I4,1X,I8,1X,7(ES24.16,1X))') IPK+1,ITT,I,XX(I),YY(I),ZZ(I),U(ITT,I),V(ITT,I),W(ITT,I),P(ITT,I)
        END DO
        END DO
        CLOSE(89)
        WRITE(*,*)'STEP_U_MIN_MAX=',IPK+1,DUMIN,DUMAXV
        WRITE(*,*)'STEP_V_MIN_MAX=',IPK+1,DVMIN,DVMAXV
        WRITE(*,*)'STEP_W_MIN_MAX=',IPK+1,DWMIN,DWMAXV
        WRITE(*,*)'STEP_P_MIN_MAX=',IPK+1,DPMIN,DPMAXV
        IF(IGAUSS_NEWTON_ONLY.NE.0) THEN
          DLM_LAMBDA_CURRENT=0.0D0
          DLM_LAMBDA_USED=0.0D0
        END IF
        WRITE(*,*)'STEP_LAMBDA_NEXT=',IPK+1,DLM_LAMBDA_CURRENT
        WRITE(*,*)'STEP_LAMBDA_USED_NEXT=',IPK+1,DLM_LAMBDA_USED,DLM_LAMBDA_CURRENT
        WRITE(*,*)'STEP_RESIDUAL_SUM=',IPK+1,WW4**2
        WRITE(*,*)'STEP_RESIDUAL_NORM=',IPK+1,WW4
        WRITE(*,*)'STEP_BLOCK_RESIDUAL_NORM CONT/V/W/U=',IPK+1,DSQRT(DRU),DSQRT(DRV),DSQRT(DRW),DSQRT(DRP)
        WRITE(*,*)'STEP_BLOCK_RESIDUAL_MAX CONT/V/W/U=',IPK+1,DRESMAX_CONT,DRESMAX_V,DRESMAX_W,DRESMAX_U
        WRITE(*,*)'STEP_ERR1_WW4=',IPK+1,ERR1,WW4
        WRITE(*,*)'STEP_ERR1_ABS_REL_L2_REL_MAX=',IPK+1,ERR1_ABS,ERR1_REL_L2,ERR1_REL_MAX
        WRITE(*,*)'STEP_SUMMARY=',IPK+1,' LAMBDA_USED=',DLM_LAMBDA_USED,' LAMBDA_NEXT=',DLM_LAMBDA_CURRENT,' ERR1=',ERR1,' WW4=',WW4
        WRITE(*,*)'STEP_SUMMARY_BLOCK_NORM CONT/V/W/U=',IPK+1,DSQRT(DRU),DSQRT(DRV),DSQRT(DRW),DSQRT(DRP)
        WRITE(*,*)'STEP_SUMMARY_UVWP_MINMAX=',IPK+1,DUMIN,DUMAXV,DVMIN,DVMAXV,DWMIN,DWMAXV,DPMIN,DPMAXV
        REWIND(2)
        DO I=1,NX
        DO ITT=1,IAT
          WRITE(2,2100) T+DT*TTC(ITT),XX(I),YY(I),ZZ(I),U(ITT,I),V(ITT,I),W(ITT,I),P(ITT,I),Re
        END DO
        END DO
        FLUSH(2)
        WRITE(*,*)'STEP_CHECKPOINT_WRITTEN=',IPK+1
        ION_WALL=0
        ION_INLET=0
        ION_POUT=0
        DO I3=1,NSB1
          IF(NSBU(I3).EQ.IUMIN) ION_WALL=1
        END DO
        DO I3=1,NSB2
          IF(NSBUI(I3).EQ.IUMIN) ION_INLET=1
        END DO
        DO I3=1,NSB3
          IF(NSBUO(I3).EQ.IUMIN) ION_POUT=1
        END DO
        WRITE(*,*)'STEP_U_MIN_NODE/TIME/XYZ/WALL/INLET/POUT=',IPK+1,IUMIN,ITUMIN,XX(IUMIN),YY(IUMIN),ZZ(IUMIN),ION_WALL,ION_INLET,ION_POUT
      END IF

      IF(IPK.EQ.0) THEN
        DUMIN=1.0D300
        DUMAXV=-1.0D300
        DVMIN=1.0D300
        DVMAXV=-1.0D300
        DWMIN=1.0D300
        DWMAXV=-1.0D300
        DPMIN=1.0D300
        DPMAXV=-1.0D300
        OPEN(88,FILE='results\newton1_uvwp_debug.dat',STATUS='REPLACE')
        WRITE(88,*)'ITT NODE X Y Z U V W P'
        DO I=1,NX
        DO ITT=1,IAT
          DUMIN=DMIN1(DUMIN,U(ITT,I))
          DUMAXV=DMAX1(DUMAXV,U(ITT,I))
          DVMIN=DMIN1(DVMIN,V(ITT,I))
          DVMAXV=DMAX1(DVMAXV,V(ITT,I))
          DWMIN=DMIN1(DWMIN,W(ITT,I))
          DWMAXV=DMAX1(DWMAXV,W(ITT,I))
          DPMIN=DMIN1(DPMIN,P(ITT,I))
          DPMAXV=DMAX1(DPMAXV,P(ITT,I))
          WRITE(88,'(I4,1X,I8,1X,7(ES24.16,1X))') ITT,I,XX(I),YY(I),ZZ(I),U(ITT,I),V(ITT,I),W(ITT,I),P(ITT,I)
        END DO
        END DO
        CLOSE(88)
        WRITE(*,*)'NEWTON1_U_MIN_MAX=',DUMIN,DUMAXV
        WRITE(*,*)'NEWTON1_V_MIN_MAX=',DVMIN,DVMAXV
        WRITE(*,*)'NEWTON1_W_MIN_MAX=',DWMIN,DWMAXV
        WRITE(*,*)'NEWTON1_P_MIN_MAX=',DPMIN,DPMAXV
      END IF
      
!     Use the residual norm as the nonlinear convergence criterion.
!     ERR1 is a diagnostic update size and can remain finite near convergence.
      IF(WW4.GT.1.0D-8.AND.IPK.LT.INEWTON_MAX) THEN
        IPK=IPK+1
        WRITE(*,*)'ERR,FX,IPK,IPK1,Re=',ERR1,WW4,IPK,IPK1,Re
        GO TO 30
      ELSE
      IF(WW4.GT.1.0D-8) WRITE(*,*)'NEWTON_MAX_REACHED IPK/MAX/WW4=',IPK,INEWTON_MAX,WW4
      IF(WW4.LE.1.0D-8) WRITE(*,*)'NEWTON_CONVERGED_BY_WW4 IPK/WW4=',IPK,WW4
      IPK1=IPK1+1
      WRITE(*,*)'ERR,FX,IPK,IPK1,Re=',ERR1,WW4,IPK,IPK1,Re
      END IF
      
        CLOSE(2)
        OPEN(92,FILE=FNAME_RESULT,STATUS='REPLACE',ACTION='WRITE')
        DO I=1,NX
        DO ITT=1,IAT
           WRITE(92,2100) T+DT*TTC(ITT),XX(I),YY(I),ZZ(I),U(ITT,I),V(ITT,I),W(ITT,I),P(ITT,I),Re
        END DO
        END DO
        FLUSH(92)
        CLOSE(92)
        WRITE(*,*)'FINAL_RESULT_WRITTEN=',TRIM(FNAME_RESULT)
        
        T=T+DT
      
      
      END DO
        
      WRITE(*,*)'SUCSESSFUL GOOD LUCK !'
      
      

        
 2100   FORMAT(8D25.16)
!           do i=1,nx*iat*4
!                     if(ssmmk(i,i)==0.0e0) then
!                 write(*,*) i,i/(nx),ssmmk(i,i)
!                   end if
!           end do
           
!     Unit 2 is closed before the final result is rewritten with unit 92.

! FLOW CHECK ------------------------
      DO J=1,NAZ
      DO I=1,NAY
        IK=I+NAY*(J-1)
        WW(IK)=W1(I)*W2(J)
      END DO
      END DO

! SECTION 1 (INFLOW)----
      IM=1
      FLOW1=0.0D0
      DO I1=1,IAN
         IF(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ) THEN
         II=IE(I1,IM)
         V1=U(IAT,II)
         FLOW1=FLOW1+V1*WW(I1-8-4*(NAX+NAY+NAZ)-2*NAX*NAY-2*NAX*NAZ-NAY*NAZ)
         END IF
      END DO
      FLOW1=FLOW1*D0*H0
      WRITE(*,*) FLOW1
! SECTION 2 (Straight Section 1)
      IM=1
      FLOW2=0.0D0
      DO I1=1,IAN
         IF(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ) THEN
         II=IE(I1,IM)
         V1=U(IAT,II)
         FLOW2=FLOW2+V1*WW(I1-8-4*(NAX+NAY+NAZ)-2*NAX*NAY-NAX*NAZ)
         END IF
      END DO
      FLOW2=FLOW2*D0*H0
      WRITE(*,*) FLOW2
! SECTION 3 (Straight Section 2)
      IM=3
      FLOW3=0.0D0
      DO I1=1,IAN
         IF(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ) THEN
         II=IE(IAT,IM)
         V1=U(IAT,II)
         FLOW3=FLOW3+V1*WW(I1-8-4*(NAX+NAY+NAZ)-2*NAX*NAY-NAX*NAZ)
         END IF
      END DO
      FLOW3=FLOW3*D0*H0
      WRITE(*,*) FLOW3
! SECTION 4 (Outlet)
      IM=5
      FLOW4=0.0D0
      DO I1=1,IAN
         IF(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ) THEN
         II=IE(I1,IM)
         V1=V(IAT,II)
         FLOW4=FLOW4+V1*WW(I1-8-4*(NAX+NAY+NAZ)-2*NAX*NAY-NAX*NAZ)
         END IF
      END DO
      FLOW4=FLOW4*D0*H0
      WRITE(*,*) FLOW4
      
      WRITE(*,*)'OUTPUT FILE NAME OF SECTION FLOW : ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)
        WRITE(2,*)'FLOW RATE AT CHANNEL SECTION --------'
        WRITE(2,*)'FLOW1(Inlet)             =',FLOW1
        WRITE(2,*)'FLOW2(Straight Section 1)=',FLOW2
        WRITE(2,*)'FLOW3(Straight Section 2)=',FLOW3
        WRITE(2,*)'FLOW4(Outlet)            =',FLOW4
     CLOSE(2)
      
!
      RETURN
      END

! *****************************************************************
!
      SUBROUTINE OZ_ADD_SLICE6(DSLICES,DVALUE)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (IOZ_SLICE_MAX=6)
      REAL(8) DSLICES(IOZ_SLICE_MAX)
      DADD=DVALUE
      DCHK=DADD-DADD
      IF(DCHK.NE.DCHK) RETURN
      DO K=1,IOZ_SLICE_MAX
        DSUM=DSLICES(K)+DADD
        DCHK=DSUM-DSUM
        IF(DCHK.NE.DCHK) RETURN
        DZ=DSUM-DSLICES(K)
        DADD=(DSLICES(K)-(DSUM-DZ))+(DADD-DZ)
        DSLICES(K)=DSUM
      END DO
      RETURN
      END

! *****************************************************************
!
      SUBROUTINE OZ_ADD_PRODUCT_SLICE6(DSLICES,DX,DY)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (DSPLIT=134217729.0D0)
      PARAMETER (IOZ_SLICE_MAX=6)
      REAL(8) DSLICES(IOZ_SLICE_MAX)
      DPROD=DX*DY
      DCHK=DPROD-DPROD
      IF(DCHK.NE.DCHK) THEN
        DPROD=0.0D0
        DERR=0.0D0
        RETURN
      END IF
      DCX=DSPLIT*DX
      DXH=DCX-(DCX-DX)
      DXL=DX-DXH
      DCY=DSPLIT*DY
      DYH=DCY-(DCY-DY)
      DYL=DY-DYH
      DERR=((DXH*DYH-DPROD)+DXH*DYL+DXL*DYH)+DXL*DYL
      DCHK=DERR-DERR
      IF(DCHK.NE.DCHK) DERR=0.0D0
      CALL OZ_ADD_SLICE6(DSLICES,DPROD)
      CALL OZ_ADD_SLICE6(DSLICES,DERR)
      RETURN
      END

! *****************************************************************
!
      SUBROUTINE OZ_CSR_MV(N,IR,IC,A,X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (IOZ_SLICE_MAX=6)
      INTEGER N,IR(N+1),IC(*)
      REAL(8) A(*),X(*),Y(*),DSLICES(IOZ_SLICE_MAX)
      !$OMP PARALLEL DO PRIVATE(I,J,K,DSLICES) IF(N.GT.1000)
      DO I=1,N
        DO K=1,IOZ_SLICE_MAX
          DSLICES(K)=0.0D0
        END DO
        DO J=IR(I),IR(I+1)-1
          CALL OZ_ADD_PRODUCT_SLICE6(DSLICES,A(J),X(IC(J)))
        END DO
        Y(I)=0.0D0
        DO K=IOZ_SLICE_MAX,1,-1
          Y(I)=Y(I)+DSLICES(K)
        END DO
      END DO
      !$OMP END PARALLEL DO
      RETURN
      END

! *****************************************************************
!
      SUBROUTINE OZ_CSR_TMV(N,IR,IC,A,X,Y,WORK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (IOZ_SLICE_MAX=6)
      INTEGER N,IR(N+1),IC(*)
      REAL(8) A(*),X(*),Y(*),WORK(IOZ_SLICE_MAX,N)
      DO I=1,N
        Y(I)=0.0D0
        DO K=1,IOZ_SLICE_MAX
          WORK(K,I)=0.0D0
        END DO
      END DO
      DO I=1,N
        DO J=IR(I),IR(I+1)-1
          ICOL=IC(J)
          CALL OZ_ADD_PRODUCT_SLICE6(WORK(1,ICOL),A(J),X(I))
        END DO
      END DO
      DO I=1,N
        DO K=IOZ_SLICE_MAX,1,-1
          Y(I)=Y(I)+WORK(K,I)
        END DO
      END DO
      RETURN
      END

! *****************************************************************
!
      SUBROUTINE OZ_CSR_TMV_OMP(N,IR,IC,A,X,Y,WORKT)
      use omp_lib
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (IOZ_SLICE_MAX=6)
      PARAMETER (KOZ_WORK_THREADS=8)
      INTEGER N,IR(N+1),IC(*)
      REAL(8) A(*),X(*),Y(*),WORKT(IOZ_SLICE_MAX,N,KOZ_WORK_THREADS),DSUM(IOZ_SLICE_MAX)
      !$OMP PARALLEL DO PRIVATE(I,K,IT) NUM_THREADS(KOZ_WORK_THREADS) IF(N.GT.1000)
      DO I=1,N
        Y(I)=0.0D0
        DO IT=1,KOZ_WORK_THREADS
          DO K=1,IOZ_SLICE_MAX
            WORKT(K,I,IT)=0.0D0
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
      !$OMP PARALLEL PRIVATE(I,J,ICOL,ITID) NUM_THREADS(KOZ_WORK_THREADS) IF(N.GT.1000)
      ITID=OMP_GET_THREAD_NUM()+1
      !$OMP DO
      DO I=1,N
        DO J=IR(I),IR(I+1)-1
          ICOL=IC(J)
          CALL OZ_ADD_PRODUCT_SLICE6(WORKT(1,ICOL,ITID),A(J),X(I))
        END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
      !$OMP PARALLEL DO PRIVATE(I,K,IT,DSUM) NUM_THREADS(KOZ_WORK_THREADS) IF(N.GT.1000)
      DO I=1,N
        DO K=1,IOZ_SLICE_MAX
          DSUM(K)=0.0D0
        END DO
        DO IT=1,KOZ_WORK_THREADS
          DO K=1,IOZ_SLICE_MAX
            CALL OZ_ADD_SLICE6(DSUM,WORKT(K,I,IT))
          END DO
        END DO
        Y(I)=0.0D0
        DO K=IOZ_SLICE_MAX,1,-1
          Y(I)=Y(I)+DSUM(K)
        END DO
      END DO
      !$OMP END PARALLEL DO
      RETURN
      END

! *****************************************************************
!
      REAL(8) FUNCTION OZ_DDOT(N,X,INCX,Y,INCY)
      use omp_lib
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N,INCX,INCY
      REAL(8) X(*),Y(*)
      PARAMETER (DSPLIT=134217729.0D0)
      PARAMETER (IOZ_SLICE_MAX=6)
      PARAMETER (KOZ_WORK_THREADS=8)
      REAL(8) DSLICES(IOZ_SLICE_MAX),DLOCAL(IOZ_SLICE_MAX),DTHREAD(IOZ_SLICE_MAX,KOZ_WORK_THREADS)
      OZ_DDOT=0.0D0
      DO K=1,IOZ_SLICE_MAX
        DSLICES(K)=0.0D0
      END DO
      DO IT=1,KOZ_WORK_THREADS
        DO K=1,IOZ_SLICE_MAX
          DTHREAD(K,IT)=0.0D0
        END DO
      END DO
      !$OMP PARALLEL PRIVATE(I,K,IX,IY,DX,DY,DPROD,DCHK,DCX,DXH,DXL,DCY,DYH,DYL,DERR,DLOCAL,ITID) NUM_THREADS(KOZ_WORK_THREADS) IF(N.GT.1000)
      ITID=OMP_GET_THREAD_NUM()+1
      DO K=1,IOZ_SLICE_MAX
        DLOCAL(K)=0.0D0
      END DO
      !$OMP DO
      DO I=1,N
        IX=1+(I-1)*INCX
        IY=1+(I-1)*INCY
        DX=X(IX)
        DY=Y(IY)
        DPROD=DX*DY
        DCHK=DPROD-DPROD
        IF(DCHK.EQ.DCHK) THEN
          DCX=DSPLIT*DX
          DXH=DCX-(DCX-DX)
          DXL=DX-DXH
          DCY=DSPLIT*DY
          DYH=DCY-(DCY-DY)
          DYL=DY-DYH
          DERR=((DXH*DYH-DPROD)+DXH*DYL+DXL*DYH)+DXL*DYL
          DCHK=DERR-DERR
          IF(DCHK.NE.DCHK) DERR=0.0D0
          CALL OZ_ADD_SLICE6(DLOCAL,DPROD)
          CALL OZ_ADD_SLICE6(DLOCAL,DERR)
        END IF
      END DO
      !$OMP END DO
      DO K=1,IOZ_SLICE_MAX
        DTHREAD(K,ITID)=DLOCAL(K)
      END DO
      !$OMP END PARALLEL
      DO IT=1,KOZ_WORK_THREADS
        DO K=1,IOZ_SLICE_MAX
          CALL OZ_ADD_SLICE6(DSLICES,DTHREAD(K,IT))
        END DO
      END DO
      DO K=IOZ_SLICE_MAX,1,-1
        OZ_DDOT=OZ_DDOT+DSLICES(K)
      END DO
      RETURN
      END

! *****************************************************************
!
      INTEGER FUNCTION GLOBAL_DOF_INDEX(ITIME,INODE,IBLOCK,NX_LOCAL,IAT_LOCAL)
      INTEGER ITIME,INODE,IBLOCK,NX_LOCAL,IAT_LOCAL
!
!     Unknown/correction block order:
!       1: U, 2: V, 3: W, 4: P
!
      GLOBAL_DOF_INDEX=ITIME+(INODE-1)*IAT_LOCAL+(IBLOCK-1)*NX_LOCAL*IAT_LOCAL
      RETURN
      END

! *****************************************************************
!
      INTEGER FUNCTION GLOBAL_EQUATION_INDEX(ITIME,INODE,IBLOCK,NX_LOCAL,IAT_LOCAL)
      INTEGER ITIME,INODE,IBLOCK,NX_LOCAL,IAT_LOCAL
!
!     Equation row block order:
!       1: continuity, 2: V momentum, 3: W momentum, 4: U momentum
!
      GLOBAL_EQUATION_INDEX=ITIME+(INODE-1)*IAT_LOCAL+(IBLOCK-1)*NX_LOCAL*IAT_LOCAL
      RETURN
      END

! *****************************************************************
!
      INTEGER FUNCTION GLOBAL_EQUATION_BLOCK(LOCAL_BLOCK)
      INTEGER LOCAL_BLOCK
!
!     FFUVP/BBUVP local equation order:
!       1: continuity, 2: U momentum, 3: V momentum, 4: W momentum
!     Global residual row order:
!       1: continuity, 2: V momentum, 3: W momentum, 4: U momentum
!
      IF(LOCAL_BLOCK.EQ.1) THEN
        GLOBAL_EQUATION_BLOCK=1
      ELSE IF(LOCAL_BLOCK.EQ.2) THEN
        GLOBAL_EQUATION_BLOCK=4
      ELSE IF(LOCAL_BLOCK.EQ.3) THEN
        GLOBAL_EQUATION_BLOCK=2
      ELSE IF(LOCAL_BLOCK.EQ.4) THEN
        GLOBAL_EQUATION_BLOCK=3
      ELSE
        GLOBAL_EQUATION_BLOCK=LOCAL_BLOCK
      END IF
      RETURN
      END

! *****************************************************************
!
       subroutine csr_S_matrix(m,n,aa,mnr,ir,mn,ic,a)
       
       real(8) aa(n,n),a(10000*n)
       integer ir(n+1),ic(10000*n)

       ii=0
       ir(1)=1
       do i=1,m
         do j=1,n
           if(aa(i,j).ne.0.0d0) then
             ii=ii+1
             a(ii)=aa(i,j)
             ic(ii)=j
           end if
         end do
         ir(i+1)=ii+1
       end do
       mn=ii
       mnr=m+1

      return
      end


!     - invgj -

      subroutine invgj(a,n)
      implicit real*8(a-h,o-z)
      dimension a(9,9),noseq(1000)  
      if(n-1) 1,2,3
    1 goto 300 
    2 a(1,1)=1.0d0/a(1,1)
    3 epsl=0.1d-70
      do nn=1,n
      noseq(nn)=nn
      end do
      do nn=1,n
      p=0.0d0
      do i=nn,n
      if(p-dabs(a(i,1))) 21,20,20
   21 p=dabs(a(i,1))
      ip=i
   20 continue
      end do
      if(p-epsl) 4,4,5
    4 WRITE(*,*)'ILL CONDITION',p
      STOP
      goto 300 
    5 nw=noseq(ip)
      noseq(ip)=noseq(nn)
      noseq(nn)=nw
      do j=1,n
      w=a(ip,j)
      a(ip,j)=a(nn,j)
      a(nn,j)=w
      end do
      w=a(nn,1)
      do j=2,n
      a(nn,j-1)=a(nn,j)/w
      end do
      a(nn,n)=1.0d0/w
      do i=1,n
      if(i-nn) 51,50,51
   51 w=a(i,1)
      do j=2,n
      a(i,j-1)=a(i,j)-w*a(nn,j-1)
      end do
      a(i,n)=-w*a(nn,n)
   50 continue
      end do
      end do
      do nn=1,n
      do j=nn,n
      if(noseq(j)-nn) 70,71,70
   70 continue
      end do
   71 noseq(j)=noseq(nn)
      do i=1,n
      w=a(i,j)
      a(i,j)=a(i,nn)
      a(i,nn)=w
      end do
      end do
  300 return
      end


