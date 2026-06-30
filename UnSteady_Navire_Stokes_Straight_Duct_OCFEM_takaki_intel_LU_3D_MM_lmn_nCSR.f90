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
      use dof_layout, only: global_dof_index, global_equation_block, block_u, block_v, block_w, block_p
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
                tbxuvp1(knaa),tbyuvp1(knaa),tbzuvp1(knaa),tcxyuvp1(knaa),tcyzuvp1(knaa),tczxuvp1(knaa),tatuvp1(KNT)



      character FNAME*100,trans*2

      INTEGER :: iparm(64)
      INTEGER(8) :: ipt(64)
      INTEGER,ALLOCATABLE :: IPIV(:),iperm(:),ir(:),ic(:),IS_CONSTRAINED(:),BC_ROW(:),IS_REPLACED_ROW(:),ISINGULAR_PIN(:)
      REAL(8),ALLOCATABLE :: FFS(:),FFS_REPLACED_ROW(:),YUVP1(:),x(:,:),a(:),FFUVP(:),FFUVP2(:),YUVWP(:),DDUVWP(:),DDUVWP_SD(:),FF1(:),BC_DELTA(:),ROW_SCALE(:),COL_SCALE(:)
      REAL(8),ALLOCATABLE :: FFS_FD_BASE(:),JCOL_FD_BASE(:),JCOL_FD_DIFF(:)
      REAL(8),ALLOCATABLE :: U_SAVE(:,:),V_SAVE(:,:),W_SAVE(:,:),P_SAVE(:,:)
      REAL(8),ALLOCATABLE :: BBUVP(:,:)
      REAL(8),ALLOCATABLE :: SSMMK(:,:),DENSE_A(:,:),DENSE_B(:,:)
      type(triplet_matrix) :: SPARSE_A, SPARSE_ELEMENT_A, SPARSE_PRESSURE_WALL_A
      N=4*NX*IAT
      IF(N.GT.20000) THEN
        WRITE(*,*)'DENSE NEWTON MATRIX IS TOO LARGE FOR THIS SOLVER VERSION.'
        WRITE(*,*)'N=',N,' REQUIRED DENSE MEMORY BYTES=',8.0D0*DBLE(N)*DBLE(N)
        WRITE(*,*)'NEXT STEP: assemble the Newton matrix directly in CSR format.'
        STOP
      END IF
      ALLOCATE (BBUVP(4*IAN*IAT,4*IAN*IAT))
      ALLOCATE (IPIV(N),iperm(N),FFS(N),YUVP1(N),x(N,1))
      ALLOCATE (FFUVP(4*IAN*IAT),FFUVP2(4*IAN*IAT),YUVWP(N),DDUVWP(N),DDUVWP_SD(N),FF1(N))
      ALLOCATE (U_SAVE(IAT,NX),V_SAVE(IAT,NX),W_SAVE(IAT,NX),P_SAVE(IAT,NX))
      ALLOCATE (IS_CONSTRAINED(N),BC_ROW(N),IS_REPLACED_ROW(N),BC_DELTA(N),FFS_REPLACED_ROW(N))
      ALLOCATE (ROW_SCALE(N),COL_SCALE(N))
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
      IOPPOSITE_TRIED=0
      IPENDING_STEP=0
      ILINESEARCH_TRIES=0
      IPRESSURE_FD_CHECK=0
      IPRESSURE_FD_MODE=0
      IPRESSURE_FD_DONE=0
      IPRESSURE_FD_COL=0
      IPRESSURE_FD_NODE=0
      IPRESSURE_FD_TIME=0
      IAPPLY_INLET_VELOCITY_BC=1
      DPRESSURE_FD_EPS=1.0D-7
      DPRESSURE_FD_OLD=0.0D0
      FFS_FD_BASE=0.0D0
      JCOL_FD_BASE=0.0D0
      JCOL_FD_DIFF=0.0D0
      call triplet_init(SPARSE_ELEMENT_A, N, N, max(1, min(1000000, 20*N)))
      call triplet_init(SPARSE_PRESSURE_WALL_A, N, N, max(1, min(1000000, 2*N)))
      
 
      !$ call omp_set_num_threads(1)

!     ---------------  JACOBIAN MATRIX FOR NEWTON RAPHSON METHOD ------
      PAI=3.1415926535897932384626433832795028841971D0

      WRITE(*,*)'OUTPUT FILE NAME OF TRANSIENT RESULT : ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)

            T=0.0D0
 
    write(*,*) 'MX,NX',mx,nx
        !$OMP parallel private(I,J)
        !$OMP do
        DO J=1,IAT
        DO I=1,NX
          U(J,I)=0.0d0
          V(J,I)=0.0d0
          W(J,I)=0.0d0
          P(J,I)=0.0d0
    !      write(*,*)'clear 005! nx,iat',i,j
        END DO
        END DO
        !$OMP end do
        !$OMP end parallel
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
      call triplet_clear(SPARSE_PRESSURE_WALL_A)
      IS_REPLACED_ROW=0
      FFS_REPLACED_ROW=0.0D0


      !$OMP parallel private(I)
      !$OMP do
      DO I=1,NX*IAT*4
        FFS(I)=0.0D0
      END DO
      !$OMP end do
      !$OMP end parallel


      INX=IAN
      write(*,*)'clear 009!'
      


      DO IM=1,MX   
        DO I=1,4*IAN*IAT
        DO J=1,4*IAN*IAT
          BBUVP(I,J)=0.0D0
        END DO
        FFUVP(I)=0.0D0
        FFUVP2(I)=0.0D0
        END DO
     
        
        CALL TRANSF(IM,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)

        
        !$OMP parallel private(I1,II,IQ1,IQ2,ITT,IBB,ISP,J1,JJ,JTT,JBB)
        !$OMP do
        
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
            TAXU=ddot(IAN, TAXUVP1, 1, UU1, 1)
            TAYU=ddot(IAN, TAYUVP1, 1, UU1, 1)
            TAZU=ddot(IAN, TAZUVP1, 1, UU1, 1)
            TAXV=ddot(IAN, TAXUVP1, 1, VV1, 1)
            TAYV=ddot(IAN, TAYUVP1, 1, VV1, 1)
            TAZV=ddot(IAN, TAZUVP1, 1, VV1, 1)
            TAXW=ddot(IAN, TAXUVP1, 1, WW1, 1)
            TAYW=ddot(IAN, TAYUVP1, 1, WW1, 1)
            TAZW=ddot(IAN, TAZUVP1, 1, WW1, 1)
            TBXU=ddot(IAN, TBXUVP1, 1, UU1, 1)
            TBYV=ddot(IAN, TBYUVP1, 1, VV1, 1)
            TBZW=ddot(IAN, TBZUVP1, 1, WW1, 1)
            TCXYU=ddot(IAN, TCXYUVP1, 1, UU1, 1)
            TCZXU=ddot(IAN, TCZXUVP1, 1, UU1, 1)
            TCYZV=ddot(IAN, TCYZUVP1, 1, VV1, 1)
            TCXYV=ddot(IAN, TCXYUVP1, 1, VV1, 1)
            TCYZW=ddot(IAN, TCYZUVP1, 1, WW1, 1)
            TCZXW=ddot(IAN, TCZXUVP1, 1, WW1, 1)
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
                TAXU=ddot(IAN, TAXUVP1, 1, UU1, 1)
                TAYU=ddot(IAN, TAYUVP1, 1, UU1, 1)
                TAZU=ddot(IAN, TAZUVP1, 1, UU1, 1)
                TAXV=ddot(IAN, TAXUVP1, 1, VV1, 1)
                TAYV=ddot(IAN, TAYUVP1, 1, VV1, 1)
                TAZV=ddot(IAN, TAZUVP1, 1, VV1, 1)
                TAXW=ddot(IAN, TAXUVP1, 1, WW1, 1)
                TAYW=ddot(IAN, TAYUVP1, 1, WW1, 1)
                TAZW=ddot(IAN, TAZUVP1, 1, WW1, 1)
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
                TAXU=ddot(IAN, TAXUVP1, 1, UU1, 1)
                TAYU=ddot(IAN, TAYUVP1, 1, UU1, 1)
                TAZU=ddot(IAN, TAZUVP1, 1, UU1, 1)
                TAXV=ddot(IAN, TAXUVP1, 1, VV1, 1)
                TAYV=ddot(IAN, TAYUVP1, 1, VV1, 1)
                TAZV=ddot(IAN, TAZUVP1, 1, VV1, 1)
                TAXW=ddot(IAN, TAXUVP1, 1, WW1, 1)
                TAYW=ddot(IAN, TAYUVP1, 1, WW1, 1)
                TAZW=ddot(IAN, TAZUVP1, 1, WW1, 1)
                TBXU=ddot(IAN, TBXUVP1, 1, UU1, 1)
                TBYV=ddot(IAN, TBYUVP1, 1, VV1, 1)
                TBZW=ddot(IAN, TBZUVP1, 1, WW1, 1)
                TCXYU=ddot(IAN, TCXYUVP1, 1, UU1, 1)
                TCZXU=ddot(IAN, TCZXUVP1, 1, UU1, 1)
                TCYZV=ddot(IAN, TCYZUVP1, 1, VV1, 1)
                TCXYV=ddot(IAN, TCXYUVP1, 1, VV1, 1)
                TCYZW=ddot(IAN, TCYZUVP1, 1, WW1, 1)
                TCZXW=ddot(IAN, TCZXUVP1, 1, WW1, 1)
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

          TAXU=ddot(IAN, TAXUVP1, 1, UU1, 1)
          TAYU=ddot(IAN, TAYUVP1, 1, UU1, 1)
          TAZU=ddot(IAN, TAZUVP1, 1, UU1, 1)
          TAXV=ddot(IAN, TAXUVP1, 1, VV1, 1)
          TAYV=ddot(IAN, TAYUVP1, 1, VV1, 1)
          TAZV=ddot(IAN, TAZUVP1, 1, VV1, 1)
          TAXW=ddot(IAN, TAXUVP1, 1, WW1, 1)
          TAYW=ddot(IAN, TAYUVP1, 1, WW1, 1)
          TAZW=ddot(IAN, TAZUVP1, 1, WW1, 1)
          TAXP=ddot(IAN, TAXUVP1, 1, PP1, 1)
          TAYP=ddot(IAN, TAYUVP1, 1, PP1, 1)
          TAZP=ddot(IAN, TAZUVP1, 1, PP1, 1)
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
          TAXU=ddot(IAN, TAXUVP1, 1, UU1, 1)
          TAYU=ddot(IAN, TAYUVP1, 1, UU1, 1)
          TAZU=ddot(IAN, TAZUVP1, 1, UU1, 1)
          TAXV=ddot(IAN, TAXUVP1, 1, VV1, 1)
          TAYV=ddot(IAN, TAYUVP1, 1, VV1, 1)
          TAZV=ddot(IAN, TAZUVP1, 1, VV1, 1)
          TAXW=ddot(IAN, TAXUVP1, 1, WW1, 1)
          TAYW=ddot(IAN, TAYUVP1, 1, WW1, 1)
          TAZW=ddot(IAN, TAZUVP1, 1, WW1, 1)
          TBXU=ddot(IAN, TBXUVP1, 1, UU1, 1)
          TBYU=ddot(IAN, TBYUVP1, 1, UU1, 1)
          TBZU=ddot(IAN, TBZUVP1, 1, UU1, 1)
          TBXV=ddot(IAN, TBXUVP1, 1, VV1, 1)
          TBYV=ddot(IAN, TBYUVP1, 1, VV1, 1)
          TBZV=ddot(IAN, TBZUVP1, 1, VV1, 1)
          TBXW=ddot(IAN, TBXUVP1, 1, WW1, 1)
          TBYW=ddot(IAN, TBYUVP1, 1, WW1, 1)
          TBZW=ddot(IAN, TBZUVP1, 1, WW1, 1)
          TAXP=ddot(IAN, TAXUVP1, 1, PP1, 1)
          TAYP=ddot(IAN, TAYUVP1, 1, PP1, 1)
          TAZP=ddot(IAN, TAZUVP1, 1, PP1, 1)
          TBXP=ddot(IAN, TBXUVP1, 1, PP1, 1)
          TBYP=ddot(IAN, TBYUVP1, 1, PP1, 1)
          TBZP=ddot(IAN, TBZUVP1, 1, PP1, 1)
          TCXYU=0.0D0
          TCYZU=0.0D0
          TCZXU=0.0D0
          TCXYV=0.0D0
          TCYZV=0.0D0
          TCZXV=0.0D0
          TCXYW=0.0D0
          TCYZW=0.0D0
          TCZXW=0.0D0
          DO IJ=1,IAN
            TCXYU=TCXYU+(DJJ(7,1)*AX(I1,IJ)+DJJ(7,2)*AY(I1,IJ)+DJJ(7,3)*AZ(I1,IJ)+DJJ(7,4)*BX(I1,IJ)+DJJ(7,5)*BY(I1,IJ)+DJJ(7,6)*BZ(I1,IJ)+&
                 DJJ(7,7)*CXY(I1,IJ)+DJJ(7,8)*CYZ(I1,IJ)+DJJ(7,9)*CZX(I1,IJ))*UU1(IJ)
            TCYZU=TCYZU+(DJJ(8,1)*AX(I1,IJ)+DJJ(8,2)*AY(I1,IJ)+DJJ(8,3)*AZ(I1,IJ)+DJJ(8,4)*BX(I1,IJ)+DJJ(8,5)*BY(I1,IJ)+DJJ(8,6)*BZ(I1,IJ)+&
                 DJJ(8,7)*CXY(I1,IJ)+DJJ(8,8)*CYZ(I1,IJ)+DJJ(8,9)*CZX(I1,IJ))*UU1(IJ)
            TCZXU=TCZXU+(DJJ(9,1)*AX(I1,IJ)+DJJ(9,2)*AY(I1,IJ)+DJJ(9,3)*AZ(I1,IJ)+DJJ(9,4)*BX(I1,IJ)+DJJ(9,5)*BY(I1,IJ)+DJJ(9,6)*BZ(I1,IJ)+&
                 DJJ(9,7)*CXY(I1,IJ)+DJJ(9,8)*CYZ(I1,IJ)+DJJ(9,9)*CZX(I1,IJ))*UU1(IJ)
            TCXYV=TCXYV+(DJJ(7,1)*AX(I1,IJ)+DJJ(7,2)*AY(I1,IJ)+DJJ(7,3)*AZ(I1,IJ)+DJJ(7,4)*BX(I1,IJ)+DJJ(7,5)*BY(I1,IJ)+DJJ(7,6)*BZ(I1,IJ)+&
                 DJJ(7,7)*CXY(I1,IJ)+DJJ(7,8)*CYZ(I1,IJ)+DJJ(7,9)*CZX(I1,IJ))*VV1(IJ)
            TCYZV=TCYZV+(DJJ(8,1)*AX(I1,IJ)+DJJ(8,2)*AY(I1,IJ)+DJJ(8,3)*AZ(I1,IJ)+DJJ(8,4)*BX(I1,IJ)+DJJ(8,5)*BY(I1,IJ)+DJJ(8,6)*BZ(I1,IJ)+&
                 DJJ(8,7)*CXY(I1,IJ)+DJJ(8,8)*CYZ(I1,IJ)+DJJ(8,9)*CZX(I1,IJ))*VV1(IJ)
            TCZXV=TCZXV+(DJJ(9,1)*AX(I1,IJ)+DJJ(9,2)*AY(I1,IJ)+DJJ(9,3)*AZ(I1,IJ)+DJJ(9,4)*BX(I1,IJ)+DJJ(9,5)*BY(I1,IJ)+DJJ(9,6)*BZ(I1,IJ)+&
                 DJJ(9,7)*CXY(I1,IJ)+DJJ(9,8)*CYZ(I1,IJ)+DJJ(9,9)*CZX(I1,IJ))*VV1(IJ)
            TCXYW=TCXYW+(DJJ(7,1)*AX(I1,IJ)+DJJ(7,2)*AY(I1,IJ)+DJJ(7,3)*AZ(I1,IJ)+DJJ(7,4)*BX(I1,IJ)+DJJ(7,5)*BY(I1,IJ)+DJJ(7,6)*BZ(I1,IJ)+&
                 DJJ(7,7)*CXY(I1,IJ)+DJJ(7,8)*CYZ(I1,IJ)+DJJ(7,9)*CZX(I1,IJ))*WW1(IJ)
            TCYZW=TCYZW+(DJJ(8,1)*AX(I1,IJ)+DJJ(8,2)*AY(I1,IJ)+DJJ(8,3)*AZ(I1,IJ)+DJJ(8,4)*BX(I1,IJ)+DJJ(8,5)*BY(I1,IJ)+DJJ(8,6)*BZ(I1,IJ)+&
                 DJJ(8,7)*CXY(I1,IJ)+DJJ(8,8)*CYZ(I1,IJ)+DJJ(8,9)*CZX(I1,IJ))*WW1(IJ)
            TCZXW=TCZXW+(DJJ(9,1)*AX(I1,IJ)+DJJ(9,2)*AY(I1,IJ)+DJJ(9,3)*AZ(I1,IJ)+DJJ(9,4)*BX(I1,IJ)+DJJ(9,5)*BY(I1,IJ)+DJJ(9,6)*BZ(I1,IJ)+&
                 DJJ(9,7)*CXY(I1,IJ)+DJJ(9,8)*CYZ(I1,IJ)+DJJ(9,9)*CZX(I1,IJ))*WW1(IJ)
          END DO
          TDIVCONV=TAXU*TAXU+U(ITT,II)*TBXU+TAXV*TAYU+V(ITT,II)*TCXYU+TAXW*TAZU+W(ITT,II)*TCZXU+&
                 TAYU*TAXV+U(ITT,II)*TCXYV+TAYV*TAYV+V(ITT,II)*TBYV+TAYW*TAZV+W(ITT,II)*TCYZV+&
                 TAZU*TAXW+U(ITT,II)*TCZXW+TAZV*TAYW+V(ITT,II)*TCYZW+TAZW*TAZW+W(ITT,II)*TBZW
          TATU=ddot(IAT, TATUVP1, 1, UUT1,1)
          TATV=ddot(IAT, TATUVP1, 1, VVT1,1)
          TATW=ddot(IAT, TATUVP1, 1, WWT1,1)

            FFUVP(IIB)=TAXU+TAYV+TAZW
            FFUVP(IIB+IAN*IAT)=TATU+U(ITT,II)*TAXU+V(ITT,II)*TAYU+W(ITT,II)*TAZU-ANYU*(TBYU+TBXU+TBZU)+1.0D0/ROU*TAXP
            FFUVP(IIB+IAN*IAT*2)=TATV+U(ITT,II)*TAXV+V(ITT,II)*TAYV+W(ITT,II)*TAZV-ANYU*(TBXV+TBYV+TBZV)+1.0D0/ROU*TAYP
            FFUVP(IIB+IAN*IAT*3)=TATW+U(ITT,II)*TAXW+V(ITT,II)*TAYW+W(ITT,II)*TAZW-ANYU*(TBXW+TBYW+TBZW)+1.0D0/ROU*TAZP
 !           WRITE(*,*)'FFUVP=',FFUVP(IIB),FFUVP(IIB+IAN*IAT),FFUVP(IIB+2*IAN*IAT),FFUVP(IIB+3*IAN*IAT)
            END IF
            
!            write(*,*)'IM,i1,itt,iib=',IM,i1,itt,iib
              
        END DO
        END DO
        !$OMP end do
        !$OMP end parallel
        

        write(*,*)'clear ohkubo 66!'

        
!       ----------------- SUPERPOSITION OF BBUVP- AND TIME -------------

        DO IK1=1,4
!       Keep this loop serial while triplet assembly uses a shared append-only buffer.
!       Parallel calls to triplet_add race on SPARSE_ELEMENT_A%nnz and corrupt CSR input.
        DO JA1=1,IAN
            IF((JA1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.JA1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ).OR.&
               (JA1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.JA1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ).OR.&
               (JA1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ.AND.JA1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ+NAX*NAY*NAZ)) THEN
        DO ITT=1,IAT 
            II1=IE(JA1,IM)
!           Source equation order in FFUVP/BBUVP is CONT,U,V,W.
!           Global residual row order is CONT,V,W,U.
            IK1G=GLOBAL_EQUATION_BLOCK(IK1)
            ISP1=GLOBAL_DOF_INDEX(ITT,II1,IK1G,NX,IAT)
            IPL1=ITT+(JA1-1)*IAT+IAN*IAT*(IK1-1)
            IIT(ISP1)=IIT(ISP1)+1
!             Unknown and Newton correction column order is always U,V,W,P.
!             Equation row order is continuity,V-momentum,W-momentum,U-momentum.
              DO IK2=1,4
              DO JA2=1,IAN
              DO JTT=1,IAT
                JJ2=IE(JA2,IM)
                ISP2=GLOBAL_DOF_INDEX(JTT,JJ2,IK2,NX,IAT)
                IPL2=JTT+(JA2-1)*IAT+IAN*IAT*(IK2-1)
                call triplet_add(SPARSE_ELEMENT_A, ISP1, ISP2, BBUVP(IPL1,IPL2))
              END DO             
              END DO
              END DO 
            FFS(ISP1)=FFS(ISP1)+FFUVP(IPL1)
        END DO
        END IF
        END DO
        END DO

      END DO
!
      MN_ELEMENT=SPARSE_ELEMENT_A%nnz
      WRITE(*,*)'ELEMENT TRIPLET CSR NNZ=',MN_ELEMENT

              write(*,*)'clear ohkubo 1TC',TC

!     ---------------- SECOND DIFFERENTIAL .EQ. 0,  U,V,W  at B3 OUT FLOW --------------------
!     SLOPE U ---------
      
      !$OMP parallel private(J1,II,ITT,ISP1,IQ1,IQ2,IJ,JTT,ISP2)
      !$OMP do
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        DO ITT=1,IAT
        ISP1=GLOBAL_DOF_INDEX(ITT,II,BLOCK_U,NX,IAT)
        
        CALL TRANSF(5,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
          DO IJ=1,IAN
            ISP2=GLOBAL_DOF_INDEX(ITT,IE(IJ,5),BLOCK_U,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, (DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ)))
                 
          END DO
          
          TBYU=0.0D0
          DO IJ=1,IAN
          JJ=IE(IJ,5)
          TBYU=TBYU+(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))*U(ITT,JJ)
          END DO
          FFS(ISP1)=TBYU
          IS_REPLACED_ROW(ISP1)=1
          FFS_REPLACED_ROW(ISP1)=TBYU
        END DO
        END IF      
      END DO
      !$OMP end do
      !$OMP end parallel

!     SLOPE V -----------------
 
      
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        DO ITT=1,IAT
        ISP1=GLOBAL_DOF_INDEX(ITT,II,BLOCK_V,NX,IAT)
        
        CALL TRANSF(5,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
          !$OMP parallel private(IJ,ISP2)
          !$OMP do
          DO IJ=1,IAN
            ISP2=GLOBAL_DOF_INDEX(ITT,IE(IJ,5),BLOCK_V,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, (DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ)))
          END DO
          !$OMP end do
          !$OMP end parallel
          
          TBYV=0.0D0
          !$OMP parallel private(IJ,JJ)
          !$OMP do reduction(+ : TBYV)
          DO IJ=1,IAN
          JJ=IE(IJ,5)
          TBYV=TBYV+(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))*V(ITT,JJ)
          END DO
          !$OMP end do
          !$OMP end parallel
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
        ISP1=GLOBAL_DOF_INDEX(ITT,II,BLOCK_W,NX,IAT)
        
        CALL TRANSF(5,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
          !$OMP parallel private(IJ,ISP2)
          !$OMP do
          DO IJ=1,IAN
            ISP2=GLOBAL_DOF_INDEX(ITT,IE(IJ,5),BLOCK_W,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, (DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ)))
          END DO
          !$OMP end do
          !$OMP end parallel
          
          TBYW=0.0D0
          !$OMP parallel private(IJ,JJ)
          !$OMP do reduction(+ : TBYW)
          DO IJ=1,IAN
          JJ=IE(IJ,5)
          TBYW=TBYW+(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))*W(ITT,JJ)
          END DO
          !$OMP end do
          !$OMP end parallel
          FFS(ISP1)=TBYW
          IS_REPLACED_ROW(ISP1)=1
          FFS_REPLACED_ROW(ISP1)=TBYW
        END DO
        END IF      
      END DO

!     ------------ SECOND DIFFERENTIAL OF PRESSURE AT INFLOW ---------------- 

!     ------- SLOPE P AT INFLOW

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
          ISP1=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_P,NX,IAT)
          I_SKIP_WALL=0
          DO I3=1,NSB1
            IF(NSBU(I3).EQ.JJ) I_SKIP_WALL=1
          END DO
          IF(I_SKIP_WALL.NE.0) CYCLE
          !$OMP parallel private(IJ,ISP2)
          !$OMP do
          DO IJ=1,IAN
            ISP2=GLOBAL_DOF_INDEX(ITT,IE(IJ,J),BLOCK_P,NX,IAT)
            call triplet_add(SPARSE_PRESSURE_WALL_A, ISP1, ISP2, (DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)))
          END DO
          !$OMP end do
          !$OMP end parallel
          TBXP=0.0D0
          !$OMP parallel private(IJ,JJ)
          !$OMP do reduction(+ : TBXP)
          DO IJ=1,IAN
          JJ=IE(IJ,J)
          TBXP=TBXP+(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ))*P(ITT,JJ)
          END DO
          !$OMP end do
          !$OMP end parallel 
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
          ISP1=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_P,NX,IAT)
          I_SKIP_POUT=0
          DO I3=1,NSB3
            IF(NSBUO(I3).EQ.JJ) I_SKIP_POUT=1
          END DO
          IF(I_SKIP_POUT.NE.0) CYCLE
          !$OMP parallel private(IJ)
          !$OMP do
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
          !$OMP end do
          !$OMP end parallel
          TAXYZP=0.0D0
          TNLAPU=0.0D0
          TNLAPV=0.0D0
          TNLAPW=0.0D0
          !$OMP parallel private(IJ,JJ4)
          !$OMP do reduction(+ : TAXYZP,TNLAPU,TNLAPV,TNLAPW)
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
          !$OMP end do
          !$OMP end parallel
          TAXYZP=TAXYZP/ROU-ANYU*(XL(I1,J1)*TNLAPU+XM(I1,J1)*TNLAPV+XN(I1,J1)*TNLAPW)
          TDXP=0.0D0
          TDYP=0.0D0
          TDZP=0.0D0
          !$OMP parallel private(IJ,JJ4)
          !$OMP do reduction(+ : TDXP,TDYP,TDZP)
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
          !$OMP end do
          !$OMP end parallel
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
          ISP1=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_U,NX,IAT)
          BCVAL=SBU(J)-U(ITT,JJ)
          IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
          IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
          IS_CONSTRAINED(ISP2)=1
          BC_ROW(ISP2)=ISP1
          BC_DELTA(ISP2)=BCVAL

          ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_V,NX,IAT)
          ISP1=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_V,NX,IAT)
          BCVAL=SBV(J)-V(ITT,JJ)
          IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
          IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
          IS_CONSTRAINED(ISP2)=1
          BC_ROW(ISP2)=ISP1
          BC_DELTA(ISP2)=BCVAL

          ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_W,NX,IAT)
          ISP1=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_W,NX,IAT)
          BCVAL=SBW(J)-W(ITT,JJ)
          IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
          IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
          IS_CONSTRAINED(ISP2)=1
          BC_ROW(ISP2)=ISP1
          BC_DELTA(ISP2)=BCVAL
        END DO
      END DO

      DO I3=1,NSB3
        DO ITT=1,IAT
          II=NSBUO(I3)
          ISP2=GLOBAL_DOF_INDEX(ITT,II,BLOCK_P,NX,IAT)
          ISP1=GLOBAL_DOF_INDEX(ITT,II,BLOCK_P,NX,IAT)
         BCVAL=SBPO(I3)-P(ITT,II)
          IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
          IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
          IS_CONSTRAINED(ISP2)=1
          BC_ROW(ISP2)=ISP1
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
            ISP1=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_U,NX,IAT)
            BCVAL=SBUI(J)-U(ITT,JJ)
            IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
            IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
            IS_CONSTRAINED(ISP2)=1
            BC_ROW(ISP2)=ISP1
            BC_DELTA(ISP2)=BCVAL

            ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_V,NX,IAT)
            ISP1=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_V,NX,IAT)
            BCVAL=SBVI(J)-V(ITT,JJ)
            IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
            IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
            IS_CONSTRAINED(ISP2)=1
            BC_ROW(ISP2)=ISP1
            BC_DELTA(ISP2)=BCVAL

            ISP2=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_W,NX,IAT)
            ISP1=GLOBAL_DOF_INDEX(ITT,JJ,BLOCK_W,NX,IAT)
            BCVAL=SBWI(J)-W(ITT,JJ)
            IF(IS_CONSTRAINED(ISP2).EQ.0) NCONSTRAINT=NCONSTRAINT+1
            IF(IS_CONSTRAINED(ISP2).NE.0.AND.(BC_ROW(ISP2).NE.ISP1.OR.DABS(BC_DELTA(ISP2)-BCVAL).GT.1.0D-10)) NCONFLICT=NCONFLICT+1
            IS_CONSTRAINED(ISP2)=1
            BC_ROW(ISP2)=ISP1
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
        IF(IS_REPLACED_ROW(I).NE.0.AND.IS_CONSTRAINED(I).NE.0) THEN
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

      DO ITARGET=1,4
        IF(ITARGET.EQ.1) ITNODE=386
        IF(ITARGET.EQ.2) ITNODE=548
        IF(ITARGET.EQ.3) ITNODE=550
        IF(ITARGET.EQ.4) ITNODE=1786
        IF(ITNODE.GE.1.AND.ITNODE.LE.NX) THEN
          IROWCHK=1+(ITNODE-1)*IAT+3*NX*IAT
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

      IF(ALLOCATED(ir)) DEALLOCATE(ir)
      IF(ALLOCATED(ic)) DEALLOCATE(ic)
      IF(ALLOCATED(a)) DEALLOCATE(a)
      call triplet_init(SPARSE_A, n, n, max(1, min(4000000, SPARSE_ELEMENT_A%nnz+SPARSE_PRESSURE_WALL_A%nnz+NCONSTRAINT)))
      DO ISP2=1,N
        IF(IS_CONSTRAINED(ISP2).NE.0) THEN
          BCVAL=BC_DELTA(ISP2)
          DO I=1,SPARSE_ELEMENT_A%nnz
            IF(SPARSE_ELEMENT_A%col(I).EQ.ISP2.AND.IS_REPLACED_ROW(SPARSE_ELEMENT_A%row(I)).EQ.0) THEN
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
        IF(IS_REPLACED_ROW(SPARSE_ELEMENT_A%row(I)).EQ.0.AND.IS_CONSTRAINED(SPARSE_ELEMENT_A%row(I)).EQ.0.AND.IS_CONSTRAINED(SPARSE_ELEMENT_A%col(I)).EQ.0) THEN
          call triplet_add(SPARSE_A, SPARSE_ELEMENT_A%row(I), SPARSE_ELEMENT_A%col(I), SPARSE_ELEMENT_A%val(I))
        END IF
      END DO
      DO I=1,SPARSE_PRESSURE_WALL_A%nnz
        IF(IS_CONSTRAINED(SPARSE_PRESSURE_WALL_A%row(I)).EQ.0.AND.IS_CONSTRAINED(SPARSE_PRESSURE_WALL_A%col(I)).EQ.0) THEN
          call triplet_add(SPARSE_A, SPARSE_PRESSURE_WALL_A%row(I), SPARSE_PRESSURE_WALL_A%col(I), SPARSE_PRESSURE_WALL_A%val(I))
        END IF
      END DO
      DO I=1,N
        IF(IS_REPLACED_ROW(I).NE.0.AND.IS_CONSTRAINED(I).EQ.0) FFS(I)=FFS_REPLACED_ROW(I)
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
          FFS_FD_BASE=FFS
          JCOL_FD_BASE=0.0D0
          DO I=1,N
            DO J=ir(I),ir(I+1)-1
              IF(ic(J).EQ.IPRESSURE_FD_COL) JCOL_FD_BASE(I)=JCOL_FD_BASE(I)+a(J)
            END DO
          END DO
          DPRESSURE_FD_OLD=P(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)
          P(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD+DPRESSURE_FD_EPS
          IPRESSURE_FD_MODE=1
          WRITE(*,*)'PRESSURE FD JACOBIAN CHECK START COL/NODE/TIME/EPS=',IPRESSURE_FD_COL,IPRESSURE_FD_NODE,IPRESSURE_FD_TIME,DPRESSURE_FD_EPS
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
                DO J=1,IPRINT-1
                  IF(I.EQ.INT(ROW_SCALE(J))) ISEEN=1
                END DO
                IF(ISEEN.EQ.0) THEN
                  DMAX_PRINT=DABS(JCOL_FD_DIFF(I))
                  IMAX_PRINT=I
                END IF
              END IF
            END DO
            IF(IMAX_PRINT.LE.0) EXIT
            ROW_SCALE(IPRINT)=DBLE(IMAX_PRINT)
            IBLOCK=(IMAX_PRINT-1)/(NX*IAT)+1
            ILOCAL=IMAX_PRINT-(IBLOCK-1)*NX*IAT
            IPNODE=(ILOCAL-1)/IAT+1
            IPTIME=ILOCAL-(IPNODE-1)*IAT
            DFD=(FFS(IMAX_PRINT)-FFS_FD_BASE(IMAX_PRINT))/DPRESSURE_FD_EPS
            WRITE(98,'(4I10,1P,6E24.14)') IMAX_PRINT,IBLOCK,IPNODE,IPTIME,JCOL_FD_BASE(IMAX_PRINT),DFD,JCOL_FD_DIFF(IMAX_PRINT),&
                 DABS(JCOL_FD_DIFF(IMAX_PRINT)),FFS_FD_BASE(IMAX_PRINT),FFS(IMAX_PRINT)
          END DO
          CLOSE(98)
          IBLOCK=(IMAX_FD-1)/(NX*IAT)+1
          ILOCAL=IMAX_FD-(IBLOCK-1)*NX*IAT
          IPNODE=(ILOCAL-1)/IAT+1
          IPTIME=ILOCAL-(IPNODE-1)*IAT
          WRITE(*,*)'PRESSURE FD JACOBIAN CHECK DONE COL/NODE/TIME=',IPRESSURE_FD_COL,IPRESSURE_FD_NODE,IPRESSURE_FD_TIME
          WRITE(*,*)'PRESSURE FD MAX ABS DIFF ROW/BLOCK/NODE/TIME=',IMAX_FD,IBLOCK,IPNODE,IPTIME
          WRITE(*,*)'PRESSURE FD MAX_ABS_DIFF SUM_ABS_DIFF MAX_REL_DIFF=',DMAX_FD,DSUM_FD,DMAX_FD_REL
          WRITE(*,*)'PRESSURE FD DETAILS FILE=results\pressure_jacobian_fd_check.dat'
          P(IPRESSURE_FD_TIME,IPRESSURE_FD_NODE)=DPRESSURE_FD_OLD
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

      WW4=0.0D0
      !$OMP parallel private(I)
      !$OMP do reduction(+ : WW4)
      DO I=1,4*NX*IAT
        WW4=WW4+FFS(I)**2
      END DO
      !$OMP end do
      !$OMP end parallel
      WW4=DSQRT(WW4)
      
      WRITE(*,*)'WW4-',WW4
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
        DRU=DRU+FFS(I)**2
        DRV=DRV+FFS(I+NX*IAT)**2
        DRW=DRW+FFS(I+2*NX*IAT)**2
        DRP=DRP+FFS(I+3*NX*IAT)**2
        IF(DABS(FFS(I)).GT.DUMAX) THEN
          DUMAX=DABS(FFS(I))
          IUMAX=I
        END IF
        IF(DABS(FFS(I+NX*IAT)).GT.DVMAX) THEN
          DVMAX=DABS(FFS(I+NX*IAT))
          IVMAX=I+NX*IAT
        END IF
        IF(DABS(FFS(I+2*NX*IAT)).GT.DWMAX) THEN
          DWMAX=DABS(FFS(I+2*NX*IAT))
          IWMAX=I+2*NX*IAT
        END IF
        IF(DABS(FFS(I+3*NX*IAT)).GT.DPMAX) THEN
          DPMAX=DABS(FFS(I+3*NX*IAT))
          IPMAX=I+3*NX*IAT
        END IF
      END DO
      WRITE(*,*)'RESIDUAL_ROW_NORM CONT/V-MOM/W-MOM/U-MOM=',DSQRT(DRU),DSQRT(DRV),DSQRT(DRW),DSQRT(DRP)
      DO IBLOCK=1,4
        IF(IBLOCK.EQ.1) ILOCAL=IUMAX
        IF(IBLOCK.EQ.2) ILOCAL=IVMAX-NX*IAT
        IF(IBLOCK.EQ.3) ILOCAL=IWMAX-2*NX*IAT
        IF(IBLOCK.EQ.4) ILOCAL=IPMAX-3*NX*IAT
        IPNODE=(ILOCAL-1)/IAT+1
        IPTIME=ILOCAL-(IPNODE-1)*IAT
        IF(IBLOCK.EQ.1) WRITE(*,*)'MAX_RESIDUAL CONT ROW/NODE/TIME/VAL=',IUMAX,IPNODE,IPTIME,FFS(IUMAX)
        IF(IBLOCK.EQ.2) WRITE(*,*)'MAX_RESIDUAL V-MOM ROW/NODE/TIME/VAL=',IVMAX,IPNODE,IPTIME,FFS(IVMAX)
        IF(IBLOCK.EQ.3) WRITE(*,*)'MAX_RESIDUAL W-MOM ROW/NODE/TIME/VAL=',IWMAX,IPNODE,IPTIME,FFS(IWMAX)
        IF(IBLOCK.EQ.4) WRITE(*,*)'MAX_RESIDUAL U-MOM ROW/NODE/TIME/VAL=',IPMAX,IPNODE,IPTIME,FFS(IPMAX)
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
        WRITE(*,*)'MAX_RESIDUAL_NODE_INFO BLOCK/NODE/XYZ/CONSTR/WALL/INLET/POUT=',IBLOCK,IPNODE,XX(IPNODE),YY(IPNODE),ZZ(IPNODE),IS_CONSTRAINED(ILOCAL+NX*IAT*(IBLOCK-1)),ION_WALL,ION_INLET,ION_POUT
      END DO

      IF(IPENDING_STEP.NE.0) THEN
        IF(WW4.LT.WW4_PREV) THEN
          WRITE(*,*)'LINESEARCH ACCEPT WW4_PREV/WW4/RELAX_U/VW/P=',WW4_PREV,WW4,RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
          IPENDING_STEP=0
          ILINESEARCH_TRIES=0
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
              DO I=1,4*IAT*NX
                DDUVWP(I)=DDUVWP_SD(I)
              END DO
              STEP_SIGN=UPDATE_SIGN
              RELAX_RETRY_U=0.05D0
              RELAX_RETRY_VW=0.05D0
              RELAX_RETRY_P=0.02D0
              DRAW_MAX_U=0.0D0
              DRAW_MAX_V=0.0D0
              DRAW_MAX_W=0.0D0
              DRAW_MAX_P=0.0D0
              DO I=1,4*IAT*NX
                IF(I.LE.IAT*NX) THEN
                  DRAW_MAX_U=DMAX1(DRAW_MAX_U,DABS(DDUVWP(I)))
                ELSE IF(I.LE.2*IAT*NX) THEN
                  DRAW_MAX_V=DMAX1(DRAW_MAX_V,DABS(DDUVWP(I)))
                ELSE IF(I.LE.3*IAT*NX) THEN
                  DRAW_MAX_W=DMAX1(DRAW_MAX_W,DABS(DDUVWP(I)))
                ELSE
                  DRAW_MAX_P=DMAX1(DRAW_MAX_P,DABS(DDUVWP(I)))
                END IF
              END DO
              IF(DRAW_MAX_U*RELAX_RETRY_U.GT.DUPDATE_CAP_U.AND.DRAW_MAX_U.GT.0.0D0) RELAX_RETRY_U=DUPDATE_CAP_U/DRAW_MAX_U
              IF(DMAX1(DRAW_MAX_V,DRAW_MAX_W)*RELAX_RETRY_VW.GT.DUPDATE_CAP_VW.AND.DMAX1(DRAW_MAX_V,DRAW_MAX_W).GT.0.0D0) THEN
                RELAX_RETRY_VW=DUPDATE_CAP_VW/DMAX1(DRAW_MAX_V,DRAW_MAX_W)
              END IF
              IF(DRAW_MAX_P*RELAX_RETRY_P.GT.DUPDATE_CAP_P.AND.DRAW_MAX_P.GT.0.0D0) RELAX_RETRY_P=DUPDATE_CAP_P/DRAW_MAX_P
              DERR1_TRIAL=0.0D0
              DO I=1,4*IAT*NX
                IF(I.LE.IAT*NX) THEN
                  DERR1_TRIAL=DERR1_TRIAL+(RELAX_RETRY_U*DDUVWP(I))**2
                ELSE IF(I.LE.3*IAT*NX) THEN
                  DERR1_TRIAL=DERR1_TRIAL+(RELAX_RETRY_VW*DDUVWP(I))**2
                ELSE
                  DERR1_TRIAL=DERR1_TRIAL+(RELAX_RETRY_P*DDUVWP(I))**2
                END IF
              END DO
              DERR1_TRIAL=DSQRT(DERR1_TRIAL)
              DERR1_CAP_SD=0.005D0
              IF(DERR1_TRIAL.GT.DERR1_CAP_SD.AND.DERR1_TRIAL.GT.0.0D0) THEN
                DRELAX_NORM_SCALE=DERR1_CAP_SD/DERR1_TRIAL
                RELAX_RETRY_U=RELAX_RETRY_U*DRELAX_NORM_SCALE
                RELAX_RETRY_VW=RELAX_RETRY_VW*DRELAX_NORM_SCALE
                RELAX_RETRY_P=RELAX_RETRY_P*DRELAX_NORM_SCALE
              END IF
              ILINESEARCH_TRIES=0
              IOPPOSITE_TRIED=1
              WRITE(*,*)'LINESEARCH TRY STEEPEST DESCENT DIRECTION'
              WRITE(*,*)'STEEPEST_RELAX U/VW/P=',RELAX_RETRY_U,RELAX_RETRY_VW,RELAX_RETRY_P
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
              GOTO 30
            END IF
            WRITE(*,*)'LINESEARCH FAILED TO DECREASE WW4. STOP NEWTON ITERATION.'
            WRITE(*,*)'WW4_PREV/WW4=',WW4_PREV,WW4
            STOP
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
      iparm(8)=10
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

      ROW_SCALE=1.0D0
      COL_SCALE=1.0D0
      DO I=1,N
        DROWSUM=0.0D0
        DO J=ir(I),ir(I+1)-1
          DROWSUM=DMAX1(DROWSUM,DABS(a(J)))
        END DO
        IF(DROWSUM.GT.0.0D0) ROW_SCALE(I)=1.0D0/DROWSUM
      END DO
      DO I=1,N
        DO J=ir(I),ir(I+1)-1
          a(J)=a(J)*ROW_SCALE(I)
        END DO
        FFS(I)=FFS(I)*ROW_SCALE(I)
      END DO
      DO I=1,N
        COL_SCALE(I)=0.0D0
      END DO
      DO I=1,N
        DO J=ir(I),ir(I+1)-1
          COL_SCALE(ic(J))=DMAX1(COL_SCALE(ic(J)),DABS(a(J)))
        END DO
      END DO
      DO I=1,N
        IF(COL_SCALE(I).GT.0.0D0) THEN
          COL_SCALE(I)=1.0D0/COL_SCALE(I)
        ELSE
          COL_SCALE(I)=1.0D0
        END IF
      END DO
      DO I=1,N
        DO J=ir(I),ir(I+1)-1
          a(J)=a(J)*COL_SCALE(ic(J))
        END DO
      END DO

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
        DRESNORM=DRESNORM+(DAXI-FF1(I)*ROW_SCALE(I))**2
        DRHSNORM=DRHSNORM+(FF1(I)*ROW_SCALE(I))**2
      END DO
      DRESNORM=DSQRT(DRESNORM)
      DRHSNORM=DSQRT(DRHSNORM)
      IF(DRHSNORM.GT.0.0D0) THEN
        DLINEAR_REL=DRESNORM/DRHSNORM
        WRITE(*,*)'LINEAR_REL_RESIDUAL CSR=',DLINEAR_REL
      ELSE
        DLINEAR_REL=0.0D0
        WRITE(*,*)'LINEAR_ABS_RESIDUAL CSR=',DRESNORM
      END IF
      DDUVWP_SD=0.0D0
      DO I=1,N
        DFSCALED=FF1(I)*ROW_SCALE(I)
        DO J=ir(I),ir(I+1)-1
          DDUVWP_SD(ic(J))=DDUVWP_SD(ic(J))+COL_SCALE(ic(J))*a(J)*DFSCALED
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
      DDIR_NEWTON=0.0D0
      DDIR_STEEPEST=0.0D0
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
      END DO
      IUSE_STEEPEST=0
      IF(DDIR_NEWTON.LE.0.0D0) THEN
        DO I=1,4*IAT*NX
          DDUVWP(I)=DDUVWP_SD(I)
        END DO
        IUSE_STEEPEST=1
      END IF
      WRITE(*,*)'DIRECTION_CHECK NEWTON/SD/LINEAR_REL/USE_SD=',DDIR_NEWTON,DDIR_STEEPEST,DLINEAR_REL,IUSE_STEEPEST
      
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
      
      RELAX_U=0.05D0
      RELAX_VW=0.05D0
      RELAX_P=0.02D0
      DUPDATE_CAP_U=0.02D0
      DUPDATE_CAP_VW=0.02D0
      DUPDATE_CAP_P=0.20D0
      DRAW_MAX_U=0.0D0
      DRAW_MAX_V=0.0D0
      DRAW_MAX_W=0.0D0
      DRAW_MAX_P=0.0D0
      DO I=1,4*IAT*NX
        IF(I.LE.IAT*NX) THEN
          DRAW_MAX_U=DMAX1(DRAW_MAX_U,DABS(DDUVWP(I)))
        ELSE IF(I.LE.2*IAT*NX) THEN
          DRAW_MAX_V=DMAX1(DRAW_MAX_V,DABS(DDUVWP(I)))
        ELSE IF(I.LE.3*IAT*NX) THEN
          DRAW_MAX_W=DMAX1(DRAW_MAX_W,DABS(DDUVWP(I)))
        ELSE
          DRAW_MAX_P=DMAX1(DRAW_MAX_P,DABS(DDUVWP(I)))
        END IF
      END DO
      IF(DRAW_MAX_U*RELAX_U.GT.DUPDATE_CAP_U.AND.DRAW_MAX_U.GT.0.0D0) RELAX_U=DUPDATE_CAP_U/DRAW_MAX_U
      IF(DMAX1(DRAW_MAX_V,DRAW_MAX_W)*RELAX_VW.GT.DUPDATE_CAP_VW.AND.DMAX1(DRAW_MAX_V,DRAW_MAX_W).GT.0.0D0) THEN
        RELAX_VW=DUPDATE_CAP_VW/DMAX1(DRAW_MAX_V,DRAW_MAX_W)
      END IF
      IF(DRAW_MAX_P*RELAX_P.GT.DUPDATE_CAP_P.AND.DRAW_MAX_P.GT.0.0D0) RELAX_P=DUPDATE_CAP_P/DRAW_MAX_P
      DERR1_CAP=1.50D0
      IF(IUSE_STEEPEST.NE.0) DERR1_CAP=0.005D0
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
      IF(DERR1_TRIAL.GT.DERR1_CAP.AND.DERR1_TRIAL.GT.0.0D0) THEN
        DRELAX_NORM_SCALE=DERR1_CAP/DERR1_TRIAL
        RELAX_U=RELAX_U*DRELAX_NORM_SCALE
        RELAX_VW=RELAX_VW*DRELAX_NORM_SCALE
        RELAX_P=RELAX_P*DRELAX_NORM_SCALE
      END IF
      WRITE(*,*)'ADAPT_RELAX U/VW/P=',RELAX_U,RELAX_VW,RELAX_P
      WRITE(*,*)'RAW_MAX_UPDATE U/V/W/P=',DRAW_MAX_U,DRAW_MAX_V,DRAW_MAX_W,DRAW_MAX_P
      WRITE(*,*)'TRIAL_ERR1_CAP/TRIAL=',DERR1_CAP,DERR1_TRIAL
      RELAX_RETRY_U=RELAX_U
      RELAX_RETRY_VW=RELAX_VW
      RELAX_RETRY_P=RELAX_P
      WW4_PREV=WW4
      IPENDING_STEP=1
      ILINESEARCH_TRIES=0
      STEP_SIGN=UPDATE_SIGN
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
           
      
      ERR1=0.0D0
      DUMAX=0.0D0
      DVMAX=0.0D0
      DWMAX=0.0D0
      DPMAX=0.0D0
      IUMAX=0
      IVMAX=0
      IWMAX=0
      IPMAX=0
      DO I=1,4*IAT*NX
        IF(IS_CONSTRAINED(I).NE.0) THEN
          IF(I.LE.IAT*NX) THEN
            ERR1=ERR1+(RELAX_RETRY_U*DDUVWP(I))**2
          ELSE IF(I.LE.3*IAT*NX) THEN
            ERR1=ERR1+(RELAX_RETRY_VW*DDUVWP(I))**2
          ELSE
            ERR1=ERR1+(RELAX_RETRY_P*DDUVWP(I))**2
          END IF
        ELSE IF(I.GE.1.AND.I.LE.IAT*NX) THEN
          ERR1=ERR1+(RELAX_RETRY_U*DDUVWP(I))**2
        ELSE IF(I.GT.IAT*NX.AND.I.LE.3*IAT*NX) THEN
          ERR1=ERR1+(RELAX_RETRY_VW*DDUVWP(I))**2
        ELSE IF(I.GT.3*IAT*NX) THEN
          ERR1=ERR1+(RELAX_RETRY_P*DDUVWP(I))**2
        END IF
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
      ERR1=DSQRT(ERR1)
      WRITE(*,*)'ERR1_CURRENT=',ERR1
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
      
      IF(ERR1.GT.1.0D-12.OR.WW4.GT.1.0D-12) THEN
        IPK=IPK+1
        WRITE(*,*)'ERR,FX,IPK,IPK1,Re=',ERR1,WW4,IPK,IPK1,Re
        GO TO 30
      ELSE
      IPK1=IPK1+1
      WRITE(*,*)'ERR,FX,IPK,IPK1,Re=',ERR1,WW4,IPK,IPK1,Re
      END IF
      
        DO I=1,NX
        DO ITT=1,IAT
           WRITE(2,2100) T+DT*TTC(ITT),XX(I),YY(I),ZZ(I),U(ITT,I),V(ITT,I),W(ITT,I),P(ITT,I),Re
          WRITE(*,*) T+DT*TTC(ITT),XX(I),YY(I),ZZ(I),U(ITT,I),V(ITT,I),W(ITT,I),P(ITT,I),Re
        END DO
        END DO
        
        T=T+DT
      
      
      END DO
        
      WRITE(*,*)'SUCSESSFUL GOOD LUCK !'
      
      

        
 2100   FORMAT(8D25.16)
!           do i=1,nx*iat*4
!                     if(ssmmk(i,i)==0.0e0) then
!                 write(*,*) i,i/(nx),ssmmk(i,i)
!                   end if
!           end do
           
      CLOSE(2)

! FLOW CHECK ------------------------
      DO J=1,NAZ
      DO I=1,NAY
        IK=I+NAY*(J-1)
        WW(IK)=W1(I)*W2(J)
      END DO
      END DO

! SECTION 1 (INFLOW)----
      FLOW1=U0*D0*H0
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


