!     NUMERICAL ANALISYS OF UNSTEADY NAVIER_STOKES EQUATION
!         BY ORTHOGONAL COLLOCATION FINITE ELEMENT METHOD
!          
!
!         PROGRAM :  (Steady_Stokes_Curve_Duct_OCFEM_takaki_intel_LU_3D_MM.f90)
!
!              CODED BY T.Okubo from 2011.3.31.(Thurs)
!              modified by T.Ohkubo  2022.7.29.(Fri)
!                                    2022.8.23.(Tues)
!                                    2025.1.9.(Thurs)
!     ************* Means of variable and parameter *******************
!
!     ------------  ALLAY PARAMETER OF DIMENSION  ----------------------
!     KNX           :MAXIMUM TOTAL COLLOCATION NUMBER
!     KMX           :MAXIMUM TOTAL ELEMENTS NUMBER
!     KNA           :MAXIMUM COLLOCATION NUMBER IN ONE DIMENSIONAL DIRECTION
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
!     YU1(KNX) :VALUE OF VEROCITY_U IN TIME & SPACE (YU1M:MODIFIED)
!     YV1(KNX) :VALUE OF VEROCITY_V IN TIME & SPACE (YV1M:MODIFIED)
!     YW1(KNX) :VALUE OF VEROCITY_V IN TIME & SPACE (YW1M:MODIFIED)
!     YP1(KNX) :VALUE OF PRESSURE_P IN TIME & SPACE (YP1M:MODIFIED)
!     YVW1(KNX) :VALUE OF VORTICITY_VW IN TIME & SPACE
!     YWU1(KNX) :VALUE OF VORTICITY_WU IN TIME & SPACE
!     YUV1(KNX) :VALUE OF VORTICITY_UV IN TIME & SPACE 
!     -----------------------------------------------------------------
!
!    -----------------------------------------------------------------
 
      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),U(1:KNX),V(1:KNX),W(1:KNX),P(1:KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA)
      COMMON/PARAM1/ Re,U0,D0,H0,w8,X1,X2,ANYU,ROU
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),JD1(1:KNX),JD2(1:KNX),NSBU(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB)

      REAL(8),ALLOCATABLE :: AX(:,:),AY(:,:),AZ(:,:),BX(:,:),BY(:,:),BZ(:,:),CXY(:,:),CYZ(:,:),CZX(:,:)
      ALLOCATE (AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA))

      
      CALL INPUTDAT(AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)

      CALL BINITI

      CALL SOLUTE(AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)

      STOP
      END


!     *****************************************************************


      SUBROUTINE INPUTDAT(AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)

      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),U(1:KNX),V(1:KNX),W(1:KNX),P(1:KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA)
      COMMON/PARAM1/ Re,U0,D0,H0,w8,X1,X2,ANYU,ROU
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),JD1(1:KNX),JD2(1:KNX),NSBU(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB)

      REAL(8)  AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA)

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


      WRITE(*,*) 'ANYU,ROU, w8= ?'
      READ(*,*) ANYU,ROU,W8 
      
      WRITE(*,*)'Velocity U0,D0,H0,X1,X2 = ?'
      READ(*,*) U0,D0,H0,X1,X2
 
      Re=U0*D0/ANYU

      RETURN
      END

!     ************************************************************      
      
      SUBROUTINE BINITI

      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),U(1:KNX),V(1:KNX),W(1:KNX),P(1:KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA)
      COMMON/PARAM1/ Re,U0,D0,H0,w8,X1,X2,ANYU,ROU
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),JD1(1:KNX),JD2(1:KNX),NSBU(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB)
      CHARACTER FNAME*100
      
      WRITE(*,*) 'FILE NAME OF INITI RESULT DATA : ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)
        DO I=1,NX
          READ(2,*) TEMP,TEMP,TEMP,U(I),V(I),W(I),P(I),TEMP
        END DO
      CLOSE(2)
      write(*,*)'clear 001!'      

      RETURN
      END




!     *************************************************************

      SUBROUTINE TRANSF(JF,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
      
      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),U(1:KNX),V(1:KNX),W(1:KNX),P(1:KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA)
      COMMON/PARAM1/ Re,U0,D0,H0,w8,X1,X2,ANYU,ROU
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),JD1(1:KNX),JD2(1:KNX),NSBU(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB)

      REAL(8)  AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA)







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
!
        CALL invgj(DJJ,9)
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

      function omp_dotp(n,a,b) result(dp)  ! dp= a dot b
         integer,intent(in) :: n
         real(8),intent(in) :: a(n),b(n)
         real(8) dp
         integer i
         dp=0.0d0
         !$OMP parallel private(i)
           !$OMP do reduction(+ : dp)
           do i=1,n
             dp=dp+a(i)*b(i)
           enddo
           !$OMP enddo
         !$OMP end parallel
      end function omp_dotp
         


!     *****************************************************************

      SUBROUTINE SOLUTE(AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)

      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),U(1:KNX),V(1:KNX),W(1:KNX),P(1:KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA)
      COMMON/PARAM1/ Re,U0,D0,H0,w8,X1,X2,ANYU,ROU
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),JD1(1:KNX),JD2(1:KNX),NSBU(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB)
      
      REAL(8) AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA)


      dimension abcxy(9),dj1(9),dj2(9),dj3(9),dj4(9),dj5(9),dj6(9),dj7(9),dj8(9),dj9(9),taxuvp1(knaa),tayuvp1(knaa),tazuvp1(knaa),uu1(knaa),vv1(knaa),ww1(knaa),pp1(knaa),&
                FFUVP(4*knaa),tbxuvp1(knaa),tbyuvp1(knaa),tbzuvp1(knaa),YUVWP(4*KNX),DDUVWP(4*KNX),FF1(4*KNX)



      character FNAME*100,trans*2

      INTEGER :: IPIV(4*nx),iperm(nx*4),ir(nx*4+1),ic(10000*nx*4),iparm(64),ipt(64)
      REAL(8) :: FFS(4*nx),YUVP1(4*nx),x(4*nx,1),a(10000*nx*4)
      REAL(8),ALLOCATABLE :: BBUVP(:,:)
      REAL(8),ALLOCATABLE :: SSMMK(:,:)
      ALLOCATE (BBUVP(KNAA*4,KNAA*4),SSMMK(4*nx,4*nx))
      
 
      !$ call omp_set_num_threads(8)

!     ------------ CALCULATION OF INITIAL PROBREM BY O.C.F.E.M ------------

      WRITE(*,*)'OUTPUT FILE NAME OF RESULTS DATA : ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)

      PAI=3.14159265358979323846D0
 
    write(*,*) 'MX,NX',mx,nx
 

        
      INX=IAN

      
      
      
        DO IB=1,NSB2
          II1=NSBUI(IB)
          U(II1)=SBUI(IB)
          V(II1)=SBVI(IB)
          W(II1)=SBWI(IB)
        END DO
          
          
        IPK=0
 30     CONTINUE
 
 
        !$OMP parallel private(I,J)
        !$OMP do
        DO I=1,NX*4
        DO J=1,NX*4
          SSMMK(I,J)=0.0D0
        END DO
        END DO
        !$OMP end do
        !$OMP end parallel
 
        !$OMP parallel private(I)
        !$OMP do 
        DO I=1,NX*4
          FFS(I)=0.0D0
        END DO
        !$OMP end do
        !$OMP end parallel
        
        
        write(*,*)'clear 001!'

      DO IM=1,MX   
        DO I=1,4*IAN
        DO J=1,4*IAN
          BBUVP(I,J)=0.0D0
        END DO
        END DO
     
        
        CALL TRANSF(IM,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)


        DO I1=1,IAN
          II=IE(I1,IM)
          IIB=I1
          ISP=II
          DO IQ1=1,9
          DO IQ2=1,9
            DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,I1)
          END DO
          END DO

          DO J1=1,IAN
            JJ=IE(J1,IM)
            JJB=J1
            IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ).OR.(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.&
                I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
              BBUVP(IIB,JJB)=(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                  DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))*XL(I1,IM)+(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+&
                  DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))*XM(I1,IM)&
                  +(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                  DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))*XN(I1,IM)
              BBUVP(IIB+INX,JJB+INX)=(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                  DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))*XL(I1,IM)+(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+&
                  DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))*XM(I1,IM)&
                  +(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                  DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))*XN(I1,IM)
              BBUVP(IIB+2*INX,JJB+2*INX)=(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                  DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))*XL(I1,IM)+(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+&
                  DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))*XM(I1,IM)&
                  +(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                  DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))*XN(I1,IM)
              BBUVP(IIB+3*INX,JJB+3*INX)=(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                  DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))*XL(I1,IM)+(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+&
                  DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))*XM(I1,IM)&
                  +(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                  DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))*XN(I1,IM)
            END IF
            IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ)+NAX*NAY*NAZ)) THEN
                BBUVP(IIB,JJB)=DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1)
                BBUVP(IIB,JJB+INX)=DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+&
                     DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1)
                BBUVP(IIB,JJB+2*INX)=DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                     DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1)
                BBUVP(IIB,JJB+3*INX)=0.0D0
                
                !$OMP parallel private(IJ)
                !$OMP do
                DO IJ=1,IAN
                TAXUVP1(IJ)=DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ)
                TAYUVP1(IJ)=DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ)
                TAZUVP1(IJ)=DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ)
                UU1(IJ)=U(IE(IJ,IM))
                VV1(IJ)=V(IE(IJ,IM))
                WW1(IJ)=W(IE(IJ,IM))
                PP1(IJ)=P(IE(IJ,IM))
                END DO
                !$OMP end do
                !$OMP end parallel
                TAXU=omp_dotp ( IAN, TAXUVP1, UU1)
                TAYU=omp_dotp ( IAN, TAYUVP1, UU1)
                TAZU=omp_dotp ( IAN, TAZUVP1, UU1)
                TAXV=omp_dotp ( IAN, TAXUVP1, VV1)
                TAYV=omp_dotp ( IAN, TAYUVP1, VV1)
                TAZV=omp_dotp ( IAN, TAZUVP1, VV1)
                TAXW=omp_dotp ( IAN, TAXUVP1, WW1)
                TAYW=omp_dotp ( IAN, TAYUVP1, WW1)
                TAZW=omp_dotp ( IAN, TAZUVP1, WW1)
                
                IF(I1.EQ.J1) THEN
                BBUVP(IIB+INX,JJB)=U(II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+TAXU+V(II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))&
                     -ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+INX,JJB+INX)=TAYU
                BBUVP(IIB+INX,JJB+2*INX)=TAZU
                END IF
                IF(I1.NE.J1) THEN
                BBUVP(IIB+INX,JJB)=U(II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))&
                     -ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+INX,JJB+INX)=0.0D0
                BBUVP(IIB+INX,JJB+2*INX)=0.0D0
                END IF
                BBUVP(IIB+INX,JJB+3*INX)=1.0D0/ROU*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))
                IF(I1.EQ.J1) THEN
                BBUVP(IIB+2*INX,JJB)=TAXV
                BBUVP(IIB+2*INX,JJB+INX)=U(II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+TAYV+W(II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))&
                     -ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+2*INX,JJB+2*INX)=TAZV
                END IF
                IF(I1.NE.J1) THEN
                BBUVP(IIB+2*INX,JJB)=0.0D0
                BBUVP(IIB+2*INX,JJB+INX)=U(II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))-&
                     ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+2*INX,JJB+2*INX)=0.0D0
                END IF                
                BBUVP(IIB+2*INX,JJB+3*INX)=1.0D0/ROU*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+&
                     DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))
                IF(I1.EQ.J1) THEN
                BBUVP(IIB+3*INX,JJB)=TAXW
                BBUVP(IIB+3*INX,JJB+INX)=TAYW
                BBUVP(IIB+3*INX,JJB+2*INX)=U(II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))+TAZW-&
                     ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                END IF
                IF(I1.NE.J1) THEN
                BBUVP(IIB+3*INX,JJB)=0.0D0
                BBUVP(IIB+3*INX,JJB+INX)=0.0D0
                BBUVP(IIB+3*INX,JJB+2*INX)=U(II)*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))+V(II)*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+&
                     DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))+W(II)*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+&
                     DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+ DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))-&
                     ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                END IF
                BBUVP(IIB+3*INX,JJB+3*INX)=1.0D0/ROU*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                     DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))
               
              END IF
              
              
          END DO
          
          
          IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ).OR.(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ&
              .AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
          !$OMP parallel private(IJ)
          !$OMP do
          DO IJ=1,IAN
            TAXUVP1(IJ)=DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ)
            TAYUVP1(IJ)=DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ)
            TAZUVP1(IJ)=DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ)
            UU1(IJ)=U(IE(IJ,IM))
            VV1(IJ)=V(IE(IJ,IM))
            WW1(IJ)=W(IE(IJ,IM))
            PP1(IJ)=P(IE(IJ,IM))
          END DO
          !$OMP end do
          !$OMP end parallel

          TAXU=omp_dotp ( IAN, TAXUVP1, UU1)
          TAYU=omp_dotp ( IAN, TAYUVP1, UU1)
          TAZU=omp_dotp ( IAN, TAZUVP1, UU1)
          TAXV=omp_dotp ( IAN, TAXUVP1, VV1)
          TAYV=omp_dotp ( IAN, TAYUVP1, VV1)
          TAZV=omp_dotp ( IAN, TAZUVP1, VV1)
          TAXW=omp_dotp ( IAN, TAXUVP1, WW1)
          TAYW=omp_dotp ( IAN, TAYUVP1, WW1)
          TAZW=omp_dotp ( IAN, TAZUVP1, WW1)
          TAXP=omp_dotp ( IAN, TAXUVP1, PP1)
          TAYP=omp_dotp ( IAN, TAYUVP1, PP1)
          TAZP=omp_dotp ( IAN, TAZUVP1, PP1)
          FFUVP(I1)=TAXU*XL(I1,IM)+TAYU*XM(I1,IM)+TAZU*XN(I1,IM)
          FFUVP(I1+IAN)=TAXV*XL(I1,IM)+TAYV*XM(I1,IM)+TAZV*XN(I1,IM)
          FFUVP(I1+IAN*2)=TAXW*XL(I1,IM)+TAYW*XM(I1,IM)+TAZW*XN(I1,IM)
          FFUVP(I1+IAN*3)=TAXP*XL(I1,IM)+TAYP*XM(I1,IM)+TAZP*XN(I1,IM)
        END IF
        
 !       write(*,*)'clear 005!'

        IF(I1.GE.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ)+NAX*NAY*NAZ) THEN
          !$OMP parallel private(IJ)
          !$OMP do
          DO IJ=1,IAN
            TAXUVP1(IJ)=DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ)
            TAYUVP1(IJ)=DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ)
            TAZUVP1(IJ)=DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ)
            TBXUVP1(IJ)=DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)
            TBYUVP1(IJ)=DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)
            TBZUVP1(IJ)=DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ)
            UU1(IJ)=U(IE(IJ,IM))
            VV1(IJ)=V(IE(IJ,IM))
            WW1(IJ)=W(IE(IJ,IM))
            PP1(IJ)=P(IE(IJ,IM))
          END DO
          !$OMP end do
          !$OMP end parallel
          TAXU=omp_dotp ( IAN, TAXUVP1, UU1)
          TAYU=omp_dotp ( IAN, TAYUVP1, UU1)
          TAZU=omp_dotp ( IAN, TAZUVP1, UU1)
          TAXV=omp_dotp ( IAN, TAXUVP1, VV1)
          TAYV=omp_dotp ( IAN, TAYUVP1, VV1)
          TAZV=omp_dotp ( IAN, TAZUVP1, VV1)
          TAXW=omp_dotp ( IAN, TAXUVP1, WW1)
          TAYW=omp_dotp ( IAN, TAYUVP1, WW1)
          TAZW=omp_dotp ( IAN, TAZUVP1, WW1)
          TBXU=omp_dotp ( IAN, TBXUVP1, UU1)
          TBYU=omp_dotp ( IAN, TBYUVP1, UU1)
          TBZU=omp_dotp ( IAN, TBZUVP1, UU1)
          TBXV=omp_dotp ( IAN, TBXUVP1, VV1)
          TBYV=omp_dotp ( IAN, TBYUVP1, VV1)
          TBZV=omp_dotp ( IAN, TBZUVP1, VV1)
          TBXW=omp_dotp ( IAN, TBXUVP1, WW1)
          TBYW=omp_dotp ( IAN, TBYUVP1, WW1)
          TBZW=omp_dotp ( IAN, TBZUVP1, WW1)
          TAXP=omp_dotp ( IAN, TAXUVP1, PP1)
          TAYP=omp_dotp ( IAN, TAYUVP1, PP1)
          TAZP=omp_dotp ( IAN, TAZUVP1, PP1)
            FFUVP(I1)=TAXU+TAYV+TAZW
            FFUVP(I1+IAN)=U(II)*TAXU+V(II)*TAYU+W(II)*TAZU-ANYU*(TBYU+TBXU+TBZU)+1.0D0/ROU*TAXP
            FFUVP(I1+IAN*2)=U(II)*TAXV+V(II)*TAYV+W(II)*TAZV-ANYU*(TBXV+TBYV+TBZV)+1.0D0/ROU*TAYP
            FFUVP(I1+IAN*3)=U(II)*TAXW+V(II)*TAYW+W(II)*TAZW-ANYU*(TBXW+TBYW+TBZW)+1.0D0/ROU*TAZP
!            WRITE(*,*)'FFUVP=',FFUVP(I1),FFUVP(I1+IAN),FFUVP(I1+2*IAN),FFUVP(I1+3*IAN)
            END IF
              
        END DO
!       ----------------- SUPERPOSITION OF BBUVP- AND TIME -------------

        DO JA1=1,IAN
            IF((JA1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.JA1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ).OR.&
               (JA1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.JA1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ).OR.&
               (JA1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ.AND.JA1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ+NAX*NAY*NAZ)) THEN
        DO IK1=1,4
            II1=IE(JA1,IM)
            ISP1=II1+NX*(IK1-1)
            IPL1=JA1+IAN*(IK1-1)
              DO JA2=1,IAN
              DO IK2=1,4
                JJ2=IE(JA2,IM)
                ISP2=JJ2+NX*(IK2-1)
                IPL2=JA2+IAN*(IK2-1)
                SSMMK(ISP1,ISP2)=SSMMK(ISP1,ISP2)+BBUVP(IPL1,IPL2)
              END DO             
              END DO 
            FFS(ISP1)=FFS(ISP1)+FFUVP(IPL1)
        END DO
            END IF
        END DO
      END DO
!
              write(*,*)'clear ohkubo 1TC',TC

!     ---------------- DIFFERENTIAL OF SLOPE .EQ. 0,  U,V,W  at B3 OUT FLOW --------------------
!     SLOPE U ---------
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II
          !$OMP parallel private(ISP2)
          !$OMP do
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          !$OMP end do
          !$OMP end parallel
          FFS(ISP1)=0.0D0
        END IF
      END DO
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II
        
        CALL TRANSF(5,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
          !$OMP parallel private(IJ,ISP2)
          !$OMP do
          DO IJ=1,IAN
            ISP2=IE(IJ,5)
            SSMMK(ISP1,ISP2)=(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))
          END DO
          !$OMP end do
          !$OMP end parallel
          
          TBYU=0.0D0
          !$OMP parallel private(IJ,JJ)
          !$OMP do reduction(+ : TBYU)
          DO IJ=1,IAN
          JJ=IE(IJ,5)
          TBYU=TBYU+(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))*U(JJ)
          END DO
          !$OMP end do
          !$OMP end parallel
          FFS(ISP1)=TBYU
        END IF      
      END DO

!     SLOPE V -----------------
 
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II+NX
          !$OMP parallel private(ISP2)
          !$OMP do
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          !$OMP end do
          !$OMP end parallel
          FFS(ISP1)=0.0D0
        END IF
      END DO

      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II+NX
        
        CALL TRANSF(5,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
          !$OMP parallel private(IJ,ISP2)
          !$OMP do
          DO IJ=1,IAN
            ISP2=IE(IJ,5)+NX
            SSMMK(ISP1,ISP2)=(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))
          END DO
          !$OMP end do
          !$OMP end parallel
          
          TBYV=0.0D0
          !$OMP parallel private(IJ,JJ)
          !$OMP do reduction(+ : TBYV)
          DO IJ=1,IAN
          JJ=IE(IJ,5)
          TBYV=TBYV+(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))*V(JJ)
          END DO
          !$OMP end do
          !$OMP end parallel
          FFS(ISP1)=TBYV
        END IF      
      END DO

!     SLOPE W -----------------

      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II+2*NX
          !$OMP parallel private(ISP2)
          !$OMP do
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          !$OMP end do
          !$OMP end parallel
          FFS(ISP1)=0.0D0
        END IF
      END DO

      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II+2*NX
        
        CALL TRANSF(5,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
          !$OMP parallel private(IJ,ISP2)
          !$OMP do
          DO IJ=1,IAN
            ISP2=IE(IJ,5)+2*NX
            SSMMK(ISP1,ISP2)=(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))
          END DO
          !$OMP end do
          !$OMP end parallel
          TBYW=0.0D0
          !$OMP parallel private(IJ,JJ)
          !$OMP do reduction(+ : TBYW)
          DO IJ=1,IAN
          JJ=IE(IJ,5)
          TBYW=TBYW+(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))*W(JJ)
          END DO
          !$OMP end do
          !$OMP end parallel
          FFS(ISP1)=TBYW
        END IF      
      END DO
      
      write(*,*)'clear 099!'

!     ------------ DIFFERENTIAL OF PRESSURE SLOPE AT INFLOW ---------------- 

!     ------- SLOPE P AT INFLOW
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ)) THEN
        II=IE(J1,1)
        ISP1=II+3*NX
          !$OMP parallel private(ISP2)
          !$OMP do
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          !$OMP end do
          !$OMP end parallel
          FFS(ISP1)=0.0D0
        END IF
      END DO

      J=1
      CALL TRANSF(J,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
      
      DO I1=1,IAN
        IF((I1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ)) THEN
                
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,I1)
        END DO
        END DO        
        
          JJ=IE(I1,J)
          ISP1=JJ+3*NX
          !$OMP parallel private(IJ,ISP2)
          !$OMP do
          DO IJ=1,IAN
            ISP2=IE(IJ,J)+3*NX
            SSMMK(ISP1,ISP2)=(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ))
          END DO
          !$OMP end do
          !$OMP end parallel
          TBXP=0.0D0
          !$OMP parallel private(IJ,JJ)
          !$OMP do reduction(+ : TBXP)
          DO IJ=1,IAN
          JJ=IE(IJ,J)
          TBXP=TBXP+(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ))*P(JJ)
          END DO
          !$OMP end do
          !$OMP end parallel 
          FFS(ISP1)=TBXP 
        END IF     
      END DO
      
      write(*,*)'clear 098!'

           
!     --------------  BOUNDARY CONDITION AT WALL --------------

      DO J=1,NSB1
          JJ=NSBU(J)
          ISP2=JJ
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*0.0D0
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=0.0D0
      END DO
      
      DO J=1,NSB1
        JJ=NSBU(J)
          ISP2=JJ+NX
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*0.0D0
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=0.0D0
      END DO

      DO J=1,NSB1
        JJ=NSBU(J)
          ISP2=JJ+2*NX
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*0.0D0
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=0.0D0
      END DO
      
      
! -------- PRESSURE AT OUTOFLOW CONDITION 
    J1=5
     DO I3=1,NSB3      
       II=NSBUO(I3)
         ISP2=II+3*NX
         DO ISP1=1,NX*4
           FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*0.0D0
           SSMMK(ISP1,ISP2)=0.0D0
           SSMMK(ISP2,ISP1)=0.0D0
         END DO
         SSMMK(ISP2,ISP2)=1.0D0
         FFS(ISP2)=0.0D0
     END DO
       

! -------- PRESSURE P BOUNDARY CONDITION AT WALL
      DO J1=1,5
        DO I1=1,IAN
          IF((I1.GE.9+4*(NAX+NAY+NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY).OR.&
            (I1.GE.1.AND.I1.LE.8).OR.(I1.GE.9.AND.I1.LE.8+4*(NAX+NAY+NAZ)).OR.&
            (I1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ).OR.&
            (I1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ)) THEN
          JJ=IE(I1,J1)
          ISP1=JJ+3*NX
          !$OMP parallel private(ISP2)
          !$OMP do
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          !$OMP end do
          !$OMP end parallel 
          FFS(ISP1)=0.0D0
          END IF
        END DO
      END DO

      DO J1=1,5
        
        CALL TRANSF(J1,AX,AY,AZ,BX,BY,BZ,CXY,CYZ,CZX)
        
        DO I1=1,IAN
        IF((I1.GE.9+4*(NAX+NAY+NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY).OR.&
            (I1.GE.1.AND.I1.LE.8).OR.(I1.GE.9.AND.I1.LE.8+4*(NAX+NAY+NAZ)).OR.&
            (I1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ).OR.&
            (I1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ)) THEN
        
          DO IQ1=1,9
          DO IQ2=1,9
            DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,I1)
          END DO
          END DO        
        
          JJ=IE(I1,J1)
          ISP1=JJ+3*NX
          DO IJ=1,IAN
            ISP2=IE(IJ,J1)+3*NX
            SSMMK(ISP1,ISP2)=SSMMK(ISP1,ISP2)+(DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+&
                 DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ))*XL(I1,J1)+(DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+&
                 DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ))*XM(I1,J1)+&
                 (DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+&
                 DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ))*XN(I1,J1)
          END DO

          TAXYZP=0.0D0
          !$OMP parallel private(IJ,JJ4)
          !$OMP do reduction(+ : TAXYZP)
          DO IJ=1,IAN
            JJ4=IE(IJ,J1)
            TAXYZP=TAXYZP+((DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+&
                 DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ))*XL(I1,J1)+(DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+&
                 DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ))*XM(I1,J1)+&
                 (DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+&
                 DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ))*XN(I1,J1))*P(JJ4)
          END DO
          !$OMP end do
          !$OMP end parallel
          FFS(ISP1)=TAXYZP
          END IF      
        END DO
      END DO

      write(*,*)'clear 097!'


!     ------------ INFLOW BOUNDARY CONDITION ------------------   
      DO J=1,NSB2
        JJ=NSBUI(J)
          ISP2=JJ
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*0.0D0
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=0.0D0
      END DO

      DO J=1,NSB2
        JJ=NSBUI(J)
          ISP2=JJ+NX
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*0.0D0
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=0.0D0
      END DO
      
      DO J=1,NSB2
        JJ=NSBUI(J)
          ISP2=JJ+2*NX
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*0.0D0
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=0.0D0
      END DO

      WRITE(*,*)'CLEAR 103 !'

      !$OMP parallel private(I)
      !$OMP do
      DO I=1,4*NX
        FF1(I)=FFS(I)
      END DO
      !$OMP end do
      !$OMP end parallel
      
      WW4=0.0D0
      !$OMP parallel private(I)
      !$OMP do reduction(+ : WW4)
      DO I=1,4*NX
        WW4=WW4+FF1(I)**2
      END DO
      !$OMP end do
      !$OMP end parallel
      WW4=DSQRT(WW4)
      
      n=nx*4
      
      call csr_S_matrix(n,n,ssmmk,mnr,ir,mn,ic,a) 
      
      
      write(*,*)'clear 99 !mnr,mn',mnr,mn 
        
      
      do i=1,64
       ipt(i)=0
       iparm(i)=0
      end do
      
      maxfct=1
      mnum=1
      mtype=11
      iphase=13 
      nrhs=1    
      iparm(3)=8 
      iparm(6)=1
      msglvl=0
      
      call pardiso (ipt, maxfct, mnum, mtype, iphase, n, a, ir, ic, iperm, nrhs, iparm, msglvl, FFS, x, ierror)

      write(*,*)'ierror=',ierror




!     -------- SOLUTION OF LINEAR SYSTEM ------------
!      trans='N'
!      n=nx*4
!      nrhs=1
!      lda=200000
!      ldb=200000
!        write(*,*)' Clear ohkubo 110!'   
!      CALL dgetrf(n,n,ssmmk,lda,ipiv,info)
!      write(*,*) '----------info=',info

!      CALL dgetrs(trans,n,nrhs,ssmmk,lda,ipiv,ffs,ldb,info)
!      write(*,*) '----------info=',info


      DO I=1,4*NX
        DDUVWP(I)=X(I,1)
      END DO
      
      DO IK=1,4
      DO I=1,NX
        ISP=I+(IK-1)*NX
        IF(ISP.GE.1.AND.ISP.LE.NX) THEN
          YUVWP(ISP)=U(I)
        END IF
        IF(ISP.GT.NX.AND.ISP.LE.2*NX) THEN
          YUVWP(ISP)=V(I)
        END IF
        IF(ISP.GT.2*NX.AND.ISP.LE.3+NX) THEN
          YUVWP(ISP)=W(I)
        END IF
        IF(ISP.GT.3*NX.AND.ISP.LE.4*NX) THEN
          YUVWP(ISP)=P(I)
        END IF
      END DO
      END DO
      
      DO I=1,4*NX
        IF(I.GE.1.AND.I.LE.3*NX) THEN
           YUVWP(I)=YUVWP(I)-0.7D0*DDUVWP(I)
        END IF
        IF(I.GT.3*NX) THEN
           YUVWP(I)=YUVWP(I)-0.3D0*DDUVWP(I)
        END IF
      END DO
           
      
      ERR1=0.0D0
      DO I=1,4*NX
        IF(I.GE.1.AND.I.LE.3*NX) THEN
          ERR1=ERR1+(0.7D0*DDUVWP(I))**2
        END IF
        IF(I.GT.3*NX) THEN
          ERR1=ERR1+(0.3D0*DDUVWP(I))**2
        END IF
      END DO
      ERR1=DSQRT(ERR1)
      
      DO IK=1,4  
      DO I=1,NX
        ISP=I+NX*(IK-1)
        IF(ISP.LE.NX) THEN
          U(I)=YUVWP(ISP)
        END IF
        IF((ISP.GT.NX).AND.(ISP.LE.2*NX)) THEN
          V(I)=YUVWP(ISP)
        END IF
        IF((ISP.GT.2*NX).AND.(ISP.LE.3*NX)) THEN
          W(I)=YUVWP(ISP)
        END IF
        IF(ISP.GT.3*NX) THEN
          P(I)=YUVWP(ISP)
        END IF
      END DO
      END DO
      
      IF(ERR1.GT.1.0D-12.OR.WW4.GT.1.0D-12) THEN
        IPK=IPK+1
        WRITE(*,*)'ERR,FX,IPK=',ERR1,WW4,IPK
        GO TO 30
      END IF
      WRITE(*,*)'ERR,FX,IPK=',ERR1,WW4,IPK
      
      
        DO I=1,NX
           WRITE(2,2100) XX(I),YY(I),ZZ(I),U(I),V(I),W(I),P(I),Re
          WRITE(*,*) XX(I),YY(I),ZZ(I),U(I),V(I),W(I),P(I),Re
        END DO
        
 2100   FORMAT(8D25.16)
!           do i=1,nx*4
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
! SECTION 2 (Straight Section)
      IM=1
      FLOW2=0.0D0
      DO I1=1,IAN
         IF(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ) THEN
         II=IE(I1,IM)
         V1=U(II)
         FLOW2=FLOW2+V1*WW(I1-8-4*(NAX+NAY+NAZ)-2*NAX*NAY-NAX*NAZ)
         END IF
      END DO
      FLOW2=FLOW2*D0*H0
      WRITE(*,*) FLOW2
! SECTION 3 (Curve Section)
      IM=3
      FLOW3=0.0D0
      DO I1=1,IAN
         IF(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ) THEN
         II=IE(I1,IM)
         V1=U(II)*DCOS(PAI/4.0D0)-V(II)*DSIN(PAI/4.0D0)
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
         V1=-V(II)
         FLOW4=FLOW4+V1*WW(I1-8-4*(NAX+NAY+NAZ)-2*NAX*NAY-NAX*NAZ)
         END IF
      END DO
      FLOW4=FLOW4*D0*H0
      WRITE(*,*) FLOW4
      
      WRITE(*,*)'OUTPUT FILE NAME OF SECTION FLOW : ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)
        WRITE(2,*)'FLOW RATE AT CHANNEL SECTION --------'
        WRITE(2,*)'FLOW1(Inlet)           =',FLOW1
        WRITE(2,*)'FLOW2(Straight Section)=',FLOW2
        WRITE(2,*)'FLOW3(Curve Section)   =',FLOW3
        WRITE(2,*)'FLOW4(Outlet)          =',FLOW4
     CLOSE(2)
      
!
      RETURN
      END
!
!
       subroutine csr_S_matrix(m,n,aa,mnr,ir,mn,ic,a)
       
       real(8) aa(n,n),a(10000*n)
       integer ia(10000*n),ir(n+1),ic(10000*n)
       
       ii=0
       do i=1,m
       do j=1,n
         if(aa(i,j).ne.0.0d0) then
           ii=ii+1
           a(ii)=aa(i,j)
           ia(ii)=i             
           ic(ii)=j
         end if
      end do
      end do
      mn=ii
      
      ik=1
      ikk=1
      ik1=1
      ir(1)=1
      do ii=1,mn
        if(ia(ii+1)-ia(ii).eq.0) then
          ik1=ik1+1
        else
          ikk=ikk+1
          ik=ik+ik1
          ir(ikk)=ik
          ik1=1
        end if
      end do
      
      mnr=ikk
      
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


