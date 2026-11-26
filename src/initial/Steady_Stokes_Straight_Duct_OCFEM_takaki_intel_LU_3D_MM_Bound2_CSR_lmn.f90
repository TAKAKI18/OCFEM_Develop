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
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),YU1(1:KNX),YV1(1:KNX),YW1(1:KNX),YP1(1:KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA)
      COMMON/PARAM1/ Re,U0,D0,H0,ANYU,ROU
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),IT(1:200000),JD1(KNX),JD2(KNX),NSBU(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB)





      CALL INPUTDAT


      CALL SOLUTE

      STOP
      END


!     *****************************************************************


      SUBROUTINE INPUTDAT

      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),YU1(1:KNX),YV1(1:KNX),YW1(1:KNX),YP1(1:KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA)
      COMMON/PARAM1/ Re,U0,D0,H0,ANYU,ROU
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),IT(1:200000),JD1(KNX),JD2(KNX),NSBU(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB)

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

    WRITE(*,*)'OUTPUT FILE NAME OF BOUNARY CONNECTION (boucone_c_.dat)'
    READ(*,*) FNAME
    OPEN(2,FILE=FNAME)
      READ(2,2000) ID
      DO I1=1,ID
        READ(2,2200) JD1(I1),JD2(I1),XL(JD1(I1),JD2(I1)),XM(JD1(I1),JD2(I1)),XN(JD1(I1),JD2(I1))
      END DO
    CLOSE(2)
2000 FORMAT(I10)
2200 FORMAT(2I10,3D25.16)


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


      WRITE(*,*) 'ANYU,ROU = ?'
      READ(*,*) ANYU,ROU 
      
      WRITE(*,*)'Dimensioless Velocity U0,D0,H0 = ?'
      READ(*,*) U0,D0,H0
 
      Re=U0*D0/ANYU

      RETURN
      END



!     *************************************************************

      SUBROUTINE TRANSF(JF)

      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),YU1(1:KNX),YV1(1:KNX),YW1(1:KNX),YP1(1:KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA)
      COMMON/PARAM1/ Re,U0,D0,H0,ANYU,ROU
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),IT(1:200000),JD1(KNX),JD2(KNX),NSBU(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB)








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
          IF (JF.EQ.5 .AND. IIu.EQ.565) THEN
            OPEN(88,FILE='fortran_debug_sum_5_565.log',STATUS='UNKNOWN',POSITION='APPEND')
            IF (JJu.EQ.1) THEN
              WRITE(88,'(A)') '# JF=5 IIu=565 JJu Ju AX*XX AX*YY AX*ZZ BX*XX BX*YY BX*ZZ CXY*XX CXY*YY CXY*ZZ CYZ*XX CYZ*YY CYZ*ZZ CZX*XX CZX*YY CZX*ZZ'
            END IF
            WRITE(88,'(2I6,1X,2I6,1P,18E20.12)') JF,IIu,JJu,Ju, &
                 AX(IIu,JJu)*XX(Ju), AX(IIu,JJu)*YY(Ju), AX(IIu,JJu)*ZZ(Ju), &
                 BX(IIu,JJu)*XX(Ju), BX(IIu,JJu)*YY(Ju), BX(IIu,JJu)*ZZ(Ju), &
                 CXY(IIu,JJu)*XX(Ju), CXY(IIu,JJu)*YY(Ju), CXY(IIu,JJu)*ZZ(Ju), &
                 CYZ(IIu,JJu)*XX(Ju), CYZ(IIu,JJu)*YY(Ju), CYZ(IIu,JJu)*ZZ(Ju), &
                 CZX(IIu,JJu)*XX(Ju), CZX(IIu,JJu)*YY(Ju), CZX(IIu,JJu)*ZZ(Ju)
            CLOSE(88)
          END IF
          IF (JF.EQ.5 .AND. IIu.EQ.565 .AND. JJu.EQ.1798) THEN
            OPEN(99,FILE='fortran_debug_5_565_1798.log',STATUS='UNKNOWN',POSITION='APPEND')
            WRITE(99,'(A,3I10,1P,12E20.12)') 'RAW: JF IIu JJu AX..CZX XX..ZZ=',JF,IIu,JJu, &
                 AX(IIu,JJu),AY(IIu,JJu),AZ(IIu,JJu),BX(IIu,JJu),BY(IIu,JJu),BZ(IIu,JJu), &
                 CXY(IIu,JJu),CYZ(IIu,JJu),CZX(IIu,JJu),XX(Ju),YY(Ju),ZZ(Ju)
            CLOSE(99)
          END IF
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

        IF (JF.EQ.5 .AND. IIu.EQ.565) THEN
          OPEN(99,FILE='fortran_debug_5_565_1798.log',STATUS='UNKNOWN',POSITION='APPEND')
          WRITE(99,'(A,1P,9E20.12)') 'TRANSF_DERIV1:', DXDG, DYDG, DZDG, DXDE, DYDE, DZDE, DXDT, DYDT, DZDT
          WRITE(99,'(A,1P,9E20.12)') 'TRANSF_DERIV2:', D2XDG2, D2YDG2, D2ZDG2, DXDG**2, DYDG**2, DZDG**2, 2*DXDG*DYDG, 2*DYDG*DZDG, 2*DZDG*DXDG
          CLOSE(99)
        END IF
          
          
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
        IF (IIu == 1 .AND. JF == 1) THEN
          WRITE(*,*) 'DEBUG: Fortran DJJ(1,1)=', DJJ(1,1), ' DJJ(1,4)=', DJJ(1,4)
        END IF
!
        DO iu1=1,9
        DO ju1=1,9
          DJJA(iu1,ju1,IIu)=DJJ(iu1,ju1)
          IF(JF==4) THEN
            write(*,*)'clear ohkubo Djj',Jf,Iu1,Ju1,iiu,djja(iu1,ju1,iiu)
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

      SUBROUTINE SOLUTE

      implicit real*8(a-h,o-z)
      parameter(KNX=50000,KMX=5,KNA=40,KNAA=10000,KNB=20000)
      COMMON/COLMAT1/ XXC(1:KNA),YYC(1:KNA),ZZC(1:KNA),AX(1:KNAA,1:KNAA),AY(1:KNAA,1:KNAA),AZ(1:KNAA,1:KNAA),BX(1:KNAA,1:KNAA),BY(1:KNAA,1:KNAA),&
                 BZ(1:KNAA,1:KNAA),CXY(1:KNAA,1:KNAA),CYZ(1:KNAA,1:KNAA),CZX(1:KNAA,1:KNAA),&
                 XX(1:KNX),YY(1:KNX),ZZ(1:KNX),YU1(1:KNX),YV1(1:KNX),YW1(1:KNX),YP1(1:KNX),W1(1:KNA),W2(1:KNA),WW(1:KNA*KNA)
      COMMON/BOUND1/ SBU(1:KNB),SBV(1:KNB),SBW(1:KNB),SBUI(1:KNB),SBVI(1:KNB),SBWI(1:KNB),SBPO(1:KNB),XL(KNAA,KMX),XM(KNAA,KMX),XN(KNAA,KMX)
      COMMON/TRANS/ DJJ1(1:9,1:9),DJJ2(1:9,1:9),DJJ(1:9,1:9),DJJA(1:9,1:9,1:KNX),GX(1:KNAA),GY(1:KNAA),GZ(1:KNAA)
      COMMON/PARAM1/ Re,U0,D0,H0,ANYU,ROU
      COMMON/PARAM2/ NSB1,NSB2,NSB3,NAX,NAY,NAZ,NNAX,NNAY,NNAZ,IAN,IAX,IAY,IAZ,NX,MX,ID,INX
      COMMON/COLLO1/ IE(1:KNAA,1:KMX),IT(1:200000),JD1(KNX),JD2(KNX),NSBU(1:KNB),JSBU(1:KNB),ISBU(1:KNB),NSBUI(1:KNB),JSBUI(1:KNB),ISBUI(1:KNB),&
                     NSBUO(1:KNB),JSBUO(1:KNB),ISBUO(1:KNB)






      character FNAME*100,trans*2

      INTEGER :: IIT(1:KNX*4),IPIV(4*nx),iperm(nx*4),ir(nx*4+1),ic(10000*nx*4),iparm(64)
      INTEGER(8) :: ipt(64)
      REAL(8) :: BBUVP(KNAA*4,KNAA*4)
      REAL(8) :: FFS(4*nx),YUVP1(4*nx),x(4*nx,1),a(10000*nx*4)
      REAL(8) :: SSMMK(4*nx,4*nx)


!     ------------ CALCULATION OF INITIAL PROBREM BY O.C.F.E.M ------------

      WRITE(*,*)'OUTPUT FILE NAME OF RESULTS DATA : ?'
      READ(*,*) FNAME
      OPEN(2,FILE=FNAME)

      PAI=3.14159265358979323846D0
 
    write(*,*) 'MX,NX',mx,nx
 
        DO I=1,NX
          YU1(I)=0.0D0
          YV1(I)=0.0D0
          YW1(I)=0.0D0
          YP1(I)=0.0D0
        END DO
        
      INX=IAN

      DO I=1,NX*4
      DO J=1,NX*4
        SSMMK(I,J)=0.0D0
      END DO
      END DO
      
      
        DO IB=1,NSB2
          II1=NSBUI(IB)
          YU1(II1)=SBUI(IB)
          YV1(II1)=SBVI(IB)
          YW1(II1)=SBWI(IB)
        END DO
          
        DO I=1,NX*4
          IT(I)=0
          IIT(I)=0
          FFS(I)=0.0D0
        END DO
        
        
 
      DO IM=1,MX   
        DO I=1,4*IAN
        DO J=1,4*IAN
          BBUVP(I,J)=0.0D0
        END DO
        END DO
     
        
        CALL TRANSF(IM)


        DO I1=1,IAN
          II=IE(I1,IM)
          IIB=I1
          ISP=II
          IT(ISP)=IT(ISP)+1
          DO IQ1=1,9
          DO IQ2=1,9
            DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,I1)
          END DO
          END DO

        DO J1=1,IAN
            JJ=IE(J1,IM)
            JJB=J1
            IF (IM.EQ.5 .AND. ((I1.EQ.565 .AND. J1.EQ.1798) .OR. (I1.EQ.2028 .AND. J1.EQ.587))) THEN
              OPEN(99,FILE='fortran_debug_5_565_1798.log',STATUS='UNKNOWN',POSITION='APPEND')
              WRITE(99,'(A,1P,9E20.12)') 'DJJ1:', DJJ(1,1),DJJ(1,2),DJJ(1,3),DJJ(1,4),DJJ(1,5),DJJ(1,6),DJJ(1,7),DJJ(1,8),DJJ(1,9)
              WRITE(99,'(A,1P,9E20.12)') 'DJJ2:', DJJ(2,1),DJJ(2,2),DJJ(2,3),DJJ(2,4),DJJ(2,5),DJJ(2,6),DJJ(2,7),DJJ(2,8),DJJ(2,9)
              WRITE(99,'(A,1P,9E20.12)') 'DJJ3:', DJJ(3,1),DJJ(3,2),DJJ(3,3),DJJ(3,4),DJJ(3,5),DJJ(3,6),DJJ(3,7),DJJ(3,8),DJJ(3,9)
              WRITE(99,'(A,1P,9E20.12)') 'AX..CZX:', AX(I1,J1),AY(I1,J1),AZ(I1,J1),BX(I1,J1),BY(I1,J1),BZ(I1,J1),CXY(I1,J1),CYZ(I1,J1),CZX(I1,J1)
              WRITE(99,'(A,1P,3E20.12)') 'XLXM XN:', XL(I1,IM),XM(I1,IM),XN(I1,IM)
              CLOSE(99)
            ENDIF
            IF (IM.EQ.5 .AND. I1.EQ.2028 .AND. J1.EQ.587) THEN
              OPEN(99,FILE='fortran_debug_5_2028_587.log',STATUS='UNKNOWN',POSITION='APPEND')
              WRITE(99,'(A,2I10,1P,9E20.12)') 'DJJ1 I1 J1:', I1, J1, DJJ(1,1),DJJ(1,2),DJJ(1,3),DJJ(1,4),DJJ(1,5),DJJ(1,6),DJJ(1,7),DJJ(1,8),DJJ(1,9)
              WRITE(99,'(A,2I10,1P,9E20.12)') 'DJJ2 I1 J1:', I1, J1, DJJ(2,1),DJJ(2,2),DJJ(2,3),DJJ(2,4),DJJ(2,5),DJJ(2,6),DJJ(2,7),DJJ(2,8),DJJ(2,9)
              WRITE(99,'(A,2I10,1P,9E20.12)') 'DJJ3 I1 J1:', I1, J1, DJJ(3,1),DJJ(3,2),DJJ(3,3),DJJ(3,4),DJJ(3,5),DJJ(3,6),DJJ(3,7),DJJ(3,8),DJJ(3,9)
              WRITE(99,'(A,2I10,1P,9E20.12)') 'AX..CZX I1 J1:', I1, J1, AX(I1,J1),AY(I1,J1),AZ(I1,J1),BX(I1,J1),BY(I1,J1),BZ(I1,J1),CXY(I1,J1),CYZ(I1,J1),CZX(I1,J1)
              WRITE(99,'(A,2I10,1P,3E20.12)') 'XLXM XN I1 J1:', I1, J1, XL(I1,IM),XM(I1,IM),XN(I1,IM)
              CLOSE(99)
            ENDIF
            IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ).OR.(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
              BBUVP(IIB,JJB)=(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                  DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))*XL(I1,IM)+(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+&
                  DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))*XM(I1,IM)&
                  +(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                  DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))*XN(I1,IM)
              IF (IM.EQ.5 .AND. ((I1.EQ.565 .AND. J1.EQ.1798) .OR. (I1.EQ.2028 .AND. J1.EQ.587))) THEN
                OPEN(99,FILE='fortran_debug_5_565_1798.log',STATUS='UNKNOWN',POSITION='APPEND')
                WRITE(99,'(A,3I10,1P,E20.12)') 'BBUVP(IIB,JJB)=',IM,I1,J1,BBUVP(IIB,JJB)
                WRITE(99,'(A,2I10)') 'GLOBAL II,JJ=',II,JJ
                CLOSE(99)
              ENDIF
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
                BBUVP(IIB+INX,JJB)=-ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+INX,JJB+INX)=0.0D0
                BBUVP(IIB+INX,JJB+2*INX)=0.0D0
                BBUVP(IIB+INX,JJB+3*INX)=1.0D0/ROU*(DJJ(1,1)*AX(I1,J1)+DJJ(1,2)*AY(I1,J1)+DJJ(1,3)*AZ(I1,J1)+DJJ(1,4)*BX(I1,J1)+DJJ(1,5)*BY(I1,J1)+DJJ(1,6)*BZ(I1,J1)+&
                     DJJ(1,7)*CXY(I1,J1)+DJJ(1,8)*CYZ(I1,J1)+DJJ(1,9)*CZX(I1,J1))
                BBUVP(IIB+2*INX,JJB)=0.0D0
                BBUVP(IIB+2*INX,JJB+INX)=-ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+2*INX,JJB+2*INX)=0.0D0
                BBUVP(IIB+2*INX,JJB+3*INX)=1.0D0/ROU*(DJJ(2,1)*AX(I1,J1)+DJJ(2,2)*AY(I1,J1)+DJJ(2,3)*AZ(I1,J1)+DJJ(2,4)*BX(I1,J1)+DJJ(2,5)*BY(I1,J1)+DJJ(2,6)*BZ(I1,J1)+&
                     DJJ(2,7)*CXY(I1,J1)+DJJ(2,8)*CYZ(I1,J1)+DJJ(2,9)*CZX(I1,J1))
                BBUVP(IIB+3*INX,JJB)=0.0D0
                BBUVP(IIB+3*INX,JJB+INX)=0.0D0
                BBUVP(IIB+3*INX,JJB+2*INX)=-ANYU*(DJJ(4,1)*AX(I1,J1)+DJJ(4,2)*AY(I1,J1)+DJJ(4,3)*AZ(I1,J1)+DJJ(4,4)*BX(I1,J1)+DJJ(4,5)*BY(I1,J1)+&
                     DJJ(4,6)*BZ(I1,J1)+DJJ(4,7)*CXY(I1,J1)+DJJ(4,8)*CYZ(I1,J1)+DJJ(4,9)*CZX(I1,J1)+DJJ(5,1)*AX(I1,J1)+DJJ(5,2)*AY(I1,J1)+DJJ(5,3)*AZ(I1,J1)+&
                     DJJ(5,4)*BX(I1,J1)+DJJ(5,5)*BY(I1,J1)+DJJ(5,6)*BZ(I1,J1)+DJJ(5,7)*CXY(I1,J1)+DJJ(5,8)*CYZ(I1,J1)+DJJ(5,9)*CZX(I1,J1)+DJJ(6,1)*AX(I1,J1)+&
                     DJJ(6,2)*AY(I1,J1)+DJJ(6,3)*AZ(I1,J1)+DJJ(6,4)*BX(I1,J1)+DJJ(6,5)*BY(I1,J1)+DJJ(6,6)*BZ(I1,J1)+DJJ(6,7)*CXY(I1,J1)+DJJ(6,8)*CYZ(I1,J1)+&
                     DJJ(6,9)*CZX(I1,J1))
                BBUVP(IIB+3*INX,JJB+3*INX)=1.0D0/ROU*(DJJ(3,1)*AX(I1,J1)+DJJ(3,2)*AY(I1,J1)+DJJ(3,3)*AZ(I1,J1)+DJJ(3,4)*BX(I1,J1)+DJJ(3,5)*BY(I1,J1)+DJJ(3,6)*BZ(I1,J1)+&
                     DJJ(3,7)*CXY(I1,J1)+DJJ(3,8)*CYZ(I1,J1)+DJJ(3,9)*CZX(I1,J1))
              END IF
              
          END DO
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
            IIT(ISP1)=IIT(ISP1)+1
            IF((IIT(ISP1).EQ.1).OR.(IIT(ISP1).EQ.2))THEN
              DO JA2=1,IAN
              DO IK2=1,4
                JJ2=IE(JA2,IM)
                ISP2=JJ2+NX*(IK2-1)
                IPL2=JA2+IAN*(IK2-1)
                IF (IM.EQ.5 .AND. ((JA1.EQ.565 .AND. JA2.EQ.1798) .OR. (JA1.EQ.2028 .AND. JA2.EQ.587)) .AND. IK1.EQ.1 .AND. IK2.EQ.1) THEN
                  OPEN(99,FILE='fortran_debug_5_565_1798.log',STATUS='UNKNOWN',POSITION='APPEND')
                  WRITE(99,'(A,4I10)') 'SSMMK add ISP1 ISP2 IPL1 IPL2=',ISP1,ISP2,IPL1,IPL2
                  WRITE(99,'(A,1P,E20.12)') 'BBUVP(IPL1,IPL2)=',BBUVP(IPL1,IPL2)
                  WRITE(99,'(A,1P,E20.12)') 'SSMMK before=',SSMMK(ISP1,ISP2)
                  CLOSE(99)
                ENDIF
                SSMMK(ISP1,ISP2)=SSMMK(ISP1,ISP2)+BBUVP(IPL1,IPL2)
                IF (IM.EQ.5 .AND. JA1.EQ.565 .AND. JA2.EQ.1798 .AND. IK1.EQ.1 .AND. IK2.EQ.1) THEN
                  OPEN(99,FILE='fortran_debug_5_565_1798.log',STATUS='UNKNOWN',POSITION='APPEND')
                  WRITE(99,'(A,1P,E20.12)') 'SSMMK after=',SSMMK(ISP1,ISP2)
                  CLOSE(99)
                ENDIF
              END DO             
              END DO 
            END IF
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
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          FFS(ISP1)=0.0D0
        END IF
      END DO
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II
        
        CALL TRANSF(5)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
     
          DO IJ=1,IAN
            ISP2=IE(IJ,5)
            SSMMK(ISP1,ISP2)=(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))
            IF (J1.EQ.565 .AND. IJ.EQ.1798) THEN
              OPEN(99,FILE='fortran_debug_5_565_1798.log',STATUS='UNKNOWN',POSITION='APPEND')
              WRITE(99,'(A,2I10,1P,E20.12)') 'SLOPE_U SSMMK(J1,IJ)=',J1,IJ,SSMMK(ISP1,ISP2)
              WRITE(99,'(A,1P,9E20.12)') 'DJJ4:', DJJ(4,1),DJJ(4,2),DJJ(4,3),DJJ(4,4),DJJ(4,5),DJJ(4,6),DJJ(4,7),DJJ(4,8),DJJ(4,9)
              CLOSE(99)
            ENDIF
          END DO
          FFS(ISP1)=0.0D0
        END IF      
      END DO

!     SLOPE V -----------------
 
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II+NX
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          FFS(ISP1)=0.0D0
        END IF
      END DO

      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II+NX
        
        CALL TRANSF(5)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
        
          DO IJ=1,IAN
            ISP2=IE(IJ,5)+NX
            SSMMK(ISP1,ISP2)=(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))
          END DO
          FFS(ISP1)=0.0D0
        END IF      
      END DO

!     SLOPE W -----------------

      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II+2*NX
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          FFS(ISP1)=0.0D0
        END IF
      END DO

      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
        II=IE(J1,5)
        ISP1=II+2*NX
        
        CALL TRANSF(5)
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,J1)
        END DO
        END DO
        
          DO IJ=1,IAN
            ISP2=IE(IJ,5)+2*NX
            SSMMK(ISP1,ISP2)=(DJJ(4,1)*AX(J1,IJ)+DJJ(4,2)*AY(J1,IJ)+DJJ(4,3)*AZ(J1,IJ)+DJJ(4,4)*BX(J1,IJ)+DJJ(4,5)*BY(J1,IJ)+DJJ(4,6)*BZ(J1,IJ)+&
                 DJJ(4,7)*CXY(J1,IJ)+DJJ(4,8)*CYZ(J1,IJ)+DJJ(4,9)*CZX(J1,IJ))
          END DO
          FFS(ISP1)=0.0D0
        END IF      
      END DO

!     ------------ DIFFERENTIAL OF PRESSURE SLOPE AT INFLOW ---------------- 

!     ------- SLOPE P AT INFLOW
      DO J1=1,IAN
        IF((J1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.J1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ)) THEN
        II=IE(J1,1)
        ISP1=II+3*NX
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          FFS(ISP1)=0.0D0
        END IF
      END DO

      J=1
      CALL TRANSF(J)
      
      DO I1=1,IAN
        IF((I1.GE.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+2*NAY*NAZ)) THEN
                
        DO IQ1=1,9
        DO IQ2=1,9
          DJJ(IQ1,IQ2)=DJJA(IQ1,IQ2,I1)
        END DO
        END DO        
        
          JJ=IE(I1,J)
          ISP1=JJ+3*NX
          DO IJ=1,IAN
            ISP2=IE(IJ,J)+3*NX
            SSMMK(ISP1,ISP2)=(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+&
                 DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ))
          END DO
          FFS(ISP1)=0.0D0 
        END IF     
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
          DO ISP2=1,NX*4
            SSMMK(ISP1,ISP2)=0.0D0
          END DO
          FFS(ISP1)=0.0D0
          END IF
        END DO
      END DO

      DO J1=1,5
        
        CALL TRANSF(J1)
        
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
            ISP2=IE(IJ,J1)
            SSMMK(ISP1,ISP2)=SSMMK(ISP1,ISP2)-ANYU*XL(I1,J1)*(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+&
                 DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)+&
                 DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+&
                 DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)+DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+&
                 DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ))
            ISP2=IE(IJ,J1)+NX
            SSMMK(ISP1,ISP2)=SSMMK(ISP1,ISP2)-ANYU*XM(I1,J1)*(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+&
                 DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)+&
                 DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+&
                 DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)+DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+&
                 DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ))
            ISP2=IE(IJ,J1)+2*NX
            SSMMK(ISP1,ISP2)=SSMMK(ISP1,ISP2)-ANYU*XN(I1,J1)*(DJJ(4,1)*AX(I1,IJ)+DJJ(4,2)*AY(I1,IJ)+DJJ(4,3)*AZ(I1,IJ)+&
                 DJJ(4,4)*BX(I1,IJ)+DJJ(4,5)*BY(I1,IJ)+DJJ(4,6)*BZ(I1,IJ)+DJJ(4,7)*CXY(I1,IJ)+DJJ(4,8)*CYZ(I1,IJ)+DJJ(4,9)*CZX(I1,IJ)+&
                 DJJ(5,1)*AX(I1,IJ)+DJJ(5,2)*AY(I1,IJ)+DJJ(5,3)*AZ(I1,IJ)+DJJ(5,4)*BX(I1,IJ)+DJJ(5,5)*BY(I1,IJ)+DJJ(5,6)*BZ(I1,IJ)+&
                 DJJ(5,7)*CXY(I1,IJ)+DJJ(5,8)*CYZ(I1,IJ)+DJJ(5,9)*CZX(I1,IJ)+DJJ(6,1)*AX(I1,IJ)+DJJ(6,2)*AY(I1,IJ)+DJJ(6,3)*AZ(I1,IJ)+&
                 DJJ(6,4)*BX(I1,IJ)+DJJ(6,5)*BY(I1,IJ)+DJJ(6,6)*BZ(I1,IJ)+DJJ(6,7)*CXY(I1,IJ)+DJJ(6,8)*CYZ(I1,IJ)+DJJ(6,9)*CZX(I1,IJ))
            ISP2=IE(IJ,J1)+3*NX
            SSMMK(ISP1,ISP2)=SSMMK(ISP1,ISP2)+((DJJ(1,1)*AX(I1,IJ)+DJJ(1,2)*AY(I1,IJ)+DJJ(1,3)*AZ(I1,IJ)+DJJ(1,4)*BX(I1,IJ)+DJJ(1,5)*BY(I1,IJ)+DJJ(1,6)*BZ(I1,IJ)+&
                 DJJ(1,7)*CXY(I1,IJ)+DJJ(1,8)*CYZ(I1,IJ)+DJJ(1,9)*CZX(I1,IJ))*XL(I1,J1)+(DJJ(2,1)*AX(I1,IJ)+DJJ(2,2)*AY(I1,IJ)+DJJ(2,3)*AZ(I1,IJ)+&
                 DJJ(2,4)*BX(I1,IJ)+DJJ(2,5)*BY(I1,IJ)+DJJ(2,6)*BZ(I1,IJ)+DJJ(2,7)*CXY(I1,IJ)+DJJ(2,8)*CYZ(I1,IJ)+DJJ(2,9)*CZX(I1,IJ))*XM(I1,J1)+&
                 (DJJ(3,1)*AX(I1,IJ)+DJJ(3,2)*AY(I1,IJ)+DJJ(3,3)*AZ(I1,IJ)+DJJ(3,4)*BX(I1,IJ)+DJJ(3,5)*BY(I1,IJ)+DJJ(3,6)*BZ(I1,IJ)+&
                 DJJ(3,7)*CXY(I1,IJ)+DJJ(3,8)*CYZ(I1,IJ)+DJJ(3,9)*CZX(I1,IJ))*XN(I1,J1))/ROU
          END DO
          FFS(ISP1)=0.0D0
          END IF      
        END DO
      END DO

!     --------------  BOUNDARY CONDITION AT WALL --------------

      DO J=1,NSB1
          JJ=NSBU(J)
          ISP2=JJ
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*SBU(J)
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=SBU(J)
      END DO
      
      DO J=1,NSB1
        JJ=NSBU(J)
          ISP2=JJ+NX
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*SBV(J)
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=SBV(J)
      END DO

      DO J=1,NSB1
        JJ=NSBU(J)
          ISP2=JJ+2*NX
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*SBW(J)
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=SBW(J)
      END DO

! -------- PRESSURE AT OUTOFLOW CONDITION 
    J1=5
     DO I3=1,NSB3      
       II=NSBUO(I3)
         ISP2=II+3*NX
         DO ISP1=1,NX*4
           FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*SBPO(I3)
           SSMMK(ISP1,ISP2)=0.0D0
           SSMMK(ISP2,ISP1)=0.0D0
         END DO
         SSMMK(ISP2,ISP2)=1.0D0
         FFS(ISP2)=SBPO(I3)
     END DO
       

      

!     ------------ INFLOW BOUNDARY CONDITION ------------------   
      DO J=1,NSB2
        JJ=NSBUI(J)
          ISP2=JJ
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*YU1(JJ)
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=YU1(JJ)
      END DO

      DO J=1,NSB2
        JJ=NSBUI(J)
          ISP2=JJ+NX
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*YV1(JJ)
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=YV1(JJ)
      END DO
      
      DO J=1,NSB2
        JJ=NSBUI(J)
          ISP2=JJ+2*NX
          DO ISP1=1,NX*4
            FFS(ISP1)=FFS(ISP1)-SSMMK(ISP1,ISP2)*YW1(JJ)
            SSMMK(ISP1,ISP2)=0.0D0
            SSMMK(ISP2,ISP1)=0.0D0
          END DO
          SSMMK(ISP2,ISP2)=1.0D0
          FFS(ISP2)=YW1(JJ)
      END DO

      WRITE(*,*)'CLEAR 103 !'
      n=nx*4
      OPEN(99,FILE='fortran_debug_5_565_1798.log',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(99,'(A,1P,E20.12)') 'SSMMK(8122,9234) before CSR=',SSMMK(8122,9234)
      CLOSE(99)

!     Output DENSE matrix to binary file for comparison with C (disabled: huge file)
      if (.false.) then
        write(*,*) 'Writing DENSE matrix to matrix_fortran_dense.bin...'
        write(*,*) '  Matrix dimension: ', n, ' x ', n
        open(unit=98, file='matrix_fortran_dense.bin', form='unformatted', &
             access='stream', status='replace')
        write(98) n           ! Matrix dimension (integer*4)
!       Write dense matrix in row-major order (like C)
        do i=1,n
          do j=1,n
            write(98) ssmmk(i,j)  ! real*8
          end do
        end do
        write(98) FFS(1:n)    ! RHS vector (real*8)
        close(98)
        write(*,*) 'Dense matrix written successfully.'
      end if

      call csr_S_matrix(n,n,ssmmk,mnr,ir,mn,ic,a) 

      write(*,*)'clear 99 !mnr,mn',mnr,mn

!     Output matrix to binary file for comparison with C
      write(*,*) 'Writing CSR matrix to matrix_fortran.bin...'
      open(unit=99, file='matrix_fortran.bin', form='unformatted', &
           access='stream', status='replace')
      write(99) n           ! Matrix dimension (integer*4)
      write(99) mn          ! Number of non-zeros (integer*4)
      write(99) ir(1:n+1)   ! Row pointers (integer*4)
      write(99) ic(1:mn)    ! Column indices (integer*4)
      write(99) a(1:mn)     ! Values (real*8)
      write(99) FFS(1:n)    ! RHS vector (real*8)
      close(99)
      write(*,*) 'Matrix written successfully.'

!     DEBUG: Print some diagonal values
      write(*,*) 'Checking diagonal values (first 20 rows):'
      do i=1,min(20,n)
        diag_val = 0.0d0
        do k=ir(i),ir(i+1)-1
          if (ic(k) .eq. i) then
            diag_val = a(k)
            exit
          end if
        end do
        write(*,'(A,I6,A,E15.6)') '  Row ', i, ' diagonal = ', diag_val
      end do

!     Count zero diagonals
      nzd = 0
      do i=1,n
        diag_val = 0.0d0
        do k=ir(i),ir(i+1)-1
          if (ic(k) .eq. i) then
            diag_val = a(k)
            exit
          end if
        end do
        if (abs(diag_val) .lt. 1.0d-30) nzd = nzd + 1
      end do
      write(*,*) 'Total rows with zero diagonal: ', nzd

      do i=1,64
        ipt(i)=0
        iparm(i)=0
      end do
      
      maxfct=1
      mnum=1
      mtype=11
      iphase=13 
      nrhs=1     
      iparm(1)=1
      iparm(2)=2
      iparm(3)=8
      iparm(4)=0
      iparm(5)=0
      iparm(6)=0
      iparm(8)=2
      iparm(10)=13
      iparm(11)=1
      iparm(13)=1
      iparm(18)=-1
      iparm(19)=-1
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

      
      DO IK=1,4  
      DO I=1,NX
        ISP=I+NX*(IK-1)
        IF(ISP.LE.NX) THEN
          YU1(I)=x(ISP,1)
        END IF
        IF((ISP.GT.NX).AND.(ISP.LE.2*NX)) THEN
          YV1(I)=x(ISP,1)
        END IF
        IF((ISP.GT.2*NX).AND.(ISP.LE.3*NX)) THEN
          YW1(I)=x(ISP,1)
        END IF
        IF(ISP.GT.3*NX) THEN
          YP1(I)=x(ISP,1)
        END IF
      END DO
      END DO
   
        DO I=1,NX
           WRITE(2,2100) XX(I),YY(I),ZZ(I),YU1(I),YV1(I),YW1(I),YP1(I),Re
 !          WRITE(*,*) XX(I),YY(I),ZZ(I),YU1(I),YV1(I),YW1(I),YP1(I),Re
        END DO
 2100   FORMAT(8D25.16)
           do i=1,nx*4
                     if(ssmmk(i,i)==0.0e0) then
!                 write(*,*) i,i/(nx),ssmmk(i,i)
                    end if
           end do
           
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
         V1=YU1(II)
         FLOW2=FLOW2+V1*WW(I1-8-4*(NAX+NAY+NAZ)-2*NAX*NAY-NAX*NAZ)
         END IF
      END DO
      FLOW2=FLOW2*D0*H0
      WRITE(*,*) FLOW2
! SECTION 3 (Straight Section)
      IM=3
      FLOW3=0.0D0
      DO I1=1,IAN
         IF(I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ) THEN
         II=IE(I1,IM)
         V1=YU1(II)
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
         V1=YU1(II)
         FLOW4=FLOW4+V1*WW(I1-8-4*(NAX+NAY+NAZ)-2*NAX*NAY-NAX*NAZ)
         END IF
      END DO
      FLOW4=FLOW4*D0*H0
      WRITE(*,*) FLOW4
      
      FNAME='out_flow_results_S_Duct_pipe_Steady_Stokes_10_11_11_3D_Re500_CSR_LMN.dat'
      OPEN(2,FILE=FNAME)
        WRITE(2,*)'FLOW RATE AT CHANNEL SECTION --------'
        WRITE(2,*)'FLOW1(Inlet)           =',FLOW1
        WRITE(2,*)'FLOW2(Straight Section)=',FLOW2
        WRITE(2,*)'FLOW3(Straight Section)=',FLOW3
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


