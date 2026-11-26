!   generation of curve pipe collocation for formulation of Stokes eq. by OCFEM
!   coded by T,Okubo  2024.4.4.(WED)
!
!
    IMPLICIT REAL*8(A-H,O-Z)
    PARAMETER(KNX=50000,KMX=5,KCN=10000,KC=40,KXLMN=20000)
    
    COMMON /PARAM1/ NX,MX,NAX,NAY,NAZ,IAX,IAY,IAZ,IAN,J1,I1,ID
    COMMON /PARAM2/ X1,Z1,Y1,D0,H0,U00,ROU
    COMMON /COLLO1/ IE(KCN,KMX),XX(KNX),YY(KNX),ZZ(KNX),XXC(KC),YYC(KC),ZZC(KC)
    COMMON /COLLO2/ GX(KCN),GY(KCN),GZ(KCN),XL(KCN,KMX),XM(KCN,KMX),XN(KCN,KMX)
    COMMON /BOUND1/ JD1(KXLMN),JD2(KXLMN),NSBU(KXLMN),SBU(KXLMN),SBV(KXLMN),SBW(KXLMN),NSBUI(KXLMN),SBUI(KXLMN),SBVI(KXLMN),SBWI(KXLMN)
    COMMON /BOUND2/ NSBUO(KXLMN),SBPO(KXLMN),JSBU(KXLMN),ISBU(KXLMN),JSBUI(KXLMN),ISBUI(KXLMN),JSBUO(KXLMN),ISBUO(KXLMN)
    COMMON /BOUND4/ NSB1,NSB2,NSB3



    CHARACTER FNAME*100
    
    MX=5
    WRITE(*,*)'INPUT FILE NAME OF COLLOCATION DATA (calcoll_ .dat)'
    READ(*,*) FNAME
    OPEN(2,FILE=FNAME)
      READ(2,1000) NAX,NAY,NAZ
      IAX=NAX+2
      IAY=NAY+2
      IAZ=NAZ+2
      IAN=(NAX+2)*(NAY+2)*(NAZ+2)
      NNA=NA
      DO I=1,IAX
        READ(2,*) XXC(I)
      END DO
      DO I=1,IAY
        READ(2,*) YYC(I)
      END DO
      DO I=1,IAZ
        READ(2,*) ZZC(I)
      END DO
    CLOSE(2)
1000 FORMAT(3I5)


     WRITE(*,*)'INPUT FILE NAME OF X1,Z1,Y1,D0,H0,U00,ROU PARAMETER:'
       READ(*,*) X1,Z1,Y1,D0,H0,U00,ROU

        DO I3=1,IAZ
        DO I2=1,IAY
        DO I1=1,IAX
!
          IF((I1.EQ.1).AND.(I2.EQ.1).AND.(I3.EQ.1)) THEN
            IK=1
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.1).AND.(I3.EQ.1)) THEN
            IK=2
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.IAY).AND.(I3.EQ.1)) THEN
            IK=3
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.1).AND.(I2.EQ.IAY).AND.(I3.EQ.1)) THEN
            IK=4
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          ENDIF
          IF((I1.EQ.1).AND.(I2.EQ.1).AND.(I3.EQ.IAZ)) THEN
            IK=5
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)            
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.1).AND.(I3.EQ.IAZ)) THEN
            IK=6
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.IAY).AND.(I3.EQ.IAZ)) THEN
            IK=7
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.1).AND.(I2.EQ.IAY).AND.(I3.EQ.IAZ)) THEN
            IK=8
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
!
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.(I2.EQ.1).AND.(I3.EQ.1)) THEN
            IK=8+(I1-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.IAX).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.(I3.EQ.1)) THEN
            IK=8+NAX+(I2-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
!            WRITE(*,*) GX(IK),GY(IK),GZ(IK)
          END IF
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.(I2.EQ.IAY).AND.(I3.EQ.1)) THEN
            IK=8+NAX+NAY+(I1-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.1).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.(I3.EQ.1)) THEN
            IK=8+2*NAX+NAY+(I2-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.(I2.EQ.1).AND.(I3.EQ.IAZ)) THEN
            IK=8+2*(NAX+NAY)+(I1-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.IAX).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.(I3.EQ.IAZ)) THEN
            IK=8+2*(NAX+NAY)+NAX+(I2-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.(I2.EQ.IAY).AND.(I3.EQ.IAZ)) THEN
            IK=8+3*(NAX+NAY)+(I1-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.1).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.(I3.EQ.IAZ)) THEN
            IK=8+3*(NAX+NAY)+NAX+(I2-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.1).AND.(I2.EQ.1).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NAX+NAY)+(I3-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.1).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NAX+NAY)+NAZ+(I3-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.IAX).AND.(I2.EQ.IAY).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NAX+NAY)+2*NAZ+(I3-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.1).AND.(I2.EQ.IAY).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NAX+NAY)+3*NAZ+(I3-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
!
          IF((((I1.GT.1).AND.(I1.LT.IAX)).AND.((I2.GT.1).AND.(I2.LT.IAY))).AND.(I3.EQ.1)) THEN
            IK=8+4*(NAX+NAY+NAZ)+(I2-2)*NAX+(I1-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((((I1.GT.1).AND.(I1.LT.IAX)).AND.((I2.GT.1).AND.(I2.LT.IAY))).AND.(I3.EQ.IAZ)) THEN
            IK=8+4*(NAX+NAY+NAZ)+NAX*NAY+(I2-2)*NAX+(I1-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I2.EQ.1).AND.(((I1.GT.1).AND.(I1.LT.IAX)).AND.((I3.GT.1).AND.(I3.LT.IAZ)))) THEN
            IK=8+4*(NAX+NAY+NAZ)+2*NAX*NAY+(I3-2)*NAX+(I1-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.IAX).AND.(((I2.GT.1).AND.(I2.LT.IAY)).AND.((I3.GT.1).AND.(I3.LT.IAZ)))) THEN
            IK=8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+(I3-2)*NAY+(I2-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I2.EQ.IAY).AND.(((I1.GT.1).AND.(I1.LT.IAX)).AND.((I3.GT.1).AND.(I3.LT.IAZ)))) THEN
            IK=8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ+(I3-2)*NAX+(I1-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
          IF((I1.EQ.1).AND.(((I2.GT.1).AND.(I2.LT.IAY)).AND.((I3.GT.1).AND.(I3.LT.IAZ)))) THEN
            IK=8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ+(I3-2)*NAY+(I2-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
!
          IF(((I1.GT.1).AND.(I1.LT.IAX)).AND.((I2.GT.1).AND.(I2.LT.IAY)).AND.((I3.GT.1).AND.(I3.LT.IAZ))) THEN
            IK=8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ)+(I3-2)*NAX*NAY+(I2-2)*NAX+(I1-1)
            GX(IK)=XXC(I1)
            GY(IK)=YYC(I2)
            GZ(IK)=ZZC(I3)
          END IF
!          write(*,*) ik,gx(ik),gy(ik),gz(ik)
        END DO
      END DO
    END DO
    
    
     DO IK=1,IAN
!         WRITE(*,*) IK,GX(IK),GY(IK),GZ(IK)
     END DO

    MX=5
    IP=0
    DO J1=1,MX
      IF(J1==1) THEN
        GO TO 10
      END IF
      J2=J1-1
      DO I1=1,IAN
          IF(I1==1) THEN
            I2=2
              IE(I1,J1)=IE(I2,J2)
              CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1==4) THEN
            I2=3
              IE(I1,J1)=IE(I2,J2)
              CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1==5) THEN
            I2=6
              IE(I1,J1)=IE(I2,J2)
              CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1==8) THEN
            I2=7
              IE(I1,J1)=IE(I2,J2)
              CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1.GT.8+2*NAX+NAY.AND.I1.LT.9+2*(NAX+NAY)) THEN
             I2=I1-NAX-NAY
              IE(I1,J1)=IE(I2,J2)
              CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1.GT.8+2*(NAX+NAY)+2*NAX+NAY.AND.I1.LT.9+4*(NAX+NAY)) THEN
             I2=I1-NAX-NAY
              IE(I1,J1)=IE(I2,J2)
              CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1.GT.8+4*(NAX+NAY).AND.I1.LT.9+4*(NAX+NAY)+NAZ) THEN
             I2=I1+NAZ
              IE(I1,J1)=IE(I2,J2)
              CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1.GT.8+4*(NAX+NAY)+3*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)) THEN
             I2=I1-NAZ
              IE(I1,J1)=IE(I2,J2)
              CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1.GT.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ)+NAY*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
             I2=I1-NAX*NAZ-NAY*NAZ
              IE(I1,J1)=IE(I2,J2)
              CALL CURVE_PIPE(I1,J1)
          END IF
          
        END DO
        
        DO I1=1,IAN
          IF((I1.EQ.2).OR.(I1.EQ.3).OR.(I1.EQ.6).OR.(I1.EQ.7)) THEN
            IP=IP+1
            IE(I1,J1)=IP
            CALL CURVE_PIPE(I1,J1)
          END IF
          IF((I1.GT.8.AND.I1.LE.8+2*NAX+NAY).OR.(I1.GT.8+2*(NAX+NAY).AND.I1.LE.8+2*(NAX+NAY)+2*NAX+NAY)) THEN
            IP=IP+1
            IE(I1,J1)=IP
            CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1.GT.8+4*(NAX+NAY)+NAZ.AND.I1.LE.8+4*(NAX+NAY)+3*NAZ) THEN
            IP=IP+1
            IE(I1,J1)=IP
            CALL CURVE_PIPE(I1,J1)
          END IF
          IF((I1.GT.8+4*(NAX+NAY+NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ)+NAY*NAZ)) THEN
            IP=IP+1
            IE(I1,J1)=IP
            CALL CURVE_PIPE(I1,J1)
          END IF
          IF(I1.GT.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ).AND.I1.LE.8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ)+NAX*NAY*NAZ) THEN
            IP=IP+1
            IE(I1,J1)=IP
            CALL CURVE_PIPE(I1,J1)
          END IF
        END DO
      GO TO 30

 10   DO I1=1,IAN
        IP=IP+1
        IE(I1,J1)=IP
        CALL CURVE_PIPE(I1,J1)
      END DO
            
 30 CONTINUE
    END DO
    
    NX=IP

    CALL NORMAL_XYZ
    
    CALL BOUNDARY_ELE_COND

    CALL OUTPUT


    END






    SUBROUTINE CURVE_PIPE(I,J)
     
    IMPLICIT REAL*8(A-H,O-Z)
    PARAMETER(KNX=50000,KMX=5,KCN=10000,KC=40,KXLMN=20000)
    
    COMMON /PARAM1/ NX,MX,NAX,NAY,NAZ,IAX,IAY,IAZ,IAN,J1,I1,ID
    COMMON /PARAM2/ X1,Z1,Y1,D0,H0,U00,ROU
    COMMON /COLLO1/ IE(KCN,KMX),XX(KNX),YY(KNX),ZZ(KNX),XXC(KC),YYC(KC),ZZC(KC)
    COMMON /COLLO2/ GX(KCN),GY(KCN),GZ(KCN),XL(KCN,KMX),XM(KCN,KMX),XN(KCN,KMX)
    COMMON /BOUND1/ JD1(KXLMN),JD2(KXLMN),NSBU(KXLMN),SBU(KXLMN),SBV(KXLMN),SBW(KXLMN),NSBUI(KXLMN),SBUI(KXLMN),SBVI(KXLMN),SBWI(KXLMN)
    COMMON /BOUND2/ NSBUO(KXLMN),SBPO(KXLMN),JSBU(KXLMN),ISBU(KXLMN),JSBUI(KXLMN),ISBUI(KXLMN),JSBUO(KXLMN),ISBUO(KXLMN)
    COMMON /BOUND4/ NSB1,NSB2,NSB3


      PAI=3.14159265358979323846D0
      
        XX(IE(I,J))=(J-1)*X1+GX(I)*X1
        YY(IE(I,J))=Y1+D0*GY(I)
        ZZ(IE(I,J))=Z1+H0*GZ(I)
    
    RETURN
    END
   
    
    SUBROUTINE NORMAL_XYZ
    
    IMPLICIT REAL*8(A-H,O-Z)
    PARAMETER(KNX=50000,KMX=5,KCN=10000,KC=40,KXLMN=20000)
    
    COMMON /PARAM1/ NX,MX,NAX,NAY,NAZ,IAX,IAY,IAZ,IAN,J1,I1,ID
    COMMON /PARAM2/ X1,Z1,Y1,D0,H0,U00,ROU
    COMMON /COLLO1/ IE(KCN,KMX),XX(KNX),YY(KNX),ZZ(KNX),XXC(KC),YYC(KC),ZZC(KC)
    COMMON /COLLO2/ GX(KCN),GY(KCN),GZ(KCN),XL(KCN,KMX),XM(KCN,KMX),XN(KCN,KMX)
    COMMON /BOUND1/ JD1(KXLMN),JD2(KXLMN),NSBU(KXLMN),SBU(KXLMN),SBV(KXLMN),SBW(KXLMN),NSBUI(KXLMN),SBUI(KXLMN),SBVI(KXLMN),SBWI(KXLMN)
    COMMON /BOUND2/ NSBUO(KXLMN),SBPO(KXLMN),JSBU(KXLMN),ISBU(KXLMN),JSBUI(KXLMN),ISBUI(KXLMN),JSBUO(KXLMN),ISBUO(KXLMN)
    COMMON /BOUND4/ NSB1,NSB2,NSB3



!   CHANNEL WALL
    PAI=3.14159265358979323846D0
    ID=0
    DO J1=1,MX
      DO I1=1,IAN
        IF((I1.GT.8+4*(NAX+NAY+NAZ).AND.I1.LT.9+4*(NAX+NAY+NAZ)+NAX*NAY).OR.&
            (I1.GE.1.AND.I1.LE.4).OR.(I1.GT.8.AND.I1.LT.9+2*(NAX+NAY))) THEN
          II=IE(I1,J1)
          IF(J1.GE.1.OR.J1.LE.5) THEN
            XL(I1,J1)=0.0D0
            XM(I1,J1)=0.0D0
            XN(I1,J1)=-1.0D0
            ID=ID+1
            JD1(ID)=I1
            JD2(ID)=J1
          END IF
        END IF
        IF((I1.GT.8+4*(NAX+NAY+NAZ)+NAX*NAY.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY).OR.&
            (I1.GE.5.AND.I1.LE.8).OR.(I1.GT.8+2*(NAX+NAY).AND.I1.LT.9+4*(NAX+NAY))) THEN
          II=IE(I1,J1)
          IF(J1.GE.1.OR.J1.LE.5) THEN
            XL(I1,J1)=0.0D0
            XM(I1,J1)=0.0D0
            XN(I1,J1)=1.0D0
            ID=ID+1
            JD1(ID)=I1
            JD2(ID)=J1
          END IF
        END IF
        IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ).OR.&
            (I1.GT.8+4*(NAX+NAY).AND.I1.LT.9+4*(NAX+NAY)+2*NAZ)) THEN
          II=IE(I1,J1)
          IF(J1.GE.1.AND.J1.LE.5) THEN
            XL(I1,J1)=0.0D0
            XM(I1,J1)=-1.0D0
            XN(I1,J1)=0.0D0
            ID=ID+1
            JD1(ID)=I1
            JD2(ID)=J1
          END IF
        END IF
        IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ).OR.&
            (I1.GT.8+4*(NAX+NAY)+2*NAZ.AND.I1.LT.9+4*(NAX+NAY)+4*NAZ)) THEN
          II=IE(I1,J1)
          IF(J1.GE.1.AND.J1.LE.5) THEN
            XL(I1,J1)=0.0D0
            XM(I1,J1)=1.0D0
            XN(I1,J1)=0.0D0
            ID=ID+1
            JD1(ID)=I1
            JD2(ID)=J1
          END IF
        END IF
      END DO
    END DO
    
! CHANNEL CONNECT
    DO J1=1,MX
      DO I1=1,IAN
      IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ)) THEN
          II=IE(I1,J1)
          IF(J1.GE.1.AND.J1.LE.5) THEN
            XL(I1,J1)=1.0D0
            XM(I1,J1)=0.0D0
            XN(I1,J1)=0.0D0
            ID=ID+1
            JD1(ID)=I1
            JD2(ID)=J1
          END IF
        END IF
        IF((I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ))) THEN
          II=IE(I1,J1)
          IF(J1.GE.1.AND.J1.LE.5) THEN
            XL(I1,J1)=-1.0D0
            XM(I1,J1)=0.0D0
            XN(I1,J1)=0.0D0
            ID=ID+1
            JD1(ID)=I1
            JD2(ID)=J1
          END IF
        END IF
        
      END DO
    END DO
    
    RETURN
    END
      
      
    SUBROUTINE BOUNDARY_ELE_COND
    
    IMPLICIT REAL*8(A-H,O-Z)
    PARAMETER(KNX=50000,KMX=5,KCN=10000,KC=40,KXLMN=20000)
    
    COMMON /PARAM1/ NX,MX,NAX,NAY,NAZ,IAX,IAY,IAZ,IAN,J1,I1,ID
    COMMON /PARAM2/ X1,Z1,Y1,D0,H0,U00,ROU
    COMMON /COLLO1/ IE(KCN,KMX),XX(KNX),YY(KNX),ZZ(KNX),XXC(KC),YYC(KC),ZZC(KC)
    COMMON /COLLO2/ GX(KCN),GY(KCN),GZ(KCN),XL(KCN,KMX),XM(KCN,KMX),XN(KCN,KMX)
    COMMON /BOUND1/ JD1(KXLMN),JD2(KXLMN),NSBU(KXLMN),SBU(KXLMN),SBV(KXLMN),SBW(KXLMN),NSBUI(KXLMN),SBUI(KXLMN),SBVI(KXLMN),SBWI(KXLMN)
    COMMON /BOUND2/ NSBUO(KXLMN),SBPO(KXLMN),JSBU(KXLMN),ISBU(KXLMN),JSBUI(KXLMN),ISBUI(KXLMN),JSBUO(KXLMN),ISBUO(KXLMN)
    COMMON /BOUND4/ NSB1,NSB2,NSB3

    
    IP=0
    DO J1=1,MX
      DO I1=1,IAN
        II=IE(I1,J1)
        IF(I1.GE.1.AND.I1.LE.8) THEN
          IP=IP+1
          NSBU(IP)=II
          JSBU(IP)=J1
          ISBU(IP)=I1
          SBU(IP)=0.0D0
          SBV(IP)=0.0D0
          SBW(IP)=0.0D0
        END IF
        IF((I1.GT.8.AND.I1.LT.9+4*(NAX+NAY+NAZ)).OR.&
          (I1.GT.8+4*(NAX+NAY+NAZ).AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ).OR.&
          (I1.GT.8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ.AND.I1.LT.9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ)) THEN
          IP=IP+1
          NSBU(IP)=II
          JSBU(IP)=J1
          ISBU(IP)=I1
          SBU(IP)=0.0D0
          SBV(IP)=0.0D0
          SBW(IP)=0.0D0
        END IF
      END DO
     END DO
     NSB1=IP
     
     
     J1=1
     IP1=0
     DO I1=9+4*(NAX+NAY+NAZ)+2*NAX*NAY+2*NAX*NAZ+NAY*NAZ,8+4*(NAX+NAY+NAZ)+2*(NAX*NAY+NAX*NAZ+NAY*NAZ)
       II=IE(I1,J1)
       IP1=IP1+1
       NSBUI(IP1)=II
       JSBUI(IP1)=J1
       ISBUI(IP1)=I1
!      Tensor-product parabolic inflow, with U00 as the centerline maximum.
       SBUI(IP1)=16.0D0*U00*GY(I1)*(1.0D0-GY(I1))*GZ(I1)*(1.0D0-GZ(I1))
       SBVI(IP1)=0.0D0
       SBWI(IP1)=0.0D0
     END DO
     NSB2=IP1
     
!    Atmospheric discharge: zero gauge pressure on the entire outlet.
     J1=5
     IP2=0       
     DO I1=9+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ,8+4*(NAX+NAY+NAZ)+2*NAX*NAY+NAX*NAZ+NAY*NAZ
       II=IE(I1,J1)
       IP2=IP2+1
       NSBUO(IP2)=II
       JSBUO(IP2)=J1
       ISBUO(IP2)=I1
       SBPO(IP2)=0.0D0
     END DO
     DO I1=9+4*(NAX+NAY)+NAZ,8+4*(NAX+NAY)+3*NAZ
       II=IE(I1,J1)
       IP2=IP2+1
       NSBUO(IP2)=II
       JSBUO(IP2)=J1
       ISBUO(IP2)=I1
       SBPO(IP2)=0.0D0
     END DO
     DO I1=9+NAX,8+NAX+NAY
       II=IE(I1,J1)
       IP2=IP2+1
       NSBUO(IP2)=II
       JSBUO(IP2)=J1
       ISBUO(IP2)=I1
       SBPO(IP2)=0.0D0
     END DO
     DO I1=9+2*(NAX+NAY)+NAX,8+2*(NAX+NAY)+NAX+NAY
       II=IE(I1,J1)
       IP2=IP2+1
       NSBUO(IP2)=II
       JSBUO(IP2)=J1
       ISBUO(IP2)=I1
       SBPO(IP2)=0.0D0
     END DO
     I1=2
       II=IE(I1,J1)
       IP2=IP2+1
       NSBUO(IP2)=II
       JSBUO(IP2)=J1
       ISBUO(IP2)=I1
       SBPO(IP2)=0.0D0
     I1=3
       II=IE(I1,J1)
       IP2=IP2+1
       NSBUO(IP2)=II
       JSBUO(IP2)=J1
       ISBUO(IP2)=I1
       SBPO(IP2)=0.0D0
     I1=6
       II=IE(I1,J1)
       IP2=IP2+1
       NSBUO(IP2)=II
       JSBUO(IP2)=J1
       ISBUO(IP2)=I1
       SBPO(IP2)=0.0D0
     I1=7
       II=IE(I1,J1)
       IP2=IP2+1
       NSBUO(IP2)=II
       JSBUO(IP2)=J1
       ISBUO(IP2)=I1
       SBPO(IP2)=0.0D0    
     NSB3=IP2
     
     RETURN
     END 
      
      
    SUBROUTINE OUTPUT
    
    IMPLICIT REAL*8(A-H,O-Z)
    PARAMETER(KNX=50000,KMX=5,KCN=10000,KC=40,KXLMN=20000)
    
    COMMON /PARAM1/ NX,MX,NAX,NAY,NAZ,IAX,IAY,IAZ,IAN,J1,I1,ID
    COMMON /PARAM2/ X1,Z1,Y1,D0,H0,U00,ROU
    COMMON /COLLO1/ IE(KCN,KMX),XX(KNX),YY(KNX),ZZ(KNX),XXC(KC),YYC(KC),ZZC(KC)
    COMMON /COLLO2/ GX(KCN),GY(KCN),GZ(KCN),XL(KCN,KMX),XM(KCN,KMX),XN(KCN,KMX)
    COMMON /BOUND1/ JD1(KXLMN),JD2(KXLMN),NSBU(KXLMN),SBU(KXLMN),SBV(KXLMN),SBW(KXLMN),NSBUI(KXLMN),SBUI(KXLMN),SBVI(KXLMN),SBWI(KXLMN)
    COMMON /BOUND2/ NSBUO(KXLMN),SBPO(KXLMN),JSBU(KXLMN),ISBU(KXLMN),JSBUI(KXLMN),ISBUI(KXLMN),JSBUO(KXLMN),ISBUO(KXLMN)
    COMMON /BOUND4/ NSB1,NSB2,NSB3

    
    CHARACTER FNAME*100
    
    WRITE(*,*) 'OUTPUT FILE NAME OF COLLOCATION ELEMNT RELATIONY.(colo_ele_pipe_c_.dat)'
    READ(*,*) FNAME
    OPEN(2,FILE=FNAME)
      WRITE(2,1000) NX,MX
      DO J=1,MX
        WRITE(2,1100) (IE(I,J),I=1,IAN)
      END DO
      DO I=1,NX
        WRITE(2,1200) I,XX(I),YY(I),ZZ(I)
      END DO
    CLOSE(2)
1000 FORMAT(2I10)
1100 FORMAT(7I10)
1200 FORMAT(I10,3D25.16)

    WRITE(*,*)'OUTPUT FILE NAME OF BOUNARY CONNECTION (boucone_c_.dat)'
    READ(*,*) FNAME
    OPEN(2,FILE=FNAME)
      WRITE(2,2000) ID
      DO I1=1,ID
        WRITE(2,2100) JD1(I1),JD2(I1),XL(JD1(I1),JD2(I1)),XM(JD1(I1),JD2(I1)),XN(JD1(I1),JD2(I1))
      END DO
    CLOSE(2)
2000 FORMAT(I10)
2100 FORMAT(2I10,3D25.16)


      
    WRITE(*,*)'OUTPUT FILE NAME OF BOUNDARY CONDITION (bound_c_.dat)'
    READ(*,*) FNAME
    OPEN(2,FILE=FNAME)
      WRITE(2,3000) NSB1
      DO I=1,NSB1
        WRITE(2,3100) NSBU(I),JSBU(I),ISBU(I),SBU(I),SBV(I),SBW(I)
      END DO
      WRITE(2,3000) NSB2
      DO I=1,NSB2
        WRITE(2,3100) NSBUI(I),JSBUI(I),ISBUI(I),SBUI(I),SBVI(I),SBWI(I)
      END DO
      WRITE(2,3000) NSB3
      DO I=1,NSB3
        WRITE(2,3200) NSBUO(I),JSBUO(I),ISBUO(I),SBPO(I)
      END DO
    CLOSE(2)
 3000 FORMAT(I10)
 3100 FORMAT(3I10,3D25.16)
 3200 FORMAT(3I10,D25.16)
    
    RETURN
    END
      
        
        
        
        
        
        
        
          
        
        
        
        
        
        
        
        
        
        
        
        
        
        
