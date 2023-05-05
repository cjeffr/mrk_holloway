      Program MRKHOL 
 
C***    FILE ISOCHORE 
            IMPLICIT DOUBLE PRECISION (A-H, O-Z)
            CHARACTER*80 ZOUT
            CHARACTER*4 NAMES(8)
            DIMENSION X(8),VOLI(5),FORMWT(8),POUT(5),DEN(5)
            DATA RBAR/83.117D0/, MIXNUM/8/
            DATA NAMES/' CO2','  CO',' CH4','  H2',' H2O',' H2S'
     $,' SO2','  N2'/
          DATA FORMWT/44.D0,28.D0,16.D0,2.D0,18.D0,34.06D0,
     164.06D0,28.0D0/
C     DATA X/0.25,0.00,0.00,0.00,0.75,0.00,0.00,0.00/
C
            WRITE(*,69)
   69 FORMAT('ENTER NAME OF OUTPUT FILE')
            READ(*,555) ZOUT
  555 FORMAT(A80)
!             OPEN (UNIT=8, FILE=ZOUT)
   75    WRITE(*,40)
   40    FORMAT('ENTER T (DEGC) AND MOLAR VOLUME (CC/MOLE)')
            READ (*,*) TSTART,VSTART
            IF(TSTART.LT.1.D-2) STOP
C       TSTART=400.0D0
C       VSTART=30.0D0
            WRITE(*,61)
   61   FORMAT('ENTER MOLE FRACTION CO2')
          READ (*,*) X(1)
            WRITE(*,62)
   62 FORMAT('ENTER MOLE FRACTION CO')
            READ(*,*) X(2)
            WRITE(*,63)      
   63 FORMAT('ENTER MOLE FRACTION CH4')
            READ(*,*) X(3)
            WRITE(*,64)
   64 FORMAT('ENTER MOLE FRACTION H2')
            READ(*,*) X(4)
            WRITE(*,65)
   65 FORMAT('ENTER MOLE FRACTION H2O')
            READ(*,*) X(5)
            WRITE(*,66)
   66 FORMAT('ENTER MOLE FRACTION H2S')
            READ(*,*) X(6)
            WRITE(*,67)
   67 FORMAT('ENTER MOLE FRACTION SO2')
            READ(*,*) X(7)
            WRITE(*,68)
   68 FORMAT('ENTER MOLE FRACTION N2')
            READ(*,*) X(8)
            XTOTAL=X(1)+X(2)+X(3)+X(4)+X(5)+X(6)+X(7)+X(8)
            IF(XTOTAL.EQ.1.00) GO TO 556
            WRITE(*,71)
   71 FORMAT('YOUR MOLE FRACTIONS DO NOT ADD UP TO 1.0! START AGAIN')
      GO TO 72
  
  556    WRITE(8,100) NAMES,X,TSTART,VSTART
  100 FORMAT(' THE MOLE FRACTIONS ARE:'/' ',8(5X,A4)/' ',8F9.3/
     1'  THE STARTING TEMPERATURE = ',F6.0,' DEG. C'/
     2'  THE STARTING MOLAR VOLUME = ',F6.3,' CM3/MOL')
            WRITE(*,100) NAMES,X,TSTART,VSTART
C       CALCULATE MEAN FORMULA WEIGHT
            FORM=0.0
            DO 3 I=1,MIXNUM
    3  FORM=FORM+X(I)*FORMWT(I)
            TK=TSTART+273.15
            PSTART=1.D3
C       CALCULATE TOTAL VOLUME FOR AN ISOCHORE
            DO 5 I=1,5
            VOLI(I)=VSTART+1.D1*I-1.D1
   5  DEN(I)=FORM/VOLI(I)
            WRITE(8,101)
  101 FORMAT('0',30x,'ISOCHORE PRESSURES IN BARS')
            WRITE(8,102) VOLI,DEN
          WRITE(*,101)
            WRITE(*,102) VOLI,DEN
  102 FORMAT('0   T DEG C  MOL VOL=',5F10.3/' ',10X,'  DENSITY=',5F10.4
     1/)
          DO 10 K=1,11
            T=TK-1.D2+1.D2*K
            TC=T-273.15
            CALL MRKMIX(TK,X,BSUM,ASUM)
C     CALCULATE ISOCHORE
            DO 20 J=1,5
            VOL=VOLI(J)
            P=0.D0
            IF(VOL.LT.(BSUM-1.0D0)) GO TO 19
            
            ATERM=((RBAR*T)/(VOL-BSUM))
            BTERM=ASUM/((DSQRT(T))*((VOL*VOL)+(BSUM*VOL)))
            P=ATERM-BTERM
  19   CONTINUE
  20   POUT(J)=P
            WRITE(8,105) TC,POUT
  105 FORMAT(' ',F10.0,10X,5F10.0)
            WRITE(*,105) TC,POUT
  10   CONTINUE
  72   WRITE(*,80)
   80 FORMAT('DO YOU WANT TO DO ANOTHER COMPOSITION? NO=0,YES=1')
            READ(*,*) INEXT
            IF(INEXT.EQ.0) GO TO 85
            IF(INEXT.EQ.1) GO TO 75
  85   STOP
      contains
      SUBROUTINE MRKMIX(T,Y,BSUM,ASUM)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION A(8),B(8),Y(8)
          INTEGER MIXNUM/8/
          DATA A/46.D6,16.98D6,31.59D6,3.56D6,35.D6,87.9D6,142.6D6,
     115.382D6/
        DATA B/2.97D1,2.738D1,2.9703D1,1.515D1,1.46D1,2.D1,3.94D1,
     12.68D1/
        DATA R/82.05D0/
          IF(T.LT.1.D-4) T=1.D0
          TCEL=T-273.15
          R2T=R*R*T**2.5
          RT=R*T**1.5
         
          AH2OM=166.8-.19308*TCEL+.1864D-3*TCEL*TCEL-.71288D-7*TCEL**3
          
        IF(TCEL.LT.6.D2) AH2OM=4.221D3-3.1227D1*TCEL+8.7485D-2*TCEL**2
     1-1.07295D-4*TCEL**3+4.86111D-8*TCEL**4
        IF(TCEL.GT.1200) AH2OM=140.-0.050*TCEL
          
          AH2OM=AH2OM*10.D5
          ACO2M=73.03-0.0714*TCEL+2.157D-5*TCEL*TCEL
          ACO2M=ACO2M*10.D5

          XK=EXP(-11.071+(5953./T)-(2.746D6/(T*T))+(4.646D8/(T*T*T)))
          CO2H2O=XK*0.5*R2T
          CO2H2O=CO2H2O+SQRT(A(1)*A(5))
          
          ASUM=0.0
          BSUM=0.0


          DO 101 I=1,MIXNUM
          BSUM=BSUM+B(I)*Y(I)
          DO 101 J=1,MIXNUM
          IF(I.EQ.J) GO TO 140
          
          IF(I.EQ.5.AND.J.EQ.1) GO TO 150
          IF(I.EQ.1.AND.J.EQ.5) GO TO 150
          
          ASUM=ASUM+Y(I)*Y(J)*SQRT(A(I)*A(J))
          GO TO 101
  140 CONTINUE
          IF(I.NE.5) GO TO 141
          ASUM=ASUM+Y(I)*Y(J)*AH2OM
          
          GO TO 101
  141 IF(I.NE.1) GO TO 142
          ASUM=ASUM+Y(I)*Y(J)*ACO2M
          
          GO TO 101
  142 ASUM=ASUM+Y(I)*Y(J)*A(I)
                
          GO TO 101
  150 ASUM=ASUM+Y(I)*Y(J)*CO2H2O
                
  101 CONTINUE
          ASUM=ASUM/1.013
C       WRITE(9,349) ASUM
C 349   FORMAT('0   THE REDLICH-KWONG "A" TERM IS ',F20.5,/)
C       WRITE(9,348) BSUM
C 348   FORMAT('0   THE REDLICH-KWONG "B" TERM IS ',F20.5,/)
C       PAUSE
C       write(9,*)tcel
C       pause
          RETURN
          END
            END
  