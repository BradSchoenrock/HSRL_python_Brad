      SUBROUTINE brill_s6(WAVELEN,TEMP,PRES,FREQS,NUMSPEC,SPECT)
C      SUBROUTINE brill_s6(WAVENUM,TEMP,PRES,CENTSHIFT,DSIGMA,
C     & NUMSPEC,SPECT,SIZESP)
C     SUBROUTINE GETRB(GRID,NDWIDTH,RBVALUE)
C     THIS SUBROUTINE DETERMINES THE VALUE OF THE RB L@HAPE FUNC'NON WHICH
C     IS CALCULATED BY TENTI'S SUBROUTINE.  VALUES OF THE RB PROFELE ARE
C     TABULATED FOR FREQUENCIES FROM ZERO TO TOPLIM, IN INCREMENTS OF GRID
C
C     WAVELEN   = wavelength in meters ewe 10/28/03
C     TEMP      = temperature in deg-Kelvin
C     PRES      = pressure in mb
C     FREQS     = frequencies in Hz

C       WIDTHRB = EFFECTIVE WIDTH OF THE RB PROFII-E (GHZ)
C       Y       = Y PARAMETER USED IN CALCULATING RB PROFILE 
C       QMAG    = MAGNITUDE OF SCATTERING WAVEVECTOR OR K 
C       V0 	= VELOCITY PARAMETER
C        


C     I have attempted to express all quantities in MKS units
C     orignal code uses a mixture of cgs and MKS (ewe 10/29/03)



      REAL PI,QMAG,Y,V0,ETALAM,ETAZET,CV,X
      REAL*8 SNHALF,PRES,TEMP,MASS
C      REAL*8 KB,DSIGMA,WAVENUM,SPDC,XP
C     new ewe 10/29/03      
      REAL*8 KB,WAVELEN,SPDC,XP

      CHARACTER*32 FNAME
      INTEGER I,NUMSPEC
      REAL*8 SPECT(NUMSPEC)
      REAL*8 FREQS(NUMSPEC)
      REAL S
   
      PI=3.1415927
      SPDC=2.998E8      
      ETALAM=0.198
      ETAZET=1.407
      CV=1.0
      FNAME="S6_model_spect.out"
C     Boltzman constant J/K
C     KB = 1.38062E-20    changed to MKS ewe 10/29/03
      KB=1.38044E-23

C      WAVE=1./WAVENUM     removed ewe 10/29/03
 
C     molecule mass for dry air 
C      MASS = 28.966/6.02297E23  changed from grams to kg ewe 10/29/03
      MASS = 28.966E-3/6.02297E23
C 
C SNHALF = SIN(theta/2) Where theta is the scattering angle (180)
      SNHALF = 1.0
C                    *****
      AXB=(TEMP+110.4)/(TEMP*TEMP)
C      ATM=PRES/1013250  now uses press in mb ewe 10/29/03
      ATM=PRES/1013.25

C     Y is supposed to be the ratio of mean-free-path to wavelength
      Y=0.2308*AXB*(ATM/SNHALF)*WAVELEN*1.0E9

C     notice that if wavelength is given in m this is just:
C     QMAG=4*pi/wavelength  with units of 1/m     ewe 10/29/03
C     QMAG=1.2566E10*SNHALF/(WAVE/1.0E-9) new version follows ewe 10/29.03
      QMAG=4*PI*SNHALF/WAVELEN

C     V0 remains dimensionless after conversions units to MKS  ewe 10/29/03
      V0 = SQRT(KB*TEMP/MASS)
      XP = 2*PI*SPDC/(SQRT(2.0)*QMAG*V0)

C        WRITE(*,*)
C	WRITE(*,1111) Y,QMAG,V0,XP,NUMSPEC
C 1111	FORMAT('Y, QMAG, V0, XP, NUMSPEC',4E14.3,I6)
C	WRITE(*,1111) FREQS(0),FREQS(1)
C 1111	FORMAT('freq0 freq1',4E14.3,I6)
      


C
C   Open file to write data
C
C      OPEN(UNIT=40,FILE=FNAME,STATUS='REPLACE')
C 109  FORMAT('% STEP (pm):  ',F10.3)
C      WRITE(40,109) DSIGMA*(WAVE**2)*1e12
C 110  FORMAT('% TEMP (K):   ',F10.3)
C      WRITE(40,110) TEMP
C 111  FORMAT('% PRES (ATM): ',F10.3)
C      WRITE(40,111) ATM
C 112  FORMAT('% Y:          ',F10.3)
C      WRITE(40,112) Y
C 113  FORMAT('% K or QMAG   ',E14.3)
C      WRITE(40,113) QMAG
C 114  FORMAT('% V0:         ',F10.3)
C      WRITE(40,114) V0
C 116  FORMAT('% SIZESP:    ',I5)
C      WRITE(40,116) NUMSPEC 
C 117  FORMAT('% DX:       ',F10.3)
C      WRITE(40,117) XP*DSIGMA
C
C 221  FORMAT('%    X       SPECTRUM(I)     I')
c 222  FORMAT(F7.3,2X,F12.3,2X,I5)
C      WRITE(40,221)
C
      S=0
      DO 300 I=1,NUMSPEC
C         X=XP*(-CENTSHIFT+DSIGMA*(I-(NUMSPEC/2.0+1))) new version follows
C        ewe 10?29/03
         X=XP*FREQS(I)/SPDC

         CALL TENT(ABS(X),CV,ETALAM,ETAZET,Y,S)        

C   FOLLOWING LINE ADDED BY J. FORKEY TO NORMALIZE S SO THAT INTEGRATING
C   ACROSS LINESHAPE IN GHZ WILL YIELD I
	 SPECT(I)=((2.*PI*3E8/SQRT(2.))/(QMAG*V0))*S
c         WRITE(40,222) X,S,I 
 300  CONTINUE
      END

C Computer Model Of Filtered Rayleigh Scattering Signals
C The following FORTRAN program, called FRSMODEL, evaluates
C equations 2.17, 2.18, 
C 3.3, and 3.4, given the appropriate input information.
C The subroutine TENT, which is 
C used to calculate the Rayleigh - Brillouin lineshape, 
C was provided by Professor G. Tenti 
C of the University of Waterloo.  Calculations in other subroutines
C  are based on expressions 
C in chapters 2 and 3, and in appendices A and B of this dissertation.

	SUBROUTINE TENT(X,CV,ETALAM,ETAZET,Y,S)

	COMPLEX* 16 INTG(10),ZERO
        REAL X,CV,ETALAM,ETAZET,Y,S
	NY=1
	NX=1
	I6=1
C       
	PI=3.1415927
C       
	ZERO=(0.D0,0.D0)
        DO 3 I=1,10
          INTG(I)=ZERO
  3     CONTINUE 
	G11=5./6.-(5./9.)*((CV/(1.5+CV))**2)*ETAZET
	GP11=SQRT(5./(2.*CV))*ETAZET*((CV/(1.5+CV))**2)/3.
	GP11=-GP11
	G10=1.5-(2./3.)*ETAZET*(CV/(1.5+CV))**2
	GP10=ETAZET*SQRT(2.*(CV**3)/3.)/(1.5+CV)**2
	GP10=-GP10
	G2=1.5-ETAZET*CV/(1.5+CV)**2
	G011=.4*(1.5+CV)**2+CV*(1.+CV/3.)*ETAZET
     &		+((CV/(1.5+CV))*ETAZET)**2/(6.*ETALAM)
	DEN=4./15.-ETALAM+(2./9.)*ETAZET*(CV/(1.5+CV))**2
	G011=G011/DEN
	G011=G011*2.*CV*ETALAM/(3.*(1.5+CV)**2)
	G011=1.5-G011
C       
	IY=1
	DTBAR=1./(ETALAM*(2.5+CV)*2.*Y)
	GBAR=(4./3.+1./ETAZET)/(2.*Y)+DTBAR/(1.5+CV)
C	Y2=1.5*Y
C	Y3=H*Y
	S6=0.
	CALL IMN(INTG,6,X,Y)
	CALL S6GET(Y,G11,GP11,G10,GP10,G2,G011,INTG,0,S6,S6E)
        S=S6
	RETURN
	END
c       Subroutines for Tenti Program
c       from G. Tenti, University of Waterloo, Sept. 7, 1989
c       from D. Seasholtz, NASA Lewis, March, 1994
c       calls from J. Forkey's programs
c       CALL TENT(x,cv,etalam,etazet,s)
c       cv=internal specific/Boltzmann's constant
c       =0 for monoatomic gases
c       =0 for diatomic gases where rotational degrees of
c	of freedom have fully kicked in (JL)
c	etalam	reciprocal of Eucker ratio = viscosity*kb/thenn cond
c	etazet	ratio of shear viscosity to bulk viscosity
c	(bulk viscosity = 0 for monoatomic gases, so
c	let etazet = 1000 for computational purposes)
c	(for diatomic gas, the bulk viscosity depends on
c	the Cv, the rotational relaxation time, and is
c	linear in pressure.  For N2 at I atm, published
c	values are between 1.367 and 1.407,( JL))
c       .. revised for MS Fortran by J. Lock
c       .. revision date 12/9/90
c       revised for MS Fortran by Jlock

C       SUB1.FOR
C       
	SUBROUTINE IMN(INTG,NMAX,X,Y)
	COMPLEX*16 Z,INTG(10),W,UNI,ZLAM,SQ2,ACPL
	REAL*8 RSQ2,XX,YY,A
	UNI=(0.0D+0,1.0D+0)
	RSQ2=DSQRT(2.0D+0)
	SQ2=DCMPLX(RSQ2,0.0D+0)
	KMAX=NMAX+1
	XX=DBLE(X)
	YY=DBLE(Y)
	Z=DCMPLX(XX,YY)
	ZLAM=UNI*SQ2*Z
	K=1
	INTG(1)=UNI*W(Z)/SQ2
	AK=1.0
	A=DBLE(1.0)
	ICG=2
	GO TO 30
 5	GO TO (10,20),ICG
 10	AK=AK*(K-2)
	A=AK
	ICG=2
	GO TO 30
 20	A=0.0
	ICG=1
 30	ACPL=DCMPLX(A,0.0D+0)
	INTG(K+1)=(ACPL-ZLAM*INTG(K))*UNI
	K=K+1
	IF(K.LE.KMAX) GO TO 5
	RETURN
	END
C       SUB2.FOR
C
	COMPLEX*16 FUNCTION W(Z)
	COMPLEX*16 Z,ZSQ,SRIES,TERM,TERM1,A0,A1,B0,B1,A,B,AN1,BN1
	COMPLEX*16 W1,W2,SUM,SQRPI
	REAL*8 PI,AA
	COMPLEX CZ,CZ1
	CZ=Z
	X=REAL(CZ)
	Y=AIMAG(CZ)
	PI=3.141592653589793
	PI=DSQRT(PI)
	SQRPI=DCMPLX(0.0D0,PI)
	IF((X.GE.0.0).AND.(Y.GE.0.0)) GO TO 1
	WRITE(*,1001)
 1001	FORMAT(' X AND Y ARE NOT IN THE FIRST QUADRANT')
	STOP
 1	IF((X.LE.5.0).AND.(Y.LE.1.5)) GO TO 50
	IF(Y.GE.1.5) GO TO 30
C       
	N=1
	ZSQ=Z*Z
	TERM=5.0D-1/ZSQ
	SRIES=DCMPLX(1.0D0,0.0D0)
 5	N=N+1
	SRIES=TERM+SRIES
	TERM1=TERM/2.0/ZSQ*(2*N-1)
        CZ=TERM
        CZ1=TERM1
	XY=ABS(REAL(CZ))+ABS(AIMAG(CZ))
	XY1=ABS(REAL(CZ1))+ABS(AIMAG(CZ1))
	IF((XY.LT.1.0E-15).OR.(XY1.GT.XY)) GO TO 6
	TERM=TERM1
	GO TO 5
 6	W1=-(SRIES)/Z
c	WRITE(*,1002) N,XY,XY1
c 1002	FORMAT('#ASYMPTOTIC SERIES--N.XY,XY1',3X,I4,2X,2E16.5)
	W=W1
	IF(Y.NE.0.0) RETURN
c       W=SQRPI*CDEXP(-ZXQ)+W1 'changed 12/9/90
	w=sqrpi*cdexp(-zsq)+w1
	RETURN
C       
 30	N=0
	A0=DCMPLX(1.0D0,0.0D0)
	A1=DCMPLX(0.0D0,0.0D0)
	B0=A1
	B1=A0
	AN1=Z
	BN1=-Z*Z+DCMPLX(5.0D-1,0.0D0)
	ICG=1
 35	A=BN1*A1+AN1*A0
	B=BN1*B1+AN1*B0
	GO TO(31,34),ICG
 31	AN1=DCMPLX(0.0D0,0.0D0)
	W1=A/B
	ICG=2
	GO TO 33
 34	W2=A/B
	W=W2
        CZ=W2-W1
	XY=ABS(REAL(CZ))+ABS(AIMAG(CZ))
        IF(XY.LE.1.0E-15) RETURN
	W1=W2
 33	A0=A1
	B0=B1
	A1=A
	B1=B
	BN1=BN1+DCMPLX(2.0D0,0.0D0)
	AA=-2.0*N-0.5
	AN1=AN1+DCMPLX(AA,0.0D0)
	N=N+1
	IF(N.LT.80) GO TO 35
	WRITE(6,1003) N,XY
 1003	FORMAT('CONTINUED FRACTION--N,XY',3X,I4,2X,E16.5)
	RETURN
C
 50	N=1
	ZSQ=Z*Z*DCMPLX(-2.0D0,0.0D0)
	SUM=DCMPLX(1.0D0,0.0D0)
	TERM=ZSQ/(3.0D0)
 55	SUM=SUM+TERM
	N=N+1
	TERM=TERM*ZSQ/(2*N+1)
	W1=TERM/SUM
        CZ=W1
	XY=ABS(REAL(CZ))+ABS(AIMAG(CZ))
	IF(XY.LT.1.0E-15) GO TO 60
	IF(N.LT.3000) GO TO 55
	WRITE(*,1004) N,XY
 1004	FORMAT('POWER SERIES--N,XY',3X,I4,2X,E16.5)
 60	SUM=SUM*(-2.0)*Z
	W=SQRPI*CDEXP(ZSQ/2.0)+SUM
	RETURN
	END
C       SUB3.FOR
c       
	SUBROUTINE S6GET(Y1,G11,GP11,G10,GP10,G2,G011,INTG,IE,S,SE)
        COMPLEX* 16 INTG(10),GAM(5,5),B(6,6),C(6),BE(4,4),CE(4)
        COMPLEX* 16 UNIT,ONE,ZERO
	REAL*8 BB(12,12),CC(12),BBE(8,8),CCE(8)
        INTEGER KS
        KS=0
	UNIT=(0.D0,1.D0)
	ONE=(1.D0,0.D0)
	ZERO=(0.D0,0.D0)
	PI=3.14159265
        do 4 I=1,5
         do 3 J=1,5
           GAM(I,J)=ZERO
  3      CONTINUE 
  4     CONTINUE
        do 6 I=1,8
         do 5 J=1,8
           BBE(I,J)=0.0
  5      CONTINUE 
           CCE(I)=0.0
  6     CONTINUE
        do 8 I=1,12
         do 7 J=1,12
           BB(I,J)=0.0
  7      CONTINUE 
           CC(I)=0.0
  8     CONTINUE
        CALL GAMMA(GAM,INTG,4)
	DO 10 I=1,4
	   IF(I.LE.3) K=I
	   IF(I.EQ.4) K=5
	   C(I)=GAM(K,1)
	   B(I,1)=-Y1*GAM(K,1)
	   B(I,2)=-Y1*GAM(K,2)
	   B(I,3)=-Y1*(G11-.5)*GAM(K,3)
	   B(I,4)=-Y1*(G10-.5)*GAM(K,5)
	   B(I,5)=-Y1*GP10*GAM(K,5)
	   B(I,6)=-Y1*GP11*GAM(K,3)
	   IF(IE.EQ.0) GO TO 10
	   CE(I)=C(I)
	   BE(I,1)=B(I,1)
	   BE(I,2)=B(I,2)
	   BE(I,3)=-Y1*GAM(K,3)/3.
	   BE(I,4)=-Y1*GAM(K,5)
 10	CONTINUE
	DO 14 I=1,4
	   DO 15 J=1,4
	      IF(I.NE.J) GO TO 15
	      B(I,J)=B(I,J)-ONE
	      IF(IE.EQ.1) BE(I,J)=BE(I,J)-ONE
 15	   CONTINUE
 14	CONTINUE
	C(5)=ZERO
	C(6)=ZERO
	B(5,1)=ZERO
	B(5,2)=ZERO
	B(5,3)=-Y1*GP11*GAM(1,2)
	B(5,4)=-Y1*GP10*GAM(1,1)
	B(5,5)=-Y1*(G2-.5)*GAM(1,1)-ONE
	B(5,6)=-Y1*(G011-.5)*GAM(1,2)
	B(6,1)=ZERO
	B(6,2)=ZERO
	B(6,3)=-Y1*GP11*GAM(2,2)
	B(6,4)=-Y1*GP10*GAM(2,1)
	B(6,5)=-Y1*(G2-.5)*GAM(2,1)
	B(6,6)=-Y1*(G011-.5)*GAM(2,2)-ONE
	DO 17 I=1,6
	   CC(I)=C(I)
	   CC(I+6)=-UNIT*C(I)
	   IF((I.GT.4).OR.(IE.EQ.0)) GO TO 19
           CCE(I)=CC(I)
c       CCE(TT+4)=-UNIT*CE(I) 'changed 12/9/90
           cce(I+4)=-unit*ce(I)
 19	   DO 18 J=1,6
	      BB(I,J)=B(I,J)
	      BB(I+6,J+6)=BB(I,J)
	      BB(I,J+6)=UNIT*B(I,J)
	      BB(I+6,J)=-BB(I,J+6)
	      IF((I.GT.4).OR.(J.GT.4).OR.(IE.EQ.0)) GO TO 18
	      BBE(I,J)=BE(I,J)
	      BBE(I+4,J+4)=BBE(I,J)
	      BBE(I,J+4)=UNIT*BE(I,J)
	      BBE(I+4,J)=-BBE(I,J+4)
 18	   CONTINUE
 17	CONTINUE
        CALL SIMQ(BB,CC,12,KS,12)
	IF(KS.EQ. 1) WRITE(1, 1005)
 1005	FORMAT(' MATRIX BB6 IS SINGULAR')
	S=CC(1)/PI
	IF(IE.EQ. 1) GO TO 32
	SE=0.
	RETURN
 32	CALL SIMQ(BBE,CCE,8,KS,8)
	IF(KS.EQ.1) WRITE(*,1006)
 1006	FORMAT('MATRIX BBE6 IS SINGULAR')
	SE=CCE(1)/PI
	RETURN
	END
C       SUB4.FOR
c       
	SUBROUTINE GAMMA(GAM,INTG,INDEX)
        COMPLEX* 16 INTG(10),INTMN,GAMA,GAM(5,5)
	REAL*8 CRL,BARNES,BIGG
C	INTEGER IRME(5)/0,0,1,0,1/
	INTEGER IRME(5)
C	INTEGER LME(5)/0,1,1,2,0/
	INTEGER LME(5)
	DO 100 IME=1,5
	 IF(IME.EQ.INDEX) GO TO 100
	 IR=IRME(IME)
	 L=LME(IME)
	 DO 200 JME=1,5
	  IF(JME.EQ.INDEX) GO TO 200
	   IRT=IRME(JME)
	   LT=LME(JME)
 10	   CRL=DBLE(1.0)/DBLE(2.0**(1.5*(L+LT)+1.0))
	   CRL=CRL*(BARNES(1.0,2*L)/BARNES(1.0,L)**2)*(BARNES(1.0,2*LT)
     &	     /BARNES(1.0,LT)**2)/DSQRT(BARNES(1.0,IR)*BARNES(1.0,IRT))
 15	   CRL=CRL*DSQRT((2*L+1)*(2*LT+1)/BARNES(0.5,L+IR+1)
     &       /BARNES(0.5,LT+IRT+1))
	   KM=L/2+1
 20	   KMT=LT/2+1
	   JM=IR+1
	   JMT=IRT+1
	   GAMA=0.0
 30	   DO 70 J= 1,JM
	    DO 70 JT=1,JMT
	     DO 70 K=1,KM
	      DO 70 KT=1,KMT
               BIGG=DBLE((-1.0)**(J+JT+K+KT))
     &	          /DBLE(2.0**(2*(K+KT)-4+JM-J+JMT-JT))
 40	       BIGG=BIGG*BARNES(JM-J+1.0,J-1)/BARNES(1.0,J-1)
     &           *BARNES(JMT-JT+1.0,JT-1)/BARNES(1.0,JT-1)
	       BIGG=BIGG*BARNES(IRT+LT+2.5-JT,JT-1)/BARNES(LT+1.5-KT,KT-1)
     &	         *BARNES(IR+L+2.5-J,J-1)/BARNES(L+1.5-K,K-1)
	       BIGG=BIGG*BARNES(L-2*K+3.0,2*(K-1))/BARNES(1.0,K-1)
     &           *BARNES(LT-2*KT+3.0,2*(KT-1))/BARNES(1.0,KT-1)
 42	       GAMA=BIGG*INTMN(IR+IRT-J-JT+K+KT,L+LT-2*(K+KT)+4,INTG)+GAMA
 70	   CONTINUE
 90	   GAMA=GAMA*CRL*DSQRT(2.0D0)
	   GAM(IME,JME)=GAMA
 200	  CONTINUE
 100	CONTINUE
       RETURN
       END
C SUB5.FOR
c
	REAL*8 FUNCTION BARNES(B,M)
	BARNES=1.0D0
	IF(M.EQ.0) RETURN
	DO 10 N=1,M
	   BARNES=BARNES*(B+DFLOAT(N-1))
 10	CONTINUE
	RETURN
	END
C SUB6.FOR
c
        COMPLEX*16 FUNCTION INTMN(M,N,INTG)
        COMPLEX* 16 INTG(10)
        REAL*8 SUBS,BARNES
 1	IP=0
	INTMN=(0.0D+0,0.0D+0)
 5	IQMX=M-IP
 7	SUBS=0.0D+0
 10	IQ=0
 15	SUBS=BARNES(0.5,IQ)/BARNES(1.0,IQ)*BARNES(0.5,M-IP-IQ)/
     &   BARNES(1.0,M-IP-IQ)+SUBS
	IQ=IQ+1
	IF((IQ-1).LT.IQMX) GO TO 15
	INTMN=SUBS*INTG(N+2*IP+1)*(2.0D+0**(M-IP))/BARNES(1.0,IP)+INTMN
	IP=IP+1
	IF((IP-1).LT.M) GO TO 5
        INTMN=INTMN*BARNES(1.0,M)
 30	RETURN
	END
C SUB7.FOR
c
	SUBROUTINE SIMQ(A,B,N,KS,M)
C       REAL*8 A(M,M),B(l),TOL,BIGA,SAVE
	REAL*8 A(M,M),B(N),TOL,BIGA,SAVE
	TOL=0.0
	KS=0
	DO 65 J=1,N
	   JY=J+1
	   BIGA=0
	   DO 30 I=J,N
c       
	      IF(DABS(BIGA)-DABS(A(I,J))) 20,30,30
 20	      BIGA=A(I,J)
	      IMAX=I
 30	   CONTINUE
c
	   IF(DABS(BIGA)-TOL) 35,35,40
 35	   KS=1
	RETURN
c
 40	DO 50 K=J,N
	   SAVE=A(J,K)
	   A(J,K)=A(IMAX,K)
	   A(IMAX,K)=SAVE
c
 50	   A(J,K)=A(J,K)/BIGA
	   SAVE=B(IMAX)
	   B(IMAX)=B(J)
	   B(J)=SAVE/BIGA
c       
	   IF(J-N) 55,70,55
 55	   DO 65 IX=JY,N
	      IT=J-IX
	   DO 60 JX=JY,N
 60	      A(IX,JX)=A(IX,JX)-(A(IX,IX+IT)*A(IX+IT,JX))
 65	      B(IX)=B(IX)-(B(J)*A(IX,IX+IT))
c       
 70	      NY=N-1
	      IT=N*N
	      DO 80 J=1,NY
		 IB=N-J
		 IC=N
		 DO 80 K=1,J
		    B(IB)=B(IB)-A(IB,IC)*B(IC)
 80		    IC=IC-1
        RETURN
	END
