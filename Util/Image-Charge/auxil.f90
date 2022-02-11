DOUBLE PRECISION FUNCTION v(u)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC        AUXILIARY FUNCTION
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC         LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! COSTH    : COS(THET)
! FDOMIN   : X*(COS(THET)+THET*SIN(THET))
! SINTH    : SIN(THET)
! THET     : ASIN(PNU/X)
! U        : ARGUMENT OF THE V FUNCTION
! UNDER    : UNDERFLOW NUMBER
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)         :: u
DOUBLE PRECISION :: d1mach,under
DOUBLE PRECISION :: thet,sinth,costh,fdomin
COMMON/parmon/thet,sinth,costh,fdomin
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+6
IF (ABS(u) < under) THEN
  v=ASIN(sinth)
ELSE
  v=ASIN(u/SINH(u)*sinth)
END IF
RETURN
END FUNCTION v

DOUBLE PRECISION FUNCTION phir(u)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC      AUXILIARY FUNCTION:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC       LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! COSTH   : COS(THET)
! COSVU   : COS(V(U))
! FDOMIN  : X*(COS(THET)+THET*SIN(THET))
! FU1     : 2*SINH(U/2)**2
! FU2     : -2*SIN(FU3/2)*SIN(0.5*(V(U)+THET))
! FU3     : ASIN(-SIN(THET)/(COS(THET)*U/SINH(U)
!           +COS(V(U))*(SINH(U)+U)*FUAC/(SINH(U)*SINH(U)))
! FUAC    : U**3/6+U**5/120+U**7/5040
! OVER    : OVERFLOW NUMBER
! PNU     : ORDER OF THE FUNCTION
! SINHU   : SINH(U)
! SINHUH  : SINH(U/2)
! SINTH   : SIN(THET)
! THET    : ASIN(PNU/X)
! U       : ARGUMENT OF THE PHIR FUNCTION
! U2      : U**2
! U3      : U**3
! U5      : U**5
! U7      : U**7
! UH      : U/2
! UNDER   : UNDERFLOW NUMBER
! V(U)    : ASIN(U/SINH(U)*SIN(THET))
! VU      : V(U)
! X       : ARGUMENT OF THE FUNCTIONS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)             :: u
DOUBLE PRECISION :: under,over,d1mach,v,uh,u2,u3,  &
    u5,u7,fu1,fu2,fu3,fuac,sinhu,sinhuh,vu,cosvu
DOUBLE PRECISION :: x,pnu,thet,sinth,costh,fdomin
COMMON/argu/x,pnu
COMMON/parmon/thet,sinth,costh,fdomin
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+6
!CC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
over=d1mach(2)*1.d-6
IF (u < 200) THEN
  IF (ABS(u) < 0.1D0) THEN
    IF (ABS(u) < under) THEN
      phir=0.d0
    ELSE
      uh=u*0.5D0
      sinhuh=SINH(uh)
      fu1=2.d0*sinhuh*sinhuh
      u2=u*u
      u3=u2*u
      u5=u3*u2
      u7=u5*u2
      fuac=u3/6.d0+u5/120.d0+u7/5040.d0
      sinhu=SINH(u)
      vu=v(u)
      cosvu=COS(vu)
      fu3=ASIN(-sinth/(costh*u/sinhu+cosvu)* (sinhu+u)*fuac/(sinhu*sinhu))
      fu2=-2.d0*SIN(fu3*0.5D0)*SIN(0.5D0*(vu+thet))
      phir=x*(fu1*cosvu+fu2+sinth*fu3)
    END IF
  ELSE
    vu=v(u)
    cosvu=COS(vu)
    phir=x*COSH(u)*cosvu+pnu*vu-fdomin
  END IF
ELSE
  phir=over
END IF
RETURN
END FUNCTION phir

DOUBLE PRECISION FUNCTION sigma(u)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC     AUXILIARY FUNCTION
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC     LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ARGU   : (COSH(MU)*U-DMUFAC)/SINH(U)
! COSCHI : COS(CHI)
! COSHM  : COSH(MU)
! D1     : COSH(MU)*U-DMUFAC
! DMU    : SOLUTION MU OF COSH(MU)=PNU/X
! DMUFAC : MU*COSH(MU)-SINH(MU)
! SINCHI : SIN(CHI)
! SINHM  : SINH(MU)
! SINHU  : SINH(U)
! U      : ARGUMENT OF THE SIGMA FUNCTION
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)             :: u
DOUBLE PRECISION :: sinhu,argu,d1,pi
DOUBLE PRECISION :: dmu,coshm,sinhm,dmufac,coschi,sinchi
COMMON/paros1/dmu,coshm,sinhm,dmufac,coschi,sinchi

pi=ACOS(-1.d0)
d1=coshm*u-dmufac
sinhu=SINH(u)
argu=d1/SINH(u)
IF (ABS(argu) > 1) THEN
  argu=1.d0
END IF
sigma=ASIN(argu)
IF (u < dmu) sigma=pi-sigma
RETURN
END FUNCTION sigma

DOUBLE PRECISION FUNCTION phib(u)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC         AUXILIARY FUNCTION
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC         LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! OVER      : OVERFLOW NUMBER
! PNU       : ORDER OF THE FUNCTIONS
! SIGMA(U)  : ASIN((COSH(MU)*U-DMUFAC)/SINH(U))
! SIGMAU    : SIGMA(U)
! U         : ARGUMENT OF THE PHIB FUNCTION
! X         : ARGUMENT OF THE FUNCTIONS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN OUT)         :: u
DOUBLE PRECISION :: sigma,sigmau,over,d1mach
DOUBLE PRECISION :: x,pnu
COMMON/argu/x,pnu
!CC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
over=d1mach(2)*1.d-8
IF (u < 200) THEN
  sigmau=sigma(u)
  phib=x*COSH(u)*COS(sigmau)+pnu*sigmau
ELSE
  phib=over
END IF
RETURN
END FUNCTION phib




SUBROUTINE aiz(ifun,ifac,x0,y0,gair,gaii,ierro)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! COMPUTATION OF THE AIRY FUNCTION AI(Z) OR ITS DERIVATIVE AI'(Z)
! THE CODE USES:
!      1. MACLAURIN SERIES FOR |Y|<3 AND -2.5<X<1.3 (Z=X+I*Y)
!      2. GAUSS-LAGUERRE QUADRATURE  FOR |Z|<15 AND  WHEN
!         MACLAURIN SERIES ARE NOT USED.
!      3. ASYMPTOTIC EXPANSION FOR |Z|>15.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  INPUTS:
!    IFUN:
!         * IFUN=1, THE CODE COMPUTES AI(Z)
!         * IFUN=2, THE CODE COMPUTES AI'(Z)
!    IFAC:
!         * IFAC=1, THE CODE COMPUTES  AI(Z) OR AI'(Z)
!         * IFAC=2, THE CODE COMPUTES NORMALIZED AI(Z) OR AI'(Z)
!    X0:   REAL PART OF THE ARGUMENT Z
!    Y0:   IMAGINARY PART OF THE ARGUMENT  Z

!  OUTPUTS:
!    GAIR: REAL PART OF AI(Z) OR AI'(Z)
!    GAII: IMAGINARY PART OF AI(Z) OR AI'(Z)

!    IERRO: ERROR FLAG
!          * IERRO=0, SUCCESSFUL COMPUTATION
!          * IERRO=1, COMPUTATION OUT OF RANGE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          MACHINE DEPENDENT CONSTANTS: FUNCTION D1MACH
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   ACCURACY:

!     1) SCALED AIRY FUNCTIONS:
!        RELATIVE ACCURACY BETTER THAN 10**(-13) EXCEPT CLOSE TO
!        THE ZEROS, WHERE 10**(-13) IS THE ABSOLUTE PRECISION.
!        GRADUAL LOSS OF PRECISION TAKES PLACE FOR |Z|>1000
!        (REACHING 10**(-8) ABSOLUTE ACCURACY FOR |Z| CLOSE
!        TO 10**(6)) IN THE CASE OF PHASE(Z) CLOSE TO PI.
!     2) UNSCALED AIRY FUNCTIONS:
!        THE FUNCTION OVERFLOWS/UNDERFLOWS FOR
!        3/2*|Z|**(3/2)>LOG(OVER).
!        FOR |Z|<30:
!        A) RELATIVE ACCURACY FOR THE MODULUS (EXCEPT AT THE
!           ZEROS) BETTER THAN 10**(-13).
!        B) ABSOLUTE ACCURACY FOR MIN(R(Z),1/R(Z)) BETTER
!           THAN 10**(-13), WHERE R(Z)=REAL(AI)/IMAG(AI)
!           OR R(Z)=REAL(AI')/IMAG(AI').
!        FOR |Z|>30, GRADUAL LOSS OF PRECISION TAKES PLACE
!        AS |Z| INCREASES.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     AUTHORS:
!        AMPARO GIL    (U. AUTONOMA DE MADRID, MADRID, SPAIN).
!                      E-MAIL: AMPARO.GIL@UAM.ES
!        JAVIER SEGURA (U. CARLOS III DE MADRID, MADRID, SPAIN).
!                      E-MAIL: JSEGURA@MATH.UC3M.ES
!        NICO M. TEMME (CWI, AMSTERDAM, THE NETHERLANDS).
!                      E-MAIL: NICO.TEMME@CWI.NL
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    REFERENCES:
!         COMPUTING AIRY FUNCTIONS BY NUMERICAL QUADRATURE.
!         NUMERICAL ALGORITHMS (2001).
!         A. GIL, J. SEGURA, N.M. TEMME
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN OUT)                  :: ifun
INTEGER, INTENT(IN OUT)                  :: ifac
DOUBLE PRECISION, INTENT(IN)             :: x0
DOUBLE PRECISION, INTENT(IN)             :: y0
DOUBLE PRECISION, INTENT(OUT)            :: gair
DOUBLE PRECISION, INTENT(OUT)            :: gaii
INTEGER, INTENT(OUT)                     :: ierro
DOUBLE PRECISION :: x,w,xd,wd
DOUBLE PRECISION :: over,under,dl1,dl2,cover,d1mach
DOUBLE PRECISION :: pi,pihal,pih3,pisr,a,alf,thet,r,th15,s1,c1,  &
    r32,facto,th025,s3,c3,f23,pi23,sqrt3,xa,ya,f23r,df1,df2,  &
    s11,c11,dex,dre,dima,gar,gai,c,s,u,v,v0,ar,ai,ar1,ai1,  &
    ro,coe1,coe2,rex,dfr,dfi,ar11,ai11,phase
INTEGER :: iexpf,iexpf2,n
DIMENSION x(25),w(25)
DIMENSION xd(25),wd(25)
COMMON/param1/pi,pihal
COMMON/param2/pih3,pisr,a,alf
COMMON/param3/thet,r,th15,s1,c1,r32
COMMON/param4/facto,th025,s3,c3
SAVE x,w
SAVE xd,wd
DATA x,w/.283891417994567679D-1,.170985378860034935D0,  &
    .435871678341770460D0,.823518257913030858D0,1.33452543254227372D0,  &
    1.96968293206435071D0,2.72998134002859938D0,3.61662161916100897D0,  &
    4.63102611052654146D0,5.77485171830547694D0,7.05000568630218682D0,  &
    8.45866437513237792D0,10.0032955242749393D0,11.6866845947722423D0,  &
    13.5119659344693551D0,15.4826596959377140D0,17.6027156808069112D0,  &
    19.8765656022785451D0,22.3091856773962780D0,24.9061720212974207D0,  &
    27.6738320739497190D0,30.6192963295084111D0,33.7506560850239946D0,  &
    37.0771349708391198D0,40.6093049694341322D0,.143720408803313866D0,  &
    .230407559241880881D0,.242253045521327626D0,.203636639103440807D0,  &
    .143760630622921410D0,.869128834706078120D-1,.454175001832915883D-1,&
    .206118031206069497D-1,.814278821268606972D-2,.280266075663377634D-2,&
     .840337441621719716D-3,.219303732907765020D-3,  &
    .497401659009257760D-4,.978508095920717661D-5,.166542824603725563D-5,&
     .244502736801316287D-6,.308537034236207072D-7,.333296072940112245D-8,&
     .306781892316295828D-9,.239331309885375719D-10,  &
    .157294707710054952D-11,.864936011664392267D-13,.394819815638647111D-14,&
    .148271173082850884D-15,.453390377327054458D-17/
DATA xd,wd/.435079659953445D-1,.205779160144678D0,  &
    .489916161318751D0,.896390483211727D0,1.42582496737580D0,  &
    2.07903190767599D0,2.85702335104978D0,3.76102058198275D0,  &
    4.79246521225895D0,5.95303247470003D0,7.24464710774066D0,  &
    8.66950223642504D0,10.2300817341775D0,11.9291866622602D0,  &
    13.7699665302828D0,15.7559563095946D0,17.8911203751898D0,  &
    20.1799048700978D0,22.6273004064466D0,25.2389175786164D0,  &
    28.0210785229929D0,30.9809287996116D0,34.1265753192057D0,  &
    37.4672580871163D0,41.0135664833476D0,.576354557898966D-1,  &
    .139560003272262D0,.187792315011311D0,.187446935256946D0,  &
    .150716717316301D0,.101069904453380D0,.575274105486025D-1,  &
    .280625783448681D-1,.117972164134041D-1,.428701743297432D-2,  &
    .134857915232883D-2,.367337337105948D-3,.865882267841931D-4,  &
    .176391622890609D-4,.309929190938078D-5,.468479653648208D-6,  &
    .607273267228907D-7,.672514812555074D-8,.633469931761606D-9,  &
    .504938861248542D-10,.338602527895834D-11,.189738532450555D-12,  &
    .881618802142698D-14,.336676636121976D-15,.104594827170761D-16/
!C CONSTANTS CCCCCCCCCCCCCCCCCCCCCCC
pi=3.1415926535897932385D0
pihal=1.5707963267948966192D0
pih3=4.71238898038469D0
f23=.6666666666666666D0
pi23=2.09439510239320D0
pisr=1.77245385090552D0
sqrt3=1.7320508075688772935D0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ya=y0
xa=x0
ierro=0
iexpf=0
iexpf2=0
IF (ya < 0.d0) ya=-ya
r=SQRT(xa*xa+ya*ya)
r32=r*SQRT(r)
thet=phase(xa,ya)
cover=2.d0/3.d0*r32*ABS(COS(1.5D0*thet))
!CC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
over=d1mach(2)*1.d-3
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+3
dl1=LOG(over)
dl2=-LOG(under)
IF (dl1 > dl2) over=1/under
IF (ifac == 1) THEN
  IF (cover >= LOG(over)) THEN
!CC OVERFLOW/UNDERFLOW PROBLEMS.
!CC   CALCULATION ABORTED
    ierro=1
    gair=0
    gaii=0
  END IF
  IF (cover >= (LOG(over)*0.2)) iexpf2=1
ELSE
  IF (cover >= (LOG(over)*0.2)) iexpf=1
END IF
IF (ierro == 0) THEN
  IF (ifun == 1) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF AI(Z) CCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCC SERIES, INTEGRALS OR EXPANSIONS CCCCCCCCCCCCCCCCCCCCCCCC
    IF ((ya < 3.d0).AND.(xa < 1.3D0).AND.(xa > -2.5D0)) THEN
!CC SERIES CCC
      CALL serai(xa,ya,gar,gai)
      IF (ifac == 2) THEN
        thet=phase(xa,ya)
        th15=1.5D0*thet
        s1=SIN(th15)
        c1=COS(th15)
        f23r=f23*r32
        df1=f23r*c1
        df2=f23r*s1
        s11=SIN(df2)
        c11=COS(df2)
        dex=EXP(df1)
        dre=dex*c11
        dima=dex*s11
        gair=dre*gar-dima*gai
        gaii=dre*gai+dima*gar
      ELSE
        gair=gar
        gaii=gai
        IF (y0 == 0.) gaii=0.d0
      END IF
    ELSE
      IF (r > 15.d0) THEN
!CC ASYMPTOTIC EXPANSIONS CCC
        thet=phase(xa,ya)
        facto=0.5D0/pisr*r**(-0.25D0)
        IF (thet > pi23) THEN
!CCCCCCCCCC CONNECTION FORMULAE CCCCCCCCCCCCCCCCCCCCCCCCCC
!CC     N= 1: TRANSFORM Z TO W= U+IV=Z EXP( 2 PI I/3)
          n=1
          c=-0.5D0
          s=n*0.5*sqrt3
          u=xa*c-ya*s
          v=xa*s+ya*c
          v0=v
          IF (v < 0.d0) v=-v
          thet=phase(u,v)
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          th025=thet*0.25D0
          s3=SIN(th025)
          c3=COS(th025)
          CALL expai(ar1,ai1)
          IF (v0 < 0.d0) ai1=-ai1
          ar=-(c*ar1-s*ai1)
          ai=-(s*ar1+c*ai1)
          IF (iexpf == 0) THEN
            IF (iexpf2 == 0) THEN
!CC     N=-1: TRANSFORM Z TO W= U+IV=Z EXP(-2 PI I/3)
              n=-1
              c=-0.5D0
              s=n*0.5*sqrt3
              u=xa*c-ya*s
              v=xa*s+ya*c
              v0=v
              IF (v < 0.d0) v=-v
              thet=phase(u,v)
              th15=1.5D0*thet
              s1=SIN(th15)
              c1=COS(th15)
              th025=thet*0.25D0
              s3=SIN(th025)
              c3=COS(th025)
              CALL expai(ar1,ai1)
              IF (v0 < 0.d0) ai1=-ai1
              thet=phase(xa,ya)
              th15=1.5D0*thet
              s1=SIN(th15)
              c1=COS(th15)
              ro=1.333333333333333D0*r32
              coe1=ro*c1
              coe2=ro*s1
              rex=EXP(coe1)
              dfr=rex*COS(coe2)
              dfi=rex*SIN(coe2)
              ar11=dfr*ar1-dfi*ai1
              ai11=dfr*ai1+dfi*ar1
              gair=ar-(c*ar11-s*ai11)
              gaii=ai-(s*ar11+c*ai11)
            ELSE
              thet=phase(xa,ya)
              th15=1.5D0*thet
              s1=SIN(th15)
              c1=COS(th15)
              gair=ar
              gaii=ai
            END IF
          ELSE
            gair=ar
            gaii=ai
          END IF
        ELSE
!CCCCCC  ASYMPTOTIC EXPANSION CCCCCCCCCCCCCCC
          thet=phase(xa,ya)
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          th025=thet*0.25D0
          s3=SIN(th025)
          c3=COS(th025)
          CALL expai(gair,gaii)
        END IF
      ELSE
!CC INTEGRALS
        a=0.1666666666666666D0
        alf=-a
        facto=0.280514117723058D0*r**(-0.25D0)
        thet=phase(xa,ya)
        IF (thet <= pihal) THEN
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          th025=thet*0.25D0
          s3=SIN(th025)
          c3=COS(th025)
          CALL airy1(x,w,gair,gaii)
        END IF
        IF ((thet > pihal).AND.(thet <= pi23)) THEN
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          th025=thet*0.25D0
          s3=SIN(th025)
          c3=COS(th025)
          CALL airy2(x,w,gair,gaii)
        END IF
        IF (thet > pi23) THEN
          n=1
          c=-0.5D0
          s=n*0.5*sqrt3
          u=xa*c-ya*s
          v=xa*s+ya*c
          v0=v
          IF (v < 0.d0) v=-v
          thet=phase(u,v)
          IF (thet <= pihal) THEN
            th15=1.5D0*thet
            s1=SIN(th15)
            c1=COS(th15)
            th025=thet*0.25D0
            s3=SIN(th025)
            c3=COS(th025)
            CALL airy1(x,w,ar1,ai1)
          END IF
          IF ((thet > pihal).AND.(thet <= pi23)) THEN
            th15=1.5D0*thet
            s1=SIN(th15)
            c1=COS(th15)
            th025=thet*0.25D0
            s3=SIN(th025)
            c3=COS(th025)
            CALL airy2(x,w,ar1,ai1)
          END IF
          IF (v0 < 0.d0) ai1=-ai1
          ar=-(c*ar1-s*ai1)
          ai=-(s*ar1+c*ai1)
          n=-1
          c=-0.5D0
          s=n*0.5*sqrt3
          u=xa*c-ya*s
          v=xa*s+ya*c
          v0=v
          IF (v < 0.d0) v=-v
          thet=phase(u,v)
          IF (thet <= pihal) THEN
            th15=1.5D0*thet
            s1=SIN(th15)
            c1=COS(th15)
            th025=thet*0.25D0
            s3=SIN(th025)
            c3=COS(th025)
            CALL airy1(x,w,ar1,ai1)
          END IF
          IF ((thet > pihal).AND.(thet <= pi23)) THEN
            th15=1.5D0*thet
            s1=SIN(th15)
            c1=COS(th15)
            th025=thet*0.25D0
            s3=SIN(th025)
            c3=COS(th025)
            CALL airy2(x,w,ar1,ai1)
          END IF
          IF (v0 < 0.d0) ai1=-ai1
          thet=phase(xa,ya)
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          ro=1.333333333333333D0*r32
          coe1=ro*c1
          coe2=ro*s1
          rex=EXP(coe1)
          dfr=rex*COS(coe2)
          dfi=rex*SIN(coe2)
          ar11=dfr*ar1-dfi*ai1
          ai11=dfr*ai1+dfi*ar1
          gair=ar-(c*ar11-s*ai11)
          gaii=ai-(s*ar11+c*ai11)
        END IF
      END IF
      IF (ifac == 1) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF THE UNSCALED AI(Z) CCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        f23r=f23*r32
        df1=f23r*c1
        df2=f23r*s1
        s11=SIN(df2)
        c11=COS(df2)
        dex=EXP(-df1)
        dre=dex*c11
        dima=-dex*s11
        gar=dre*gair-dima*gaii
        gai=dre*gaii+dima*gair
        gair=gar
        gaii=gai
        IF (y0 == 0.) gaii=0.d0
      END IF
    END IF
  ELSE
!CCC CALCULATION OF AI�(Z) CCCCCCCCCCC
    alf=0.1666666666666666D0
    facto=-0.270898621247918D0*r**0.25D0
!CCCCCCCCCCCCCC SERIES OR INTEGRALS CCCCCCCCCCCCCCCCCCCCCCCCCC
    IF ((ya < 3.d0).AND.(xa < 1.3D0).AND.(xa > -2.5D0)) THEN
!CC SERIES
      CALL seraid(xa,ya,gar,gai)
      IF (ifac == 2) THEN
        thet=phase(xa,ya)
        th15=1.5D0*thet
        s1=SIN(th15)
        c1=COS(th15)
        f23r=f23*r32
        df1=f23r*c1
        df2=f23r*s1
        s11=SIN(df2)
        c11=COS(df2)
        dex=EXP(df1)
        dre=dex*c11
        dima=dex*s11
        gair=dre*gar-dima*gai
        gaii=dre*gai+dima*gar
      ELSE
        gair=gar
        gaii=gai
        IF (y0 == 0.) gaii=0.d0
      END IF
    ELSE
      IF (r > 15.d0) THEN
!CC  ASYMPTOTIC EXPANSIONS CCCCCCCCCCCCC
        thet=phase(xa,ya)
        facto=0.5D0/pisr*r**0.25D0
        IF (thet > pi23) THEN
!CCCCCC CONNECTION FORMULAE CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC     N= 1: TRANSFORM Z TO W= U+IV=Z EXP( 2 PI I/3)
          n=1
          c=-0.5D0
          s=n*0.5*sqrt3
          u=xa*c-ya*s
          v=xa*s+ya*c
          v0=v
          IF (v < 0.d0) v=-v
          thet=phase(u,v)
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          th025=thet*0.25D0
          s3=SIN(th025)
          c3=COS(th025)
          CALL expaid(ar1,ai1)
          IF (v0 < 0.d0) ai1=-ai1
          ar=-(c*ar1+s*ai1)
          ai=-(-s*ar1+c*ai1)
          IF (iexpf == 0) THEN
            IF (iexpf2 == 0) THEN
!CC     N=-1: TRANSFORM Z TO W= U+IV=Z EXP(-2 PI I/3)
              n=-1
              c=-0.5D0
              s=n*0.5*sqrt3
              u=xa*c-ya*s
              v=xa*s+ya*c
              v0=v
              IF (v < 0.d0) v=-v
              thet=phase(u,v)
              th15=1.5D0*thet
              s1=SIN(th15)
              c1=COS(th15)
              th025=thet*0.25D0
              s3=SIN(th025)
              c3=COS(th025)
              CALL expaid(ar1,ai1)
              IF (v0 < 0.d0) ai1=-ai1
              thet=phase(xa,ya)
              th15=1.5D0*thet
              s1=SIN(th15)
              c1=COS(th15)
              ro=1.333333333333333D0*r32
              coe1=ro*c1
              coe2=ro*s1
              rex=EXP(coe1)
              dfr=rex*COS(coe2)
              dfi=rex*SIN(coe2)
              ar11=dfr*ar1-dfi*ai1
              ai11=dfr*ai1+dfi*ar1
              gair=ar-(c*ar11+s*ai11)
              gaii=ai-(-s*ar11+c*ai11)
            ELSE
              thet=phase(xa,ya)
              th15=1.5D0*thet
              s1=SIN(th15)
              c1=COS(th15)
              gair=ar
              gaii=ai
            END IF
          ELSE
            gair=ar
            gaii=ai
          END IF
        ELSE
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          th025=thet*0.25D0
          s3=SIN(th025)
          c3=COS(th025)
          CALL expaid(gair,gaii)
        END IF
      ELSE
!CC INTEGRALS CCCCCCCCCCCCCCCC
        thet=phase(xa,ya)
        IF (thet <= pihal) THEN
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          th025=thet*0.25D0
          s3=SIN(th025)
          c3=COS(th025)
          CALL airy1d(xd,wd,gair,gaii)
        END IF
        IF ((thet > pihal).AND.(thet <= pi23)) THEN
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          th025=thet*0.25D0
          s3=SIN(th025)
          c3=COS(th025)
          CALL airy2d(xd,wd,gair,gaii)
        END IF
        IF (thet > pi23) THEN
          n=1
          c=-0.5D0
          s=n*0.5*sqrt3
          u=xa*c-ya*s
          v=xa*s+ya*c
          v0=v
          IF (v < 0.d0) v=-v
          thet=phase(u,v)
          IF (thet <= pihal) THEN
            th15=1.5D0*thet
            s1=SIN(th15)
            c1=COS(th15)
            th025=thet*0.25D0
            s3=SIN(th025)
            c3=COS(th025)
            CALL airy1d(xd,wd,ar1,ai1)
          END IF
          IF ((thet > pihal).AND.(thet <= pi23)) THEN
            th15=1.5D0*thet
            s1=SIN(th15)
            c1=COS(th15)
            th025=thet*0.25D0
            s3=SIN(th025)
            c3=COS(th025)
            CALL airy2d(xd,wd,ar1,ai1)
          END IF
          IF (v0 < 0.d0) ai1=-ai1
          ar=-(c*ar1+s*ai1)
          ai=-(-s*ar1+c*ai1)
          n=-1
          c=-0.5D0
          s=n*0.5*sqrt3
          u=xa*c-ya*s
          v=xa*s+ya*c
          v0=v
          IF (v < 0.d0) v=-v
          thet=phase(u,v)
          IF (thet <= pihal) THEN
            th15=1.5D0*thet
            s1=SIN(th15)
            c1=COS(th15)
            th025=thet*0.25D0
            s3=SIN(th025)
            c3=COS(th025)
            CALL airy1d(xd,wd,ar1,ai1)
          END IF
          IF ((thet > pihal).AND.(thet <= pi23)) THEN
            th15=1.5D0*thet
            s1=SIN(th15)
            c1=COS(th15)
            th025=thet*0.25D0
            s3=SIN(th025)
            c3=COS(th025)
            CALL airy2d(xd,wd,ar1,ai1)
          END IF
          IF (v0 < 0.d0) ai1=-ai1
          thet=phase(xa,ya)
          th15=1.5D0*thet
          s1=SIN(th15)
          c1=COS(th15)
          ro=1.333333333333333D0*r32
          coe1=ro*c1
          coe2=ro*s1
          rex=EXP(coe1)
          dfr=rex*COS(coe2)
          dfi=rex*SIN(coe2)
          ar11=dfr*ar1-dfi*ai1
          ai11=dfr*ai1+dfi*ar1
          gair=ar-(c*ar11+s*ai11)
          gaii=ai-(-s*ar11+c*ai11)
        END IF
      END IF
      IF (ifac == 1) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF THE UNSCALED AI'(z) CCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        f23r=f23*r32
        df1=f23r*c1
        df2=f23r*s1
        s11=SIN(df2)
        c11=COS(df2)
        dex=EXP(-df1)
        dre=dex*c11
        dima=-dex*s11
        gar=dre*gair-dima*gaii
        gai=dre*gaii+dima*gair
        gair=gar
        gaii=gai
        IF (y0 == 0) gaii=0.d0
      END IF
    END IF
  END IF
END IF
IF (y0 < 0.d0) gaii=-gaii
RETURN
END SUBROUTINE aiz

SUBROUTINE  airy1(x,w,gair,gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
!CC              0 <= PHASE(Z) <= PI/2                        C
!CC                                                           C
!CC INPUTS:                                                   C
!CC      X,W,      NODES AND WEIGHTS FOR THE GAUSSIAN         C
!CC                QUADRATURE                                 C
!CC OUTPUTS:                                                  C
!CC      GAIR, GAII,  REAL AND IMAGINARY PARTS OF AI(Z)       C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)             :: x(25)
DOUBLE PRECISION, INTENT(IN)             :: w(25)
DOUBLE PRECISION, INTENT(OUT)            :: gair
DOUBLE PRECISION, INTENT(OUT)            :: gaii

DOUBLE PRECISION :: pih3,pisr,a,alf,thet,r,th15,s1,c1,  &
    r32,facto,th025,s3,c3,sumar,sumai,df1,df1c1,phi,phi6,  &
    s2,c2,dmodu,dmodu2,funr,funi,fac1,fac2,phase
INTEGER :: i

COMMON/param2/pih3,pisr,a,alf
COMMON/param3/thet,r,th15,s1,c1,r32
COMMON/param4/facto,th025,s3,c3

sumar=0.d0
sumai=0.d0
DO  i=1,25
  df1=1.5D0*x(i)/r32
  df1c1=df1*c1
  phi=phase(2.d0+df1c1,df1*s1)
  phi6=phi/6.d0
  s2=SIN(phi6)
  c2=COS(phi6)
  dmodu=SQRT(4.d0+df1*df1+4.d0*df1c1)
  dmodu2=dmodu**alf
  funr=dmodu2*c2
  funi=dmodu2*s2
  sumar=sumar+w(i)*funr
  sumai=sumai+w(i)*funi
END DO
fac1=facto*c3
fac2=facto*s3
gair=fac1*sumar+fac2*sumai
gaii=fac1*sumai-fac2*sumar
RETURN
END SUBROUTINE  airy1

SUBROUTINE  airy2(x,w,gair,gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
!CC              PI/2 < PHASE(Z) <= 2PI/3                     C
!CC                                                           C
!CC INPUTS:                                                   C
!CC      X,W,        NODES AND WEIGHTS FOR THE GAUSSIAN       C
!CC                  QUADRATURE                               C
!CC OUTPUTS:                                                  C
!CC      GAIR, GAII, REAL AND IMAGINARY PARTS OF AI(Z)        C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)             :: x(25)
DOUBLE PRECISION, INTENT(IN)             :: w(25)
DOUBLE PRECISION, INTENT(OUT)            :: gair
DOUBLE PRECISION, INTENT(OUT)            :: gaii

DOUBLE PRECISION :: pih3,pisr,a,alf,thet,r,th15,s1,c1,  &
    r32,facto,th025,s3,c3,sumar,sumai,df1,df1c1,phi,phi6,  &
    s2,c2,dmodu,dmodu2,funr,funi,fac1,fac2,phase
DOUBLE PRECISION :: sqr2,sqr2i,tau,tgtau,b,ang,ctau,cfac,ct,st,  &
    sumr,sumi,ttau,beta
INTEGER :: i

COMMON/param2/pih3,pisr,a,alf
COMMON/param3/thet,r,th15,s1,c1,r32
COMMON/param4/facto,th025,s3,c3

sqr2=1.41421356237310D0
sqr2i=0.707106781186548D0
tau=th15-pih3*0.5D0
tgtau=DTAN(tau)
b=5.d0*a
ang=tau*b
ctau=COS(tau)
cfac=ctau**(-b)
ct=COS(ang)
st=SIN(ang)
sumr=0.d0
sumi=0.d0
DO  i=1,25
  df1=3.d0*x(i)/(ctau*r32)
  df1c1=df1*sqr2i*0.5D0
  phi=phase(2.d0-df1c1,df1c1)
  phi6=phi/6.d0
  ttau=x(i)*tgtau
  beta=phi6-ttau
  s2=SIN(beta)
  c2=COS(beta)
  dmodu=SQRT(4.d0+df1*df1*0.25D0-sqr2*df1)
  dmodu2=dmodu**alf
  funr=dmodu2*c2
  funi=dmodu2*s2
  sumr=sumr+w(i)*funr
  sumi=sumi+w(i)*funi
END DO
sumar=cfac*(ct*sumr-st*sumi)
sumai=cfac*(ct*sumi+st*sumr)
fac1=facto*c3
fac2=facto*s3
gair=fac1*sumar+fac2*sumai
gaii=fac1*sumai-fac2*sumar
RETURN
END SUBROUTINE  airy2

SUBROUTINE airy1d(x,w,gair,gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI'(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
!CC              0 <= PHASE(Z) <= PI/2                         C
!CC                                                            C
!CC INPUTS:                                                    C
!CC       X,W,      NODES AND WEIGHTS FOR THE GAUSSIAN         C
!CC                 QUADRATURE                                 C
!CC OUTPUTS:                                                   C
!CC       GAIR,GAII, REAL AND IMAGINARY PARTS OF AI'(Z)        C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)             :: x(25)
DOUBLE PRECISION, INTENT(IN)             :: w(25)
DOUBLE PRECISION, INTENT(OUT)            :: gair
DOUBLE PRECISION, INTENT(OUT)            :: gaii

DOUBLE PRECISION :: pih3,pisr,a,alf,thet,r,th15,s1,c1,  &
    r32,facto,th025,s3,c3,sumar,sumai,df1,df1c1,phi,phi6,  &
    s2,c2,dmodu,dmodu2,funr,funi,fac1,fac2,phase
INTEGER :: i

COMMON/param2/pih3,pisr,a,alf
COMMON/param3/thet,r,th15,s1,c1,r32
COMMON/param4/facto,th025,s3,c3

sumar=0.d0
sumai=0.d0
DO  i=1,25
  df1=1.5D0*x(i)/r32
  df1c1=df1*c1
  phi=phase(2.d0+df1c1,df1*s1)
  phi6=-phi*alf
  s2=SIN(phi6)
  c2=COS(phi6)
  dmodu=SQRT(4.d0+df1*df1+4.d0*df1c1)
  dmodu2=dmodu**alf
  funr=dmodu2*c2
  funi=dmodu2*s2
  sumar=sumar+w(i)*funr
  sumai=sumai+w(i)*funi
END DO
fac1=facto*c3
fac2=facto*s3
gair=fac1*sumar-fac2*sumai
gaii=fac1*sumai+fac2*sumar
RETURN
END SUBROUTINE airy1d

SUBROUTINE  airy2d(x,w,gair,gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI'(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
!CC              PI/2 < PHASE(Z) <= 3PI/2                      C
!CC                                                            C
!CC INPUTS:                                                    C
!CC      X,W,   NODES AND WEIGHTS FOR THE GAUSSIAN             C
!CC                   QUADRATURE                               C
!CC OUTPUTS:                                                   C
!CC      GAIR,GAII, REAL AND IMAGINARY PARTS OF AI'(Z)         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)             :: x(25)
DOUBLE PRECISION, INTENT(IN)             :: w(25)
DOUBLE PRECISION, INTENT(OUT)            :: gair
DOUBLE PRECISION, INTENT(OUT)            :: gaii

DOUBLE PRECISION :: pih3,pisr,a,alf,thet,r,th15,s1,c1,  &
    r32,facto,th025,s3,c3,sumar,sumai,df1,df1c1,phi,phi6,  &
    s2,c2,dmodu,dmodu2,funr,funi,fac1,fac2,phase
DOUBLE PRECISION :: sqr2,sqr2i,tau,tgtau,b,ang,ctau,cfac,ct,st,  &
    sumr,sumi,ttau,beta
INTEGER :: i

COMMON/param2/pih3,pisr,a,alf
COMMON/param3/thet,r,th15,s1,c1,r32
COMMON/param4/facto,th025,s3,c3

sqr2=1.41421356237310D0
sqr2i=0.707106781186548D0
tau=th15-pih3*0.5D0
tgtau=DTAN(tau)
b=7.d0*alf
ang=tau*b
ctau=COS(tau)
cfac=ctau**(-b)
ct=COS(ang)
st=SIN(ang)
sumr=0.d0
sumi=0.d0
DO  i=1,25
  df1=3.d0*x(i)/(ctau*r32)
  df1c1=df1*sqr2i*0.5D0
  phi=phase(2.d0-df1c1,df1c1)
  phi6=-phi/6.d0
  ttau=x(i)*tgtau
  beta=phi6-ttau
  s2=SIN(beta)
  c2=COS(beta)
  dmodu=SQRT(4.d0+df1*df1*0.25D0-sqr2*df1)
  dmodu2=dmodu**alf
  funr=dmodu2*c2
  funi=dmodu2*s2
  sumr=sumr+w(i)*funr
  sumi=sumi+w(i)*funi
END DO
sumar=cfac*(ct*sumr-st*sumi)
sumai=cfac*(ct*sumi+st*sumr)
fac1=facto*c3
fac2=facto*s3
gair=fac1*sumar-fac2*sumai
gaii=fac1*sumai+fac2*sumar
RETURN
END SUBROUTINE  airy2d

DOUBLE PRECISION FUNCTION phase(x,y)

DOUBLE PRECISION, INTENT(IN)         :: x
DOUBLE PRECISION, INTENT(IN)         :: y
DOUBLE PRECISION :: pi,pihal, ay,p
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  COMPUTES THE PHASE OF Z = X + IY, IN (-PI,PI]
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
COMMON/param1/pi,pihal

IF ((x == 0).AND.(y == 0)) THEN
  p=0.d0
ELSE
  ay=ABS(y)
  IF (x >= ay) THEN
    p=ATAN(ay/x)
  ELSE IF ((x+ay) >= 0.d0) THEN
    p=pihal-ATAN(x/ay)
  ELSE
    p=pi+ATAN(ay/x)
  END IF
  IF (y < 0.d0) p=-p
END IF
phase=p
END FUNCTION phase

SUBROUTINE fgp(x,y,eps,fr,fi,gr,gi)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    COMPUTES THE FUNCTIONS F AND G FOR THE SERIES  C
!    OF AI'(Z).                                     C
!    THIS ROUTINE IS CALLED BY SERAID.              C
!                                                   C
!    INPUTS:                                        C
!         X,Y,  REAL AND IMAGINARY PARTS OF Z       C
!         EPS,  PRECISION FOR THE COMPUTATION OF    C
!               THE SERIES                          C
!    OUTPUTS:                                       C
!         FR,FI, REAL AND IMAGINARY PARTS OF F      C
!         GR,GI, REAL AND IMAGINARY PARTS OF G      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)             :: x
DOUBLE PRECISION, INTENT(IN)             :: y
DOUBLE PRECISION, INTENT(IN OUT)         :: eps
DOUBLE PRECISION, INTENT(OUT)            :: fr
DOUBLE PRECISION, INTENT(OUT)            :: fi
DOUBLE PRECISION, INTENT(OUT)            :: gr
DOUBLE PRECISION, INTENT(OUT)            :: gi

INTEGER :: a,b,k3
DOUBLE PRECISION :: x2,y2,u,v,p,q,cr,ci,dr,di

x2=x*x
y2=y*y
k3=0
u=x*(x2-3*y2)
v=y*(3*x2-y2)
cr=0.5D0
ci=0.d0
dr=1.d0
di=0.d0
fr=0.5D0
fi=0.d0
gr=1.d0
gi=0.d0
70     a=(k3+5)*(k3+3)
b=(k3+1)*(k3+3)
p=(u*cr-v*ci)/a
q=(v*cr+u*ci)/a
cr=p
ci=q
p=(u*dr-v*di)/b
q=(v*dr+u*di)/b
dr=p
di=q
fr=fr+cr
fi=fi+ci
gr=gr+dr
gi=gi+di
k3=k3+3
IF ((ABS(cr)+ABS(dr)+ABS(ci)+ABS(di)) >= eps) GO TO 70
u=x2-y2
v=2.d0*x*y
p=u*fr-v*fi
q=u*fi+v*fr
fr=p
fi=q
RETURN
END SUBROUTINE fgp

SUBROUTINE fg(x,y,eps,fr,fi,gr,gi)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    COMPUTES THE FUNCTIONS F AND G IN EXPRESSION   C
!    10.4.2 OF ABRAMOWITZ & STEGUN FOR THE SERIES   C
!    OF AI(Z).                                      C
!    THIS ROUTINE IS CALLED BY SERAI.               C
!                                                   C
!    INPUTS:                                        C
!          X,Y,  REAL AND IMAGINARY PARTS OF Z      C
!          EPS,  PRECISION FOR THE COMPUTATION      C
!                OF THE SERIES.                     C
!    OUTPUTS:                                       C
!          FR,FI, REAL AND IMAGINARY PARTS OF F     C
!          GR,GI, REAL AND IMAGINARY PARTS OF G     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)             :: x
DOUBLE PRECISION, INTENT(IN)             :: y
DOUBLE PRECISION, INTENT(IN OUT)         :: eps
DOUBLE PRECISION, INTENT(OUT)            :: fr
DOUBLE PRECISION, INTENT(OUT)            :: fi
DOUBLE PRECISION, INTENT(OUT)            :: gr
DOUBLE PRECISION, INTENT(OUT)            :: gi
INTEGER :: a,b,k3
DOUBLE PRECISION :: x2,y2,u,v,p,q,cr,ci,dr,di


x2=x*x
y2=y*y
k3=0
u=x*(x2-3.d0*y2)
v=y*(3.d0*x2-y2)
cr=1.d0
ci=0.d0
dr=1.d0
di=0.d0
fr=1.d0
fi=0.d0
gr=1.d0
gi=0.d0
71     a=(k3+2)*(k3+3)
b=(k3+4)*(k3+3)
p=(u*cr-v*ci)/a
q=(v*cr+u*ci)/a
cr=p
ci=q
p=(u*dr-v*di)/b
q=(v*dr+u*di)/b
dr=p
di=q
fr=fr+cr
fi=fi+ci
gr=gr+dr
gi=gi+di
k3=k3+3
IF ((ABS(cr)+ABS(dr)+ABS(ci)+ABS(di)) >= eps) GO TO 71
p=x*gr-y*gi
q=x*gi+y*gr
gr=p
gi=q
RETURN
END SUBROUTINE fg

SUBROUTINE serai(x,y,air,aii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI(Z), TAYLOR, COMPLEX Z      CCC
!CC                                      CCC
!CC   INPUTS:                            CCC
!CC        X,Y,    REAL AND IMAGINARY    CCC
!CC                PARTS OF Z            CCC
!CC   OUTPUTS:                           CCC
!CC        AIR,AII, REAL AND IMAGINARY   CCC
!CC                 PARTS OF AI(Z)       CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN OUT)         :: x
DOUBLE PRECISION, INTENT(IN OUT)         :: y
DOUBLE PRECISION, INTENT(OUT)            :: air
DOUBLE PRECISION, INTENT(OUT)            :: aii
DOUBLE PRECISION :: eps
DOUBLE PRECISION :: fzr,fzi,gzr,gzi,cons1,cons2
DOUBLE PRECISION :: d1mach

eps=d1mach(3)
IF (eps < 1.d-15) eps=1.d-15
cons1=0.355028053887817239260D0
cons2=0.258819403792806798405D0
CALL fg(x,y,eps,fzr,fzi,gzr,gzi)
air=cons1*fzr-cons2*gzr
aii=cons1*fzi-cons2*gzi
RETURN
END SUBROUTINE serai

SUBROUTINE seraid(x,y,air,aii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI'(Z), TAYLOR, COMPLEX Z     CCC
!CC                                      CCC
!CC   INPUTS:                            CCC
!CC        X,Y,   REAL AND IMAGINARY     CCC
!CC               PARTS OF Z             CCC
!CC   OUTPUTS:                           CCC
!CC        AIR,AII, REAL AND IMAGINARY   CCC
!CC                 PARTS OF AI'(Z)      CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN OUT)         :: x
DOUBLE PRECISION, INTENT(IN OUT)         :: y
DOUBLE PRECISION, INTENT(OUT)            :: air
DOUBLE PRECISION, INTENT(OUT)            :: aii
DOUBLE PRECISION :: eps
DOUBLE PRECISION :: fzr,fzi,gzr,gzi,cons1,cons2
DOUBLE PRECISION :: d1mach

eps=d1mach(3)
IF (eps < 1.d-15) eps=1.d-15
cons1=0.355028053887817239260D0
cons2=0.258819403792806798405D0
CALL fgp(x,y,eps,fzr,fzi,gzr,gzi)
air=cons1*fzr-cons2*gzr
aii=cons1*fzi-cons2*gzi
RETURN
END SUBROUTINE seraid

SUBROUTINE expai(gair,gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI(Z), ASYMPTOTIC EXPANSION, COMPLEX Z CCC
!CC                                               CCC
!CC   OUTPUTS:                                    CCC
!CC        GAIR, GAII,  REAL AND IMAGINARY        CCC
!CC                     PARTS OF AI(Z)            CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(OUT)            :: gair
DOUBLE PRECISION, INTENT(OUT)            :: gaii
DOUBLE PRECISION :: eps
DOUBLE PRECISION :: thet,r,th15,s1, c1,r32,facto,th025,s3,c3
DOUBLE PRECISION :: df1,psiir,psiii,ck,dfrr,dfii,sumar,sumai,  &
    dfr,dfi,deltar,deltai,fac1,fac2
DOUBLE PRECISION :: co,df
DOUBLE PRECISION :: d1mach
INTEGER :: k
DIMENSION co(20)
COMMON/param3/thet,r,th15,s1,c1,r32
COMMON/param4/facto,th025,s3,c3
SAVE co
DATA co/-6.944444444444445D-2,3.713348765432099D-2,  &
    -3.799305912780064D-2,5.764919041266972D-2,-0.116099064025515D0,  &
    0.291591399230751D0,-0.877666969510017D0,3.07945303017317D0,  &
    -12.3415733323452D0,55.6227853659171D0,-278.465080777603D0,  &
    1533.16943201280D0,-9207.20659972641D0,59892.5135658791D0,  &
    -419524.875116551D0,3148257.41786683D0,-25198919.8716024D0,  &
    214288036.963680D0,-1929375549.18249D0,18335766937.8906D0/

eps=d1mach(3)
IF (eps < 1.d-15) eps=1.d-15
df1=1.5D0/r32
psiir=c1
psiii=-s1
k=0
ck=1.d0
df=1.d0
dfrr=1.d0
dfii=0.d0
sumar=1.d0
sumai=0.d0
80      df=df*df1
ck=co(k+1)*df
dfr=dfrr
dfi=dfii
dfrr=dfr*psiir-dfi*psiii
dfii=dfr*psiii+dfi*psiir
deltar=dfrr*ck
deltai=dfii*ck
sumar=sumar+deltar
sumai=sumai+deltai
k=k+1
IF (sumar /= 0) THEN
  IF (ABS(deltar/sumar) > eps)  GO TO 80
END IF
IF (sumai /= 0) THEN
  IF (ABS(deltai/sumai) > eps) GO TO 80
END IF
fac1=facto*c3
fac2=facto*s3
gair=fac1*sumar+fac2*sumai
gaii=fac1*sumai-fac2*sumar
RETURN
END SUBROUTINE expai

SUBROUTINE expaid(gair,gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI'(Z), ASYMPTOTIC EXPANSION, COMPLEX Z CCC
!CC                                                CCC
!CC   OUTPUTS:                                     CCC
!CC        GAIR, GAII,  REAL AND IMAGINARY         CCC
!CC                     PARTS OF AI'(Z)            CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(OUT)            :: gair
DOUBLE PRECISION, INTENT(OUT)            :: gaii
DOUBLE PRECISION :: eps
DOUBLE PRECISION :: thet,r,th15,s1, c1,r32,facto,th025,s3,c3
DOUBLE PRECISION :: df1,psiir,psiii,vk,dfrr,dfii,sumar,sumai,  &
    dfr,dfi,deltar,deltai,fac1,fac2
DOUBLE PRECISION :: co,df
DOUBLE PRECISION :: d1mach
INTEGER :: k
DIMENSION co(20)
COMMON/param3/thet,r,th15,s1,c1,r32
COMMON/param4/facto,th025,s3,c3
SAVE co
DATA co/9.722222222222222D-2,-4.388503086419753D-2,  &
    4.246283078989484D-2,-6.266216349203230D-2,  &
    0.124105896027275D0,-0.308253764901079D0,  &
    0.920479992412945D0,-3.21049358464862D0,  &
    12.8072930807356D0,-57.5083035139143D0,  &
    287.033237109221D0,-1576.35730333710D0,  &
    9446.35482309593D0,-61335.7066638521D0,  &
    428952.400400069D0,-3214536.52140086D0,  &
    25697908.3839113D0,-218293420.832160D0,  &
    1963523788.99103D0,-18643931088.1072D0/

eps=d1mach(3)
IF (eps < 1.d-15) eps=1.d-15
df1=1.5D0/r32
psiir=c1
psiii=-s1
k=0
df=1.d0
dfrr=1.d0
dfii=0.d0
sumar=1.d0
sumai=0.d0
81     df=df*df1
vk=co(k+1)*df
dfr=dfrr
dfi=dfii
dfrr=dfr*psiir-dfi*psiii
dfii=dfr*psiii+dfi*psiir
deltar=dfrr*vk
deltai=dfii*vk
sumar=sumar+deltar
sumai=sumai+deltai
k=k+1
IF (sumar /= 0) THEN
  IF (ABS(deltar/sumar) > eps)  GO TO 81
END IF
IF (sumai /= 0) THEN
  IF (ABS(deltai/sumai) > eps) GO TO 81
END IF
fac1=facto*c3
fac2=facto*s3
gair=-(fac1*sumar-fac2*sumai)
gaii=-(fac1*sumai+fac2*sumar)
RETURN
END SUBROUTINE expaid

SUBROUTINE biz(ifun,ifac,x0,y0,gbir,gbii,ierro)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! COMPUTATION OF THE AIRY FUNCTION BI(Z) OR ITS DERIVATIVE BI'(Z)
! THE CODE USES THE CONNECTION OF BI(Z) WITH AI(Z).
!                BIZ CALLS THE ROUTINE AIZ
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  INPUTS:
!    IFUN:
!         * IFUN=1, THE CODE COMPUTES BI(Z)
!         * IFUN=2, THE CODE COMPUTES BI'(Z)
!    IFAC:
!         * IFAC=1, THE CODE COMPUTES  BI(Z) OR BI'(Z)
!         * IFAC=2, THE CODE COMPUTES NORMALIZED BI(Z) OR BI'(Z)
!    X0:   REAL PART OF THE ARGUMENT Z
!    Y0:   IMAGINARY PART OF THE ARGUMENT  Z

!  OUTPUTS:
!    GBIR: REAL PART OF BI(Z) OR BI'(Z)
!    GBII: IMAGINARY PART OF BI(Z) OR BI'(Z)

!    IERRO: ERROR FLAG
!          * IERRO=0, SUCCESSFUL COMPUTATION
!          * IERRO=1, COMPUTATION OUT OF RANGE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          MACHINE DEPENDENT CONSTANTS: FUNCTION D1MACH
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   ACCURACY:

!     1) SCALED AIRY FUNCTIONS:
!        RELATIVE ACCURACY BETTER THAN 10**(-13) EXCEPT CLOSE TO
!        THE ZEROS, WHERE 10**(-13) IS THE ABSOLUTE PRECISION.
!        GRADUAL LOSS OF PRECISION TAKES PLACE FOR |Z|>1000
!        IN THE CASE OF PHASE(Z) CLOSE TO +3*PI/2 OR -3*PI/2.
!     2) UNSCALED AIRY FUNCTIONS:
!        THE FUNCTION OVERFLOWS/UNDERFLOWS FOR
!        3/2*|Z|**(3/2)>LOG(OVER).
!        FOR |Z|<30:
!        A) RELATIVE ACCURACY FOR THE MODULUS (EXCEPT AT THE
!           ZEROS) BETTER THAN 10**(-13).
!        B) ABSOLUTE ACCURACY FOR MIN(R(Z),1/R(Z)) BETTER
!           THAN 10**(-13), WHERE R(Z)=REAL(BI)/IMAG(BI)
!           OR R(Z)=REAL(BI')/IMAG(BI').
!        FOR |Z|>30, GRADUAL LOSS OF PRECISION TAKES PLACE
!        AS |Z| INCREASES.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     AUTHORS:
!        AMPARO GIL    (U. AUTONOMA DE MADRID, MADRID, SPAIN).
!                      E-MAIL: AMPARO.GIL@UAM.ES
!        JAVIER SEGURA (U. CARLOS III DE MADRID, MADRID, SPAIN).
!                      E-MAIL: JSEGURA@MATH.UC3M.ES
!        NICO M. TEMME (CWI, AMSTERDAM, THE NETHERLANDS).
!                      E-MAIL: NICO.TEMME@CWI.NL
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    REFERENCES:
!         COMPUTING AIRY FUNCTIONS BY NUMERICAL QUADRATURE.
!         NUMERICAL ALGORITHMS (2001).
!         A. GIL, J. SEGURA, N.M. TEMME
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN OUT)                  :: ifun
INTEGER, INTENT(IN OUT)                  :: ifac
DOUBLE PRECISION, INTENT(IN)             :: x0
DOUBLE PRECISION, INTENT(IN)             :: y0
DOUBLE PRECISION, INTENT(OUT)            :: gbir
DOUBLE PRECISION, INTENT(OUT)            :: gbii
INTEGER, INTENT(OUT)                     :: ierro

DOUBLE PRECISION :: over,under,dl1,dl2,cover,d1mach
DOUBLE PRECISION :: pi,pi3,pi23,sqrt3,c,s,c1,s1,u,v,x,y,ar,ai,  &
    apr,api,br,bi,bpr,bpi,bbr,bbi,bbpr,bbpi,phase
DOUBLE PRECISION :: thet,r,r32,thet32,a1,b1,df1,expo,expoi,  &
    etar,etai,etagr,etagi,pihal
INTEGER :: iexpf,ierr
COMMON/param1/pi,pihal

sqrt3=1.7320508075688772935D0
pi=3.1415926535897932385D0
pihal=1.5707963267948966192D0
pi3=pi/3.d0
pi23=2.d0*pi3
x=x0
c=0.5D0*sqrt3
s=0.5D0
ierro=0
iexpf=0
IF (y0 < 0.d0) THEN
  y=-y0
ELSE
  y=y0
END IF
r=SQRT(x*x+y*y)
r32=r*SQRT(r)
thet=phase(x,y)
cover=2.d0/3.d0*r32*ABS(COS(1.5D0*thet))
!CC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
over=d1mach(2)*1.d-3
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+3
dl1=LOG(over)
dl2=-LOG(under)
IF (dl1 > dl2) over=1/under
IF (ifac == 1) THEN
  IF (cover >= LOG(over)) THEN
!CC OVERFLOW/UNDERFLOW PROBLEMS.
!CC   CALCULATION ABORTED
    ierro=1
    gbir=0
    gbii=0
  END IF
ELSE
  IF (cover >= (LOG(over)*0.2)) iexpf=1
END IF
IF (ierro == 0) THEN
  IF (ifac == 1) THEN
    IF (y == 0.d0) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  TAKE TWICE THE REAL PART OF EXP(-PI I/6) AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      c1=-0.5D0
      s1=-0.5D0*sqrt3
      u=x*c1-y*s1
      v=x*s1+y*c1
      CALL aiz(ifun,ifac,u,v,ar,ai,ierr)
      IF (ifun == 1) THEN
        br=sqrt3*ar+ai
        bi=0.d0
      ELSE
        u=ar*c1-ai*s1
        v=ar*s1+ai*c1
        apr=u
        api=v
        bpr=sqrt3*apr+api
        bpi=0.d0
      END IF
    ELSE
      IF ((x < 0.d0).AND.(y < -x*sqrt3)) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC      2 PI/3  < PHASE(Z) < PI
!CC      BI(Z)=EXP(I PI/6) AI_(-1)(Z) + EXP(-I PI/6) AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        c1=-0.5D0
        s1=0.5D0*sqrt3
        u=x*c1-y*s1
        v=x*s1+y*c1
        CALL aiz(ifun,ifac,u,v,ar,ai,ierr)
        IF (ifun == 1) THEN
          br=c*ar-s*ai
          bi=c*ai+s*ar
        ELSE
          u=ar*c1-ai*s1
          v=ar*s1+ai*c1
          apr=u
          api=v
          bpr=c*apr-s*api
          bpi=c*api+s*apr
        END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   WE NEED ALSO AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        c1=-0.5D0
        s1=-0.5D0*sqrt3
        u=x*c1-y*s1
        v=x*s1+y*c1
        CALL aiz(ifun,ifac,u,v,ar,ai,ierr)
        IF (ifun == 1) THEN
          s=-s
          br=br+c*ar-s*ai
          bi=bi+c*ai+s*ar
        ELSE
          u=ar*c1-ai*s1
          v=ar*s1+ai*c1
          apr=u
          api=v
          s=-s
          bpr=bpr+c*apr-s*api
          bpi=bpi+c*api+s*apr
        END IF
      ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   BI(Z) = I AI(Z) + 2 EXP(-I PI/6) AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        c1=-0.5D0
        s1=-0.5D0*sqrt3
        u=x*c1-y*s1
        v=x*s1+y*c1
        CALL aiz(ifun,ifac,u,v,ar,ai,ierr)
        IF (ifun == 1) THEN
          br=sqrt3*ar+ai
          bi=-ar+sqrt3*ai
        ELSE
          u=ar*c1-ai*s1
          v=ar*s1+ai*c1
          apr=u
          api=v
          bpr=sqrt3*apr+api
          bpi=-apr+sqrt3*api
        END IF
        CALL aiz(ifun,ifac,x,y,ar,ai,ierr)
        IF (ifun == 1) THEN
          br=br-ai
          bi=bi+ar
        ELSE
          bpr=bpr-ai
          bpi=bpi+ar
        END IF
      END IF
    END IF
  ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC         SCALED BI AIRY FUNCTIONS         C
!CC   WE USE THE FOLLOWING NORMALIZATION:    C
!CC   LET ARGZ=ARG(Z), THEN:                 C
!CC   A) IF  0 <= ARGZ <= PI/3               C
!CC      BI=EXP(-2/3Z^3/2)BI                 C
!CC   B) IF  PI/3 <= ARGZ <= PI              C
!CC      BI=EXP(2/3Z^3/2)BI                  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    thet=phase(x,y)
    IF (thet <= pi3) THEN
      c1=-0.5D0
      s1=-0.5D0*sqrt3
      u=x*c1-y*s1
      v=x*s1+y*c1
      CALL aiz(ifun,ifac,u,v,ar,ai,ierr)
      IF (ifun == 1) THEN
        br=sqrt3*ar+ai
        bi=-ar+sqrt3*ai
      ELSE
        u=ar*c1-ai*s1
        v=ar*s1+ai*c1
        apr=u
        api=v
        bpr=sqrt3*apr+api
        bpi=-apr+sqrt3*api
      END IF
      IF (iexpf == 0) THEN
        r=SQRT(x*x+y*y)
        r32=r*SQRT(r)
        thet32=thet*1.5D0
        a1=COS(thet32)
        b1=SIN(thet32)
        df1=4.d0/3.d0*r32
        expo=EXP(df1*a1)
        expoi=1.d0/expo
        etar=expo*COS(df1*b1)
        etai=expo*SIN(df1*b1)
        etagr=expoi*COS(-df1*b1)
        etagi=expoi*SIN(-df1*b1)
        CALL aiz(ifun,ifac,x,y,ar,ai,ierr)
        IF (ifun == 1) THEN
          br=br-ar*etagi-etagr*ai
          bi=bi+ar*etagr-etagi*ai
        ELSE
          bpr=bpr-ar*etagi-etagr*ai
          bpi=bpi+ar*etagr-etagi*ai
        END IF
      END IF
    END IF
    IF ((thet > pi3).AND.(thet <= pi23)) THEN
      IF (iexpf == 0) THEN
        r=SQRT(x*x+y*y)
        r32=r*SQRT(r)
        thet32=thet*1.5D0
        a1=COS(thet32)
        b1=SIN(thet32)
        df1=4.d0/3.d0*r32
        expo=EXP(df1*a1)
        expoi=1.d0/expo
        etar=expo*COS(df1*b1)
        etai=expo*SIN(df1*b1)
        etagr=expoi*COS(-df1*b1)
        etagi=expoi*SIN(-df1*b1)
        c1=-0.5D0
        s1=-0.5D0*sqrt3
        u=x*c1-y*s1
        v=x*s1+y*c1
        CALL aiz(ifun,ifac,u,v,ar,ai,ierr)
        IF (ifun == 1) THEN
          bbr=sqrt3*ar+ai
          bbi=-ar+sqrt3*ai
          br=bbr*etar-bbi*etai
          bi=bbi*etar+bbr*etai
        ELSE
          u=ar*c1-ai*s1
          v=ar*s1+ai*c1
          apr=u
          api=v
          bbpr=sqrt3*apr+api
          bbpi=-apr+sqrt3*api
          bpr=bbpr*etar-bbpi*etai
          bpi=bbpi*etar+bbpr*etai
        END IF
      ELSE
        IF (ifun == 1) THEN
          br=0.d0
          bi=0.d0
        ELSE
          bpr=0.d0
          bpi=0.d0
        END IF
      END IF
      CALL aiz(ifun,ifac,x,y,ar,ai,ierr)
      IF (ifun == 1) THEN
        br=br-ai
        bi=bi+ar
      ELSE
        bpr=bpr-ai
        bpi=bpi+ar
      END IF
    END IF
    IF (thet > pi23) THEN
      c1=-0.5D0
      s1=0.5D0*sqrt3
      u=x*c1-y*s1
      v=x*s1+y*c1
      CALL aiz(ifun,ifac,u,v,ar,ai,ierr)
      IF (ifun == 1) THEN
        br=c*ar-s*ai
        bi=c*ai+s*ar
      ELSE
        u=ar*c1-ai*s1
        v=ar*s1+ai*c1
        apr=u
        api=v
        bpr=c*apr-s*api
        bpi=c*api+s*apr
      END IF
      IF (iexpf == 0) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   WE NEED ALSO AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        r=SQRT(x*x+y*y)
        r32=r*SQRT(r)
        thet32=thet*1.5D0
        a1=COS(thet32)
        b1=SIN(thet32)
        df1=4.d0/3.d0*r32
        expo=EXP(df1*a1)
        expoi=1.d0/expo
        etar=expo*COS(df1*b1)
        etai=expo*SIN(df1*b1)
        etagr=expoi*COS(-df1*b1)
        etagi=expoi*SIN(-df1*b1)
        c1=-0.5D0
        s1=-0.5D0*sqrt3
        u=x*c1-y*s1
        v=x*s1+y*c1
        CALL aiz(ifun,ifac,u,v,ar,ai,ierr)
        IF (ifun == 1) THEN
          s=-s
          bbr=c*ar-s*ai
          bbi=c*ai+s*ar
          br=br+etar*bbr-etai*bbi
          bi=bi+bbi*etar+etai*bbr
        ELSE
          u=ar*c1-ai*s1
          v=ar*s1+ai*c1
          apr=u
          api=v
          s=-s
          bbpr=c*apr-s*api
          bbpi=c*api+s*apr
          bpr=bpr+etar*bbpr-etai*bbpi
          bpi=bpi+bbpi*etar+etai*bbpr
        END IF
      END IF
    END IF
  END IF
  IF (y0 < 0) THEN
    bi=-bi
    bpi=-bpi
  END IF
  IF (ifun == 1) THEN
    gbir=br
    gbii=bi
  ELSE
    gbir=bpr
    gbii=bpi
  END IF
END IF
RETURN
END SUBROUTINE biz




DOUBLE PRECISION FUNCTION d1mach(i)
INTEGER, INTENT(IN)                      :: i
!  DOUBLE-PRECISION MACHINE CONSTANTS
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  D1MACH( 5) = LOG10(B)
INTEGER :: small(2)
INTEGER :: large(2)
INTEGER :: right(2)
INTEGER :: diver(2)
INTEGER :: LOG10(2)
INTEGER :: sc, cray1(38), j
COMMON /d9mach/ cray1
SAVE small, large, right, diver, LOG10, sc
DOUBLE PRECISION :: dmach(5)
EQUIVALENCE (dmach(1),small(1))
EQUIVALENCE (dmach(2),large(1))
EQUIVALENCE (dmach(3),right(1))
EQUIVALENCE (dmach(4),diver(1))
EQUIVALENCE (dmach(5),LOG10(1))
!  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
!  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
!  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
!  MANY MACHINES YET.
!  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
!  ON THE NEXT LINE
DATA sc/0/
!  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
!  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
!          mail netlib@research.bell-labs.com
!          send old1mach from blas
!  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.

!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
!      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
!      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
!      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
!      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/

!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGERS.
!      DATA SMALL(1),SMALL(2) /    8388608,           0 /
!      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
!      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
!      DATA DIVER(1),DIVER(2) /  620756992,           0 /
!      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/

!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
!      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
!      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
!      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
!      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/

!     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
IF (sc /= 987) THEN
  dmach(1) = 1.d13
  IF (      small(1) == 1117925532 .AND. small(2) == -448790528) THEN
!           *** IEEE BIG ENDIAN ***
    small(1) = 1048576
    small(2) = 0
    large(1) = 2146435071
    large(2) = -1
    right(1) = 1017118720
    right(2) = 0
    diver(1) = 1018167296
    diver(2) = 0
    LOG10(1) = 1070810131
    LOG10(2) = 1352628735
  ELSE IF ( small(2) == 1117925532  &
        .AND. small(1) == -448790528) THEN
!           *** IEEE LITTLE ENDIAN ***
    small(2) = 1048576
    small(1) = 0
    large(2) = 2146435071
    large(1) = -1
    right(2) = 1017118720
    right(1) = 0
    diver(2) = 1018167296
    diver(1) = 0
    LOG10(2) = 1070810131
    LOG10(1) = 1352628735
  ELSE IF ( small(1) == -2065213935  &
        .AND. small(2) == 10752) THEN
!               *** VAX WITH D_FLOATING ***
    small(1) = 128
    small(2) = 0
    large(1) = -32769
    large(2) = -1
    right(1) = 9344
    right(2) = 0
    diver(1) = 9472
    diver(2) = 0
    LOG10(1) = 546979738
    LOG10(2) = -805796613
  ELSE IF ( small(1) == 1267827943  &
        .AND. small(2) == 704643072) THEN
!               *** IBM MAINFRAME ***
    small(1) = 1048576
    small(2) = 0
    large(1) = 2147483647
    large(2) = -1
    right(1) = 856686592
    right(2) = 0
    diver(1) = 873463808
    diver(2) = 0
    LOG10(1) = 1091781651
    LOG10(2) = 1352628735
  ELSE IF ( small(1) == 1120022684  &
        .AND. small(2) == -448790528) THEN
!           *** CONVEX C-1 ***
    small(1) = 1048576
    small(2) = 0
    large(1) = 2147483647
    large(2) = -1
    right(1) = 1019215872
    right(2) = 0
    diver(1) = 1020264448
    diver(2) = 0
    LOG10(1) = 1072907283
    LOG10(2) = 1352628735
  ELSE IF ( small(1) == 815547074  &
        .AND. small(2) == 58688) THEN
!           *** VAX G-FLOATING ***
    small(1) = 16
    small(2) = 0
    large(1) = -32769
    large(2) = -1
    right(1) = 15552
    right(2) = 0
    diver(1) = 15568
    diver(2) = 0
    LOG10(1) = 1142112243
    LOG10(2) = 2046775455
  ELSE
    dmach(2) = 1.d27 + 1
    dmach(3) = 1.d27
    large(2) = large(2) - right(2)
    IF (large(2) == 64 .AND. small(2) == 0) THEN
      cray1(1) = 67291416
      DO  j = 1, 20
        cray1(j+1) = cray1(j) + cray1(j)
      END DO
      cray1(22) = cray1(21) + 321322
      DO  j = 22, 37
        cray1(j+1) = cray1(j) + cray1(j)
      END DO
      IF (cray1(38) == small(1)) THEN
!                  *** CRAY ***
        CALL i1mcry(small(1), j, 8285, 8388608, 0)
        small(2) = 0
        CALL i1mcry(large(1), j, 24574, 16777215, 16777215)
        CALL i1mcry(large(2), j, 0, 16777215, 16777214)
        CALL i1mcry(right(1), j, 16291, 8388608, 0)
        right(2) = 0
        CALL i1mcry(diver(1), j, 16292, 8388608, 0)
        diver(2) = 0
        CALL i1mcry(LOG10(1), j, 16383, 10100890, 8715215)
        CALL i1mcry(LOG10(2), j, 0, 16226447, 9001388)
      ELSE
        WRITE(*,9000)
        STOP 779
      END IF
    ELSE
      WRITE(*,9000)
      STOP 779
    END IF
  END IF
  sc = 987
END IF
!    SANITY CHECK
IF (dmach(4) >= 1.0D0) STOP 778
IF (i < 1 .OR. i > 5) THEN
  WRITE(*,*) 'D1MACH(I): I =',i,' is out of bounds.'
  STOP
END IF
d1mach = dmach(i)
RETURN
9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/  &
    ' appropriate for your machine.')
! /* Standard C source for D1MACH -- remove the * in column 1 */
!#include <stdio.h>
!#include <float.h>
!#include <math.h>
!double d1mach_(long *i)
!{
! switch(*i){
!   case 1: return DBL_MIN;
!   case 2: return DBL_MAX;
!   case 3: return DBL_EPSILON/FLT_RADIX;
!   case 4: return DBL_EPSILON;
!   case 5: return log10((double)FLT_RADIX);
!   }
! fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
! exit(1); return 0; /* some compilers demand return values */
!}
END FUNCTION d1mach

SUBROUTINE i1mcry(a, a1, b, c, d)
!*** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****


INTEGER, INTENT(OUT)                     :: a
INTEGER, INTENT(OUT)                     :: a1
INTEGER, INTENT(IN)                      :: b
INTEGER, INTENT(IN)                      :: c
INTEGER, INTENT(IN)                      :: d

a1 = 16777216*b + c
a = 16777216*a1 + d
END SUBROUTINE i1mcry




SUBROUTINE dkia(ifac,x,a,dki,dkid,ierro)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC CALCULATION OF K_IA(X) AND K'_IA(X), WHERE A = PNU.
!CCC THE RANGE OF THE PARAMETERS (X,A) FOR THE COMPUTATION IS:
!CCC        0 < X <= 1500,  -1500 <= A <= 1500
!CCC WE USE:
!CCC    1) SERIES,  IF  ABS(A) > 0.044*ABS(X-3.1)**1.9+(X-3.1)
!CCC    2) CONTINUED FRACTION,
!CCC        IF X>3 AND ABS(A) < 380*(ABS((X-3)/2300))**0.572
!CCC    3) AIRY-TYPE ASYMPTOTIC EXPANSION,
!CCC        IF ABS(A) > 0.4*X+7.5 AND ABS(A) < 3.7*X-10
!CCC    4) QUADRATURES,
!CCC        IN THE REST OF CASES
!CCC IFAC:
!CCC     * IFAC=1, THE CODE COMPUTES K_IA(X), K'_IA(X)
!CCC     * IFAC=2, THE CODE COMPUTES NORMALIZED K_IA(X), K'_IA(X)
!CCC IERRO: ERROR FLAG
!CCC     * IF IERRO=0, COMPUTATION SUCCESSFUL.
!CCC     * IF IERRO=1, COMPUTATION OUT OF RANGE.
!CCC     * IF IERRO=2, ARGUMENT X AND/OR ORDER A,  OUT OF RANGE.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC     MACHINE DEPENDENT CONSTANTS: FUNCTION D1MACH
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC ACCURACY:
!CCC THE RELATIVE ACCURACY IS:
!CCC  BETTER THAN 10**(-13)   FOR (X,A) IN (0,200]X[-200,200];
!CCC  BETTER THAN 5.10**(-13) FOR (X,A) IN (0,500]X[-500,500];
!CCC  CLOSE TO    10**(-12)   FOR (X,A) IN (0,1500]X[-1500,1500].
!CCC NEAR THE ZEROS OF THE FUNCTIONS (THERE ARE INFINITELY
!CCC MANY OF THEM IN THE OSCILLATORY REGION) RELATIVE PRECISION
!CCC LOOSES MEANING AND ONLY ABSOLUTE PRECISION MAKES SENSE.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC  VARIABLES:
! DKI    :  KIA(X) FUNCTION
! DKID   :  DERIVATIVE OF THE KIA(X) FUNCTION
! IFAC   :
!          * IF IFAC=1, COMPUTATION OF KIA(X), KIA'(X)
!          * IF IFAC=2, COMPUTATION OF SCALED KIA(X), KIA'(X)
! PNU    :  ORDER OF THE FUNCTIONS
! X      :  ARGUMENT OF THE FUNCTIONS
! Z      :  RATIO X/PNU
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     AUTHORS:
!        AMPARO GIL    (U. AUTONOMA DE MADRID, MADRID, SPAIN).
!                      E-MAIL: AMPARO.GIL@UAM.ES
!        JAVIER SEGURA (U. CANTABRIA, SANTANDER, SPAIN).
!                      E-MAIL: SEGURAJJ@UNICAN.ES
!        NICO M. TEMME (CWI, AMSTERDAM, THE NETHERLANDS).
!                      E-MAIL: NICO.TEMME@CWI.NL
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    REFERENCES:
!     1) 'COMPUTATION OF THE MODIFIED BESSEL FUNCTION OF THE
!         THIRD KIND FOR IMAGINARY ORDERS'
!         A. GIL, J. SEGURA, N.M. TEMME
!         J. COMPUT. PHYS. 175 (2002) 398-411
!     2) 'COMPUTATION OF THE MODIFIED BESSEL FUNCTIONS OF THE
!         THIRD KIND OF IMAGINARY ORDERS:
!         UNIFORM AIRY-TYPE ASYMPTOTIC EXPANSION'
!         A. GIL, J. SEGURA, N.M. TEMME, PROCEEDINGS OPSFA 2001
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN)                      :: ifac
DOUBLE PRECISION, INTENT(IN)             :: x
DOUBLE PRECISION, INTENT(IN)             :: a
DOUBLE PRECISION, INTENT(OUT)            :: dki
DOUBLE PRECISION, INTENT(OUT)            :: dkid
INTEGER, INTENT(OUT)                     :: ierro
DOUBLE PRECISION :: pnu
DOUBLE PRECISION :: d,df1,df2,df3,df4


ierro=0
pnu=a
IF ((x > 1500).OR.(x <= 0.d0)) THEN
  ierro=2
  dki=0.d0
  dkid=0.d0
END IF
IF (ABS(pnu) > 1500) THEN
  ierro=2
  dki=0.d0
  dkid=0.d0
END IF
IF (ierro == 0) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC THESE FUNCTIONS ARE EVEN FUNCTIONS IN THE  C
!CC PARAMETER A (=PNU)                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  IF (pnu < 0.d0) pnu=-pnu
  d=x-3.1D0
  df1=0.044D0*ABS(d)**1.9D0+d
  df2=380.d0*(ABS((x-3.d0)/2300.d0))**0.572D0
  df3=0.4D0*x+7.5D0
  df4=3.7D0*x-10.d0
  IF (pnu > df1) THEN
    CALL serkia(ifac,x,pnu,dki,dkid,ierro)
  ELSE IF ((x > 3.d0).AND.(pnu < df2)) THEN
    CALL frakia(ifac,x,pnu,dki,dkid,ierro)
  ELSE IF ((pnu > df3).AND.(pnu < df4)) THEN
    CALL aiexki(ifac,x,pnu,dki,dkid,ierro)
  ELSE
    CALL dkint(ifac,x,pnu,dki,dkid,ierro)
  END IF
END IF
RETURN
END SUBROUTINE dkia

SUBROUTINE dkint(ifac,xx,pnua,dkinf,dkind,ierro)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  CALCULATION OF KIA, KIA' USING NON-OSCILLATING INTEGRAL   C
!CC  REPRESENTATIONS                                           C
!CC                                                            C
!CC  INPUT:                                                    C
!CC         XX:  ARGUMENT OF THE KIA, KIA' FUNCTIONS           C
!CC         PNUA: ORDER OF THE KIA, KIA' FUNCTIONS             C
!CC         IFAC:                                              C
!CC               =1,  NON SCALED KIA, KIA'                    C
!CC               =2,  SCALED KIA, KIA'                        C
!CC  OUTPUT:                                                   C
!CC         DKINF: KIA FUNCTION                                C
!CC         DKIND: KIA' FUNCTION                               C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  METHOD OF COMPUTATION:
!CC  * IF XX>=PNUA, THE NON-OSCILLATING INTEGRAL REPRESENTATIONS
!CC                 FOR THE MONOTONIC REGION ARE USED
!CC  * IF XX<PNUA,  THE NON-OSCILLATING INTEGRAL REPRESENTATIONS
!CC                 FOR THE OSCILLATORY REGION ARE USED
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  VARIABLES:
! CHI      : X*SINH(MU)-PNU*MU
! CONTR1   : CONTRIBUTION OF THE SEMI-INFINITE INTEGRAL
!            IN THE OSCILLATORY CASE (INCLUDING ADDITIONAL
!            FACTORS: CONTR1=FOS1*FACTORS).
! COSCHI   : COS(CHI)
! COSH2M   : COSH(2*MU)
! COSHM    : COSH(MU)
! COSTH    : COS(THET), THET=ASIN(PNU/X)
! DF1      : FACTOR (DEPENDING ON IFAC)
! DKIND    : DERIVATIVE OF THE KIA(X) FUNCTION
! DKINF    : KIA(X) FUNCTION
! DMAIN    : PI*PNU*0.5
! DMU      : SOLUTION MU OF COSH(MU)=PNU/X
! DMU2     : 2*MU
! DMU3     : 3*MU
! DMU5     : 5*MU
! DMU7     : 7*MU
! DMU9     : 9*MU
! DMUFAC   : MU*COSH(MU)-SINH(MU)
! DMUTAN   : MU-TANH(MU)
! FDOMIN   : X*(COS(THET)+THET*SIN(THET))
! FOS1     : CONTRIBUTION OF THE SEMI-INFINITE INTEGRAL
!            IN THE OSCILLATORY CASE (KIA(X)).
! FOSD1    : CONTRIBUTION OF THE SEMI-INFINITE INTEGRAL
!            IN THE OSCILLATORY CASE (KIA'(X)).
! HIR      : MONOTONIC CASE, CONTRIBUTION OF THE INTEGRAL
!            (KIA(X)).
! HIRD     : MONOTONIC CASE, CONTRIBUTION OF THE INTEGRAL
!            (KIA'(X)).
! IERRO    : FLAG FOR UNDERFLOW/OVERFLOW ERRORS
! IFAC     : INDICATES THE CHOICE FOR SCALING
! PI       : 3.1415...
! PINU     : PI*PNU
! PNU      : ORDER OF THE FUNCTION
! PNUA     : SAME AS PNU
! SINCHI   : SIN(CHI)
! SINH2M   : SINH(2*MU)
! SINHM    : SINH(MU)
! SINTH    : SIN(THET)
! THET     : ASIN(PNU/X)
! UNDER    : UNDERFLOW NUMBER
! X        : ARGUMENT OF THE FUNCTION
! XX       : SAME AS X
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN)                      :: ifac
DOUBLE PRECISION, INTENT(IN)             :: xx
DOUBLE PRECISION, INTENT(IN)             :: pnua
DOUBLE PRECISION, INTENT(OUT)            :: dkinf
DOUBLE PRECISION, INTENT(OUT)            :: dkind
INTEGER, INTENT(OUT)                     :: ierro
DOUBLE PRECISION :: x,pnu,sinth,thet,  &
    costh,fdomin,hir,hird,pi,coshm,dmu,cosh2m,sinhm,  &
    sinh2m,dmutan,chi,dmufac,dmu2,dmu3,dmu5,dmu7,dmu9,coschi,  &
    sinchi,pinu,dmain,fos1,df1,contr1,fosd1,d1mach,under

COMMON/argu/x,pnu
COMMON/parmon/thet,sinth,costh,fdomin
COMMON/paros1/dmu,coshm,sinhm,dmufac,coschi,sinchi
COMMON/paros2/cosh2m,sinh2m,dmain
COMMON/paros3/dmutan
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+8
x=xx
pnu=pnua
ierro=0
IF (x >= pnu) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCC
!CC  MONOTONIC REGION   CC
!CCCCCCCCCCCCCCCCCCCCCCCCC
  sinth=pnu/x
  thet=ASIN(sinth)
  costh=COS(thet)
  fdomin=x*(costh+thet*sinth)
  IF (ifac == 1) THEN
    IF (-fdomin <= LOG(under)) ierro=1
  END IF
  IF (ierro == 0) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF KIA:
!CC    CALL TO THE TRAPEZOIDAL ROUTINE  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    CALL trapre(1,hir)
    IF (ifac == 1) THEN
      dkinf=0.5D0*hir*EXP(-fdomin)
    ELSE
      dkinf=0.5D0*hir
    END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF KIA':
!CC    CALL TO THE TRAPEZOIDAL ROUTINE  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    CALL trapre(10,hird)
    IF (ifac == 1) THEN
      dkind=-hird*EXP(-fdomin)
    ELSE
      dkind=-hird
    END IF
  ELSE
    dkinf=0.d0
    dkind=0.d0
  END IF
ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  OSCILLATORY REGION   CC
!CCCCCCCCCCCCCCCCCCCCCCCCCCC
  pi=ACOS(-1.d0)
  IF (ifac == 1) THEN
    IF ((-pi*pnu*0.5D0) <= LOG(under)) ierro=1
  END IF
  IF (ierro == 0) THEN
    coshm=pnu/x
    dmu=LOG(coshm+SQRT((coshm-1.d0)*(coshm+1.d0)))
    dmu2=2.d0*dmu
    cosh2m=COSH(dmu2)
    sinhm=SINH(dmu)
    sinh2m=SINH(dmu2)
    IF (dmu > 0.1D0) THEN
      chi=x*sinhm-pnu*dmu
      dmufac=dmu*coshm-sinhm
    ELSE
      dmu2=dmu*dmu
      dmu3=dmu2*dmu
      dmu5=dmu3*dmu2
      dmu7=dmu5*dmu2
      dmu9=dmu7*dmu2
      chi=-2.d0*x*(1.d0/6.d0*dmu3+1.d0/60.d0*dmu5+  &
          3.d0/5040.d0*dmu7+4.d0/362880.d0*dmu9)
      dmufac=dmu3/3.d0+dmu5/30.d0+dmu7/840.d0+dmu9/45360.d0
    END IF
    dmutan=dmu-TANH(dmu)
    coschi=COS(chi)
    sinchi=SIN(chi)
    pinu=pi*pnu
    dmain=pinu*0.5D0
!CCCCCCCCCCCCCCCCC
!CC  INTEGRALS  CC
!CCCCCCCCCCCCCCCCC
    CALL trapre(2,fos1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC THEN, THE KIA FUNCTION IS GIVEN BY ... CC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    IF (ifac == 1) THEN
      df1=EXP(-dmain)
      contr1=df1*fos1
    ELSE
      df1=1.d0
      contr1=df1*fos1
    END IF
    dkinf=contr1
!CCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF KIA' CC
!CCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCC
!CC  INTEGRALS  CC
!CCCCCCCCCCCCCCCCC
    CALL trapre(20,fosd1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC THEN, KIA' IS GIVEN BY ...
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    IF (ifac == 1) THEN
      df1=EXP(-dmain)
      contr1=df1*fosd1
    ELSE
      df1=1.d0
      contr1=df1*fosd1
    END IF
    dkind=contr1
  ELSE
    dkinf=0.d0
    dkind=0.d0
  END IF
END IF
RETURN
END SUBROUTINE dkint

DOUBLE PRECISION FUNCTION fa(u)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC     MONOTONIC PART: KIA(X) FUNCTION
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC        LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PHIR(U) : X*COSH(U)*COS(V(U))+PNU*V(U)
!           -FDOMIN, WHERE
!           V(U)=ASIN(U/SINH(U)*SIN(THET))
!           AND FDOMIN=X*(COS(THET)+
!                      THET*SIN(THET))
! U       : ARGUMENT OF THE FA FUNCTION
! UNDER   : UNDERFLOW NUMBER
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


DOUBLE PRECISION, INTENT(IN)         :: u
DOUBLE PRECISION :: phir,d1mach,under
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+8
IF ((-phir(u)) <= LOG(under)) THEN
  fa=0.d0
ELSE
  fa=EXP(-phir(u))
END IF
RETURN
END FUNCTION fa

DOUBLE PRECISION FUNCTION fad(t)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC MONOTONIC PART: (KIA'(X) FUNCTION)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC        LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ANG      :  THET-V(U)
! ANGH     :  ANG/2
! COSTH    :  COS(THET)
! COSVU    :  COS(V(U))
! DFAC     :  COS(THET)+(FU1+2.D0*SIN(ANG/2)*SIN(ANG/2))/COS(VU)
! DJACO    :  COSH(T)*EXP(S)/(1.D0+EXP(S))
! EXPO     :  EXP(-PHIR(U))
! FDOMIN   :  X*(COS(THET)+THET*SIN(THET))
! FU1      :  2*SINH(U/2)**2
! FUAC     :  U**3/6+U**5/120+U**7/5040
! PHIR(U)  :  X*COSH(U)*COS(V(U))+PNU*V(U)-FDOMIN
! S        :  SINH(T)
! SINANH   :  SIN(ANG/2)
! SINHU    :  SINH(U)
! SINHUH   :  SINH(U/2)
! SINTH    :  SIN(THET)
! T        :  INTEGRATION ABCISSA
! THET     :  ASIN(PNU/X)
! U        :  LOG(1+EXP(S))
! U2       :  U**2
! U3       :  U**3
! U5       :  U**5
! U7       :  U**7
! UH       :  U/2
! UNDER    :  UNDERFLOW NUMBER
! V(U)     :  ASIN(U/SINH(U)*SIN(THET))
! VU       :  V(U)
! Y        :  EXP(S)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)         :: t
DOUBLE PRECISION :: d1mach,under,v,phir,s,y,u,  &
    djaco,dfac,uh,sinhuh,fu1,u2,u3,u5,u7,fuac, sinhu,vu,cosvu,ang,angh,sinanh,expo
DOUBLE PRECISION :: thet,sinth,costh,fdomin
COMMON/parmon/thet,sinth,costh,fdomin
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+8
s=SINH(t)
y=EXP(s)
u=LOG(1.d0+y)
djaco=COSH(t)*y/(1.d0+y)
IF (ABS(u) < 1.d-1) THEN
  IF (ABS(u) < under) THEN
    dfac=costh
  ELSE
    uh=u*0.5D0
    sinhuh=SINH(uh)
    fu1=2.d0*sinhuh*sinhuh
    u2=u*u
    u3=u2*u
    u5=u2*u3
    u7=u5*u2
    fuac=u3/6.d0+u5/120.d0+u7/5040.d0
    sinhu=SINH(u)
    vu=v(u)
    cosvu=COS(vu)
    ang=ASIN(-sinth/(costh*u/sinhu+cosvu)* (sinhu+u)*fuac/(sinhu*sinhu))
    angh=0.5D0*ang
    sinanh=SIN(angh)
    dfac=costh+(fu1+2.d0*sinanh*sinanh)/cosvu
  END IF
ELSE
  vu=v(u)
  ang=thet-vu
  fu1=COSH(u)-1.d0
  angh=0.5D0*ang
  sinanh=SIN(angh)
  dfac=costh+(fu1+2.d0*sinanh*sinanh)/COS(vu)
END IF
IF ((-phir(u)) <= LOG(under)) THEN
  expo=0.d0
ELSE
  expo=EXP(-phir(u))
END IF
fad=expo*dfac*djaco
RETURN
END FUNCTION fad

DOUBLE PRECISION FUNCTION fdtau2(t)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   OSCILLATORY PART: TAU CONTRIBUTION ROUTINE
!CC     WE USE THIS INTEGRAL FOR LARGE PNU
!CC   (KIA'(X) FUNCTION). SEMI-INFINITE INTEGRAL.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC         LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ARGU    : (COSH(MU)*U-DMUFAC)/SINH(U)
! ARGU2   : ARGU**2
! COSCHI  : COSH(CHI), CHI=X*SINH(MU)-PNU*MU
! COSH2M  : COSH(2*MU)
! COSHM   : COSH(MU)
! COSHU   : COSH(U)
! COSS    : COS(SIGMA(U))
! D1      : COSH(MU)-DMUFAC
! DELTA   : -SIN(CHI)*COS(SIGMA(U))*COSH(U)
!           +COS(CHI)*SIN(SIGMA(U))*SINH(U)
! DERI    : (COSH(MU)/SINH(U)-D1*COSH(U)/SINH(U)**2)
!           /SQRT(1-ARGU2)
! DJACO   : COSH(T)*EXP(S)/SQRT(1+EXP(S)**2)
! DMAIN   : PI*PNU*0.5
! DMU     : SOLUTION MU OF COSH(MU)=PNU/X
! DMU2    : 2*MU
! DMUFAC  : MU*COSH(MU)-SINH(MU)
! EXPON   : EXP(-(PHIB(U)-PI*PNU*0.5))
! F1      : SINH(U-MU)/(U-MU)
! G1      : Z/6+Z3/120+Z5/5040+Z7/362880, Z=2*Y
! GAMMA   : SIN(CHI)*SIN(SIGMA(U))*SINH(U)+
!           COS(CHI)*COS(SIGMA(U))*COSH(U)
! PHIB(U) : X*COSH(U)*COS(SIGMA(U))+PNU*SIGMA(U),
!           WHERE X=ARGUMENT OF THE FUNCTION,
!           SIGMA(U)=ARCSIN((COSH(MU)*U-DMUFAC)/SINH(U)
! RESTO   : -SIN(CHI)*COS(SIGMA(U))*COSH(U)
!           +COS(CHI)*SIN(SIGMA(U))*SINH(U)+
!           (SIN(CHI)*SIN(SIGMA(U))*SINH(U)+
!            COS(CHI)*COS(SIGMA(U))*COSH(U))*DERI
! S       : SINH(T)
! SIGMA(U): ASIN((COSH(MU)*U-DMUFAC)/SINH(U))
! SIGMAU  : SIGMA(U)
! SINCHI  : SIN(CHI)
! SINH2M  : SINH(2*MU)
! SINHM   : SINH(MU)
! SINHU   : SINH(U)
! SINHU2  : 2*SINH(U)
! SINS    : SIN(SIGMA(U))
! T       : INTEGRATION ABCISSA
! U       : MU+LOG(X+SQRT(X**2+1))
! UNDER   : UNDERFLOW NUMBER
! X       : EXP(S)
! Y       : U-MU
! Z       : 2*Y
! Z2      : Z**2
! Z3      : Z**3
! Z5      : Z**5
! Z7      : Z**7
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)         :: t
DOUBLE PRECISION :: phib,d1mach,under,sigma,d1,sigmau
DOUBLE PRECISION :: s,x,u,y,coshu,sinhu,f1,z,z2,z3,  &
    z5,z7,g1,dmu2,deri,argu,argu2,sinhu2,djaco, expon,sins,coss,gamma,delta,resto
DOUBLE PRECISION :: dmu,coshm,sinhm,dmufac,coschi,sinchi
DOUBLE PRECISION :: cosh2m,sinh2m,dmain,dmutan
COMMON/paros1/dmu,coshm,sinhm,dmufac,coschi,sinchi
COMMON/paros2/cosh2m,sinh2m,dmain
COMMON/paros3/dmutan
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+8
s=SINH(t)
x=EXP(s)
u=dmutan+LOG(x+SQRT(x*x+1.d0))
y=u-dmu
coshu=COSH(u)
sinhu=SINH(u)
IF (ABS(y) <= 1.d-1) THEN
  IF (ABS(y) > under) THEN
    f1=SINH(y)/y
  ELSE
    f1=1.d0
  END IF
  z=y*2.d0
  z2=z*z
  z3=z2*z
  z5=z3*z2
  z7=z5*z2
  g1=z/6.d0+z3/120.d0+z5/5040.d0+z7/362880.d0
  dmu2=2.d0*dmu
  deri=sinhu/(f1-coshm*coshu)*SQRT(COSH(dmu2)*f1*f1+  &
      2.d0*SINH(dmu2)*g1-coshm*coshm)
  deri=1.d0/deri
ELSE
  d1=coshm*u-dmufac
  argu=d1/sinhu
  argu2=argu*argu
  sinhu2=sinhu*sinhu
  deri=1.d0/SQRT(1.d0-argu2)*(coshm/sinhu- d1*coshu/sinhu2)
  IF (u < dmu) deri=-deri
END IF
djaco=COSH(t)*x/SQRT(1.d0+x*x)
IF ((-(phib(u)-dmain)) <= LOG(under)) THEN
  expon=0.d0
  fdtau2=0.d0
ELSE
  expon=EXP(-(phib(u)-dmain))
  sigmau=sigma(u)
  sins=SIN(sigmau)
  coss=COS(sigmau)
  gamma=coschi*sins*sinhu-sinchi*coss*coshu
  delta=-coschi*coss*coshu-sinchi*sins*sinhu
  resto=delta+gamma*deri
  fdtau2=expon*resto*djaco
END IF
RETURN
END FUNCTION fdtau2

DOUBLE PRECISION FUNCTION fstau2(t)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC OSCILLATORY PART: TAU CONTRIBUTION ROUTINE
!CC   WE USE THIS INTEGRAL FOR LARGE PNU
!CC (KIA(X) FUNCTION). SEMI-INFINITE INTEGRAL.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC        LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ARGU    : (COSH(MU)*U-DMUFAC)/SINH(U)
! ARGU2   : ARGU**2
! COSCHI  : COSH(CHI), CHI=X*SINH(MU)-PNU*MU
! COSH2M  : COSH(2*MU)
! COSHM   : COSH(MU)
! D1      : COSH(MU)-DMUFAC
! DERI    : (COSH(MU)/SINH(U)-D1*COSH(U)/SINH(U)**2)
!           /SQRT(1-ARGU2)
! DJACO   : COSH(T)*EXP(S)/SQRT(1+EXP(S)**2)
! DMAIN   : PI*PNU*0.5
! DMU     : SOLUTION MU OF COSH(MU)=PNU/X
! DMU2    : 2*MU
! DMUFAC  : MU*COSH(MU)-SINH(MU)
! EXPON   : EXP(-(PHIB(U)-PI*PNU*0.5))
! F1      : SINH(U-MU)/(U-MU)
! G1      : Z/6+Z3/120+Z5/5040+Z7/362880, Z=2*Y
! PHIB(U) : X*COSH(U)*COS(SIGMA(U))+PNU*SIGMA(U),
!           WHERE X=ARGUMENT OF THE FUNCTION,
!           SIGMA(U)=ARCSIN((COSH(MU)*U-DMUFAC)/SINH(U)
! RESTO   : COS(CHI)+SIN(CHI)*DERI
! S       : SINH(T)
! SINCHI  : SIN(CHI)
! SINH2M  : SINH(2*MU)
! SINHM   : SINH(MU)
! SINHU   : SINH(U)
! SINHU2  : 2*SINH(U)
! T       : INTEGRATION ABCISSA
! U       : MU+LOG(X+SQRT(X**2+1))
! UNDER   : UNDERFLOW NUMBER
! X       : EXP(S)
! Y       : U-MU
! Z       : 2*Y
! Z2      : Z**2
! Z3      : Z**3
! Z5      : Z**5
! Z7      : Z**7
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

DOUBLE PRECISION, INTENT(IN)         :: t
DOUBLE PRECISION :: d1mach,under, s,x,u,y,f1,z,z2,  &
    z3,z5,z7,g1,dmu2,deri,d1,sinhu,argu,argu2,sinhu2, djaco,phib,expon,resto
DOUBLE PRECISION :: dmu,coshm,sinhm,dmufac,coschi,sinchi
DOUBLE PRECISION :: cosh2m,sinh2m,dmain,dmutan
COMMON/paros1/dmu,coshm,sinhm,dmufac,coschi,sinchi
COMMON/paros2/cosh2m,sinh2m,dmain
COMMON/paros3/dmutan
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+8
s=SINH(t)
x=EXP(s)
u=dmutan+LOG(x+SQRT(x*x+1.d0))
y=u-dmu
IF (ABS(y) > under) THEN
  f1=SINH(y)/y
ELSE
  f1=1.d0
END IF
z=y*2.d0
z2=z*z
IF (ABS(y) <= 0.1D0) THEN
  z3=z2*z
  z5=z3*z2
  z7=z5*z2
  g1=z/6.d0+z3/120.d0+z5/5040.d0+z7/362880.d0
  dmu2=2.d0*dmu
  deri=SINH(u)/(f1-coshm*COSH(u))*  &
      SQRT(COSH(dmu2)*f1*f1+2.d0*SINH(dmu2)*g1- coshm*coshm)
  deri=1.d0/deri
ELSE
  d1=coshm*u-dmufac
  sinhu=SINH(u)
  argu=d1/sinhu
  argu2=argu*argu
  sinhu2=sinhu*sinhu
  deri=1.d0/SQRT(1.d0-argu2)* (coshm/sinhu-d1*COSH(u)/sinhu2)
  IF (u < dmu) deri=-deri
END IF
djaco=COSH(t)*x/SQRT(1.d0+x*x)
IF ((-(phib(u)-dmain)) <= LOG(under)) THEN
  expon=0.d0
  fstau2=0.d0
ELSE
  expon=EXP(-(phib(u)-dmain))
  resto=(coschi+sinchi*deri)
  fstau2=expon*resto*djaco
END IF
RETURN
END FUNCTION fstau2

SUBROUTINE trapre(ic,ti)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC IMPLEMENTATION OF AN ADAPTIVE TRAPEZOIDAL RULE
!CC
!CC   INPUT:
!CC          IC,   CHOICE OF THE INTEGRAND
!CC   OUTPUT:
!CC          TI,  CONTRIBUTION TO KIA(X)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC            LIST OF VARIABLES
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! A        :  LOWER INTEGRATION LIMIT
! B        :  UPPER INTEGRATION LIMIT
! DELTA    :  CALCULATES THE RELATIVE PRECISION
! EPS      :  RELATIVE PRECISION PARAMETER USED IN
!             THE CALCULATION
! H        :  INTEGRATION STEP
! IC       :  CHOICE OF THE INTEGRAND
! PNU      :  ORDER OF THE FUNCTION
! SUM      :  ACCUMULATES THE ELEMENTARY
!             CONTRIBUTIONS
! TI, TIN  :  EVALUATED INTEGRAL
! X        :  ARGUMENT
! XAC      :  INTEGRATION ABCISSA
! Z        :  X/PNU
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN)                      :: ic
DOUBLE PRECISION, INTENT(OUT)            :: ti
DOUBLE PRECISION :: eps,a,b,h,delta,sum,xac,  tin,finte,d1mach
DOUBLE PRECISION :: x,pnu,z
INTEGER :: i,ifi, n
COMMON/argu/x,pnu

eps=d1mach(3)
IF (eps < 1.d-14) eps=1.d-14
n=0
!CCCC   INTEGRATION LIMITS: A,B
z=x/pnu
IF ((z > 0.999D0).AND.(z < 1.001D0)) THEN
  IF ((ic == 1).OR.(ic == 10)) THEN
    a=-3.5D0
  ELSE
    a=-4.5D0
  END IF
ELSE
  IF ((ic == 2).OR.(ic == 20)) THEN
    a=-2.5D0
  ELSE
    a=-4.5D0
  END IF
END IF
b=-a
h=b-a
ti=0.5D0*h*(finte(ic,a)+finte(ic,b))
delta=1.d0+eps
11     n=n+1
h=0.5D0*h
IF (n == 1) THEN
  ifi=1
ELSE
  ifi=2*ifi
END IF
sum=0.d0
DO  i=1,ifi
  xac=a+(2*i-1)*h
  sum=sum+finte(ic,xac)
END DO
tin=0.5D0*ti+h*sum
IF ((tin /= 0).AND.(n > 4)) THEN
  delta=ABS(1.d0-ti/tin)
END IF
ti=tin
IF ((delta > eps).AND.(n < 9)) GO TO 11
RETURN
END SUBROUTINE trapre

FUNCTION finte(ic,t)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   INPUT:                                      C
!CC      T,      INTEGRATION ABSCISSA             C
!CC   OUTPUT :                                    C
!CC      FINTE,  FUNCTION                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC       LIST OF VARIABLES
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  FA      : MONOTONIC PART CONTRIBUTION
!            (KIA(X) FUNCTION).
!  FAD     : MONOTONIC PART CONTRIBUTION
!            (KIA'(X) FUNCTION).
!  FDTAU2  : OSCILLATORY PART: TAU CONTRIBUTION ROUTINE
!            (KIA'(X) FUNCTION). SEMI-INFINITE INTEGRAL.
!                USED FOR LARGE PNU
!  FINTE   : ELEMENTARY CONTRIBUTION
!  FSTAU2  : OSCILLATORY PART: TAU CONTRIBUTION ROUTINE
!            (KIA(X) FUNCTION). SEMI-INFINITE INTEGRAL.
!                  USED FOR LARGE PNU
!  IC      : CHOICE OF FUNCTION
!  T       : INTEGRATION ABCISSA
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN)                  :: ic
DOUBLE PRECISION, INTENT(IN)         :: t
DOUBLE PRECISION :: finte
DOUBLE PRECISION :: fa,fad,fstau2,fdtau2


IF (ic == 1) THEN
  finte=fa(t)
END IF
IF (ic == 2) THEN
  finte=fstau2(t)
END IF
IF (ic == 10) THEN
  finte=fad(t)
END IF
IF (ic == 20) THEN
  finte=fdtau2(t)
END IF
END FUNCTION finte

SUBROUTINE frakia(ifac,x,pnu,pser,pserd,ierro)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   IMPLEMENTATION OF THE CONTINUED FRACTION     C
!CC   METHOD FOR THE CALCULATION OF KIA AND KIA'   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC           LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! A      : PNU, ORDER OF THE FUNCTIONS
! A2     : A**2
! AA     : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! AB     : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! B      : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! CC     : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! COSTH  : COS(THET)
! D      : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! DELS   : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! DELTA  : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! FC     : CONTINUED FRACTION FOR K(PNU+1)/K(PNU)
! FDOMIN : X*(COS(THET)+THET*SIN(THET))
! IERRO  : FLAG FOR UNDERFLOW/OVERFLOW ERRORS
! K      : KIA(X) FUNCTION
! KP     : KIA'(X) FUNCTION
! PI     : 3.1415..
! PIA    : PI*PNU
! PISQ   : SQRT(PI)
! PNU    : ORDER OF THE FUNCTIONS
! PRECI  : RELATIVE PRECISION USED IN THE CALCULATION
! PSER   : KIA(X) FUNCTION
! PSERD  : KIA'(X) FUNCTION
! Q0B    : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! Q1B    : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! QB     : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! QQB    : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! S      : AUXILIARY VARIABLE IN THE IMPLEMENTATION
!          OF THE LENZ-THOMPSON ALGORITHM
! SINTH  : SIN(THET)
! THET   : ASIN(PNU/X)
! UNDER  : UNDERFLOW NUMBER
! X      : ARGUMENT OF THE FUNCTIONS
! Z0     : 1/(1+S)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN)                      :: ifac
DOUBLE PRECISION, INTENT(IN)             :: x
DOUBLE PRECISION, INTENT(IN)             :: pnu
DOUBLE PRECISION, INTENT(OUT)            :: pser
DOUBLE PRECISION, INTENT(OUT)            :: pserd
INTEGER, INTENT(OUT)                     :: ierro
DOUBLE PRECISION :: k,kp
DOUBLE PRECISION :: pi,d1mach,preci,a,a2,pisq,ab,aa,cc,q0b,  &
    q1b,qb,qqb,b,d,delta,fc,s,z0,dels,pia,sinth,thet, costh,fdomin,under
INTEGER :: mm
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+8
pi=ACOS(-1.d0)
ierro=0
IF (ifac == 1) THEN
  IF (x >= pnu) THEN
    sinth=pnu/x
    thet=ASIN(sinth)
    costh=COS(thet)
    fdomin=x*(costh+thet*sinth)
    IF (-fdomin <= LOG(under)) ierro=1
  ELSE
    IF ((-pi*pnu*0.5D0) <= LOG(under)) ierro=1
  END IF
END IF
IF (ierro == 0) THEN
  preci=d1mach(3)*10
  a=pnu
  a2=a*a
  pisq=pi**0.5D0
!CC CF FOR K(NU+1)/K(NU)
  mm=2
  ab=-(0.25D0+a2)
  aa=ab
  cc=-ab
  q0b=0.d0
  q1b=1.d0
  qqb=cc
  b=2.d0*(1.d0+x)
  d=1.d0/b
  delta=d
  fc=delta
  s=qqb*delta
  ab=-2.d0+ab
  b=b+2.d0
  91      d=1.d0/(b+ab*d)
  delta=(b*d-1.d0)*delta
  fc=fc+delta
  cc=-ab*cc/mm
  qb=(q0b-(b-2.d0)*q1b)/ab
  q0b=q1b
  q1b=qb
  qqb=qqb+cc*q1b
  dels=qqb*delta
  s=s+dels
  b=b+2.d0
  ab=-2.d0*mm+ab
  mm=mm+1
  IF (mm < 10000) THEN
    IF (ABS(dels/s) > preci) GO TO 91
  END IF
  z0=1.d0/(1.d0+s)
  IF (ifac == 1) THEN
    k=pisq*(2.d0*x)**(-0.5D0)*EXP(-x)*z0
    kp=-k/x*(0.5D0+x-(a2+0.25D0)*fc)
  ELSE
    IF (x < a) THEN
      pia=pi*a
      k=pisq*(2.d0*x)**(-0.5D0)*EXP(-x+pia*0.5D0)*z0
      kp=-k/x*(0.5D0+x-(a2+0.25D0)*fc)
    ELSE
      sinth=a/x
      thet=ASIN(sinth)
      costh=COS(thet)
      fdomin=x*(costh+thet*sinth)
      k=pisq*(2.d0*x)**(-0.5D0)*EXP(-x+fdomin)*z0
      kp=-k/x*(0.5D0+x-(a2+0.25D0)*fc)
    END IF
  END IF
  pser=k
  pserd=kp
ELSE
  pser=0.d0
  pserd=0.d0
END IF
RETURN
END SUBROUTINE frakia

SUBROUTINE serkia(ifac,x,pnu,pser,pserd,ierro)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   CALCULATION OF SERIES FOR K, K'
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC         LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! A       : ORDER OF THE FUNCTIONS
! A2      : A**2
! A2H     : A**2/2
! A2N     : A**(2*N)
! ACCP    : ACCUMULATES THE P COEFFICIENTS
! ACCQ    : ACCUMULATES THE Q COEFFICIENTS
! ARGU    : SIG(0)-A*LOG(X/2)
! C       : X**(2*K)/K!
! COCI    : 1/(K**2+A**2)
! COSTH   : COS(THET)
! DELTAK  : ACCUMULATES THE CONTRIBUTION FOR THE
!           KIA(X) FUNCTION
! DELTKP  : ACCUMULATES THE CONTRIBUTION FOR THE
!           KIA'(X) FUNCTION
! DF1     : FACTOR (DEPENDING ON IFAC)
! ETA0    : PARAMETER FOR THE CALCULATION OF THE
!           COULOMB PHASE SHIFT
! ETA02   : ETA0**2
! F(K)    : SIN(PHI(A,K)-A*LOG(X/2))
!           /(A**2*(1+A**2)...(K**2+A**2))**1/2,
!           WHERE PHI(A,K)=PHASE(GAMMA(1+K+IA))
! FDOMIN  : X*(COS(THET)+THET*SIN(THET))
! IERRO   : FLAG FOR UNDERFLOW/OVERFLOW ERRORS
! K       : CONTRIBUTION TO THE KIA(X) FUNCTION
! KP      : CONTRIBUTION TO THE KIA'(X) FUNCTION
! OVER    : OVERFLOW NUMBER
! P0      : PARAMETERS FOR THE CALCULATION OF THE COULOMB
!           PHASE SHIFT
! P1      : PARAMETERS FOR THE CALCULATION OF THE COULOMB
!           PHASE SHIFT
! P2      : PARAMETERS FOR THE CALCULATION OF THE COULOMB
!           PHASE SHIFT
! PI      : 3.1415...
! PIA     : PI*A
! PIA2    : 2*PI*A
! PNU     : ORDER OF THE FUNCTIONS
! PRECI   : RELATIVE PRECISION USED IN THE CALCULATION
! PSER    : KIA(X) FUNCTION
! PSERD   : KIA'(X) FUNCTION
! Q0      : PARAMETERS FOR THE CALCULATION OF THE COULOMB
!           PHASE SHIFT
! Q1      : PARAMETERS FOR THE CALCULATION OF THE COULOMB
!           PHASE SHIFT
! Q2      : PARAMETERS FOR THE CALCULATION OF THE COULOMB
!           PHASE SHIFT
! R(K)    : F(K)*A/TAN(PHI(A,K)-A*LOG(X/2))
! SIG0    : COULOMB PHASE SHIFT
! SINTH   : SIN(THET)
! THET    : ASIN(A/X)
! UNDER   : UNDERFLOW NUMBER
! X       : ARGUMENT OF THE FUNCTIONS
! X2      : X*X
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN)                      :: ifac
DOUBLE PRECISION, INTENT(IN)             :: x
DOUBLE PRECISION, INTENT(IN)             :: pnu
DOUBLE PRECISION, INTENT(OUT)            :: pser
DOUBLE PRECISION, INTENT(OUT)            :: pserd
INTEGER, INTENT(OUT)                     :: ierro

DOUBLE PRECISION :: p0(0:9),q0(0:9),p1(0:9),q1(0:9),  &
    p2(0:6),q2(0:6),r(0:500),f(0:500)
DOUBLE PRECISION :: pi,eta0,preci,d1mach,a,eta02,a2,  &
    accp,accq,sig0,a2n,pia,pia2,argu,c,k,kp,x2,deltak,dee,  &
    deltkp,coci,df1,sinth,thet,costh,fdomin,over,under
DOUBLE PRECISION :: dds,dsmall,euler
INTEGER :: n,m
INTEGER :: l
SAVE p0,q0,p1,q1,p2,q2
DATA p0/1.08871504904797411683D5,3.64707573081160914640D5,  &
    4.88801471582878013158D5,3.36275736298197324009D5,  &
    1.26899226277838479804D5,  &
    2.60795543527084582682D4,2.73352480554497990544D3,  &
    1.26447543569902963184D2, 1.85446022125533909390D0,1.90716219990037648146D-3/
DATA q0/6.14884786346071135090D5,2.29801588515708014282D6,  &
    3.50310844128424021934D6,  &
    2.81194990286041080264D6,1.28236441994358406742D6,  &
    3.35209348711803753154D5,  &
    4.84319580247948701171D4,3.54877039006873206531D3,  &
    1.11207201299804390166D2,1.d0/
DATA p1/-1.044100987526487618670D10,-1.508574107180079913696D10,  &
    -5.582652833355901160542D9,4.052529174369477275446D8,  &
    5.461712273118594275192D8,  &
    9.510404403068169395714D7,6.281126609997342119416D6,  &
    1.651178048950518520416D5,  &
    1.498824421329341285521D3,2.974686506595477984776D0/
DATA q1/1.808868161493543887787D10,3.869142051704700267785D10,  &
    3.003264575147162634046D10,1.075554651494601843525D10,  &
    1.901298501823290694245D9,  &
    1.665999832151229472632D8,6.952188089169487375936D6,  &
    1.253235080625688652718D5,7.904420414560291396996D2,1.d0/
DATA p2/7.08638611024520906826D-3,-6.54026368947801591128D-2,  &
    2.92684143106158043933D-1,4.66821392319665609167D0,  &
    -3.43943790382690949054D0,  &
    -7.72786486869252994370D0,-9.88841771200290647461D-01/
DATA q2/-7.08638611024520908189D-3,6.59931690706339630254D-2,  &
    -2.98754421632058618922D-1,-4.63752355513412248006D0,  &
    3.79700454098863541593D0,7.06184065426336718524D0,1.d0/

pi=ACOS(-1.d0)
eta0=1.8055470716051069198764D0
euler=0.577215664901532860606512D0
preci=d1mach(3)*10
!CC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
over=d1mach(2)*1.d-8
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+8
ierro=0
IF (ifac == 1) THEN
  IF (x >= pnu) THEN
    sinth=pnu/x
    thet=ASIN(sinth)
    costh=COS(thet)
    fdomin=x*(costh+thet*sinth)
    IF (-fdomin <= LOG(under)) ierro=1
  ELSE
    IF ((-pi*pnu*0.5D0) <= LOG(under)) ierro=1
  END IF
END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COEFFICIENTS FOR THE CALCULATION OF THE COULOMB PHASE SHIFT
!CC FROM CODY & HILLSTROM, MATH. COMPUT. 24(111) 1970
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
IF (ierro == 0) THEN
  a=pnu
  eta02=eta0*eta0
  a2=a*a
  n=0
  accp=0.d0
  accq=0.d0
  IF (a <= 2.d0) THEN
    33        a2n=a2**n
    accp=accp+p0(n)*a2n
    accq=accq+q0(n)*a2n
    n=n+1
    IF (n <= 9) GO TO 33
    sig0=a*(a2-eta02)*accp/accq
  ELSE
    IF ((a > 2.d0).AND.(a <= 4.d0)) THEN
      44          a2n=a2**n
      accp=accp+p1(n)*a2n
      accq=accq+q1(n)*a2n
      n=n+1
      IF (n <= 9) GO TO 44
      sig0=a*accp/accq
    ELSE
      55          a2n=a2**n
      accp=accp+p2(n)/a2n
      accq=accq+q2(n)/a2n
      n=n+1
      IF (n <= 6) GO TO 55
      sig0=ATAN(a)*0.5D0+a*(LOG(1.d0+a2)*0.5D0+accp/accq)
    END IF
  END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC EVALUATION OF F(0), R(0), R(1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  pia=pi*a
  argu=sig0-a*LOG(x*0.5D0)
  r(0)=COS(argu)
  r(1)=1.d0/(1.d0+a2)*(r(0)-a*SIN(argu))
  IF (a < under) THEN
    f(0)=-(euler+LOG(x*0.5D0))
  ELSE
    f(0)=1.d0/a*SIN(argu)
  END IF
  c=1.d0
  k=f(0)
  kp=-0.5D0*r(0)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC RECURSION FOR F(K), R(K), C(K)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  x2=0.25D0*x*x
  m=1
  coci=1.d0/(m*m+a2)
  deltak=k*10
  deltkp=kp*10
  66      f(m)=(m*f(m-1)+r(m-1))*coci
  c=x2*c/m
  deltak=f(m)*c
  k=k+deltak
  deltkp=(m*f(m)-0.5D0*r(m))*c
  kp=kp+deltkp
  m=m+1
  coci=1.d0/(m*m+a2)
  r(m)=((2.d0*m-1.d0)*r(m-1)-r(m-2))*coci
  IF (k > over) m=500
  IF (kp > over) m=500
  IF (m < 500) THEN
    IF ((ABS(deltak/k) > preci).OR. (ABS(deltkp/kp) > preci)) GO TO 66
  END IF
  pia2=2.d0*pia
  IF (ifac == 1) THEN
    IF (-pia2 <= LOG(under)) THEN
      dee=0.d0
    ELSE
      dee=EXP(-pia2)
    END IF
    IF (a < 1.d-1) THEN
      IF (a < under) THEN
        df1=1.d0
      ELSE
        l=0
        dds=1.d0
        dsmall=1.d0
        47            l=l+1
        dds=-dds*pia2/(l+1.d0)
        dsmall=dsmall+dds
        IF (ABS(dds/dsmall) > preci) GO TO 47
        df1=EXP(pia*0.5D0)*SQRT(dsmall)
      END IF
    ELSE
      df1=EXP(pia*0.5D0)*SQRT((1.d0-dee)/pia2)
    END IF
    pser=k/df1
    pserd=kp*2.d0/x/df1
  ELSE
    IF (x < a) THEN
      IF (-pia2 <= LOG(under)) THEN
        dee=0.d0
      ELSE
        dee=EXP(-pia2)
      END IF
      IF (a < 1.d-1) THEN
        IF (a < under) THEN
          df1=1.d0
        ELSE
          l=0
          dds=1.d0
          dsmall=1.d0
          48              l=l+1
          dds=-dds*pia2/(l+1.d0)
          dsmall=dsmall+dds
          IF (ABS(dds/dsmall) > preci) GO TO 48
          df1=SQRT(dsmall)
        END IF
      ELSE
        df1=SQRT((1.d0-dee)/pia2)
      END IF
    ELSE
      sinth=a/x
      thet=ASIN(sinth)
      costh=COS(thet)
      fdomin=x*(costh+thet*sinth)
      IF (-pia2 <= LOG(under)) THEN
        dee=0.d0
      ELSE
        dee=EXP(-pia2)
      END IF
      IF (a < 1.d-1) THEN
        IF (a < under) THEN
          df1=1.d0
        ELSE
          l=0
          dds=1.d0
          dsmall=1.d0
          49              l=l+1
          dds=-dds*pia2/(l+1.d0)
          dsmall=dsmall+dds
          IF (ABS(dds/dsmall) > preci) GO TO 49
          df1=EXP(pia*0.5D0-fdomin)*SQRT(dsmall)
        END IF
      ELSE
        df1=EXP(pia*0.5D0-fdomin)*SQRT((1.d0-dee)/pia2)
      END IF
    END IF
    pser=k/df1
    pserd=kp*2.d0/x/df1
  END IF
ELSE
  pser=0.d0
  pserd=0.d0
END IF
RETURN
END SUBROUTINE serkia

SUBROUTINE aiexki(ifac,x,a,dkai,dkaid,ierrok)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC AIRY-TYPE ASYMPTOTIC EXPANSION FOR THE KIA
!CC AND KIA' FUNCTIONS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC           LIST OF VARIABLES:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! A       : ORDER OF THE FUNCTIONS
! A0EX    : EXACT VALUE OF THE COEFFICIENT A_0
! A13     : A**(1/3)
! A2      : A**2
! A23     : A**(2/3)
! A2K     : A**(2*K)
! AC      : COEFFICIENTS A_S(PSI) USED IN THE EXPANSIONS
! AII     : IMAGINARY PART OF THE AIRY FUNCTION AI(Z)
! AIR     : REAL PART OF THE AIRY FUNCTION AI(Z)
! APIHAL  : A*PI/2
! ARG     : -A**(2/3)*PSI
! AS      : ACCUMULATES THE CONTRIBUTION OF THE AC COEFFICIENTS
!           (FOR THE KIA(X) FUNCTION)
! ASP     : ACCUMULATES THE CONTRIBUTION OF THE AC COEFFICIENTS
!           (FOR THE KIA'(X) FUNCTION)
! B0EX    : EXACT VALUE OF THE COEFFICIENT B_0
! B0PEX   : EXACT VALUE OF THE COEFFICIENT B'_0
! BC      : COEFFICIENTS B_S(PSI) USED IN THE EXPANSIONS
! BS      : ACCUMULATES THE CONTRIBUTION OF THE BC COEFFICIENTS
!           (FOR THE KIA(X) FUNCTION)
! BSO     : ACCUMULATES THE OLD VALUES OF BS
! BSP     : ACCUMULATES THE CONTRIBUTION OF THE BC COEFFICIENTS
!           (FOR THE KIA'(X) FUNCTION)
! BSPO    : ACCUMULATES THE OLD VALUES OF BSP
! C0EX    : EXACT VALUE OF THE COEFFICIENT C_0
! CHI     : ACCUMULATES THE CONTRIBUTION OF THE CHIN COEFFICIENTS
! CHIEX   : EXACT VALUE OF THE FUNCTION CHI(PSI)
! CHIN    : COEFFICIENTS FOR THE EXPANSION OF THE FUNCTION
!           CHI(PSI)
! COSTH   : COS(THET)
! D0EX    : EXACT VALUE OF THE COEFFICIENT D_0
! DAII    : IMAGINARY PART OF THE AIRY FUNCTION AI'(Z)
! DAIR    : REAL PART OF THE AIRY FUNCTION AI'(Z)
! DF      : 2**(1/3)
! DKAI    : KIA(X) FUNCTION
! DKAID   : KIA'(X) FUNCTION
! DZZ     : (ABS(1-Z**2))**(1/2)
! ETA     : PSI/2**(1/3)
! ETAJ    : ETA**J
! ETAK    : ETA**K
! ETAL    : ETA**L
! EXPAM   : EXP(A*PI/2-FDOMIN)
! EXPAPI  : EXP(A*PI/2)
! F2      : (-)**K/A**(2*K)
! F4      : (A**(1/3))**4
! FAC     : FACTOR FOR THE KIA(X) FUNCTION
! FACD    : FACTOR FOR THE KIA'(X) FUNCTION
! FDOMIN  : X*(COS(THET)+THET*SIN(THET))
! IERRO   : ERROR FLAG FOR THE AIRY FUNCTIONS
! IFAC    :
!           *IF IFAC=1, COMPUTATION OF UNSCALED KIA(X),KIA'(X)
!           *IF IFAC=2, COMPUTATION OF SCALED KIA(X),KIA'(X)
! IFACA   :
!           *IF IFACA=1, COMPUTATION OF UNSCALED AIRY AI(Z),AI'(Z)
!           *IF IFACA=2, COMPUTATION OF SCALED AIRY AI(Z),AI'(Z)
! IFUN    :
!           *IF IFUN=1, COMPUTATION OF AIRY AI(Z)
!           *IF IFUN=2, COMPUTATION OF AIRY AI'(Z)
! PHI     : COEFFICIENTS FOR THE EXPANSION OF THE FUNCTION
!           PHI(PSI)
! PHIEX   : EXACT VALUE OF THE FUNCTION PHI(PSI)
! PHIEX2  : PHI(PSI)**2
! PHIS    : VALUE OF THE FUNCTION PHI(PSI)
! PI      : 3.1415...
! PIHALF  : PI/2
! PSI     :
!          *IF 0<=Z<=1, 2/3*PSI**3/2=LOG((1+SQRT(1-Z**2))/Z)
!                       -SQRT(1-Z**2)
!          *IF Z>=1, 2/3*(-PSI)**3/2=SQRT(Z**2-1)-ARCCOS(1/Z)
! PSI12   : SQRT(PSI)
! PSI2    : PSI*PSI
! PSI3    : PSI**3
! SAS     : ACCUMULATES THE CONTRIBUTION OF A_S(PSI)
! SBS     : ACCUMULATES THE CONTRIBUTION OF B_S(PSI)
! SCS     : ACCUMULATES THE CONTRIBUTION OF C_S(PSI)
! SDS     : ACCUMULATES THE CONTRIBUTION OF D_S(PSI)
! SIG     : (-)**K
! SINTH   : SIN(THET)
! THET    : ASIN(A/X)
! UNDER   : UNDERFLOW NUMBER
! X       : ARGUMENT OF THE FUNCTIONS
! Y       : Z-1
! Z       : X/A
! Z2      : Z**2
! Z21M    : 1-Z**2
! ZDS     : SQRT(1-Z**2)
! ZMASF   : EXPANSION IN TERMS OF (Z-1) FOR THE
!           CALCULATION OF PSI(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN)                      :: ifac
DOUBLE PRECISION, INTENT(IN)             :: x
DOUBLE PRECISION, INTENT(IN)             :: a
DOUBLE PRECISION, INTENT(OUT)            :: dkai
DOUBLE PRECISION, INTENT(OUT)            :: dkaid
INTEGER, INTENT(OUT)                     :: ierrok

DOUBLE PRECISION :: cozmas(0:12),phi(0:35),chin(0:26),  &
    ac(0:5,0:20),bc(0:5,0:20)
DOUBLE PRECISION :: pi,pihalf,df,fac,facd,sinth,thet,costh,  &
    fdomin,f4,z,a2,a13,a23,apihal,expapi,zds,psi,psi2,psi3,  &
    y,zmasf,eta,expam,arg,phiex,phiex2,a0ex,b0ex,chiex,d0ex,b0pex,  &
    c0ex,phis,sas,sbs,scs,sds,bspo,as,bs,bso,asp,bsp,etak,a2k,sig,  &
    f2,etaj,air,aii,dair,daii,chi,etal,z2,z21m,psi12,dzz,d1mach, under
INTEGER :: inda(0:5),indb(0:5), ifaca,j,l,k,ifun,ierro
SAVE phi,chin,ac,bc
!CCC VALUES OF THE COEFFICIENTS
DATA phi/1.d0,0.2D0,.25714285714285714286D-1,  &
    -.56507936507936507937D-2,-.39367965367965367965D-2,  &
    -.5209362066504924D-3,.3708541264731741D-3,  &
    .2123827840293627D-3,.2150629788753145D-4,  &
    -.2636904062981799D-4,-.1405469826493129D-4,  &
    -.1149328899029441D-5,.1972641193938624D-5,  &
    .1014324305593329D-5,.7034331100192200D-7,  &
    -.1525044777392676D-6,-.7677866256900572D-7,  &
    -.4669842638693018D-8,.1206673645965329D-7,  &
    .5990877668092344D-8,.3269102150077715D-9,  &
    -.97138350244428335095D-9,-.47745934295232233834D-9,  &
    -.23750318642839155779D-10,.79244598109106655567D-10,  &
    .38653584230817865528D-10,.17734105846426807873D-11,  &
    -.65332110030864751956D-11,-.31673512094772575686D-11,  &
    -.13524195177121030660D-12,.54324103217903338951D-12,  &
    .26204918647967626464D-12,.10488134973664584898D-13,  &
    -.45490420121539001435D-13,-.21851238232690420878D-13,  &
    -.82456080260379042800D-15/
DATA chin/.2D0,.1142857142857142857142857D-1,  &
    -.2438095238095238095238095D-1,-.1003471449185734900020614D-1,  &
    .8811404468547325690182833D-3,.2318339655809043564145605D-2,  &
    .7794494413564441575646057D-3,-.1373952504077228744949558D-3,  &
    -.2162301322540308393022684D-3,-.6585004634375583559042795D-4,  &
    .1502851367097217578058824D-4,.1991904617871647675455262D-4,  &
    .5731972706719818525615427D-5,-.1496427612747891044606885D-5,  &
    -.1821231830428939670133992D-5,-.5054254875930882534538217D-6,  &
    .1433283859497625931203415D-6,.1657681319078678321113361D-6,  &
    .4485592642702540575627044D-7,-.1346138242826094098161839D-7,  &
    -.1504601862773535117708677D-7,-.3995660198654805921651406D-8,  &
    .1250124931952495738882300D-8,.1363187174221864073749532D-8,  &
    .3567608459777506132029204D-9,-.1152749290139859167119863D-9,  &
    -.1233547289799408517691696D-9/
DATA cozmas/.9428090415820634D0,-.4242640687119285D0,  &
    .2904188565587606D0,-.2234261009999160D0,  &
    .1821754153944745D0,-.1539625154198624D0,  &
    .1333737583085222D0,-.1176617834148007D0,  &
    .1052687194772381D0,-.9524025714638822D-1,  &
    .8695738197073783D-1,-.8000034897653656D-1, .7407421797273086D-1/
DATA ac(0,0)/1.d0/
DATA ac(1,0)/-.44444444444444444445D-2/
DATA ac(1,1)/-.18441558441558441558D-2/
DATA ac(1,2)/.11213675213675213675D-2/
DATA ac(1,3)/.13457752124418791086D-2/
DATA ac(1,4)/.0003880626562979504D0/
DATA ac(1,5)/-.0001830686723781799D0/
DATA ac(1,6)/-.0001995460887806733D0/
DATA ac(1,7)/-.00005256191234041587D0/
DATA ac(1,8)/.00002460619652459158D0/
DATA ac(1,9)/.00002519246780924541D0/
DATA ac(1,10)/.6333157376533242D-5/
DATA ac(1,11)/-.2957485733830202D-5/
DATA ac(1,12)/-.2925255920564838D-5/
DATA ac(1,13)/-.7159702610502009D-6/
DATA ac(1,14)/.3331510720390949D-6/
DATA ac(1,15)/.3227670475692310D-6/
DATA ac(1,16)/.7767729381664199D-7/
DATA ac(1,17)/-.3600954237921120D-7/
DATA ac(1,18)/-.3441724449034226D-7/
DATA ac(1,19)/-.8188194356398772D-8/
DATA ac(1,20)/.3783148485152038D-8/
DATA ac(2,0)/.69373554135458897365D-3/
DATA ac(2,1)/.46448349036584330703D-3/
DATA ac(2,2)/-.42838130171535112460D-3/
DATA ac(2,3)/-.0007026702868771135D0/
DATA ac(2,4)/-.0002632580046778811D0/
DATA ac(2,5)/.0001663853666288703D0/
DATA ac(2,6)/.0002212087687818584D0/
DATA ac(2,7)/.00007020345615329662D0/
DATA ac(2,8)/-.00004000421782540614D0/
DATA ac(2,9)/-.00004786324966453962D0/
DATA ac(2,10)/-.00001394600741473631D0/
DATA ac(2,11)/.7536186591273727D-5/
DATA ac(2,12)/.8478502161067667D-5/
DATA ac(2,13)/.2345355228453912D-5/
DATA ac(2,14)/-.1225943294710883D-5/
DATA ac(2,15)/-.1325082343401027D-5/
DATA ac(2,16)/-.3539954776569997D-6/
DATA ac(2,17)/.1808291719376674D-6/
DATA ac(2,18)/.1900383515233655D-6/
DATA ac(3,0)/-.35421197145774384076D-3/
DATA ac(3,1)/-.31232252789031883276D-3/
DATA ac(3,2)/.3716442237502298D-3/
DATA ac(3,3)/.0007539269155977733D0/
DATA ac(3,4)/.0003408300059444739D0/
DATA ac(3,5)/-.0002634968172069594D0/
DATA ac(3,6)/-.0004089275726648432D0/
DATA ac(3,7)/-.0001501108759563460D0/
DATA ac(3,8)/.00009964015205538056D0/
DATA ac(3,9)/.0001352492955751283D0/
DATA ac(3,10)/.00004443117087272903D0/
DATA ac(3,11)/-.00002713205071914117D0/
DATA ac(3,12)/-.00003396796969771860D0/
DATA ac(3,13)/-.00001040708865273043D0/
DATA ac(3,14)/.6024639065414414D-5/
DATA ac(3,15)/.7143919607846883D-5/
DATA ac(4,0)/.378194199201773D-3/
DATA ac(4,1)/.000404943905523634D0/
DATA ac(4,2)/-.000579130526946453D0/
DATA ac(4,3)/-.00138017901171011D0/
DATA ac(4,4)/-.000722520056780091D0/
DATA ac(4,5)/.000651265924036825D0/
DATA ac(4,6)/.00114674563328389D0/
DATA ac(4,7)/.000474423189340405D0/
DATA ac(4,8)/-.000356495172735468D0/
DATA ac(4,9)/-.000538157791035111D0/
DATA ac(4,10)/-.000195687390661225D0/
DATA ac(4,11)/.000132563525210293D0/
DATA ac(4,12)/.000181949256267291D0/
DATA ac(5,0)/-.69114139728829416760D-3/
DATA ac(5,1)/-.00085995326611774383285D0/
DATA ac(5,2)/.0014202335568143511489D0/
DATA ac(5,3)/.0038535426995603052443D0/
DATA ac(5,4)/.0022752811642901374595D0/
DATA ac(5,5)/-.0023219572034556988366D0/
DATA ac(5,6)/-.0045478643394434635622D0/
DATA ac(5,7)/-.0020824431758272449919D0/
DATA ac(5,8)/.0017370443573195808719D0/
DATA bc(0,0)/.14285714285714285714D-1/
DATA bc(0,1)/.88888888888888888889D-2/
DATA bc(0,2)/.20482374768089053803D-2/
DATA bc(0,3)/-.57826617826617826618D-3/
DATA bc(0,4)/-.60412089799844901886D-3/
DATA bc(0,5)/-.0001472685745626922D0/
DATA bc(0,6)/.00005324102148009784D0/
DATA bc(0,7)/.00005206561006583416D0/
DATA bc(0,8)/.00001233115050894939D0/
DATA bc(0,9)/-.4905932728531366D-5/
DATA bc(0,10)/-.4632230987136350D-5/
DATA bc(0,11)/-.1077174523455235D-5/
DATA bc(0,12)/.4475963978932822D-6/
DATA bc(0,13)/.4152586188464624D-6/
DATA bc(0,14)/.9555819293589234D-7/
DATA bc(0,15)/-.4060599208403059D-7/
DATA bc(0,16)/-.3731367187988482D-7/
DATA bc(0,17)/-.8532670645553778D-8/
DATA bc(0,18)/.3673017245573624D-8/
DATA bc(0,19)/.3355960460784536D-8/
DATA bc(0,20)/.7643107095110475D-9/
DATA bc(1,0)/-.11848595848595848596D-2/
DATA bc(1,1)/-.13940630797773654917D-2/
DATA bc(1,2)/-.48141005586383737645D-3/
DATA bc(1,3)/.26841705366016142958D-3/
DATA bc(1,4)/.0003419706982709903D0/
DATA bc(1,5)/.0001034548234902078D0/
DATA bc(1,6)/-.00005418191982095504D0/
DATA bc(1,7)/-.00006202184830690167D0/
DATA bc(1,8)/-.00001724885886056087D0/
DATA bc(1,9)/.8744675992887053D-5/
DATA bc(1,10)/.9420684216180929D-5/
DATA bc(1,11)/.2494922112085850D-5/
DATA bc(1,12)/-.1238458608836357D-5/
DATA bc(1,13)/-.1285461713809769D-5/
DATA bc(1,14)/-.3299710862537507D-6/
DATA bc(1,15)/.1613441105788315D-6/
DATA bc(1,16)/.1633623194402374D-6/
DATA bc(1,17)/.4104252949605779D-7/
DATA bc(1,18)/-.1984317042326989D-7/
DATA bc(1,19)/-.1973948142769706D-7/
DATA bc(1,20)/-.4882194808588752D-8/
DATA bc(2,0)/.43829180944898810994D-3/
DATA bc(2,1)/.71104865116708668943D-3/
DATA bc(2,2)/.31858383945387580576D-3/
DATA bc(2,3)/-.0002404809426804458D0/
DATA bc(2,4)/-.0003722966038621536D0/
DATA bc(2,5)/-.0001352752059595618D0/
DATA bc(2,6)/.00008691694372704142D0/
DATA bc(2,7)/.0001158750753591377D0/
DATA bc(2,8)/.00003724965927846447D0/
DATA bc(2,9)/-.00002198334949606935D0/
DATA bc(2,10)/-.00002686449633870452D0/
DATA bc(2,11)/-.8023061612032524D-5/
DATA bc(2,12)/.4494756592180126D-5/
DATA bc(2,13)/.5193504763856015D-5/
DATA bc(2,14)/.1477156191529617D-5/
DATA bc(2,15)/-.7988793826096784D-6/
DATA bc(3,0)/-.37670439477105454219D-3/
DATA bc(3,1)/-.75856271658798642365D-3/
DATA bc(3,2)/-.0004103253968775048D0/
DATA bc(3,3)/.0003791263310429010D0/
DATA bc(3,4)/.0006850981673903450D0/
DATA bc(3,5)/.0002878310571932216D0/
DATA bc(3,6)/-.0002157010636115705D0/
DATA bc(3,7)/-.0003260863991373500D0/
DATA bc(3,8)/-.0001181317008748678D0/
DATA bc(3,9)/.00007887526841582582D0/
DATA bc(3,10)/.0001072081833420685D0/
DATA bc(3,11)/.00003544595251288735D0/
DATA bc(3,12)/-.00002201447920733824D0/
DATA bc(3,13)/-.00002789336359620813D0/
DATA bc(4,0)/.00058453330122076187255D0/
DATA bc(4,1)/.0013854690422372401251D0/
DATA bc(4,2)/.00086830374184946900245D0/
DATA bc(4,3)/-.00093502904801345951693D0/
DATA bc(4,4)/-.0019175486005525492059D0/
DATA bc(4,5)/-.00090795047113308137941D0/
DATA bc(4,6)/.00077050429806392235104D0/
DATA bc(4,7)/.0012953100128255983488D0/
DATA bc(4,8)/.00051933869471899550762D0/
DATA bc(4,9)/-.00038482631948759834653D0/
DATA bc(4,10)/-.00057335393099012476502D0/
DATA bc(5,0)/-.0014301070053470410656D0/
DATA bc(5,1)/-.0038637811942002539408D0/
DATA bc(5,2)/-.0027322816261168328245D0/
DATA bc(5,3)/.0033294980346743452748D0/
DATA bc(5,4)/.0075972237660887795911D0/
DATA bc(5,5)/.0039816655673062060620D0/
DATA bc(5,6)/-.0037510180460986006595D0/
DATA inda/0,20,18,15,12,8/
DATA indb/20,20,15,13,10,6/
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under=d1mach(1)*1.d+8
!CCC CONSTANTS CCCCCCCCCCCCCC
pi=ACOS(-1.d0)
pihalf=pi*0.5D0
ierrok=0
IF (ifac == 1) THEN
  IF (x >= a) THEN
    sinth=a/x
    thet=ASIN(sinth)
    costh=COS(thet)
    fdomin=x*(costh+thet*sinth)
    IF (-fdomin <= LOG(under)) ierrok=1
  ELSE
    IF ((-pi*a*0.5D0) <= LOG(under)) ierrok=1
  END IF
END IF
IF (ierrok == 0) THEN
  df=2.d0**(1.d0/3.d0)
!CCC VARIABLES CCCCCCCCCCCCCCC
  z=x/a
  a2=a*a
  a13=a**(1.d0/3.d0)
  a23=a13*a13
  IF (ifac == 1) THEN
    apihal=a*pihalf
    expapi=EXP(-apihal)
    fac=pi*expapi/a13
    facd=2.d0*pi*expapi/z/a23
  ELSE
    IF (x < a) THEN
      fac=pi/a13
      facd=2.d0*pi/z/a23
    ELSE
      sinth=a/x
      thet=ASIN(sinth)
      costh=COS(thet)
      fdomin=x*(costh+thet*sinth)
      expam=EXP(-a*pihalf+fdomin)
      fac=pi*expam/a13
      facd=2.d0*pi*expam/z/a23
    END IF
  END IF
  f4=a13**4
  IF (z <= 0.9D0) THEN
    zds=SQRT((1.d0-z)*(1.d0+z))
    psi=(1.5D0*(LOG((1.d0+zds)/z)-zds))**(2.d0/3.d0)
  ELSE IF (z > 1.1D0) THEN
    zds=SQRT((z-1.d0)*(1.d0+z))
    psi=-(1.5D0*(zds-ACOS(1.d0/z)))**(2.d0/3.d0)
  ELSE
    y=z-1.d0
    zmasf=cozmas(12)
    DO  k=0,11
      j=11-k
      zmasf=cozmas(j)+y*zmasf
    END DO
    zmasf=ABS(y)**(1.5D0)*zmasf
    IF (z < 1.d0) THEN
      psi=(1.5D0*zmasf)**(2.d0/3.d0)
    ELSE
      psi=-(1.5D0*zmasf)**(2.d0/3.d0)
    END IF
  END IF
  eta=psi/df
  arg=-a23*psi
  psi2=psi*psi
  psi3=psi2*psi
  IF ((z > 0.8D0).AND.(z < 1.2D0)) THEN
    phis=0.d0
    chi=0.d0
    sas=0.d0
    sbs=0.d0
    sds=0.d0
    scs=1.d0
    bs=0.d0
    bspo=0.d0
    DO  l=0,20
      IF (l == 0) THEN
        etal=1.d0
      ELSE
        etal=etal*eta
      END IF
      chi=chi+chin(l)*etal
    END DO
    chi=chi/df
    DO  k=0,20
      bso=bs
      bspo=bsp
      as=0.d0
      bs=0.d0
      asp=0.d0
      bsp=0.d0
      IF (k == 0) THEN
        etak=1.d0
        a2k=1.d0
        sig=1.d0
      ELSE
        etak=etak*eta
        a2k=a2k*a2
        sig=-1.d0*sig
      END IF
      phis=phis+phi(k)*etak
      f2=sig/a2k
      IF (k <= 5) THEN
        DO  j=0,20
          IF (j == 0) THEN
            etaj=1.d0
          ELSE
            etaj=etaj*eta
          END IF
          IF (j <= inda(k)) THEN
            as=as+ac(k,j)*etaj
          END IF
          IF (j <= indb(k)) THEN
            bs=bs+bc(k,j)*etaj
          END IF
          IF ((j+1) <= inda(k)) THEN
            asp=asp+(j+1)*ac(k,j+1)*etaj
          END IF
          IF ((j+1) <= indb(k)) THEN
            bsp=bsp+(j+1)*bc(k,j+1)*etaj
          END IF
        END DO
        asp=1.d0/df*asp
        bs=bs*df
        sas=sas+as*f2
        sbs=sbs+bs*f2/df
        sds=sds-(chi*as+asp+psi*bs)*f2
      END IF
      IF ((k > 0).AND.(k <= 6)) THEN
        scs=scs+(as+chi*bso+bspo)*f2
      END IF
    END DO
    phis=df*phis
    sbs=df*sbs
  ELSE
!CCC EXACT VALUES CCCCCCCCCCCCCCCCCCCCC
    z2=z*z
    z21m=1.d0-z2
    phiex=(4.d0*psi/z21m)**0.25D0
    phiex2=phiex*phiex
    a0ex=1.d0
    b0ex=-5.d0/48.d0/(psi*psi)+phiex2/48.d0/psi* (5.d0/z21m-3.d0)
    chiex=0.25D0/psi*(1.d0-z2*phiex2**3*0.25D0)
    chi=chiex
    d0ex=-(7.d0/48.d0/psi+phiex2/48.d0* (9.d0-7.d0/z21m))
    IF (psi > 0.d0) THEN
      psi12=SQRT(psi)
      dzz=SQRT(z21m)
    ELSE
      psi12=SQRT(-psi)
      dzz=SQRT(ABS(z21m))
    END IF
    b0pex=5.d0/24.d0/psi3+phiex2/48.d0*((2.d0*chi*psi-1.d0)/  &
        psi2*(5.d0/z21m-3.d0)-10.d0*z2*psi12/ z21m**2/dzz/psi)
    c0ex=1.d0
    sas=a0ex
    sbs=b0ex
    scs=c0ex
    sds=d0ex
    bs=0.d0
    bsp=0.d0
    a2k=1.d0
    sig=1.d0
    DO  k=1,6
      bso=bs
      bspo=bsp
      as=0.d0
      bs=0.d0
      asp=0.d0
      bsp=0.d0
      a2k=a2k*a2
      sig=-1.d0*sig
      f2=sig/a2k
      IF (k <= 5) THEN
        DO  j=0,20
          IF (j == 0) THEN
            etaj=1.d0
          ELSE
            etaj=etaj*eta
          END IF
          IF (j <= inda(k)) THEN
            as=as+ac(k,j)*etaj
          END IF
          IF (j <= indb(k)) THEN
            bs=bs+bc(k,j)*etaj
          END IF
          IF ((j+1) <= inda(k)) THEN
            asp=asp+(j+1)*ac(k,j+1)*etaj
          END IF
          IF ((j+1) <= indb(k)) THEN
            bsp=bsp+(j+1)*bc(k,j+1)*etaj
          END IF
        END DO
        bs=df*bs
        asp=asp/df
        sas=sas+as*f2
        sbs=sbs+bs*f2
        sds=sds-(chi*as+asp+psi*bs)*f2
      END IF
      IF (k >= 2) THEN
        scs=scs+(as+chi*bso+bspo)*f2
      ELSE
        scs=scs+(as+chi*b0ex+b0pex)*f2
      END IF
    END DO
    phis=phiex
  END IF
!CCCC CALL THE AIRY ROUTINE CCCCCCCCC
  ifun=1
  ifaca=1
  CALL aiz(ifun,ifaca,arg,0.d0,air,aii,ierro)
  ifun=2
  ifaca=1
  CALL aiz(ifun,ifaca,arg,0.d0,dair,daii,ierro)
  dkai=fac*phis*(air*sas+dair*sbs/f4)
  dkaid=facd/phis*(dair*scs+air*sds/a23)
ELSE
  dkai=0.d0
  dkaid=0.d0
END IF
RETURN
END SUBROUTINE aiexki
