  MODULE Conical
  !------------------------------------------------
  ! Authors:
  !  T.Mark Dunster (SDSU, San Diego, USA)
  !                 e-mail: mdunster@mail.sdsu.edu   
  !  Amparo Gil     (U. Cantabria, Santander, Spain)
  !                  e-mail: amparo.gil@unican.es
  !  Javier Segura  (U. Cantabria, Santander, Spain)
  !                  e-mail: javier.segura@unican.es
  !  Nico M. Temme  (CWI, The Netherlands)
  !                  e-mail: nico.temme@cwi.nl
  ! -------------------------------------------------------------
  !  References:
  !     1. T.M. Dunster, A. Gil, J. Segura, N.M. Temme.
  !           Accompanying paper in Computer Physics Commun.  
  !     2. A. Gil, J. Segura, N.M. Temme. Comput Phys Commun 
  !                                       183 (2012) 794-799
  ! -------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, PARAMETER  :: r8 = KIND(0.0d0)
  PRIVATE
  PUBLIC  :: conicpr, conicp, conicr
CONTAINS
    SUBROUTINE conicpr(x,mu,tau,pm,pmd,rm,rmd,ierr)
    ! ------------------------------------------------------------
    ! Calculation of the conical functions P^(mu)_(-1/2+i*tau)(x),
    ! R^(mu)_(-1/2+i*tau)(x) and their first order derivatives for
    ! 1<x, 0<tau, 0 <=mu.
    !   
    !  In order to avoid, overflow/underflow problems in IEEE double
    !  precision arithmetic, the admissible parameter ranges 
    !  for computation are:
    !        1<x<=100, 0<tau<=100, 0<=mu<=100
    ! ------------------------------------------------------------
    ! Inputs:
    !   x ,    argument of the functions (real value, x>1).
    !   mu,    function parameter (integer value, m=0,1,2,...).
    !   tau,   function parameter (real value, tau>0).
    !
    ! Outputs:
    !   pm ,  function P^(mu)_(-1/2+i*tau)(x)
    !   pmd , function d(P^(mu)_(-1/2+i*tau)(x))/dx
    !   rm ,  function R^(mu)_(-1/2+i*tau)(x)
    !   rmd , function d(R^(mu)_(-1/2+i*tau)(x))/dx 
    !   ierr , error flag
    !          ierr=0, computation succesful. 
    !          ierr=1, an error occurred during computation.
    !                  The function values and the derivatives
    !                  are set to 0. 
    ! -----------------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: tau
    REAL(r8), INTENT(OUT) :: pm
    REAL(r8), INTENT(OUT) :: pmd
    REAL(r8), INTENT(OUT) :: rm
    REAL(r8), INTENT(OUT) :: rmd
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: p0, p1, r0, r1
    INTEGER :: ierr1, ierr2, ierr3, ierr4
    ierr=0
    CALL conicp(x,mu,tau,p0,ierr1)  
    CALL conicp(x,mu+1,tau,p1,ierr2)  
    CALL conicr(x,mu,tau,r0,ierr3)
    CALL conicr(x,mu+1,tau,r1,ierr4)
    pm=p0
    rm=r0
    IF (((ierr1==0).AND.(ierr2==0)).AND.&
         ((ierr3==0).AND.(ierr4==0))) THEN
       !Computation of the derivatives
      pmd=-p1/sqrt(x*x-1.0_r8)+mu*x/(x*x-1.0_r8)*p0
      rmd=-r1/sqrt(x*x-1.0_r8)+mu*x/(x*x-1.0_r8)*r0
    ELSE
      ierr=1
      pm=0.0_r8
      pmd=0.0_r8
      rm=0.0_r8
      rmd=0.0_r8
    ENDIF  
    END SUBROUTINE conicpr
  
    SUBROUTINE conicp(x,mu,tau,pmtau,ierr)
    ! ------------------------------------------------------------
    ! Calculation of the conical functions P^(mu)_(-1/2+i*tau)(x),
    !  -1<x, 0<tau, 0 <=mu.
    !   
    !  In order to avoid, overflow/underflow problems in IEEE double
    !  precision arithmetic, the admissible parameter ranges 
    !  for computation are:
    !       -1<x<1,   0<tau<=100, 0<=mu<=40
    !        1<x<=100, 0<tau<=100, 0<=mu<=100
    ! ------------------------------------------------------------
    ! Inputs:
    !   x ,    argument of the functions (real value greater than -1)
    !   mu,    function parameter (integer positive value)
    !   tau,   function parameter (real positive value)
    !
    ! Outputs:
    !   pmtau
    !   ierr , error flag
    !          ierr=0, computation succesful. 
    !          ierr=1, overflow/underflow problems. 
    !          ierr=2, any of the arguments (x, mu, tau) is out of 
    !                  range. 
    ! -----------------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: tau
    REAL(r8), INTENT(OUT) :: pmtau   
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: pmu1, pmu2, beta, q, argu, xargu, xfac, tau2, pmm, pm1, pm2
    REAL(r8) :: x2, x3, x4, fxt
    INTEGER :: mm, mm2, mu2, itau, im
    ierr=0
    itau=0
    IF ((x<-1.0_r8).OR.(mu<0).OR.(tau<0).OR.&
        (mu>100).OR.(tau>100).OR.(x>100)) THEN
      ierr=2
      pmtau=0
    ELSEIF ((abs(x)<1.0_r8).AND.(mu>40)) THEN
      ierr=2
      pmtau=0
    ELSE
      IF (x<=0) THEN
        IF (mu<=60) THEN        
          CALL pmuint(x,mu,tau,pmm,ierr)
          pmtau=pmm
        ELSE 
          CALL expaelem(x,mu,tau,pmu1,ierr)
          pmtau=pmu1
        ENDIF
      ELSE
        IF (x<1.0_r8) THEN
          IF (mu<=80) THEN
            mm=90
            CALL expaelem(x,mm,tau,pmu1,ierr)
            mm=mm-1
            CALL expaelem(x,mm,tau,pmu2,ierr)
            xargu=1.0_r8-x*x
            xfac=x/sqrt(xargu)
            tau2=tau*tau
            DO im=mm,mu+1,-1
              pmm=(2.0_r8*im*xfac*pmu2+pmu1)/(tau2+0.25_r8+im*(im-1.0_r8))
              pmu1=pmu2
              pmu2=pmm
            ENDDO
            pmtau=pmm
          ELSE
            CALL expaelem(x,mu,tau,pmu1,ierr)
            pmtau=pmu1
          ENDIF
        ELSE
          IF (abs(x-1.0_r8)<dwarf) THEN
            pmtau=0.0_r8
          ELSE     
            mm=mu
            mm2=mm*mm
            IF (mu==0) THEN
              beta=1.0_r8  
            ELSE   
              beta=sqrt(1.0_r8+tau*tau/mm2)/(tau/mm);
            ENDIF
            IF (x>19.7_r8) THEN             
              x2=x*x
              x3=x2*x
              x4=x3*x
              fxt=2.266e-07_r8*x4-7.042e-05_r8*x3 +0.008317_r8*x2 -0.4637_r8*x &
                  +12.15_r8
              IF (tau<fxt+0.2_r8) THEN
                itau=1
              ENDIF
           ENDIF
            IF (x<beta) THEN
              ! Monotonic Region
              IF (itau==1) THEN
                IF ((x>84.0_r8).AND.(tau<1.7_r8)) THEN
                  mm=140
                ELSE   
                  mm=160
                ENDIF   
              ELSE                
                IF (x<1.01_r8) THEN
                  mm=160
                ELSE 
                  mm=120
                ENDIF
              ENDIF   
              mm2=mm*mm
              CALL expakia(x,mm,tau,pmu1,ierr)
              mm=mm-1
              CALL expakia(x,mm,tau,pmu2,ierr)
              xargu=x*x-1.0_r8
              xfac=x/sqrt(xargu)
              tau2=tau*tau
              DO im=mm,mu+1,-1
                pmm=(+2.0_r8*im*xfac*pmu2-pmu1)/(tau2+0.25_r8+im*(im-1.0_r8))
                pmu1=pmu2
                pmu2=pmm           
              ENDDO          
              pmtau=pmm
            ELSE
              ! Oscillatory region
              IF (mu==0) THEN
                q=1.0_r8
              ELSE   
                q=x/sqrt((tau/mu)*(tau/mu)*(x*x-1.0_r8)-1.0_r8)
              ENDIF   
              IF ((q<0.25_r8).AND.(mu>80)) THEN
                CALL expaosc(x,mu,tau,pmtau,ierr)
              ELSE
                IF (((tau<34.0_r8).AND.((x>=1.07).AND.(x<1.6_r8))).OR.&
                     (x<1.04_r8).OR.((tau<45.0_r8).AND.((x>=1.04_r8).AND.&
                     (x<1.07_r8))).OR.((tau<15.0_r8).AND.((x>=1.6_r8).AND.&
                     (x<5.0_r8))).OR.((tau<10.0_r8).AND.((x>=5.0_r8).AND.&
                     (x<8.0_r8))).OR.((tau<5.0_r8).AND.(x>=8.0_r8))) THEN
                  IF ((mu>95).AND.(x<20.0_r8)) THEN
                    mm=mu
                    mm2=mm*mm
                    argu=(tau/mm)*(tau/mm)*(x*x-1.0_r8)-1.0_r8
                    IF (argu>0) THEN
                      q=x/sqrt(argu)
                      IF (q<0.25_r8) THEN
                        CALL expaosc(x,mm,tau,pmu1,ierr)
                      ELSE
                        CALL expakia(x,mm,tau,pmu1,ierr)
                      ENDIF
                    ELSE
                      CALL expakia(x,mm,tau,pmu1,ierr)
                    ENDIF
                    pmtau=pmu1    
                  ELSE
                    IF (x<1.6_r8) THEN
                      IF (mu<5) THEN
                        IF (x<1.01_r8) THEN
                          mm=160
                        ELSE   
                         mm=140
                        ENDIF       
                      ENDIF   
                      IF (x<1.01_r8) THEN
                        mm=160
                      ELSE
                        mm=120  
                      ENDIF
                    ELSEIF (x<2.45_r8) THEN    
                      mm=95
                    ELSE
                      IF (itau==1) THEN
                        IF ((x>84.0_r8).AND.(tau<1.7_r8)) THEN
                          IF ((x>99).AND.(tau<1.4_r8)) THEN                              
                            mm=153  
                          ELSE   
                            mm=140
                          ENDIF
                        ELSE   
                          mm=160
                        ENDIF    
                      ELSE
                        mm=100
                      ENDIF   
                    ENDIF
                    mm2=mm*mm
                    argu=(tau/mm)*(tau/mm)*(x*x-1.0_r8)-1.0_r8
                    IF (argu>0) THEN
                       q=x/sqrt(argu)
                      IF (q<0.25_r8) THEN
                        CALL expaosc(x,mm,tau,pmu1,ierr)
                      ELSE
                        CALL expakia(x,mm,tau,pmu1,ierr)
                      ENDIF
                    ELSE
                      CALL expakia(x,mm,tau,pmu1,ierr)
                    ENDIF
                    mm=mm-1
                    argu=(tau/mm)*(tau/mm)*(x*x-1.0_r8)-1.0_r8
                    IF (argu>0) THEN
                      q=x/sqrt(argu)
                      IF (q<0.25_r8) THEN
                        CALL expaosc(x,mm,tau,pmu2,ierr)
                      ELSE
                        CALL expakia(x,mm,tau,pmu2,ierr)
                      ENDIF
                    ELSE
                      CALL expakia(x,mm,tau,pmu2,ierr)
                    ENDIF
                    xargu=x*x-1.0_r8
                    xfac=x/sqrt(xargu)
                    tau2=tau*tau
                    DO im=mm,mu+1,-1
                      pmm=(+2.0_r8*im*xfac*pmu2-pmu1)/(tau2+0.25_r8+im*(im-1.0_r8))
                      pmu1=pmu2
                      pmu2=pmm           
                    ENDDO          
                    pmtau=pmm
                  ENDIF  
                ELSE
                  IF (mu>1) THEN
                  ! Forward recurrence
                    mm=0   
                    CALL tauex1(x,mm,tau,pmu1,ierr)
                    mm=1
                    CALL tauex1(x,mm,tau,pmu2,ierr)
                    xargu=x*x-1.0_r8
                    xfac=x/sqrt(xargu)
                    tau2=tau*tau
                    DO im=1,mu-1
                      pmm=2.0_r8*im*xfac*pmu2-pmu1*(tau2+0.25_r8+im*(im-1.0_r8))
                      pmu1=pmu2
                      pmu2=pmm           
                    ENDDO          
                    pmtau=pmm
                  ELSE
                    CALL tauex1(x,mu,tau,pmtau,ierr)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF        
      ENDIF
    ENDIF
    END SUBROUTINE conicp
  
    SUBROUTINE conicr(x,mu,tau,rm,ierr)
    ! ------------------------------------------------------------
    ! Calculation of the conical functions R^(mu)_(-1/2+i*tau)(x),
    !  1<x, 0<tau, 0 <=mu. This function is a real valued
    !  numerically satisfactory companion of the function
    !  P^(mu)_(-1/2+i*tau)(x) for x>1.   
    !   
    !  In order to avoid, overflow/underflow problems in IEEE double
    !  precision arithmetic, the admissible parameter ranges 
    !  for computation are:
    !        1<x<=100, 0<tau<=100, 0<=mu<=100
    ! ------------------------------------------------------------
    ! Inputs:
    !   x ,    argument of the functions (real value greater than -1)
    !   mu,    function parameter (integer positive value)
    !   tau,   function parameter (real positive value)
    !
    ! Outputs:
    !   rm,    function R^(mu)_(-1/2+i*tau)(x)
    !   rmd,   first order derivative of
    !                 R^(mu)_(-1/2+i*tau)(x)
    !   ierr , error flag
    !          ierr=0, computation succesful. 
    !          ierr=1, overflow/underflow problems. 
    !          ierr=2, any of the arguments (x, mu, tau) is out of 
    !                  range. The function values are set to 0. 
    ! -----------------------------------------------------------
 
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: tau
    REAL(r8), INTENT(OUT) :: rm
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: rmu1, rmu2, xargu, xfac, tau2, rmm, rmd
    INTEGER :: itau, im, mm
    ierr=0
    rmd=0.0_r8
    IF ((x<1.001_r8).OR.(mu<0).OR.(tau<0).OR.&
        (mu>100).OR.(tau>100).OR.(x>100)) THEN
      ierr=2
      rm=0.0_r8
      rmd=0.0_r8
    ELSE
      IF ((x<1.5_r8).OR.((x<1.65_r8).AND.(tau>40))) THEN
         IF (x>1.15_r8) THEN
           itau=0 
           IF (x<1.38) THEN
             IF (tau>20.0_r8-30*(x-1.226_r8)) THEN
                itau=1
             ENDIF   
           ELSEIF (x<1.5) THEN
             IF (tau>16.0_r8) THEN
               itau=1
             ENDIF 
           ELSE    
             IF (tau>20.0_r8-30*(x-1.226_r8)) THEN
               itau=1
             ENDIF 
           ENDIF   
           IF (itau==1) THEN
             IF (mu==0) THEN
               CALL conirkum(0, tau, x, rm, ierr)
             ELSEIF (mu==1) THEN
               CALL conirkum(1, tau, x, rm, ierr)
             ELSE
               ! Kummer expansions to compute the starting values  
               CALL conirkum(0, tau, x, rmu1, ierr)
               CALL conirkum(1, tau, x, rmu2, ierr)
               xargu=x*x-1.0_r8
               xfac=x/sqrt(xargu)
               tau2=tau*tau
               DO im=1,mu-1
                 rmm=2.0_r8*im*xfac*rmu2-rmu1*(tau2+0.25_r8+im*(im-1.0_r8))
                 rmu1=rmu2
                 rmu2=rmm           
               ENDDO
               rm=rmm
             ENDIF  
           ELSE
             CALL conirsmall(tau, x, rmu1, rmu2)
             IF (mu==0) THEN
               rm=rmu1
             ELSEIF (mu==1) THEN
               rm=rmu2
             ELSE   
               xargu=x*x-1.0_r8
               xfac=x/sqrt(xargu)
               tau2=tau*tau
               DO im=1,mu-1
                 rmm=2.0_r8*im*xfac*rmu2-rmu1*(tau2+0.25_r8+im*(im-1.0_r8))
                 rmu1=rmu2
                 rmu2=rmm           
               ENDDO
               rm=rmm
             ENDIF  
           ENDIF
        ELSE
          IF (tau<23.0_r8) THEN
            CALL conirsmall(tau, x, rmu1, rmu2)
            IF (mu==0) THEN
              rm=rmu1
            ELSEIF (mu==1) THEN
              rm=rmu2
            ELSE   
              xargu=x*x-1.0_r8
              xfac=x/sqrt(xargu)
              tau2=tau*tau
              DO im=1,mu-1
                rmm=2.0_r8*im*xfac*rmu2-rmu1*(tau2+0.25_r8+im*(im-1.0_r8))
                rmu1=rmu2
                rmu2=rmm           
              ENDDO
              rm=rmm
            ENDIF  
          ELSE
            IF (mu<2) THEN
            ! Direct computation  
              CALL conirkum(mu, tau, x, rm, ierr)
            ELSE   
              IF ((mu<5).AND.(tau>50.0_r8)) THEN
               ! Direct computation  
                CALL conirkum(mu, tau, x, rm, ierr)
              ELSE
               ! Use of the recurrence relation
                mm=0   
                CALL conirkum(mm, tau, x, rmu1, ierr)
                mm=1
                CALL conirkum(mm, tau, x, rmu2, ierr)
                xargu=x*x-1.0_r8
                xfac=x/sqrt(xargu)
                tau2=tau*tau
                DO im=1,mu-1
                  rmm=2.0_r8*im*xfac*rmu2-rmu1*(tau2+0.25_r8+im*(im-1.0_r8))
                  rmu1=rmu2
                  rmu2=rmm           
                ENDDO
                rm=rmm
              ENDIF
            ENDIF 
          ENDIF
        ENDIF
     ELSE
       CALL conirlar(mu, tau, x, rm, ierr)
      ENDIF
    ENDIF
    END SUBROUTINE conicr

    SUBROUTINE conirsmall(tau, x, r0, r1)
    !-------------------------------------------------  
    !Computation of r0, r1 using the power series
    !expansion valid for small x-1 
    !--------------------------------------------------
    USE gammaCHI
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(IN) :: tau
    REAL(r8), INTENT(OUT) :: r0
    REAL(r8), INTENT(OUT) :: r1
    REAL(r8) :: w, w2, sw, xxi, tau2, psizt, psilnw, s, pocha, sn, sn2, toln, &
                logw, z, zz, psik, gammk1, psifac, psifac2, s2, ss, ss2    
    INTEGER :: k
    tau2=tau*tau
    toln=1e-20_r8;
    w2=(x-1.0_r8)/(x+1.0_r8)
    w=sqrt(w2)
    xxi=1.0_r8/(x*x-1.0_r8)
    sw=0.5_r8*(-1.0_r8+w2)
    logw=log(w);
    z=(1.0_r8-x)*0.5_r8;
    psizt=psizr(0.5_r8,tau);
    psilnw=psizt+logw;
    s=0.0_r8
    s2=0.0_r8
    k=0;
    psik=-0.577215664901532860606512090082_r8
    sn=1.0_r8;
    sn2=1.0_r8
    zz=1.0_r8
    pocha=1.0_r8
    DO WHILE ((sn>toln).AND.(sn2>toln))  
      psifac=psik-psilnw;
      psifac2=k*psifac+sw
      gammk1=gamma(1.0_r8+k)
      ss=(pocha/gammk1)*zz*psifac/gammk1;
      ss2=(pocha/gammk1)*zz*psifac2/gammk1;
      s=s+ss;
      s2=s2+ss2;
      sn=abs(ss/s);
      sn2=abs(ss2/s2);
      psik=psik+1.0_r8/(k+1.0_r8);
      pocha=pocha*((0.5_r8+k)*(0.5_r8+k)+tau2)
      k=k+1;
      zz=zz*z
    ENDDO
    r0=s
    r1=-s2/w
    END SUBROUTINE conirsmall   

    SUBROUTINE conirkum(mu, tau, x, ras, ierr)
    !--------------------------------------------------------
    ! Computation of rmu using the expansion in terms of the 
    ! Kummer functions with real computations  
    !--------------------------------------------------------
    USE Someconstants
    USE BesselJY  
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: tau
    REAL(r8), INTENT(OUT) :: ras   
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: fk(0:7)  
    REAL(r8) :: a, b, d, alpha, eta, q, nu, z, w, f0, h
    REAL(r8) :: z2, z3, z4, z5, z6, z7, z8, d2, d3, d4, d5, d6, &
                d7, d8, b2, b3, b4, b5, b6, b7, argu, mum, nu1, nu2
    REAL(r8) :: fac2,beta,cosb,sinb,rephi0,imphi0,&
                fac3,j2,y2,rephi1,imphi1,s,g,sing,cosg,&
                rphik(0:7),iphik(0:7),taui,atau,coef1,coef2       
    REAL(r8) :: jax1, jax2, jap1, jap2, yax1, yax2, yap1, yap2
    INTEGER :: k
    ierr=0
    mum=mu
    w=sqrt((x-1.0_r8)*(x+1.0_r8));
    z=0.5_r8/(w*(w+x));
    z2=z*z
    z3=z2*z
    z4=z3*z
    z5=z4*z
    z6=z5*z
    z7=z6*z
    z8=z7*z
    alpha=log(1.0_r8+1.0_r8/z);
    eta=z/(1.0_r8+z); 
    f0=1.0_r8/((z+1.0_r8)*alpha);
    a=alpha;
    d=z*alpha;
    d2=d*d
    d3=d2*d
    d4=d3*d
    d5=d4*d
    d6=d5*d
    d7=d6*d
    d8=d7*d
    b=-mu-0.5_r8;
    b2=b*b
    b3=b2*b
    b4=b3*b
    b5=b4*b
    b6=b5*b
    b7=b6*b
    f0=f0**b;
    fk(0)=1.0_r8;
    fk(1)=0.5_r8*(2.0_r8*d*z+d-2.0_r8*z)*b/d;
    fk(2)=(1.0_r8/24.0_r8)*(12.0_r8*z2+12.0_r8*b*z2+d2-12.0_r8*d2*z&
          -12.0_r8*d2*z2-24.0_r8*b*d*z2+12.0_r8*b*d2*z+12.0_r8*b*d2*z2&
           +3.0_r8*b*d2-12.0_r8*b*d*z)*b/d2;
    fk(3)=(1.0_r8/48.0_r8)*(-16.0_r8*z3-24.0_r8*b*z3+8.0_r8*d3*z+&
          24.0_r8*d3*z2+24.0_r8*b*d*z3+16.0_r8*d3*z3+24.0_r8*b2*d*z3&
          -24.0_r8*b2*d2*z2-24.0_r8*b2*d2*z3-10.0_r8*b*d3*z-36.0_r8*b*d3*z2&
          -24.0_r8*b*d3*z3+6.0_r8*b2*d3*z+12.0_r8*b2*d3*z2+8.0_r8*b2*d3*z3&
          +24.0_r8*b*d2*z2+24.0_r8*b*d2*z3+b2*d3-6.0_r8*b2*d2*z+12.0_r8*b2*d*z2&
          +b*d3-2.0_r8*b*d2*z+12.0_r8*b*d*z2-8.0_r8*b2*z3)*b/d3;
    fk(4)= (1.0_r8/5760.0_r8)*(-1440.0_r8*b*d2*z4-2880.0_r8*b2*d*z4&
          -1920.0_r8*b*d3*z4+2880.0_r8*b2*d3*z4-2880.0_r8*d4*z3&
          -240.0_r8*d4*z-1680.0_r8*d4*z2+2640.0_r8*b*d4*z4+240.0_r8*b3*d4*z4&
          -960.0_r8*b3*d*z4+1440.0_r8*b3*d2*z4-960.0_r8*b3*d3*z4-1440.0_r8*b2*d4*z4&
          -1920.0_r8*b*d*z4-2.0_r8*d4-960.0_r8*b*d*z3-1440.0_r8*b2*d*z3&
          -2880.0_r8*b*d3*z3+4320.0_r8*b2*d3*z3-1440.0_r8*b*d2*z3&
          +120.0_r8*b*d2*z2-1440.0_r8*d4*z4-2880.0_r8*b2*d4*z3+5.0_r8*b*d4&
          +15.0_r8*b3*d4+30.0_r8*b2*d4-120.0_r8*b3*d3*z+360.0_r8*b3*d2*z2&
          -480.0_r8*b3*d*z3+1440.0_r8*b3*d2*z3-1440.0_r8*b3*d3*z3+360.0_r8*b*d4*z&
          +3000.0_r8*b*d4*z2+5280.0_r8*b*d4*z3-720.0_r8*b3*d3*z2+120.0_r8*b3*d4*z&
          +360.0_r8*b3*d4*z2+480.0_r8*b3*d4*z3-240.0_r8*b2*d4*z-1680.0_r8*b2*d4*z2&
          +480.0_r8*b2*d2*z2-960.0_r8*b*d3*z2-120.0_r8*b2*d3*z+1200.0_r8*b2*d3*z2&
          +240.0_r8*b3*z4+2640.0_r8*b*z4+1440.0_r8*b2*z4+1440.0_r8*z4)*b/d4;
    fk(5)=(1.0_r8/11520.0_r8)*(2880.0_r8*b*d4*z5+2880.0_r8*b3*d4*z5&
          +2880.0_r8*b3*d*z5-1920.0_r8*b3*d2*z5-1920.0_r8*b3*d3*z5&
          -5280.0_r8*b2*d4*z5+2880.0_r8*b*d*z5+1920.0_r8*b*d2*z5&
          +5280.0_r8*b2*d*z5+960.0_r8*b2*d2*z5+1920.0_r8*b*d3*z5&
          -960.0_r8*b2*d3*z5+5760.0_r8*d5*z4+1440.0_r8*d5*z2+96.0_r8*d5*z&
          +4800.0_r8*d5*z3-4800.0_r8*b*d5*z5+96.0_r8*b4*d5*z5-480.0_r8*b4*d4*z5&
          -960.0_r8*b3*d5*z5+3360.0_r8*b2*d5*z5+480.0_r8*b4*d*z5-960.0_r8*b4*d2*z5&
          +960.0_r8*b4*d3*z5+1920.0_r8*b*d2*z4+2640.0_r8*b2*d*z4+2880.0_r8*b*d3*z4&
          -1440.0_r8*b2*d3*z4+5760.0_r8*b*d4*z4+5760.0_r8*b3*d4*z4+1440.0_r8*b3*d*z4&
          -1920.0_r8*b3*d2*z4-2880.0_r8*b3*d3*z4-10560.0_r8*b2*d4*z4+1440.0_r8*b*d*z4&
          -720.0_r8*b2*d2*z3+960.0_r8*b*d3*z3-240.0_r8*b2*d3*z3-160.0_r8*b*d2*z3&
          +2304.0_r8*d5*z5+960.0_r8*b2*d2*z4-30.0_r8*b4*d4*z+120.0_r8*b4*d3*z2&
          -240.0_r8*b4*d2*z3+240.0_r8*b4*d*z4+130.0_r8*b2*d5*z+2040.0_r8*b2*d5*z2&
          +6960.0_r8*b2*d5*z3+8400.0_r8*b2*d5*z4-164.0_r8*b*d5*z-2880.0_r8*b*d5*z2&
          -9920.0_r8*b*d5*z3-12000.0_r8*b*d5*z4+30.0_r8*b4*d5*z+120.0_r8*b4*d5*z2&
          +240.0_r8*b4*d5*z3+240.0_r8*b4*d5*z4-60.0_r8*b3*d5*z-720.0_r8*b3*d5*z2&
          -2080.0_r8*b3*d5*z3-2400.0_r8*b3*d5*z4+720.0_r8*b4*d3*z3&
          -240.0_r8*b4*d4*z2-720.0_r8*b4*d4*z3-960.0_r8*b4*d4*z4&
          -960.0_r8*b4*d2*z4+1440.0_r8*b4*d3*z4-6000.0_r8*b2*d4*z3-2.0_r8*b*d5&
          +3.0_r8*b4*d5+10.0_r8*b3*d5+5.0_r8*b2*d5-800.0_r8*b3*d2*z3&
          -480.0_r8*b3*d3*z3+4.0_r8*b*d4*z+480.0_r8*b*d4*z2+3360.0_r8*b*d4*z3&
          +240.0_r8*b3*d3*z2-60.0_r8*b3*d4*z+480.0_r8*b3*d4*z2+3360.0_r8*b3*d4*z3&
          -10.0_r8*b2*d4*z-720.0_r8*b2*d4*z2+120.0_r8*b2*d3*z2-96.0_r8*b4*z5&
          -960.0_r8*b3*z5-3360.0_r8*b2*z5-4800.0_r8*b*z5-2304.0_r8*z5)*b/d5;
    fk(6)=(1.0_r8/2903040.0_r8)*(302400.0_r8*b2*d4*z6-580608.0_r8*b*d*z6&
          -362880.0_r8*b*d2*z6-1209600.0_r8*b2*d*z6-302400.0_r8*b2*d2*z6&
          -322560.0_r8*b*d3*z6+60480.0_r8*b5*d4*z6+4032.0_r8*b5*d6*z6&
          -60480.0_r8*b4*d6*z6-907200.0_r8*b2*d6*z6-24192.0_r8*b5*d*z6&
          +342720.0_r8*b3*d6*z6-580608.0_r8*b*d5*z6+241920.0_r8*b4*d5*z6&
          -302400.0_r8*b4*d4*z6-846720.0_r8*b3*d5*z6+1209600.0_r8*b2*d5*z6&
          -241920.0_r8*b4*d*z6+302400.0_r8*b4*d2*z6-362880.0_r8*b*d4*z6&
          +302400.0_r8*b3*d4*z6-846720.0_r8*b3*d*z6+302400.0_r8*b3*d2*z6&
          +403200.0_r8*b3*d3*z6-4032.0_r8*d6*z-1572480.0_r8*d6*z4-1451520.0_r8*d6*z5&
          -124992.0_r8*d6*z2-725760.0_r8*d6*z3+1104768.0_r8*b*d6*z6-80640.0_r8*b5*d3*z6&
          +60480.0_r8*b5*d2*z6-24192.0_r8*b5*d5*z6+16.0_r8*d6-725760.0_r8*b*d4*z5&
          +604800.0_r8*b3*d4*z5-423360.0_r8*b3*d*z5+302400.0_r8*b3*d2*z5&
          +604800.0_r8*b3*d3*z5+604800.0_r8*b2*d4*z5-290304.0_r8*b*d*z5&
          -362880.0_r8*b*d2*z5-604800.0_r8*b2*d*z5-302400.0_r8*b2*d2*z5&
          -483840.0_r8*b*d3*z5-1451520.0_r8*b*d5*z5+604800.0_r8*b4*d5*z5&
          -604800.0_r8*b4*d4*z5-2116800.0_r8*b3*d5*z5+3024000.0_r8*b2*d5*z5&
          -120960.0_r8*b4*d*z5+302400.0_r8*b4*d2*z5+30240.0_r8*b*d2*z4&
          -161280.0_r8*b*d3*z4-40320.0_r8*b2*d3*z4-423360.0_r8*b*d4*z4&
          +332640.0_r8*b3*d4*z4+196560.0_r8*b3*d2*z4+100800.0_r8*b3*d3*z4&
          +332640.0_r8*b2*d4*z4-20160.0_r8*b2*d3*z3-483840.0_r8*d6*z6&
          +3578400.0_r8*b*d6*z4+3314304.0_r8*b*d6*z5+271656.0_r8*b*d6*z2&
          +1632960.0_r8*b*d6*z3+756.0_r8*b5*d6*z+15120.0_r8*b5*d6*z4+&
          12096.0_r8*b5*d6*z5+3780.0_r8*b5*d6*z2+10080.0_r8*b5*d6*z3&
          -1260.0_r8*b4*d6*z-206640.0_r8*b4*d6*z4-181440.0_r8*b4*d6*z5&
          -26460.0_r8*b4*d6*z2-110880.0_r8*b4*d6*z3+3780.0_r8*b3*d6*z&
          -756.0_r8*b5*d5*z+3780.0_r8*b5*d4*z2-10080.0_r8*b5*d3*z3&
          +15120.0_r8*b5*d2*z4-12096.0_r8*b5*d*z5+1123920.0_r8*b3*d6*z4&
          +1028160.0_r8*b3*d6*z5+99540.0_r8*b3*d6*z2+534240.0_r8*b3*d6*z3&
          -60480.0_r8*b5*d5*z4-60480.0_r8*b5*d3*z4+30240.0_r8*b5*d4*z3&
          +90720.0_r8*b5*d4*z4+120960.0_r8*b5*d4*z5+60480.0_r8*b5*d2*z5&
          -120960.0_r8*b5*d3*z5-7560.0_r8*b5*d5*z2-30240.0_r8*b5*d5*z3&
          -60480.0_r8*b5*d5*z5-6804.0_r8*b2*d6*z-2938320.0_r8*b2*d6*z4&
          -2721600.0_r8*b2*d6*z5-223524.0_r8*b2*d6*z2-1340640.0_r8*b2*d6*z3&
          +7560.0_r8*b*d6*z-42.0_r8*b*d6+146160.0_r8*b2*d2*z4+63.0_r8*b5*d6&
          +315.0_r8*b4*d6-91.0_r8*b2*d6+315.0_r8*b3*d6+504.0_r8*b2*d5*z&
          +41328.0_r8*b2*d5*z2+725760.0_r8*b2*d5*z3+2499840.0_r8*b2*d5*z4&
          -24192.0_r8*b*d5*z2-362880.0_r8*b*d5*z3-1209600.0_r8*b*d5*z4&
          -2520.0_r8*b4*d5*z+15120.0_r8*b4*d5*z2+181440.0_r8*b4*d5*z3&
          +524160.0_r8*b4*d5*z4-1260.0_r8*b3*d5*z-32760.0_r8*b3*d5*z2&
          -514080.0_r8*b3*d5*z3-1753920.0_r8*b3*d5*z4-40320.0_r8*b4*d3*z3&
          +11340.0_r8*b4*d4*z2-30240.0_r8*b4*d4*z3-332640.0_r8*b4*d4*z4&
          +95760.0_r8*b4*d2*z4-80640.0_r8*b4*d3*z4+30240.0_r8*b2*d4*z3&
          -50400.0_r8*b3*d3*z3-504.0_r8*b*d4*z2-60480.0_r8*b*d4*z3&
          +8820.0_r8*b3*d4*z2+30240.0_r8*b3*d4*z3+756.0_r8*b2*d4*z2&
          +60480.0_r8*b4*z6+907200.0_r8*b2*z6+342720.0_r8*b3*z6&
          +1104768.0_r8*b*z6+4032.0_r8*b5*z6+483840.0_r8*z6)*b/d6;
    fk(7)=0.717638521457965902410346854791e-9_r8*(2247840.0_r8*b3*d4*z4&
         -120960.0_r8*b*d4*z4+80640.0_r8*b2*d4*z4-4838400.0_r8*b2*d4*z5&
         +29030400.0_r8*b2*d5*z5-58060800.0_r8*b*d5*z5-17902080.0_r8*b3*d3*z5&
         +2741760.0_r8*b3*d5*z4+4515840.0_r8*b4*d4*z4-80640.0_r8*b3*d5*z3&
         +62899200.0_r8*b3*d5*z5+10886400.0_r8*b3*d4*z5-14515200.0_r8*b*d4*z5&
         -20563200.0_r8*b4*d3*z5-5806080.0_r8*b2*d3*z5-36288000.0_r8*b4*d5*z5&
         -3870720.0_r8*b*d5*z4+80640.0_r8*b2*d5*z3-2136960.0_r8*b4*d5*z4&
         +806400.0_r8*b2*d5*z4-665280.0_r8*b4*d5*z3-193536000.0_r8*b*d5*z6&
         -137088000.0_r8*b4*d5*z6+222566400.0_r8*b3*d5*z6+109670400.0_r8*b2*d5*z6&
         +9676800.0_r8*d2*b*z6-46448640.0_r8*d3*b*z6+51125760.0_r8*d2*b2*z6&
         +129427200.0_r8*d4*b3*z6-4838400.0_r8*d4*b2*z6-101606400.0_r8*d4*b*z6&
         +18385920.0_r8*d3*b3*z6+84430080.0_r8*d2*b3*z6+61286400.0_r8*d2*b4*z6&
         +220147200.0_r8*d6*b5*z6-435456000.0_r8*d6*b4*z6+153619200.0_r8*d6*b3*z6&
         -24192000.0_r8*d5*b5*z6+481420800.0_r8*d6*b2*z6-377395200.0_r8*d6*b*z6&
         -31449600.0_r8*d4*b5*z6-28546560.0_r8*d3*b5*z6+21772800.0_r8*d2*b5*z6&
         +34560.0_r8*z8*b7+967680.0_r8*z8*b6+11128320.0_r8*z8*b5+67737600.0_r8*z8*b4&
         +233936640.0_r8*z8*b3+453841920.0_r8*z8*b2+451630080.0_r8*z8*b&
         +404.0_r8*d8*b+7862400.0_r8*d7*b6*z4-151200000.0_r8*d7*b5*z5&
         +950745600.0_r8*d7*b4*z6-45964800.0_r8*d6*b6*z6+1391040.0_r8*d7*b6*z3&
         -48585600.0_r8*d7*b5*z4+598752000.0_r8*d7*b4*z5-2091156480.0_r8*d7*b3*z6&
         -24192000.0_r8*d6*b6*z5+30240.0_r8*d7*b6*z2-6350400.0_r8*d7*b5*z3&
         +174867840.0_r8*d7*b4*z4-1299110400.0_r8*d7*b3*z5+2272112640.0_r8*d7*b2*z6&
         -5443200.0_r8*d6*b6*z4+27417600.0_r8*d5*b6*z6-15120.0_r8*d7*b6*z-151200.0_r8*d7*b5*z2&
         +19172160.0_r8*d7*b4*z3-366912000.0_r8*d7*b3*z4+1412812800.0_r8*d7*b2*z5&
         -928972800.0_r8*d7*b*z6-120960.0_r8*d6*b6*z3+135.0_r8*d8*b7+1260.0_r8*d8*b6&
         +3150.0_r8*d8*b5+840.0_r8*d8*b4-2345.0_r8*d8*b3+540.0_r8*d8*b2-199065600.0_r8*d7*b*z8&
         -48384000.0_r8*d*b5*z8-870912000.0_r8*d8*z5-352719360.0_r8*d8*z4-66769920.0_r8*d8*z3&
         -4389120.0_r8*d8*z2-696729600.0_r8*d8*z7+13547520.0_r8*d2*b6*z8-5806080.0_r8*d*b6*z8&
         -13547520.0_r8*d3*b6*z8+5806080.0_r8*d7*b6*z8-48384000.0_r8*d7*b5*z8+203212800.0_r8*d7*b4*z8&
         -13547520.0_r8*d6*b6*z8-449003520.0_r8*d7*b3*z8+487710720.0_r8*d7*b2*z8+13547520.0_r8*d5*b6*z8&
         +34560.0_r8*d8*b7*z8-967680.0_r8*d8*b6*z8+11128320.0_r8*d8*b5*z8-276480.0_r8*d7*b7*z8&
         -67737600.0_r8*d8*b4*z8+233936640.0_r8*d8*b3*z8+967680.0_r8*d6*b7*z8-453841920.0_r8*d8*b2*z8&
         -34560.0_r8*d8*z-1103155200.0_r8*d8*z6+967680.0_r8*d2*b7*z8+2419200.0_r8*d4*b7*z8&
         -1935360.0_r8*d3*b7*z8+451630080.0_r8*d8*b*z8-1935360.0_r8*d5*b7*z8-276480.0_r8*d*b7*z8&
         -116121600.0_r8*d6*b*z8-33868800.0_r8*d4*b5*z8-13547520.0_r8*d3*b5*z8+67737600.0_r8*d2*b5*z8&
         +1209600.0_r8*d4*b6*z4-12096000.0_r8*d3*b6*z6-3840.0_r8*d7*b2*z-276480.0_r8*d7*b*z2&
         -2661120.0_r8*d3*b6*z5+3709440.0_r8*d2*b6*z6+7257600.0_r8*d5*b6*z5-25200.0_r8*d7*b5*z&
         +346080.0_r8*d7*b4*z2-37013760.0_r8*d7*b3*z3+399813120.0_r8*d7*b2*z4-580608000.0_r8*d7*b*z5&
         +90720.0_r8*d6*b6*z2-604800.0_r8*d5*b6*z4+4838400.0_r8*d4*b6*z6+1680.0_r8*d7*b4*z&
         -544320.0_r8*d7*b3*z2+40400640.0_r8*d7*b2*z3-166440960.0_r8*d7*b*z4-383040.0_r8*d5*b6*z3&
         +4838400.0_r8*d4*b6*z5+10080.0_r8*d7*b3*z+556800.0_r8*d7*b2*z2-17418240.0_r8*d7*b*z3&
         +28546560.0_r8*d7*b6*z6+20563200.0_r8*d7*b6*z5-229824000.0_r8*d7*b5*z6+241920.0_r8*d2*b7*z6&
         -138240.0_r8*d*b7*z7-2160.0_r8*d7*b7*z+15120.0_r8*d6*b7*z2-60480.0_r8*d5*b7*z3+151200.0_r8*d4*b7*z4&
         -92897280.0_r8*b*d5*z8-67737600.0_r8*b4*d5*z8+108380160.0_r8*b3*d5*z8+&
         54190080.0_r8*b2*d5*z8-116121600.0_r8*d2*b*z8-199065600.0_r8*d*b*z8-54190080.0_r8*d3*b2*z8&
         -92897280.0_r8*d3*b*z8-149022720.0_r8*d2*b2*z8-487710720.0_r8*d*b2*z8+118540800.0_r8*d4*b3*z8&
         -87091200.0_r8*d4*b*z8+108380160.0_r8*d3*b3*z8+47416320.0_r8*d2*b3*z8-449003520.0_r8*d*b3*z8&
         +67737600.0_r8*d3*b4*z8+135475200.0_r8*d2*b4*z8-203212800.0_r8*d*b4*z8+67737600.0_r8*d6*b5*z8&
         -135475200.0_r8*d6*b4*z8+47416320.0_r8*d6*b3*z8-13547520.0_r8*d5*b5*z8+149022720.0_r8*d6*b2*z8&
         +52080.0_r8*d8*b3*z-10680960.0_r8*d8*b2*z2+169102080.0_r8*d8*b*z3+2903040.0_r8*d6*b7*z7+&
         15120.0_r8*d8*b7*z2-887040.0_r8*d8*b6*z3+27720000.0_r8*d8*b5*z4-348364800.0_r8*d8*b4*z5&
         +1482727680.0_r8*d8*b3*z6-1815367680.0_r8*d8*b2*z7-604800.0_r8*d7*b7*z4+3628800.0_r8*d6*b7*z6&
         +2160.0_r8*d8*b7*z-120960.0_r8*d8*b6*z2+6572160.0_r8*d8*b5*z3-147490560.0_r8*d8*b4*z4&
         +1173070080.0_r8*d8*b3*z5-2870784000.0_r8*d8*b2*z6+1806520320.0_r8*d8*b*z7&
         -181440.0_r8*d7*b7*z3+2419200.0_r8*d6*b7*z5-4838400.0_r8*d5*b7*z7+&
         635040.0_r8*d8*b5*z2-30481920.0_r8*d8*b4*z3+477378720.0_r8*d8*b3*z4&
         -2258565120.0_r8*d8*b2*z5+2857559040.0_r8*d8*b*z6-30240.0_r8*d7*b7*z2+907200.0_r8*d6*b7*z4&
         -4838400.0_r8*d5*b7*z6+138240.0_r8*d8*b7*z7+241920.0_r8*d8*b7*z6-3870720.0_r8*d8*b6*z7&
         +241920.0_r8*d8*b7*z5-6451200.0_r8*d8*b6*z6+44513280.0_r8*d8*b5*z7-967680.0_r8*d7*b7*z7&
         +151200.0_r8*d8*b7*z4-5806080.0_r8*d8*b6*z5+72092160.0_r8*d8*b5*z6-270950400.0_r8*d8*b4*z7&
         -1451520.0_r8*d7*b7*z6+60480.0_r8*d8*b7*z3-3024000.0_r8*d8*b6*z4+60480000.0_r8*d8*b5*z5&
         -432230400.0_r8*d8*b4*z6+935746560.0_r8*d8*b3*z7-1209600.0_r8*d7*b7*z5&
         -241920.0_r8*d3*b7*z5-7257600.0_r8*d3*b4*z6-144.0_r8*d8-38707200.0_r8*d3*b2*z6&
         +101606400.0_r8*d6*b5*z5+17539200.0_r8*d6*b5*z4-193536000.0_r8*d6*b4*z5+604800.0_r8*d6*b5*z3&
         -29756160.0_r8*d6*b4*z4+70156800.0_r8*d6*b3*z5-2419200.0_r8*d5*b5*z5+151200.0_r8*d6*b5*z2&
         -725760.0_r8*d6*b4*z3+11551680.0_r8*d6*b3*z4+217728000.0_r8*d6*b2*z5-201600.0_r8*d5*b5*z4&
         +53760.0_r8*d6*b4*z2+181440.0_r8*d6*b3*z3+35199360.0_r8*d6*b2*z4-174182400.0_r8*d6*b*z5&
         -826560.0_r8*d5*b5*z3+2419200.0_r8*d4*b5*z5-31920.0_r8*d6*b3*z2+846720.0_r8*d6*b2*z3&
         -29998080.0_r8*d6*b*z4+3528000.0_r8*d4*b5*z4-6240.0_r8*d6*b2*z2-967680.0_r8*d6*b*z3&
         -10886400.0_r8*d3*b5*z5+3840.0_r8*d6*b*z2+91344960.0_r8*d8*b3*z3-907643520.0_r8*d8*b2*z4&
         +2249856000.0_r8*d8*b*z5+181440.0_r8*d6*b7*z3-2419200.0_r8*d5*b7*z5+4838400.0_r8*d4*b7*z7&
         -26880.0_r8*d8*b4*z+6170640.0_r8*d8*b3*z2-168940800.0_r8*d8*b2*z3+905627520.0_r8*d8*b*z4&
         -604800.0_r8*d5*b7*z4+3628800.0_r8*d4*b7*z6+1209600.0_r8*d4*b7*z5-2903040.0_r8*d3*b7*z7&
         -76800.0_r8*d8*b2*z+10735680.0_r8*d8*b*z2-1451520.0_r8*d3*b7*z6+73920.0_r8*d8*b*z+&
         967680.0_r8*d2*b7*z7+10080.0_r8*d8*b5*z-2365440.0_r8*d8*b4*z2-174182400.0_r8*d8*z8&
         -20321280.0_r8*d3*b6*z7+20321280.0_r8*d7*b6*z7-169344000.0_r8*d7*b5*z7+711244800.0_r8*d7*b4*z7&
         -40642560.0_r8*d6*b6*z7-1571512320.0_r8*d7*b3*z7+1706987520.0_r8*d7*b2*z7+33868800.0_r8*d5*b6*z7&
         -696729600.0_r8*d7*b*z7-24192000.0_r8*d*b5*z7+13547520.0_r8*d2*b6*z7-2903040.0_r8*d*b6*z7&
         +101606400.0_r8*d3*b4*z7+135475200.0_r8*d2*b4*z7-101606400.0_r8*d*b4*z7&
         +203212800.0_r8*d6*b5*z7-406425600.0_r8*d6*b4*z7+142248960.0_r8*d6*b3*z7&
         -33868800.0_r8*d5*b5*z7+447068160.0_r8*d6*b2*z7-348364800.0_r8*d6*b*z7&
         -67737600.0_r8*d4*b5*z7-20321280.0_r8*d3*b5*z7+67737600.0_r8*d2*b5*z7&
         -232243200.0_r8*b*d5*z7-169344000.0_r8*b4*d5*z7+270950400.0_r8*b3*d5*z7&
         +135475200.0_r8*b2*d5*z7-116121600.0_r8*d2*b*z7-99532800.0_r8*d*b*z7&
         -81285120.0_r8*d3*b2*z7-139345920.0_r8*d3*b*z7-149022720.0_r8*d2*b2*z7&
         -243855360.0_r8*d*b2*z7+237081600.0_r8*d4*b3*z7-174182400.0_r8*d4*b*z7&
         +162570240.0_r8*d3*b3*z7+47416320.0_r8*d2*b3*z7-224501760.0_r8*d*b3*z7+174182400.0_r8*z8)*b/d8
    argu=alpha*tau*0.5_r8
    nu1=mum
    nu2=mum-1.0_r8
    IF (nu2<0) nu2=-nu2     
    CALL bessJYN(nu2, argu, jax2, jap2, yax2, yap2, ierr)
    CALL bessJYN(nu1, argu, jax1, jap1, yax1, yap1, ierr)
    IF (mum-1.0_r8<0) THEN
      jax2=-jax2
      yax2=-yax2 
    ENDIF
    IF (mu==0) THEN
      fac2=-0.5_r8*sqrt(pi)
    ELSE   
      fac2=-0.5_r8*sqrt(pi)*(tau/alpha)**mu
    ENDIF   
    beta=0.5_r8*(alpha*tau+pi)
    cosb=cos(beta)
    sinb=sin(beta)
    rephi0=fac2*(cosb*jax1+sinb*yax1)
    imphi0=fac2*(sinb*jax1-cosb*yax1)
    fac3=-fac2*0.5_r8*alpha
    j2=yax1+jax2
    y2=yax2-jax1
    imphi1=-fac3*(cosb*j2+sinb*y2)
    rephi1=fac3*(sinb*j2-cosb*y2)
    g=tau*log(x+w)
    sing=sin(g)
    cosg=cos(g)
    rphik(0)=rephi0
    iphik(0)=imphi0
    rphik(1)=rephi1
    iphik(1)=imphi1     
    s=rphik(0)*cosg+iphik(0)*sing+fk(1)*(rphik(1)*cosg+iphik(1)*sing);
    taui=1.0_r8/tau
    atau=alpha*taui 
    DO k=1,6
      IF (mu==0) THEN
        coef1=k*taui
        coef2=(k-0.5_r8)*atau
      ELSE   
        coef1=(k-2.0_r8*mu)*taui
        coef2=(k-mu-0.5_r8)*atau
      ENDIF  
      rphik(k+1)=coef1*iphik(k)-alpha*rphik(k)+coef2*iphik(k-1);
      iphik(k+1)=-coef1*rphik(k)-alpha*iphik(k)-coef2*rphik(k-1);
      s=s+fk(k+1)*(rphik(k+1)*cosg+iphik(k+1)*sing)
    ENDDO 
    ras=sqrt(pi*0.5_r8)*alpha**(mu+0.5_r8)/sqrt(w)*s    
    END SUBROUTINE conirkum
        
    SUBROUTINE conirlar(mu, tau, x, rmu, ierr)
    !--------------------------------------------------------  
    ! Computation of rmu for moderate or large values of x
    !--------------------------------------------------------  
    USE Someconstants      
    USE gammaCHI
    REAL(r8), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: tau
    REAL(r8), INTENT(OUT) :: rmu 
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: w,z,phi,s,zz,tauk,factk,rn,an,taun,cn,gamn,gamk
    REAL(r8) :: tau2, gre1, gre2, gim1, gim2, gmod2, r1, r2, rmod, rpha, phi2,&
                uk(0:100), vk(0:100), wk(0:100), fac, t, uvmod, uvpha, ang, &
                tkm1, ss
    REAL(r8) :: r1a, r2a
    INTEGER :: mu2,k,n
    ierr=0
    w=sqrt(x*x-1.0_r8);
    z=(-1.0_r8/(2.0_r8*w*(w+x)));
    phi=tau*log(x+w);
    tau2=tau*tau
    mu2=mu*mu 
    ! Computation of the ratio of gamma functions. 
    ! The gamma functions are computed using our algorithm
    ! implemented in gammac
    CALL gammac(0.5_r8+mu,tau,gre1,gim1)
    CALL gammac(1.0_r8,tau,gre2,gim2)
    gmod2=gre2*gre2+gim2*gim2
    r1=(gre1*gre2+gim1*gim2)/gmod2
    r2=(gim1*gre2-gim2*gre1)/gmod2  
    rmod=sqrt(r1*r1+r2*r2)
    rpha=phase(r1,r2)
    phi2=phi-rpha
    uk(0)=1.0_r8;
    vk(0)=0.0_r8;
    wk(0)=1.0_r8
    k=1;
    s=cos(phi2)
    zz=1.0_r8
    fac=1.0_r8
    t=1.0_r8
    DO WHILE (t>1.e-20_r8)
      zz=z*zz 
      uk(k)=k*uk(k-1)+tau*vk(k-1);
      vk(k)=k*vk(k-1)-tau*uk(k-1);
      wk(k)=(k*k+tau2)*wk(k-1)
      uvmod=sqrt(uk(k)*uk(k)+vk(k)*vk(k))
      uvpha=atan(vk(k)/uk(k))
      uvpha=phase(uk(k),vk(k))       
      ang=phi2-uvpha
      tkm1=2*k-1.0_r8
      fac=(tkm1*tkm1*0.25_r8-mu2)*fac/k
      ss=fac*zz*cos(ang)*uvmod/wk(k)
      s=s+ss
      t=abs(ss/s)
      k=k+1;
    ENDDO
    rmu=sqrt(pi*0.5_r8)*rmod*s/sqrt(w);
    END SUBROUTINE conirlar
  
    RECURSIVE SUBROUTINE gammac(zre,zim,gre,gim)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN)  :: zre, zim
    REAL(r8), INTENT(OUT) :: gre, gim
    COMPLEX(r8) :: uyi, z, gamm, gam
    REAL(r8) :: zre2
    INTEGER :: k, k1, n
    k=nint(zre)
    k1=k-1
    uyi=dcmplx(0.0_r8,1.0_r8)
    IF (zre<0.45_r8) THEN
      z=dcmplx(zre,zim)  
      CALL gammac(1.0_r8-zre,zim,gre,gim)
      gamm=dcmplx(gre,gim)
      gam=pi/(sin(pi*z)*gamm)
      gre=real(gam)
      gim=aimag(gam)
    ELSEIF (zre<3.0_r8) THEN
      IF (k>zre) THEN
        k=k1
      ENDIF
      k1=3-k;
      zre2=k1+zre;
      z=dcmplx(zre2,zim)
      CALL gammac(zre2,zim,gre,gim)
      gamm=dcmplx(gre,gim)
      DO n=1,k1 
        gamm=gamm/(z-n)
      ENDDO
      gre=real(gamm)
      gim=aimag(gamm)
    ELSE
      z=dcmplx(zre,zim)   
      gam=sqrttwopi*exp(-z+(z-0.5_r8)*log(z)+stirlingc(z)) 
      gre=real(gam)
      gim=aimag(gam) 
    ENDIF       
    END SUBROUTINE gammac

    FUNCTION stirlingc(x)            
    !Stirling series, function corresponding with}
    !asymptotic series for log(gamma(z))}
    !that is:  1/(12z)-1/(360z**3)...; |z|>= 3}
    !Complex argument
    !----------------------------------------  
    USE Someconstants
    IMPLICIT NONE
    COMPLEX(r8) :: stirlingc, x, z
    REAL(r8) a(0:17), c(0:6)
    IF (abs(x)<12.0_r8) THEN
      a(0)=1.996379051590076518221_r8;
      a(1)=-0.17971032528832887213e-2_r8;
      a(2)=0.131292857963846713e-4_r8;
      a(3)=-0.2340875228178749e-6_r8;
      a(4)=0.72291210671127e-8_r8;
      a(5)=-0.3280997607821e-9_r8;
      a(6)=0.198750709010e-10_r8;
      a(7)=-0.15092141830e-11_r8;
      a(8)=0.1375340084e-12_r8;
      a(9)=-0.145728923e-13_r8;
      a(10)=0.17532367e-14_r8;
      a(11)=-0.2351465e-15_r8;
      a(12)=0.346551e-16_r8;
      a(13)=-0.55471e-17_r8;
      a(14)=0.9548e-18_r8;
      a(15)=-0.1748e-18_r8;
      a(16)=0.332e-19_r8;
      a(17)=-0.58e-20_r8;
      z=18.0_r8/(x*x)-1.0_r8;
      stirlingc=chepolsumc(17,z,a)/(12.0_r8*x);
    ELSE
      z=1.0_r8/(x*x);
      IF (abs(x)<1000.0_r8) THEN
        c(0)=0.25721014990011306473e-1_r8;
        c(1)=0.82475966166999631057e-1_r8;
        c(2)=-0.25328157302663562668e-2_r8;
        c(3)=0.60992926669463371e-3_r8;
        c(4)=-0.33543297638406e-3_r8;
        c(5)=0.250505279903e-3_r8;
        c(6)=0.30865217988013567769_r8;
        stirlingc=((((((c(5)*z+c(4))*z+c(3))*z+c(2))*z+c(1))*z+c(0))/(c(6)+z)/x)
      ELSE
        stirlingc=(((-z*0.000595238095238095238095238095238_r8+&
                 0.000793650793650793650793650793651_r8)*z&
                -0.00277777777777777777777777777778_r8)*z+&
                 0.0833333333333333333333333333333_r8)/x
      ENDIF
    ENDIF 
    END FUNCTION stirlingc

    FUNCTION chepolsumc(n,t,ak)
    USE Someconstants
    IMPLICIT NONE
    COMPLEX(r8), INTENT(IN) :: t
    REAL(r8), INTENT(IN) :: ak(0:8)
    INTEGER, INTENT(IN) :: n
    COMPLEX(r8) :: chepolsumc, u0, u1, u2, s, tt;
    INTEGER :: k
    u0=0; u1=0; k=n; tt=t+t;
    DO WHILE (k>=0)
      u2=u1; 
      u1=u0; 
      u0=tt*u1-u2+ ak(k); 
      k= k-1 
    ENDDO
    s=(u0-u2)/2.0_r8
    chepolsumc=s
    END FUNCTION chepolsumc
    
    FUNCTION phase(x,y)
    !-----------------------------------------  
    ! Computes the phase of z=x+iy in (-pi,pi)
    !-----------------------------------------  
    USE Someconstants 
    IMPLICIT NONE
    REAL(r8) :: phase, x, y, p, ay
    IF ((x==0.0_r8).AND.(y==0.0_r8)) THEN
      p=0.0_r8
    ELSE
      ay=abs(y) 
      IF (x>=ay) THEN
        p=atan(ay/x)
      ELSEIF ((x+ay)>=0.0_r8) THEN
        p=pihalf-atan(x/ay)
      ELSE
        p=pi+atan(ay/x)
      ENDIF
      IF (y<0.0_r8) p=-p
    ENDIF
    phase=p
    END FUNCTION phase

   FUNCTION psizr(x,y)
   !Computes psi(x+I*y) with real computations, x > 0
   USE Someconstants 
   IMPLICIT NONE
   COMPLEX(r8) :: psizr, uyi
   REAL(r8) :: x, y      
   REAL(r8) :: y2, x0, x02, x1, x2, s, t, ck, uk(0:9), vk(0:9), u, v, sre, sim, tre,&
               tim, r, r2, ber(0:9) 
   INTEGER :: k, m
   uyi=dcmplx(0.0_r8,1.0_r8)
   ber(0)=1.0_r8;
   ber(1)=.166666666666666666666666666667_r8;
   ber(2)=-0.333333333333333333333333333333e-1_r8;
   ber(3)=0.238095238095238095238095238095e-1_r8
   ber(4)=-0.333333333333333333333333333333e-1_r8
   ber(5)=0.757575757575757575757575757576e-1_r8
   ber(6)=-.253113553113553113553113553114_r8
   ber(7)=1.16666666666666666666666666667_r8
   ber(8)=-7.09215686274509803921568627451_r8
   ber(9)=54.9711779448621553884711779449_r8
   y2=y*y;
   r=x*x+y2;
   m=0;
   IF (r<144.0_r8) THEN
     x1=sqrt(144.0_r8-y2)-x;
     m=1+int(x1)
  ENDIF
   x0=x+m;
   x02=x0*x0;
   r=x02+y2;
   r2=r*r;
   u=x02-y2;
   v=2.0_r8*x0*y;
   uk(0)=1.0_r8;
   vk(0)=0.0_r8;
   sre=0.0_r8;
   sim=0.0_r8;
   k=1;
   t=1.0_r8;
   DO WHILE ((k < 9).AND.(t>1.0e-15_r8))
     ck=0.5_r8*ber(k)/k; 
     uk(k)=(u*uk(k-1)+v*vk(k-1))/r2;
     vk(k)=(u*vk(k-1)-v*uk(k-1))/r2;
     tre=ck*uk(k);
     tim=ck*vk(k);
     sre=sre+tre;
     sim=sim+tim;
     k=k+1;
     t=sqrt(tre*tre+tim*tim);
   ENDDO
   sre=log(r)*0.5_r8-0.5_r8*x0/r-sre;
   sim=atan(y/x0)+0.5_r8*y/r-sim;
   DO k=0,m-1 
     r=(x+k)*(x+k)+y2;
     sre=sre-(x+k)/r;
     sim= sim+y/r
   ENDDO
   psizr=sre+uyi*sim
   psizr=sre
   END FUNCTION psizr

    FUNCTION pochhamc(a,n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    COMPLEX(r8), INTENT(IN) :: a
    COMPLEX(r8) :: pochhamc, p
    INTEGER :: k
    p= 1.0_r8;
    DO k=1,n 
      p= p*(a+k-1.0_r8) 
    ENDDO 
    pochhamc=p
    END FUNCTION pochhamc 
   
    FUNCTION pochham(a,n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(r8), INTENT(IN) :: a
    REAL(r8) :: pochham, p
    INTEGER :: k
    p= 1.0_r8;
    DO k=1,n 
      p= p*(a+k-1.0_r8) 
    ENDDO 
    pochham=p
    END FUNCTION pochham

    FUNCTION gammah(n)
    USE Someconstants
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(r8) :: gammah
    INTEGER :: i,j
    i=n
    j=2*i-1
    gammah=1.0_r8
    DO WHILE ((j>=1).AND.(gammah<giant)) 
       gammah=gammah*j/2.0_r8
       i=i-1
       j=2*i-1
    ENDDO  
    gammah=gammah*sqrtpi
    IF (j>1) THEN 
       gammah=0.0_r8
    ENDIF   
    END FUNCTION gammah
                                  
    FUNCTION xpowy(x,y)
    IMPLICIT NONE
    REAL(r8) :: x,y,xpowy
    xpowy=x**y
    END FUNCTION xpowy       

    SUBROUTINE expaelem(x,mu,tau,pmmu,ierr)
    ! --------------------------------------------------
    ! Calculation of P^(mu)_(-1/2+i*tau)(x) by using 
    ! asymptotic expansions in terms of elementary 
    ! functions, -1<x<1 (also x>1, in the monotonic 
    ! region)
    ! --------------------------------------------------
    ! Inputs:
    !   x ,    argument of the functions
    !   mu,    function parameter (integer value)
    !   tau,   function parameter (real value)
    ! Outputs:
    !   pmmu , P^(mu)_(-1/2+i*tau)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! ---------------------------------------------------------
    USE Someconstants
    REAL(r8), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: tau
    REAL(r8), INTENT(OUT) :: pmmu
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: beta, argu2, p, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13,&
                p14, p15, p16, p17, p18, p19, p20, p21, beta2, beta4, beta6, beta8, beta10,&
                beta12, beta14, mu2, mu3, mu4, mu5, mu6, mu7, tau2, be2p1, be2p12, be2p13, be2p14,&
                be2p15, be2p16, be2p17, p1pbe,  phit0, S, fac2, Pas 
    REAL(r8), DIMENSION(1:15) :: uk
    INTEGER n
    ierr=0
    beta=tau/mu
    beta2=beta*beta
    beta4=beta2*beta2
    beta6=beta2*beta4
    beta8=beta2*beta6
    beta10=beta2*beta8    
    beta12=beta2*beta10   
    beta14=beta2*beta12 
    mu2=mu*mu
    mu3=mu2*mu
    mu4=mu3*mu 
    mu5=mu4*mu
    mu6=mu5*mu
    mu7=mu6*mu
    tau2=tau*tau  
    IF (x<1) THEN
      argu2=sqrt(1.0_r8-x*x)
    ELSE
      argu2=sqrt(x*x-1.0_r8)
    ENDIF  
    p=x/sqrt(1.0_r8+beta2*(1.0_r8-x*x))
    p2=p*p
    p3=p2*p
    p4=p3*p
    p5=p4*p
    p6=p5*p
    p7=p6*p
    p8=p7*p
    p9=p8*p
    p10=p9*p
    p11=p10*p
    p12=p11*p
    p13=p12*p   
    p14=p13*p
    p15=p14*p
    p16=p14*p
    p17=p16*p
    p18=p17*p
    p19=p18*p
    p20=p19*p
    p21=p20*p
    be2p1=1.0_r8+beta2
    be2p12=be2p1*be2p1
    be2p13=be2p12*be2p1
    be2p14=be2p13*be2p1 
    be2p15=be2p14*be2p1 
    be2p16=be2p15*be2p1
    be2p17=be2p16*be2p1
    p1pbe=p*be2p1
    phit0=log(x*(1.0_r8+p)/p1pbe)+beta*acos(x*(1.0_r8-p*beta2)/p1pbe)
    uk(1)= -1.0_r8/24.0_r8*(-beta2+5*beta2*p3-3*beta2*p+3*p)/be2p1;
    uk(2)= 1.0_r8/1152.0_r8*(462*beta2*p4+385*beta4*p6-522*beta2*p2+6*beta4*p+ &
          81*beta4*p2+81*p2-462*beta4*p4+beta4+72*beta2-6*beta2*p-10*beta4*p3-72)&
          /be2p12;
    uk(3)= -1.0_r8/414720.0_r8*(-32400*p+194400*beta2*p+7128*beta2-1300806*beta4*p5-6930*&
           beta4*p4-564165*beta2*p3-32355*beta4*p+7830*beta4*p2-1215*beta2*p2&
           +6930*beta6*p4+564165*beta4*p3-5775*beta6*p6+30375*p3+1944*beta4&
           +1003*beta6-30300*beta6*p3+369603*beta2*p5-765765*beta6*p7-1215*beta6*p2&
           +425425*beta6*p9+369603*beta6*p5-45*beta6*p+765765*beta4*p7)/be2p13;
    uk(4)= 1.0_r8/39813120.0_r8*(1632960+1606608*beta4-1701700*beta8*p9-446185740*beta8*p10&
            -12036*beta8*p-1478412*beta8*p5-94110126*beta8*p6+2430*beta8*p2+446185740*beta6*p10&
            -4027*beta8+349922430*beta4*p8+94121676*beta2*p6-748417428*beta4*p6+490305150*beta4*p4&
            -1478412*beta4*p5+2377620*beta4*p3+4465125*p4+141560*beta8*p3+4451265*beta8*p4+3063060*beta8*p7&
            +349922430*beta8*p8-178143300*beta2*p4+185910725*beta8*p12+105348*beta6*p-93484530*beta4*p2&
            +748417428*beta6*p6-1022624460*beta6*p8+5203224*beta6*p5-3063060*beta6*p7-178129440*beta6*p4&
            -2196180*beta6*p3+6107940*beta6*p2-121500*beta2*p3-9936*beta6+202176*beta2*p-9486720*beta2&
            -813888*beta4*p+93486960*beta2*p2-6123600*p2)/be2p14;
     uk(5)= -1.0_r8/6688604160.0_r8*(284499769554.0_r8*beta4*p9-2571912000.0_r8*p3+43742592.0_r8*beta4&
              +1686284061462.0_r8*beta8*p9-3123300180.0_r8*beta8*p10-284151*beta8*p+110699860005.0_r8*beta8*p5&
              -5222622636.0_r8*beta8*p6-46765908*beta8*p2+614135872350.0_r8*beta10*p11+188699385875.0_r8*beta10*p15&
              -1301375075.0_r8*beta10*p12-566098157625.0_r8*beta10*p13-566195*beta10*p3&
              -37743279*beta10*p4+49276227897.0_r8*beta10*p7-2449457010.0_r8*beta10*p8&
              -284493813604.0_r8*beta10*p9+3123300180.0_r8*beta10*p10+84567*beta10*p&
              -1513861083.0_r8*beta10*p5+664257902.0_r8*beta10*p6+1137402*beta10*p2-5128423*beta10&
              -24659280*beta8+49286948607.0_r8*beta2*p7-724843496988.0_r8*beta4*p7-110718071289.0_r8*beta2*p5&
              -658851732.0_r8*beta4*p6+1286121564.0_r8*beta4*p4+622277581842.0_r8*beta4*p5-197794576350.0_r8*beta4*p3&
              -2564437050.0_r8*beta8*p3+1233931272.0_r8*beta8*p4-724832776278.0_r8*beta8*p7+7158371220.0_r8*beta8*p8&
              -31255875.0_r8*beta2*p4+566098157625.0_r8*beta8*p13-1060081344.0_r8*beta6*p&
              -1635658724700.0_r8*beta8*p11-695178288*beta4*p2+1063056960.0_r8*p+5271520716.0_r8*beta6*p6&
              -1686284061462.0_r8*beta6*p9-2449457010.0_r8*beta6*p8-622272407400.0_r8*beta6*p5+&
              1618196016762.0_r8*beta6*p7-3451695282.0_r8*beta6*p4+77336405370.0_r8*beta6*p3+&
              640305162.0_r8*beta6*p2+1519035525.0_r8*p5+77345150400.0_r8*beta2*p3-60734160*beta6&
              -15939322560.0_r8*beta2*p-94665024*beta2+15938360928.0_r8*beta4*p&
              +49723632*beta2*p2+614135872350.0_r8*beta6*p11)/be2p15; 
     uk(6)= 1.0_r8/4815794995200.0_r8*(-5817276719424.0_r8*beta4+5104696716244125.0_r8*beta8*p14&
             +50660565015060.0_r8*beta8*p9+21585355164186606.0_r8*beta8*p10+36918192960.0_r8*beta8*p&
             +18631155590712.0_r8*beta8*p5+3246246809762985.0_r8*beta8*p6+3198809877165.0_r8*beta8*p2&
             -12854717138513250.0_r8*beta10*p14+49069761741000.0_r8*beta10*p11+3685299006138750.0_r8*beta10&
             *p16+17102195452922580.0_r8*beta10*p12-16982944728750.0_r8*beta10*p13+86578912260.0_r8*beta10*p3&
             +5551483305690.0_r8*beta10*p4+21701971799820.0_r8*beta10*p7+3086709988029210.0_r8*beta10*p8&
             -50552500258260.0_r8*beta10*p9-10694231441748504.0_r8*beta10*p10+938847690.0_r8*beta10*p&
             -3326778949158.0_r8*beta10*p5-330805922867910.0_r8*beta10*p6+760821390.0_r8*beta10*p2&
             +825861096.0_r8*beta10+168359651.0_r8*beta12+5104696716244125.0_r8*beta12*p14&
             -18424076170500.0_r8*beta12*p11-5660981576250.0_r8*beta12*p15-3685299006138750.0_r8*beta12*p16&
             -3369012547635735.0_r8*beta12*p12+16982944728750.0_r8*beta12*p13-1622314950.0_r8*beta12*p3&
             +1023694168371875.0_r8*beta12*p18+ 664187895.0_r8*beta12*p4-1500114202470.0_r8*beta12*p7&
             -127540556499600.0_r8*beta12*p8+8546940722320.0_r8*beta12*p9+1050713924955201.0_r8*beta12*p10&
             +461558070.0_r8*beta12*p+55950996402.0_r8*beta12*p5+2747003910420.0_r8*beta12*p6&
             -34249635.0_r8*beta12*p2-10694278291251204.0_r8*beta6*p10+1951093656.0_r8*beta8&
             -3086817363597510.0_r8*beta4*p8+127577298354750.0_r8*beta2*p8-1478608458210.0_r8*beta4*p7&
             -330884017717050.0_r8*beta2*p6-45571065750.0_r8*beta2*p5+3246326371554525.0_r8*beta4*p6&
             -1430972525967075.0_r8*beta4*p4+3384132189102.0_r8*beta4*p5-2424891245760.0_r8*beta4*p3&
             -5569796925000.0_r8*p4-2301274458180.0_r8*beta8*p3-287851315212090.0_r8*beta8*p4&
             -48610719356940.0_r8*beta8*p7-12549445147755630.0_r8*beta8*p8+287903384029800.0_r8*beta2*p4&
             +1050760774457901.0_r8*beta4*p10-17102195452922580.0_r8*beta8*p12-470422624320.0_r8*beta6*p&
             -18424076170500.0_r8*beta8*p11-394457817600.0_r8+226528604270640.0_r8*beta4*p2&
             -6565164928398720.0_r8*beta6*p6-8534993086620.0_r8*beta6*p9+12549481889610780.0_r8*beta6*p8&
             -18857316121308.0_r8*beta6*p5+21874982617800.0_r8*beta6*p7+1430953116204780.0_r8*beta6*p4&
             +5978615182380.0_r8*beta6*p3-90407686919760.0_r8*beta6*p2+2757049477875.0_r8*p6&
             +82301184000.0_r8*beta2*p3+394653741696.0_r8*beta6-44320867200.0_r8*beta2*p+5820457305600.0_r8&
             *beta2+512985052800.0_r8*beta4*p-90418726137600.0_r8*beta2*p2+3208203028800.0_r8*p2+&
             3369032068261860.0_r8*beta6*p12)/be2p16;
      uk(7)= -1.0_r8/115579079884800.0_r8*(-460119768008171850.0_r8*beta4*p9+138799253740521843.0_r8*beta4*p11&
             +347555328120000.0_r8*p3-7198694208000.0_r8*beta4-5491739740694062590.0_r8*beta8*p9&
             +10713168010739844.0_r8*beta8*p10-77949544996080.0_r8*beta8*p-319704673216628322.0_r8*beta8*p5&
             +6581529649081800.0_r8*beta8*p6+89694981527760.0_r8*beta8*p2-5104696716244125.0_r8*beta10*p14&
             -7566929675484345180.0_r8*beta10*p11+1570320948552481125.0_r8*beta10*p17&
             -6208169160222604575.0_r8*beta10*p15+17110066169376180.0_r8*beta10*p12+&
              9708066889428193650.0_r8*beta10*p13-346412669599215.0_r8*beta10*p3+287941809852270.0_r8&
              *beta10*p4-569838665576990925.0_r8*beta10*p7+12545106096419010.0_r8*beta10*p8&
              +3028730899409633295.0_r8*beta10*p9-21594800023930926.0_r8*beta10*p10-21635411040.0_r8*beta10*p&
             +38163770321137083.0_r8*beta10*p5-3240019737486945.0_r8*beta10*p6-3602562002505.0_r8*beta10*p2&
             +1421661369384.0_r8*beta10+471807529704.0_r8*beta12+12854717138513250.0_r8*beta12*p14&
             +2198267726616662715.0_r8*beta12*p11-3808508329237862250.0_r8*beta12*p17+&
             6208169160222604575.0_r8*beta12*p15-3685299006138750.0_r8*beta12*p16-17098260094695780.0_r8*&
             beta12*p12-5105814396159838725.0_r8*beta12*p13-47628489105.0_r8*beta12*p3+&
             931766432052080625.0_r8*beta12*p19-6583668032970.0_r8*beta12*p4+36681022640663100.0_r8*beta12*p7&
             -3086590051932930.0_r8*beta12*p8-460094509768835520.0_r8*beta12*p9+&
             10687966101587424.0_r8*beta12*p10-1150513095.0_r8*beta12*p-467533548497547.0_r8*beta12*p5&
             +333998042835870.0_r8*beta12*p6+51959928330.0_r8*beta12*p2-1050760774457901.0_r8*beta6*p10&
             -3447494407896.0_r8*beta14*p6-931766432052080625.0_r8*beta14*p19+1933966502784.0_r8*beta8&
             -1023694168371875.0_r8*beta14*p18+12493049053044375.0_r8*beta2*p9-127577298354750.0_r8*beta4*p8&
             -36691852040413425.0_r8*beta2*p7+569863003356096435.0_r8*beta4*p7-2757049477875.0_r8*beta2*p6&
             +38173067351808480.0_r8*beta2*p5+334868752992186.0_r8*beta4*p6-295885306798200.0_r8*beta4*p4&
             -319714173583489125.0_r8*beta4*p5+76639466793153600.0_r8*beta4*p3+613221795981706275.0_r8*beta6*p13&
             +16158998818811925.0_r8*beta8*p3-1429564436066340.0_r8*beta8*p4+1987125347761092225.0_r8*beta8*p7&
             -12585368560751100.0_r8*beta8*p8+1698039130.0_r8*beta14*p3-33204271842.0_r8*beta14*p5&
             +199689155040375.0_r8*p7+221849150488590625.0_r8*beta14*p21-505078953.0_r8*beta14*p+&
              5758832457000.0_r8*beta2*p4+613213304509341900.0_r8*beta14*p13+3370337347462085.0_r8*beta14*&
              p12+3685299006138750.0_r8*beta14*p16-5105822887632203100.0_r8*beta8*p13&
              -3369032068261860.0_r8*beta8*p12+2188420277457600.0_r8*beta6*p+7566938887522430430.0_r8*beta8*p11&
              +1347119637570231525.0_r8*beta8*p15+94819714767360.0_r8*beta4*p2-78180246144000.0_r8*p&
              +68168266699.0_r8*beta14-3276464475636765.0_r8*beta6*p6+3028756265713726425.0_r8*beta6*p9&
              +3101631679593990.0_r8*beta6*p8+617897733369567273.0_r8*beta6*p5-1987135610787026100.0_r8*beta6*p7&
              +1448347343094225.0_r8*beta6*p4-76638243241186560.0_r8*beta6*p3-228506150206800.0_r8*beta6*p2&
              -1053893444538441.0_r8*beta14*p10-1347116807079443400.0_r8*beta14*p15+12488769564195740.0_r8*beta14*p9&
              -469199692962000.0_r8*p5-16162007951678400.0_r8*beta2*p3-6231033945.0_r8*beta14*p2+&
               8593152421824.0_r8*beta6+1570320948552481125.0_r8*beta14*p17+2188694172672000.0_r8*beta2*p&
              +2458168404480.0_r8*beta2-5104696716244125.0_r8*beta14*p14-138790041702436593.0_r8*beta14*p11&
              +130034103735780.0_r8*beta14*p8-198928264661685.0_r8*beta14*p7-5471599704570432.0_r8*beta4*p&
              -3561174331200.0_r8*beta2*p2-2198292261497533215.0_r8*beta6*p11+66889614015.0_r8*beta14*p4)/be2p17;
    S= 1.0_r8+uk(1)/mu+uk(2)/mu2+1*uk(3)/mu3+uk(4)/mu4+uk(5)/mu5+uk(6)/mu6+uk(7)/mu7  
    fac2=1.0_r8  
    DO n=1,mu
      fac2=argu2*fac2
    ENDDO
    fac2=sqrt(p/(x*mu))*fac2
    pmmu=exp(-mu*phit0)*gammah(mu)*fac2*cosh(pi*tau)/pi*S
    END SUBROUTINE expaelem

    SUBROUTINE expaosc(x,mu,tau,pmmu,ierr)
    ! --------------------------------------------------
    ! Calculation of P^(mu)_(-1/2+i*tau)(x) by using 
    ! asymptotic expansions in terms of elementary 
    ! functions, (x>1, in the oscillatory region)
    ! --------------------------------------------------
    ! Inputs:
    !   x ,    argument of the functions
    !   mu,    function parameter (integer value)
    !   tau,   function parameter (real value)
    ! Outputs:
    !   pmmu , P^(mu)_(-1/2+i*tau)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! ---------------------------------------------------------
    USE Someconstants
    REAL(r8), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: tau
    REAL(r8), INTENT(OUT) :: pmmu
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: beta, beta2, beta4, beta6, beta8, x2, q, q2, q4, q6, q8, q10, q12,&
                mu2, mu3, mu4, be2p1, be2p12, be2p13, be2p14, xc, chxi, xi, shxi,&
                U, V, chi, S, argu2, fac2, fac
    REAL(r8), DIMENSION(0:4) :: uk, vk
    INTEGER n
    ierr=0
    beta=tau/mu
    beta2=beta*beta
    beta4=beta2*beta2
    beta6=beta2*beta4
    beta8=beta2*beta6
    x2=x*x
    q=x/sqrt(beta2*(x2-1.0_r8)-1.0_r8)
    q2=q*q
    q4=q2*q2
    q6=q4*q2
    q8=q6*q2
    q10=q8*q2
    q12=q10*q2
    mu2=mu*mu
    mu3=mu2*mu
    mu4=mu3*mu   
    be2p1=1.0_r8+beta2
    be2p12=be2p1*be2p1
    be2p13=be2p12*be2p1
    be2p14=be2p13*be2p1
    uk(0)=1.0_r8;
    vk(0)=0.0_r8;
    uk(1)=(1.0_r8/24.0_r8)*beta2/be2p1;
    vk(1)=-(1.0_r8/24.0_r8)*q*(5.0_r8*beta2*q2+3.0_r8*beta2-3.0_r8)/be2p1;
    uk(2)=-(1.0_r8/1152.0_r8)*(385.0_r8*beta4*q6+81.0_r8*beta4*q2+462.0_r8*beta4*q4-beta4&
             -72*beta2-522*beta2*q2-462*beta2*q4+72+81*q2)/be2p12;
    vk(2)=-(1.0_r8/576.0_r8)*beta2*q*(5*beta2*q2+3*beta2-3.0_r8)/be2p12;
    uk(3)=-(1.0_r8/414720.0_r8)*beta2*(5775.0_r8*beta4*q6+1215.0_r8*beta4*q2+&
          1003.0_r8*beta4+6930.0_r8*beta4*q4-7830.0_r8*beta2*q2+1944.0_r8*beta2&
          -6930.0_r8*beta2*q4+1215.0_r8*q2+7128.0_r8)/be2p13;
    vk(3)=(1.0_r8/414720.0_r8)*q*(-45.0_r8*beta6+765765.0_r8*beta6*q6+425425.0_r8*beta6*q8+&
          30300.0_r8*beta6*q2+369603.0_r8*beta6*q4-564165.0_r8*beta4*q2-1300806.0_r8*&      
          beta4*q4-32355.0_r8*beta4-765765.0_r8*beta4*q6+369603.0_r8*beta2*q4+194400.0_r8*beta2+&
          564165.0_r8*beta2*q2-32400.0_r8-30375.0_r8*q2)/be2p13;
    uk(4)= (1.0_r8/39813120.0_r8)*(1632960.0_r8-9486720.0_r8*beta2-93486960.0_r8*beta2*q2+&
         93484530.0_r8*beta4*q2+6123600.0_r8*q2-178143300.0_r8*beta2*q4+490305150.0_r8*beta4*q4&
         +748417428.0_r8*beta4*q6-94121676.0_r8*beta2*q6+1606608.0_r8*beta4+4465125.0_r8*q4&
         +349922430.0_r8*beta4*q8-9936.0_r8*beta6-6107940.0_r8*beta6*q2-178129440.0_r8*beta6*q4&
         -748417428.0_r8*beta6*q6-1022624460.0_r8*beta6*q8-446185740.0_r8*beta6*q10&
         -4027.0_r8*beta8-2430.0_r8*beta8*q2+4451265.0_r8*beta8*q4+94110126.0_r8*beta8*q6&
         +349922430.0_r8*beta8*q8+446185740.0_r8*beta8*q10+185910725.0_r8*beta8*q12)/be2p14;
    vk(4)=(1.0_r8/9953280.0_r8)*beta2*q*(765765.0_r8*beta6*q6+425425.0_r8*beta6*q8+&
          3009.0_r8*beta6+35390.0_r8*beta6*q2+369603.0_r8*beta6*q4-549045.0_r8*beta4*q2-&
          1300806.0_r8*beta4*q4-26337.0_r8*beta4-765765.0_r8*beta4*q6+594405.0_r8*beta2*q2+&
          369603.0_r8*beta2*q4+203472.0_r8*beta2-50544.0_r8-30375.0_r8*q2)/be2p14;
    xc=sqrt(be2p1)/beta;
    chxi= x/xc; xi= log(chxi+sqrt(chxi*chxi-1.0_r8)); 
    shxi= sqrt(chxi*chxi-1.0_r8);
    U=uk(0)+uk(1)/mu+uk(2)/mu2+uk(3)/mu3+uk(4)/mu4
    V=vk(0)+vk(1)/mu+vk(2)/mu2+vk(3)/mu3+vk(4)/mu4
    chi= mu*(beta*xi-atan(1.0_r8/q))-pi/4.0_r8;
    S=cos(chi)*U-sin(chi)*V;
    argu2=sqrt(be2p1)
    fac2=1.0_r8
    DO n=1,mu
      fac2=argu2*fac2
    ENDDO
    fac2=2.0_r8*sqrt(q/(x*mu))*fac2
    pmmu=exp(-tau*(pi-atan(1.0_r8/beta)))*gammah(mu)*fac2*cosh(pi*tau)/pi*S
    END SUBROUTINE expaosc

    SUBROUTINE expakia(x,mu,tau,pmmu,ierr)
    ! --------------------------------------------------
    ! Calculation of P^(mu)_(-1/2+i*tau)(x) by using 
    ! asymptotic expansions in terms of the Kia functions 
    ! functions
    ! --------------------------------------------------
    ! Inputs:
    !   x ,    argument of the functions
    !   mu,    function parameter (integer value)
    !   tau,   function parameter (real value)
    ! Outputs:
    !   pmmu , P^(mu)_(-1/2+i*tau)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! ---------------------------------------------------------
    USE Someconstants
    REAL(r8), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: tau
    REAL(r8), INTENT(OUT) :: pmmu
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: conv,toln,x2,tau2,mu2,mu3,mu4,beta,beta2,beta3,beta4,beta5,beta6,beta7,beta8,&
                beta9,beta10,beta11,beta12,beta13,beta14,beta15,beta16,be2p1,be2p12,be2p13,&
                be2p14,zc,p,sqq,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,rhsi,f,fp,zeta,zeta2,znew,&
                lambda,boz,omega,c,W,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,W14,W16,W18,W20,&
                c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,&
                c18,c19,c20,c21,cdenom,cdenom2,cdenom3,cdenom4,cdenom5,cdenom6,&
                cdenom7,cdenom8,cdenom9,dki,dkid,A0,Ac,Bc,fac2
    REAL(r8), DIMENSION(0:5) :: ak, bk
    INTEGER ierro, j, iter, n
    conv=1.0_r8
    toln=epss
    ak(4)=0;
    bk(4)=0;
    x2=x*x
    ierr=0
    tau2=tau*tau
    beta=tau/mu
    mu2=mu*mu
    mu3=mu2*mu
    mu4=mu3*mu
    beta2=beta*beta
    beta3=beta2*beta
    beta4=beta3*beta
    beta5=beta4*beta
    beta6=beta5*beta
    beta7=beta6*beta
    beta8=beta7*beta
    beta9=beta8*beta
    beta10=beta9*beta
    beta11=beta10*beta
    beta12=beta11*beta
    beta13=beta12*beta
    beta14=beta13*beta
    beta15=beta14*beta
    beta16=beta15*beta
    be2p1=1.0_r8+beta2
    be2p12=be2p1*be2p1
    be2p13=be2p12*be2p1
    be2p14=be2p13*be2p1
    lambda=0.5_r8*(log((x2-1.0_r8)/be2p1)+beta*acos((1.0_r8-beta2)/be2p1));
    DO j=0,4
      ak(j)=0; bk(j)=0;
    ENDDO
    ak(0)=1.0_r8
    ak(1)=0.041666666666666666666666666666666666666666666666666666666666666_r8*beta2/be2p1
    zc=sqrt(be2p1)/beta
    IF (x<zc) THEN
      p=x/sqrt(1.0_r8+beta2*(1.0_r8-x2));
    ELSE 
      p=x/sqrt(beta2*(x2-1.0_r8)-1.0_r8); 
    ENDIF
    sqq=sqrt((x-1.0_r8)*(x+1.0_r8))
    p2=p*p
    p3=p2*p
    p4=p3*p
    p5=p4*p
    p6=p5*p
    p7=p6*p
    p8=p7*p
    p9=p8*p
    p10=p9*p
    p11=p10*p
    p12=p11*p
! Determination of zeta: Newton method
    IF (x<zc) THEN
      rhsi=log((1.0_r8+p)/(p-1.0_r8))-beta*acos((beta2*p2-1.0_r8)/(beta2*p2+1.0_r8))
      omega=rhsi*0.5_r8/beta
      IF (omega>1.0_r8) THEN
        zeta=beta*(omega+pi*0.5_r8-0.5_r8/omega+pi*0.25_r8/(omega*omega))
      ELSE
        zeta=beta*(1.0_r8+xpowy(1.5_r8*omega/sqrt(2.0_r8),2.0_r8/3.0_r8))
      ENDIF      
      zeta2=zeta*zeta
      iter=0
      DO WHILE ((conv>toln).AND.(iter<100))
        iter=iter+1
        f=2.0_r8*(sqrt(zeta2-beta2)-beta*acos(beta/zeta))-rhsi;
        fp=2.0_r8*sqrt(1.0_r8-beta2/zeta2)
        znew=zeta-f/fp
        conv=abs(znew-zeta)
        zeta=znew
        zeta2=zeta*zeta
      ENDDO
      W=sqrt(zeta2-beta2)
      W2=W*W
      W3=W2*W
      W4=W2*W2
      W5=W4*W 
      W6=W4*W2
      W7=W6*W
      W8=W6*W2
      W9=W8*W
      W10=W8*W2
      W11=W10*W
      W12=W10*W2
      W13=W12*W
      W14=W12*W2
      W16=W14*W2
      W18=W16*W2
      W20=W18*W2
      IF (W<0.2_r8) THEN
        c=xpowy(1.0_r8/(1.0_r8+beta2),1.0_r8/3.0_r8)
        c2=c*c; c3=c2*c; c4=c3*c; c5=c4*c; c6=c5*c; c7=c6*c; c8=c7*c;
        c9=c8*c; c10=c9*c; c11=c10*c; c12=c11*c; c13=c12*c; c14=c13*c;
        c15=c14*c; c16=c15*c; c17=c16*c; c18=c17*c; c19=c18*c; c20=c19*c;
        c21=c20*c
        cdenom=c2+c+1.0_r8
        cdenom2=cdenom*cdenom
        cdenom3=cdenom2*cdenom
        cdenom4=cdenom3*cdenom         
        cdenom5=cdenom4*cdenom
        cdenom6=cdenom5*cdenom
        cdenom7=cdenom6*cdenom  
        cdenom8=cdenom7*cdenom
        cdenom9=cdenom8*cdenom  
        ak(2)=(-1.0_r8/201600.0_r8*(1393.0_r8*c8+1393*c7+97*c6-1730*c5&
             -1730*c4-434*c3+145*c2+145*c+721.0_r8)/cdenom*W6+1.0_r8/11088000.0_r8&
             *(61641.0_r8*c10+123282.0_r8*c9+167675.0_r8*c8+88786.0_r8*c7&
             +17025.0_r8*c6-28864.0_r8*c5-17504.0_r8*c4-13272.0_r8*c3&
             -14128.0_r8*c2-10592.0_r8*c-10224.0_r8)*c2/cdenom2*W8-1.0_r8/2522520000.0_r8&
             *(7150374.0_r8*c12+21451122.0_r8*c11+36020889.0_r8*c10+32983740.0_r8*c9&
              +14301635.0_r8*c8-6262716.0_r8*c7-14919053.0_r8*c6-14610316.0_r8*c5&
              -11834100.0_r8*c4-9657260.0_r8*c3-6688048.0_r8*c2-3310224.0_r8*c-707168.0_r8)*c4&
              /cdenom3*W10+1.0_r8/17463600000.0_r8*(22888390.0_r8*c14+91553560.0_r8*c13&
              +186272716.0_r8*c12+218657894.0_r8*c11+132101885.0_r8*c10-23474560.0_r8*c9&
              -135459915.0_r8*c8-154804164.0_r8*c7-115901691.0_r8*c6-68398480.0_r8*c5&
              -31760320.0_r8*c4-6014080.0_r8*c3+7134928.0_r8*c2+7581952.0_r8*c+2937760.0_r8)*c6&
              /cdenom4*W12-1.0_r8/1500899400000.0_r8*(860675886.0_r8*c16+4303379430.0_r8*c15&
              +10459577580.0_r8*c14+14858486859.0_r8*c13+11534759793.0_r8*c12+380167041.0_r8*c11&
              -11018380556.0_r8*c10-15148048045.0_r8*c9-11672657236.0_r8*c8-5518835692.0_r8*c7&
              -755351417.0_r8*c6+1624555844.0_r8*c5+2194550780.0_r8*c4+1646476012.0_r8*c3+&
               723229184.0_r8*c2+132301600.0_r8*c-36402688.0_r8)*c8/cdenom5*W14+1.0_r8&
              /141159588570000000.0_r8*(34172732883195.0_r8*c18+205036397299170.0_r8*c17&
              +585142109710263.0_r8*c16+982070424905148.0_r8*c15+943000129848528.0_r8*c14&
              +222975520459740.0_r8*c13-789007393738082.0_r8*c12-1362631504977174.0_r8*c11&
              -1163145597443907.0_r8*c10-509586706230712.0_r8*c9+74713644479970.0_r8*c8+&
              339651442452438.0_r8*c7+337238315362017.0_r8*c6+213724589653944.0_r8*c5+&
              80569717917264.0_r8*c4-945754950760.0_r8*c3-19890984887376.0_r8*c2-10453159775136.0_r8*c&
              -807559327280.0_r8)*c10/cdenom6*W16)/W6;
        ak(3)=(1.0_r8/14515200.0_r8*(c-1.0_r8)*(145649.0_r8*c8+145649.0_r8*c7+141761.0_r8*c6&
              +29390.0_r8*c5+29390.0_r8*c4+33278.0_r8*c3+36065.0_r8*c2+36065.0_r8*c+37793.0_r8)*W6&
              -1.0_r8/266112000.0_r8*(61641.0_r8*c10+123282.0_r8*c9+167675.0_r8*c8+88786.0_r8*c7+&
               17025.0_r8*c6-28864.0_r8*c5-17504.0_r8*c4-13272.0_r8*c3-14128.0_r8*c2-10592.0_r8*c-10224.0_r8)&
               *(c-1.0_r8)*c2/cdenom*W8+1.0_r8/60540480000.0_r8*(7150374.0_r8*c12+21451122.0_r8*c11&
               +36020889.0_r8*c10+32983740.0_r8*c9+14301635.0_r8*c8-6262716.0_r8*c7-14919053.0_r8*c6&
               -14610316.0_r8*c5-11834100.0_r8*c4-9657260.0_r8*c3-6688048.0_r8*c2-3310224.0_r8*c-707168.0_r8)&
               *(c-1.0_r8)*c4/cdenom2*W10-1.0_r8/419126400000.0_r8*(22888390.0_r8*c14+91553560.0_r8*c13&
               +186272716.0_r8*c12+218657894.0_r8*c11+132101885.0_r8*c10-23474560.0_r8*c9-&
               135459915.0_r8*c8-154804164.0_r8*c7-115901691.0_r8*c6-68398480.0_r8*c5-31760320.0_r8*c4&
               -6014080.0_r8*c3+7134928.0_r8*c2+7581952.0_r8*c+2937760.0_r8)*(c-1.0_r8)*c6/cdenom3*W12&
               +1.0_r8/36021585600000.0_r8*(860675886.0_r8*c16+4303379430.0_r8*c15+10459577580.0_r8*c14&
               +14858486859.0_r8*c13+11534759793.0_r8*c12+380167041.0_r8*c11-11018380556.0_r8*c10&
               -15148048045.0_r8*c9-11672657236.0_r8*c8-5518835692.0_r8*c7-755351417.0_r8*c6+&
                1624555844.0_r8*c5+2194550780.0_r8*c4+1646476012.0_r8*c3+723229184.0_r8*c2+132301600.0_r8*c&
               -36402688.0_r8)*(c-1.0_r8)*c8/cdenom4*W14-1.0_r8/3387830125680000000.0_r8*(34172732883195.0_r8*c18&
               +205036397299170.0_r8*c17+585142109710263.0_r8*c16+982070424905148.0_r8*c15+943000129848528.0_r8*c14&
               +222975520459740.0_r8*c13-789007393738082.0_r8*c12-1362631504977174.0_r8*c11-1163145597443907.0_r8*c10&
               -509586706230712.0_r8*c9+74713644479970.0_r8*c8+339651442452438.0_r8*c7+337238315362017.0_r8*c6&
               +213724589653944.0_r8*c5+80569717917264.0_r8*c4-945754950760.0_r8*c3-19890984887376.0_r8*c2&
               -10453159775136.0_r8*c-807559327280.0_r8)*(c-1.0_r8)*c10/cdenom5*W16)/W6;
        bk(1)=(-1.0_r8/280.0_r8*c*zeta*(9*c5+9*c4+9*c3+4*c+4)/cdenom*W4&
              +1.0_r8/12600.0_r8*(98*c7+196*c6+213*c5+83*c4-47*c3-96*c2&
              -76*c-56)*zeta*c3/cdenom2*W6-1.0_r8/1617000.0_r8*&
              (4077*c9+12231*c8+16916*c7+9978*c6-3237*c5-11410*c4&
              -10413*c3-5592*c2-1828*c+828)*zeta*c5/cdenom3*W8+&
              1.0_r8/315315000_r8*(271656.0_r8*c11+1086624.0_r8*c10+1921545*c9+&
              1573920.0_r8*c8-107700.0_r8*c7-1684734.0_r8*c6-1875571.0_r8*c5-975175*c4-83470*c3&
              +273160*c2+252628.0_r8*c+22792)*zeta*c7/cdenom4*W10-1.0_r8/33108075000.0_r8*&
              (9926995.0_r8*c13+49634975.0_r8*c12+108971493.0_r8*c11+117996680.0_r8*c10+&
               23304685.0_r8*c9-112142745.0_r8*c8&
              -164160315.0_r8*c7-99982062.0_r8*c6-2998230*c5+44249930.0_r8*c4+37346320.0_r8*c3+13501080.0_r8*c2&
              -2899976.0_r8*c-1250080.0_r8)*zeta*c9/cdenom5*W12+1.0_r8/938062125000.0_r8*&
               (98709912.0_r8*c15+592259472.0_r8*c14+1566631407.0_r8*c13+2144609910.0_r8*c12+993484611.0_r8*c11&
              -1648415217.0_r8*c10-3459671744.0_r8*c9-2675618409.0_r8*c8-303627696.0_r8*c7+1332420910.0_r8*c6&
              +1322416896.0_r8*c5+514721244.0_r8*c4-73283846.0_r8*c3-186087792.0_r8*c2-37851504.0_r8*c+&
              4317096)*zeta*c11/cdenom6*W14-1.0_r8/61757319999375000.0_r8*(2290875403815.0_r8*c17&
              +16036127826705.0_r8*c16+49847663782476.0_r8*c15+83192099634792.0_r8*c14+60433694450568.0_r8*c13&
              -42447122909352.0_r8*c12-150313744710496.0_r8*c11-151930198189654.0_r8*c10&
              -40647228995323.0_r8*c9+71949889801323.0_r8*c8+94022483075388.0_r8*c7+44189908944534.0_r8*c6&
              -5389203336702.0_r8*c5-19139960492979.0_r8*c4-9728875988691.0_r8*c3+299051082664.0_r8*c2&
              +984902727172.0_r8*c+51375356260.0_r8)*zeta*c13/cdenom7*W16+&
               1.0_r8/308786599996875000.0_r8*(4048556587768.0_r8*c19&
               +32388452702144.0_r8*c18+115966656910053.0_r8*c17+229357243890444.0_r8*c16&
               +226393507645548.0_r8*c15-38068980041772.0_r8*c14-439693813180464.0_r8*c13&
               -587811603844524.0_r8*c12-270394076368940.0_r8*c11+226230760410240.0_r8*c10&
               +439866482466809.0_r8*c9+268458673161834.0_r8*c8-5466610953768.0_r8*c7-122882117855892.0_r8*c6&
               -81933726464274.0_r8*c5-10976641917615.0_r8*c4+15049928649936.0_r8*c3+5555561805944.0_r8*c2&
               -336943037932.0_r8*c-125602833664.0_r8)*zeta*c15/cdenom8*W18&
               -1.0_r8/3231451768967296875000.0_r8*(15000147517319766.0_r8*c21+&
                135001327655877894.0_r8*c20+547902204206413110.0_r8*c19+1256083696114375320.0_r8*c18&
               +1564644811187692800.0_r8*c17+338396649780494412.0_r8*c16-2373267185032712572.0_r8*c15&
               -4279605957210423960.0_r8*c14-2915791867689908985.0_r8*c13+973459027265184775.0_r8*c12&
               +3730399203798277107.0_r8*c11+3042903732882256788.0_r8*c10+396117195396233010.0_r8*c9&
               -1328396918245677270.0_r8*c8-1187627541412476150.0_r8*c7-291199621343628524.0_r8*c6+&
                211393950537241899.0_r8*c5+181498437006912315.0_r8*c4+19874016923212235.0_r8*c3&
               -16503501533468100.0_r8*c2-2542411653780036.0_r8*c+155654367677916.0_r8)*zeta*c17/cdenom9*W20)&
                /W4       
        bk(2)=(1.0_r8/6720.0_r8*c*(9*c5+9*c4+9*c3+4*c+4.0_r8)*(c-1.0_r8)*zeta*W4&
                 -1.0_r8/302400.0_r8*(98*c7+196*c6+213*c5+83*c4-47*c3-96*c2-76*c-56.0_r8)*c3*(c-1.0_r8)&
                 *zeta/cdenom*W6+1.0_r8/38808000.0_r8*(4077*c9+12231.0_r8*c8+16916.0_r8*c7&
                  +9978.0_r8*c6-3237*c5-11410.0_r8*c4-10413.0_r8*c3-5592*c2&
                  -1828.0_r8*c+828)*c5*(c-1.0_r8)*zeta/cdenom2*W8-1.0_r8/7567560000.0_r8*&
                  (271656.0_r8*c11+1086624.0_r8*c10+1921545.0_r8*c9+1573920.0_r8*c8-107700.0_r8*&
                   c7-1684734.0_r8*c6-1875571.0_r8*c5-975175.0_r8*c4-83470.0_r8*c3+273160.0_r8*c2&
                   +252628.0_r8*c+22792.0_r8)*c7*(c-1.0_r8)*zeta/cdenom3*W10+&
                  1.0_r8/794593800000.0_r8*(9926995.0_r8*c13+49634975.0_r8*c12&
                   +108971493.0_r8*c11+117996680.0_r8*c10+23304685.0_r8*c9&
                   -112142745.0_r8*c8-164160315.0_r8*c7-99982062.0_r8*c6&
                   -2998230.0_r8*c5+44249930.0_r8*c4+37346320.0_r8*c3+13501080.0_r8*c2&
                   -2899976.0_r8*c-1250080.0_r8)*c9*(c-1.0_r8)*zeta/cdenom4*W12&
                   -1.0_r8/22513491000000.0_r8*(98709912.0_r8*c15+592259472.0_r8*c14+&
                   1566631407.0_r8*c13+2144609910.0_r8*c12+993484611.0_r8*c11&
                   -1648415217.0_r8*c10-3459671744.0_r8*c9-2675618409.0_r8*c8-303627696.0_r8*c7&
                   +1332420910.0_r8*c6+1322416896.0_r8*c5+514721244.0_r8*c4-73283846.0_r8*c3&
                   -186087792.0_r8*c2-37851504.0_r8*c+4317096.0_r8)*c11*(c-1.0_r8)*zeta/cdenom5*W14&
                   +1.0_r8/1482175679985000000.0_r8*(2290875403815.0_r8*c17+&
                   16036127826705.0_r8*c16+49847663782476.0_r8*c15+83192099634792.0_r8*c14&
                   +60433694450568.0_r8*c13-42447122909352.0_r8*c12-150313744710496.0_r8*c11&
                   -151930198189654.0_r8*c10-40647228995323.0_r8*c9+71949889801323.0_r8*c8+&
                    94022483075388.0_r8*c7+44189908944534.0_r8*c6-5389203336702.0_r8*c5-19139960492979.0_r8*c4&
                   -9728875988691.0_r8*c3+299051082664.0_r8*c2+984902727172.0_r8*c+51375356260.0_r8)*c13*&
                    (c-1.0_r8)*zeta/cdenom6*W16-1.0_r8/7410878399925000000.0_r8*(4048556587768.0_r8*c19&
                   +32388452702144.0_r8*c18+115966656910053.0_r8*c17+229357243890444.0_r8*c16&
                   +226393507645548.0_r8*c15-38068980041772.0_r8*c14-439693813180464.0_r8*c13&
                   -587811603844524.0_r8*c12-270394076368940.0_r8*c11+226230760410240.0_r8*c10&
                   +439866482466809.0_r8*c9+268458673161834.0_r8*c8-5466610953768.0_r8*c7&
                   -122882117855892.0_r8*c6-81933726464274.0_r8*c5-10976641917615.0_r8*c4&
                   +15049928649936.0_r8*c3+5555561805944.0_r8*c2-336943037932.0_r8*c-&
                    125602833664.0_r8)*c15*(c-1.0_r8)*zeta/cdenom7*W18+1.0_r8/77554842455215125000000.0_r8&
                    *(15000147517319766.0_r8*c21+135001327655877894.0_r8*c20+547902204206413110.0_r8*c19&
                    +1256083696114375320.0_r8*c18+1564644811187692800.0_r8*c17+338396649780494412.0_r8*c16&
                    -2373267185032712572.0_r8*c15-4279605957210423960.0_r8*c14-2915791867689908985.0_r8*c13&
                    +973459027265184775.0_r8*c12+3730399203798277107.0_r8*c11+3042903732882256788.0_r8*c10&
                    +396117195396233010.0_r8*c9-1328396918245677270.0_r8*c8-1187627541412476150.0_r8*c7&
                    -291199621343628524.0_r8*c6+211393950537241899.0_r8*c5+181498437006912315.0_r8*c4&
                    +19874016923212235.0_r8*c3-16503501533468100.0_r8*c2-2542411653780036.0_r8*c&
                   +155654367677916.0_r8)*c17*(c-1.0_r8)*zeta/cdenom8*W20)/W4;
        bk(3)=(1.0_r8/524160000.0_r8*c*zeta*(6059583.0_r8*c13+12119166.0_r8*c12+18178749.0_r8*c11&
               +9111312.0_r8*c10+95615.0_r8*c9-8920082.0_r8*c8-5698629.0_r8*c7-2545036.0_r8*c6+&
               608557.0_r8*c5+546510.0_r8*c4+1181699.0_r8*c3+1816888.0_r8*c2+1229112.0_r8*c+614556.0_r8)&
               /cdenom2*W10-1.0_r8/1816214400000.0_r8*(9727357705.0_r8*c15+29182073115.0_r8*c14&
               +54164855211.0_r8*c13+55512024253.0_r8*c12+33367820181.0_r8*c11-1784732145.0_r8*c10&
               -19542687620.0_r8*c9-20194526124.0_r8*c8-12274421967.0_r8*c7-8011196494.0_r8*c6-7270887045.0_r8*c5&
               -7771225740.0_r8*c4-7435063897.0_r8*c3-6252124236.0_r8*c2-3350967852.0_r8*c-1258951720.0_r8)*&
               zeta*c3/cdenom3*W12+1.0_r8/51459408000000.0_r8*(144811913265.0_r8*c17+579247653060.0_r8*c16&
               +1282754051665.0_r8*c15+1728105185620.0_r8*c14+1451348623879.0_r8*c13+484314490976.0_r8*c12&
               -489620336905.0_r8*c11-945041189740.0_r8*c10-873657521970.0_r8*c9-631216670541.0_r8*c8&
               -445274929539.0_r8*c7-335069642665.0_r8*c6-245861617660.0_r8*c5-156878390545.0_r8*c4&
               -77977426823.0_r8*c3-17492885312.0_r8*c2+3279403240.0_r8*c+6170395620.0_r8)*zeta*c5/cdenom4*W14&
               -1.0_r8/20532303792000000.0_r8*(29048679511047.0_r8*c19+145243397555235.0_r8*c18+&
                377950239272970.0_r8*c17+611632330911882.0_r8*c16+633363321305990.0_r8*c15+334989811415962.0_r8*c14&
               -124793259618935.0_r8*c13-462893943016455.0_r8*c12-532059257514943.0_r8*c11-&
                409102997160979.0_r8*c10-246328251031409.0_r8*c9-125738641603410.0_r8*c8-49813792932765.0_r8*c7&
               -1352246462804.0_r8*c6+26758233574976.0_r8*c5+35762086322945.0_r8*c4+30037846475850.0_r8*c3&
               +16475189210860.0_r8*c2+5900059068820.0_r8*c+687614243288.0_r8)*zeta*c7/cdenom5*W16)/W10;
        IF (abs(zeta-beta)<1.e-10_r8) THEN
          A0=beta2*c+2.0_r8*beta*(3.0_r8+2.0_r8*beta2+2.0_r8/c2)*(zeta-beta)&
             /5.0_r8/be2p1
        ELSE 
          A0=xpowy((zeta2-beta2)/(1.0_r8+beta2*(1.0_r8-x2)),0.25_r8)
        ENDIF
      ELSE 
         c=xpowy(1.0_r8/be2p1,1.0_r8/3.0_r8);
         ak(2)=1.0_r8/1152.0_r8*(-594*W2*beta6+W6*beta4-594*W2*beta2+42*W3*p*beta2&
              -270*W4*beta2-455*beta8+72*W6*beta2+81*W6*p2-72*W6-54*W5*p*beta4+54*W5*p&
              +70*W3*beta6*p3+70*W3*beta4*p3-910*beta6-455*beta4-135*W4-42*W3*beta6*p&
              -1188*W2*beta4+385*W6*beta4*p6-462.0_r8*W6*beta4*p4+462*W6*beta2*p4&
              -522*W6*beta2*p2+90*W5*p3*beta4+90*W5*p3*beta2-135*W4*beta4+81*W6*beta4*p2)&
               /(W6*be2p12);
        ak(3)=1.0_r8/414720.0_r8*(-6930.0_r8*W6*beta4*p4-7128*W6+1215*W6*p2&
              -1944*W6*beta2-2025.0_r8*W4*beta4-13650.0_r8*beta6-6825*beta4+1215.0_r8*&
               W6*beta4*p2-4050.0_r8*W4*beta2+6930*W6*beta2*p4-1003.0_r8*W6*beta4+&
               810*W5*p+1350*W5*p3*beta2+1350*W5*p3*beta4-7830*W6*beta2*p2+&
               5775*W6*beta4*p6+1050*beta6*W3*p3-8910*W2*beta6-630*W3*beta6*p&
              -17820*W2*beta4+1050*W3*p3*beta4-6825*beta8+630*W3*p*beta2-&
               8910*W2*beta2-810*W5*p*beta4-2025*W4)*beta2/(W6*(beta6+1.0_r8+3*beta2+3*beta4));
         ak(4)=1.0_r8/39813120.0_r8*(4451265.0_r8*beta8*W12*p4+2430.0_r8*beta8*&
               W12*p2-748417428.0_r8*beta4*W12*p6+490305150.0_r8*beta4*W12*p4&
               -1022624460.0_r8*beta6*W12*p8-5740875.0_r8*W8+6107940.0_r8*beta6*W12*p2&
               -93484530.0_r8*beta4*W12*p2-178129440.0_r8*beta6*W12*p4+748417428.0_r8*&
                beta6*W12*p6-94110126.0_r8*beta8*W12*p6+349922430.0_r8*beta8*W12*p8&
               +349922430.0_r8*beta4*W12*p8+1632960.0_r8*W12+4914000.0_r8*W6*beta10*p2&
               -202076875.0_r8*beta8+11411400.0_r8*beta8*W3*p-10510500.0_r8*W6*beta10*p6&
               -808307500.0_r8*beta10+1871100.0_r8*beta8*W10*p4-328050.0_r8*beta8*W10*p2&
               -178143300.0_r8*W12*p4*beta2-111234708.0_r8*beta2*W6-1212461250.0_r8*beta12&
               -808307500.0_r8*beta14-202076875.0_r8*beta16-295650.0_r8*W10*beta4-299700.0_r8*W10&
               *beta6+11911900.0_r8*beta8*W9*p9-5255250.0_r8*W6*beta12*p6+6415200.0_r8*W8*beta8*p2&
               +28528500.0_r8*beta10*W3*p3+93486960.0_r8*beta2*W12*p2-6123600.0_r8*W12*p2&
               -1443420.0_r8*W8*beta10*p2+291600.0_r8*W10+113400.0_r8*beta2*W9*p&
               -1620.0_r8*W11*beta8*p+15315300.0_r8*W11*beta8*p9-27567540.0_r8*W11*beta8*p7&
               -493152660.0_r8*beta6*W2-443956032.0_r8*beta4*W6-396578754.0_r8*beta4*W4&
               -21680460.0_r8*beta2*W8+13305708.0_r8*W11*beta8*p5-6306300.0_r8*W6*beta8*p4&
               -1090800.0_r8*W11*beta8*p3+9509500.0_r8*beta8*W3*p3+5420844.0_r8*beta2*W7*p&
               -33162210.0_r8*beta4*W8-1559250.0_r8*W10*beta8*p6+28528500.0_r8*beta12*W3*p3&
               -4050.0_r8*beta8*W10-4027.0_r8*beta8*W12-9936.0_r8*beta6*W12+185910725.0_r8*beta8&
               *W12*p12-666425448.0_r8*beta6*W6-1586315016.0_r8*beta6*W4-1972610640.0_r8*beta8*W2&
               -446185740.0_r8*beta8*W12*p10-9486720.0_r8*beta2*W12+291600.0_r8*W10*beta2+&
                6306300.0_r8*W6*beta10*p4+4465125.0_r8*W12*p4-5255250.0_r8*W6*beta8*p6&
               +1458000.0_r8*W10*p2*beta2+5705700.0_r8*beta6*W3*p-19216440.0_r8*W11*beta2*p3&
               -1105650.0_r8*W6*beta12*p2+1606608.0_r8*beta4*W12+94121676.0_r8*beta2*W12*p6&
               -445935282.0_r8*W6*beta8-112244808.0_r8*W6*beta10-13650.0_r8*W6*beta12&
               -2379472524.0_r8*W4*beta8-1586315016.0_r8*W4*beta10-396578754.0_r8*W4*beta12&
               -2958915960.0_r8*beta10*W2-1972610640.0_r8*beta12*W2-493152660.0_r8*beta14*W2&
               -6306300.0_r8*beta6*W6*p4+4914000.0_r8*beta6*W6*p2-1105650.0_r8*beta4*W6*p2&
               +12039300.0_r8*W6*beta8*p2+6306300.0_r8*W6*beta12*p4+10602900.0_r8*W5*p*beta4&
               +21205800.0_r8*W5*beta6*p-21205800.0_r8*W5*p*beta10-5705700.0_r8*beta14*W3*p+&
                17671500.0_r8*beta6*W5*p3+53014500.0_r8*beta8*W5*p3+53014500.0_r8*beta10*W5*p3&
               -11411400.0_r8*beta12*W3*p+9509500.0_r8*beta14*W3*p3+17671500.0_r8*beta12*W5*p3&
               -10602900.0_r8*beta12*W5*p+13305708.0_r8*beta2*W11*p5-1559250.0_r8*W10*beta4*p6&
               +4536000.0_r8*beta4*W9*p+3516660.0_r8*beta6*W9*p-1417500.0_r8*beta8*W9*p&
               -1260.0_r8*beta10*W9*p+8232840.0_r8*W8*beta8*p4+8232840.0_r8*W8*beta10*p4&
               -24264360.0_r8*W8*beta6-7059555.0_r8*W8*beta8-17820.0_r8*W8*beta10+10841688.0_r8*beta4*W7*p&
               -10841688.0_r8*beta8*W7*p-5420844.0_r8*beta10*W7*p+11911900.0_r8*beta10*W9*p9&
               -21441420.0_r8*beta10*W9*p7+10348884.0_r8*beta10*W9*p5-848400.0_r8*beta10*W9*p3&
               -26073684.0_r8*beta8*W9*p5+15798720.0_r8*beta8*W9*p3+21441420.0_r8*beta6*W9*p7&
               -1871100.0_r8*W10*beta4*p4-26073684.0_r8*beta6*W9*p5+2551500.0_r8*beta6*W9*p3+&
                6415200.0_r8*W8*p2*beta4+15717240.0_r8*W8*p2*beta6+510300.0_r8*W9*p&
               -1443420.0_r8*W8*p2*beta2+1701000.0_r8*W9*beta2*p3-12394620.0_r8*W9*beta4*p3&
               -1871100.0_r8*W10*beta2*p4+5832000.0_r8*W11*p*beta2-1166400.0_r8*W11*beta6*p+&
                19219140.0_r8*W11*beta6*p3+5833620.0_r8*W11*p*beta4+9034740.0_r8*beta4*W7*p3+27104220.0_r8*&
                beta8*W7*p3+27104220.0_r8*beta6*W7*p3+9034740.0_r8*beta10*W7*p3-1166400.0_r8*W11*p+&
                1093500.0_r8*W11*p3-8232840.0_r8*beta6*W8*p4-8232840.0_r8*beta4*W8*p4+&
                10348884.0*beta4*W9*p5-6860700.0_r8*W8*beta10*p6-6860700.0_r8*W8*beta6*p6&
               -13721400.0_r8*W8*beta8*p6+27567540.0_r8*W11*beta4*p7+446185740.0_r8*beta6*W12*p10&
               -33523308.0_r8*W11*beta4*p5-33523308.0_r8*W11*beta6*p5-328050.0_r8*W10*p2&
               +1871100.0_r8*W10*beta6*p4-3118500.0_r8*W10*beta6*p6+15315300.0_r8*W11*beta6*p9&
               +3572100.0_r8*W10*p2*beta4+1458000.0_r8*W10*p2*beta6)/(W12*be2p14);
         bk(1)=-1.0_r8/24.0_r8*(-5*beta4-3*W3*p*beta2-3*W2*beta2+5*W3*p3*beta2&
               -5*beta2-3*W2+3*W3*p)*zeta/(W4*be2p1);
         bk(2)=-1.0_r8/576.0_r8*zeta*beta2*(-5*beta4-3*W3*p*beta2-3*W2*beta2+&
              5*W3*p3*beta2-5*beta2-3*W2+3*W3*p)/(W4*be2p12);
         bk(3)=-1.0_r8/414720.0_r8*(-20790.0_r8*W5*p*beta8-34650.0_r8*W6*beta4*p4&
               -30375.0_r8*W6-32355.0_r8*W9*beta4*p-1276275.0_r8*beta10-2297295.0_r8*W2*beta8&
               -765765.0_r8*W2*beta10+194400.0_r8*beta2*W9*p-28875.0_r8*W6*beta6*p6&
               +3240*W8-20790.0_r8*W8*p4*beta2-3645*W7*p*beta4+6075*W7*p3*beta2+&
                12150*W7*p3*beta4-85725.0_r8*W6*beta2-30300.0_r8*W9*beta6*p3+369603.0_r8*W9*beta6*p5&
               +765765.0_r8*W9*beta4*p7+57750.0_r8*W3*p3*beta8-75*W6*beta8-1108809.0_r8*W4*beta4&
               -1108809*W4*beta6-20790.0_r8*W5*p*beta6-425425.0_r8*beta6+33075.0_r8*W6*beta4*p2&
               +28875.0_r8*W3*beta10*p3-6075*W6*beta8*p2-35850.0_r8*W6*beta6-17325.0_r8*W3*beta10*p&
               -369603.0_r8*W4*beta2-32400*W9*p+3645*W7*p*beta2+3645*W7*p+6075*W7*p3*beta6+&
                564165*W9*beta4*p3-91125.0_r8*W6*beta4-17325.0_r8*W3*beta8*p-17325.0_r8*W8*beta6*p6&
               -1300806.0_r8*W9*beta4*p5+33075*W6*beta6*p2-425425.0_r8*beta12+34650.0_r8*W5*p3*beta4&
               -28875.0_r8*W6*beta8*p6+34650.0_r8*W6*beta8*p4+20790.0_r8*W8*beta6*p4&
               -17325.0_r8*W8*beta4*p6-3645*W8*p2+425425.0_r8*W9*beta6*p9-6075*W6*beta2*p2&
               -369603.0_r8*W4*beta8-765765.0_r8*W9*beta6*p7+28875.0_r8*beta6*W3*p3-2297295.0_r8*&
                W2*beta6+17325*W3*beta6*p-765765.0_r8*W2*beta4-1276275.0_r8*beta8+17325.0_r8*W3*beta4*p&
               -45*W8*beta6-3645*W8*beta6*p2+30375.0_r8*W9*p3-3645.0_r8*W7*p*beta6+20790.0_r8*W5*p*beta4&
               +20790.0_r8*W5*p*beta2-564165.0_r8*W9*beta2*p3+369603.0_r8*W9*beta2*p5&
               +19845.0_r8*W8*beta2*p2+19845.0_r8*W8*beta4*p2-3285*W8*beta4+34650.0_r8*W5*p3*beta8&
               -45*W9*beta6*p+69300*W5*p3*beta6)*zeta/(be2p13*W10);
         bk(4)=-1.0_r8/9953280.0_r8*zeta*(21384.0_r8*W8-30375.0_r8*W6-425425.0_r8*beta6&
               -1276275.0_r8*beta8-17325.0_r8*beta8*W3*p+20790.0_r8*W5*p*beta2-&
                1276275.0_r8*beta10-55485.0_r8*beta2*W6-369603.0_r8*beta2*W4-425425.0_r8*beta12&
               -765765.0_r8*beta4*W2+28875.0_r8*beta10*W3*p3+425425.0_r8*beta6*W9*p9&
               +28875.0_r8*beta6*W3*p3+203472.0_r8*beta2*W9*p-2297295.0_r8*beta6*W2&
               -45765.0_r8*beta4*W6-1108809.0_r8*beta4*W4+27216.0_r8*beta2*W8+34650.0_r8*W6*beta8*p4&
               +57750.0_r8*beta8*W3*p3+3645.0_r8*beta2*W7*p+8841.0_r8*beta4*W8-28875.0_r8*&
                beta6*W6*p6-17325.0_r8*beta10*W3*p-15640.0_r8*beta6*W6-1108809.0_r8*beta6*W4&
               -2297295.0_r8*beta8*W2-28875.0_r8*W6*beta8*p6+17325.0_r8*beta6*W3*p+&
                17325.0_r8*beta4*W3*p+5015.0_r8*W6*beta8-369603.0_r8*W4*beta8-765765.0_r8*beta10*W2&
               +33075.0_r8*beta6*W6*p2-34650.0_r8*beta4*W6*p4+33075.0_r8*beta4*W6*p2&
               -6075.0_r8*W6*beta8*p2+20790.0_r8*W5*p*beta4-20790.0_r8*W5*beta6*p&
               -20790.0_r8*W5*p*beta8-6075.0_r8*W6*p2*beta2+34650.0_r8*beta4*W5*p3&
               +69300.0_r8*beta6*W5*p3+34650.0_r8*beta8*W5*p3-26337.0_r8*beta4*W9*p+&
                3009.0_r8*beta6*W9*p+3009.0_r8*W8*beta6-3645.0_r8*beta4*W7*p-3645.0_r8*beta6*&
                W7*p-765765.0_r8*beta6*W9*p7-3645.0_r8*W8*p2+3645.0_r8*W7*p+369603.0_r8*beta6*W9*p5&
               -35390.0_r8*beta6*W9*p3+19845.0_r8*W8*p2*beta4-3645.0_r8*W8*p2*beta6&
               -50544.0_r8*W9*p+19845.0_r8*W8*p2*beta2-594405.0_r8*W9*beta2*p3&
               +549045.0_r8*W9*beta4*p3+30375.0_r8*W9*p3+12150.0_r8*beta4*W7*p3+6075.0_r8*beta2*W7*p3&
               +6075.0_r8*beta6*W7*p3-20790.0_r8*beta2*W8*p4+20790.0_r8*beta6*W8*p4&
               -1300806.0_r8*beta4*W9*p5+369603.0_r8*beta2*W9*p5+765765.0_r8*beta4*W9*p7&
               -17325.0_r8*W8*beta6*p6-17325.0_r8*beta4*W8*p6)*beta2/(W8*be2p14);
          A0=xpowy((zeta2-beta2)/(1.0_r8+beta2*(1.0_r8-x2)),0.25_r8)
      ENDIF 
    ELSE
      rhsi=2.0_r8*(pi*0.5_r8-atan(p))-beta*log((beta*p+1.0_r8)/(beta*p-1.0_r8))
      omega=rhsi*0.5_r8/beta
      IF ((omega>-0.5_r8).AND.(omega<0.0_r8)) THEN
        zeta=beta*(1.0_r8-xpowy(-1.5_r8*omega/sqrt(2.0_r8),2.0_r8/3.0_r8))
      ELSE
        zeta=beta*2.0_r8*exp(omega-1.0_r8)
      ENDIF      
      zeta2=zeta*zeta
      iter=0
      DO WHILE ((conv>toln).AND.(iter<100))
        iter=iter+1
        boz=beta/zeta
        f=2.0_r8*(sqrt(beta2-zeta2)-beta*log(boz+sqrt(boz*boz-1.0_r8)))-rhsi;
        fp=2.0_r8*sqrt(beta2/zeta2-1.0_r8)
        znew=zeta-f/fp
        conv=abs(znew-zeta)
        zeta=znew
        zeta2=zeta*zeta
      ENDDO
      W=sqrt(beta2-zeta2)
      W2=W*W
      W3=W2*W
      W4=W2*W2
      W5=W4*W
      W6=W4*W2
      W7=W6*W
      W8=W6*W2
      W9=W8*W
      W10=W8*W2
      W11=W10*W
      W12=W10*W2
      W13=W12*W
      W14=W12*W2
      W16=W14*W2
      W18=W16*W2
      W20=W18*W2
      IF (W<0.2) THEN
        c=xpowy(1.0_r8/(1.0_r8+beta2),1.0_r8/3.0_r8)
        c2=c*c; c3=c2*c; c4=c3*c; c5=c4*c; c6=c5*c; c7=c6*c; c8=c7*c;
        c9=c8*c; c10=c9*c; c11=c10*c; c12=c11*c; c13=c12*c; c14=c13*c;
        c15=c14*c; c16=c15*c; c17=c16*c; c18=c17*c; c19=c18*c; c20=c19*c;
        c21=c20*c
        cdenom=c2+c+1.0_r8
        cdenom2=cdenom*cdenom
        cdenom3=cdenom2*cdenom
        cdenom4=cdenom3*cdenom         
        cdenom5=cdenom4*cdenom
        cdenom6=cdenom5*cdenom
        cdenom7=cdenom6*cdenom  
        cdenom8=cdenom7*cdenom
        cdenom9=cdenom8*cdenom
        ak(2)=(1.0_r8/201600.0_r8*(1393*c8+1393.0_r8*c7+97*c6-1730.0_r8*c5&
              -1730*c4-434*c3+145*c2+145*c+721.0_r8)/cdenom*W6+1.0_r8/11088000.0_r8&
              *(61641*c10+123282.0_r8*c9+167675.0_r8*c8+88786.0_r8*c7+17025.0_r8*c6&
               -28864.0_r8*c5-17504.0_r8*c4-13272.0_r8*c3-14128.0_r8*c2-&
              10592.0_r8*c-10224.0_r8)*c2/cdenom2*W8+1.0_r8/2522520000.0_r8*&
               (7150374.0_r8*c12+21451122.0_r8*c11+36020889.0_r8*c10+32983740.0_r8*c9&
              +14301635.0_r8*c8-6262716.0_r8*c7-14919053.0_r8*c6-14610316.0_r8*c5&
              -11834100.0_r8*c4-9657260.0_r8*c3-6688048.0_r8*c2-3310224.0_r8*c-707168.0_r8)*c4&
              /cdenom3*W10+1.0_r8/17463600000.0_r8*(22888390.0_r8*c14+91553560.0_r8*c13&
              +186272716.0_r8*c12+218657894.0_r8*c11+132101885.0_r8*c10-23474560.0_r8*c9&
              -135459915.0_r8*c8-154804164.0_r8*c7-115901691.0_r8*c6-68398480.0_r8*c5&
              -31760320.0_r8*c4-6014080.0_r8*c3+7134928.0_r8*c2+7581952.0_r8*c+2937760.0_r8)*c6&
              /cdenom4*W12+1.0_r8/1500899400000.0_r8*(860675886.0_r8*c16+4303379430.0_r8*c15&
              +10459577580.0_r8*c14+14858486859.0_r8*c13+11534759793.0_r8*c12+380167041.0_r8*c11&
              -11018380556.0_r8*c10-15148048045.0_r8*c9-11672657236.0_r8*c8-5518835692.0_r8*c7&
              -755351417.0_r8*c6+1624555844.0_r8*c5+2194550780.0_r8*c4+1646476012.0_r8*c3+&
               723229184.0_r8*c2+132301600.0_r8*c-36402688.0_r8)*c8/cdenom5*W14&
               -1.0_r8/141159588570000000.0_r8*(34172732883195.0_r8*c18+205036397299170.0_r8*c17&
               +585142109710263.0_r8*c16+982070424905148.0_r8*c15+943000129848528.0_r8*c14&
               +222975520459740.0_r8*c13-789007393738082.0_r8*c12-1362631504977174.0_r8*c11&
               -1163145597443907.0_r8*c10-509586706230712.0_r8*c9+74713644479970.0_r8*c8+&
               339651442452438.0_r8*c7+337238315362017.0_r8*c6+213724589653944.0_r8*c5+&
               80569717917264.0_r8*c4-945754950760.0_r8*c3-19890984887376.0_r8*c2-10453159775136.0_r8*c&
              -807559327280.0_r8)*c10/cdenom6*W16)/(-W6);
         ak(3)=(-1.0_r8/14515200.0_r8*(c-1.0_r8)*(145649.0_r8*c8+145649.0_r8*c7+141761.0_r8*c6&
               +29390.0_r8*c5+29390.0_r8*c4+33278.0_r8*c3+36065.0_r8*c2+36065.0_r8*c+37793)*W6&
               -1.0_r8/266112000.0_r8*(61641*c10+123282.0_r8*c9+167675.0_r8*c8+88786.0_r8*c7&
               +17025.0_r8*c6-28864.0_r8*c5-17504.0_r8*c4-13272.0_r8*c3-14128.0_r8*c2&
               -10592.0_r8*c-10224.0_r8)*(c-1.0_r8)*c2/cdenom*W8-1.0_r8/60540480000.0_r8*&
               (7150374.0_r8*c12+21451122.0_r8*c11+36020889.0_r8*c10+32983740.0_r8*c9+&
                14301635.0_r8*c8-6262716.0_r8*c7-14919053.0_r8*c6-14610316.0_r8*c5&
                -11834100.0_r8*c4-9657260.0_r8*c3-6688048.0_r8*c2-3310224.0_r8*c-707168.0_r8)*&
                (c-1.0_r8)*c4/cdenom2*W10-1.0_r8/419126400000.0_r8*(22888390.0_r8*c14+&
                91553560.0_r8*c13+186272716.0_r8*c12+218657894.0_r8*c11+132101885.0_r8*c10&
                -23474560.0_r8*c9-135459915.0_r8*c8-154804164.0_r8*c7-115901691.0_r8*c6&
                -68398480.0_r8*c5-31760320.0_r8*c4-6014080.0_r8*c3+7134928.0_r8*c2+7581952.0_r8*c&
                +2937760)*(c-1.0_r8)*c6/cdenom3*W12-1.0_r8/36021585600000.0_r8*(860675886.0_r8*c16&
                +4303379430.0_r8*c15+10459577580.0_r8*c14+14858486859.0_r8*c13+11534759793.0_r8*c12&
                +380167041.0_r8*c11-11018380556.0_r8*c10-15148048045.0_r8*c9-&
                11672657236.0_r8*c8-5518835692.0_r8*c7-755351417.0_r8*c6+1624555844.0_r8*c5&
                +2194550780.0_r8*c4+1646476012.0_r8*c3+723229184.0_r8*c2+132301600.0_r8*c-36402688.0_r8)*&
                (c-1.0_r8)*c8/cdenom4*W14-1.0_r8/3387830125680000000.0_r8*(34172732883195.0_r8*c18&
                +205036397299170.0_r8*c17+585142109710263.0_r8*c16+982070424905148.0_r8*c15&
                +943000129848528.0_r8*c14+222975520459740.0_r8*c13-789007393738082.0_r8*c12&
                -1362631504977174.0_r8*c11-1163145597443907.0_r8*c10-509586706230712.0_r8*c9&
                 +74713644479970.0_r8*c8+339651442452438.0_r8*c7+337238315362017.0_r8*c6&
                 +213724589653944.0_r8*c5+80569717917264.0_r8*c4-945754950760.0_r8*c3&
                -19890984887376.0_r8*c2-10453159775136.0_r8*c-807559327280.0_r8)*(c-1.0_r8)*c10/cdenom5*&
                 W16)/(-W6);
        bk(1)=(-1.0_r8/280.0_r8*c*zeta*(9*c5+9*c4+9*c3+4*c+4)/cdenom*W4&
              -1/12600.0_r8*(98*c7+196*c6+213*c5+83*c4-47*c3-96*c2-76*c-56)*zeta*c3&
               /cdenom2*W6-1.0_r8/1617000.0_r8*(4077*c9+12231*c8+16916*c7+9978*c6-3237*c5&
               -11410*c4-10413*c3-5592*c2-1828*c+828)*zeta*c5/cdenom3*W8-1.0_r8&
               /315315000_r8*(271656*c11+1086624*c10+1921545*c9+1573920*c8-107700*c7&
               -1684734*c6-1875571*c5-975175*c4-83470*c3+273160*c2+252628*c+22792)*zeta*&
                c7/cdenom4*W10-1.0_r8/33108075000.0_r8*(9926995.0_r8*c13+49634975.0_r8*c12+&
                108971493.0_r8*c11+117996680.0_r8*c10+23304685.0_r8*c9-112142745.0_r8*c8&
                -164160315.0_r8*c7-99982062.0_r8*c6-2998230*c5+44249930.0_r8*c4+37346320.0_r8*c3&
                +13501080.0_r8*c2-2899976.0_r8*c-1250080.0_r8)*zeta*c9/cdenom5*W12&
                -1.0_r8/938062125000.0_r8*(98709912.0_r8*c15+592259472.0_r8*c14+&
                 1566631407.0_r8*c13+2144609910.0_r8*c12+993484611.0_r8*c11-1648415217.0_r8*c10&
                -3459671744.0_r8*c9-2675618409.0_r8*c8-303627696.0_r8*c7+1332420910.0_r8*c6&
                +1322416896.0_r8*c5+514721244.0_r8*c4-73283846.0_r8*c3-186087792.0_r8*c2-&
                 37851504.0_r8*c+4317096)*zeta*c11/cdenom6*W14-1.0_r8/61757319999375000.0_r8&
                 *(2290875403815.0_r8*c17+16036127826705.0_r8*c16+49847663782476.0_r8*c15&
                 +83192099634792.0_r8*c14+60433694450568.0_r8*c13-42447122909352.0_r8*c12&
                 -150313744710496.0_r8*c11-151930198189654.0_r8*c10-40647228995323.0_r8*c9&
                 +71949889801323.0_r8*c8+94022483075388.0_r8*c7+44189908944534.0_r8*c6-5389203336702.0_r8*c5&
                 -19139960492979.0_r8*c4-9728875988691.0_r8*c3+299051082664.0_r8*c2+984902727172.0_r8*c&
                 +51375356260.0_r8)*zeta*c13/cdenom7*W16-1.0_r8/308786599996875000.0_r8*(4048556587768.0_r8*&
                  c19+32388452702144.0_r8*c18+115966656910053.0_r8*c17+229357243890444.0_r8*c16+&
                  226393507645548.0_r8*c15-38068980041772.0_r8*c14-439693813180464.0_r8*c13&
                 -587811603844524.0_r8*c12-270394076368940.0_r8*c11+226230760410240.0_r8*c10&
                 +439866482466809.0_r8*c9+268458673161834.0_r8*c8-5466610953768.0_r8*c7&
                 -122882117855892.0_r8*c6-81933726464274.0_r8*c5-10976641917615.0_r8*c4&
                 +15049928649936.0_r8*c3+5555561805944.0_r8*c2-336943037932.0_r8*c&
                 -125602833664.0_r8)*zeta*c15/cdenom8*W18-1.0_r8/3231451768967296875000.0_r8*&
                 (15000147517319766.0_r8*c21+135001327655877894.0_r8*c20+547902204206413110.0_r8*c19&
                 +1256083696114375320.0_r8*c18+1564644811187692800.0_r8*c17+338396649780494412.0_r8*c16&
                 -2373267185032712572.0_r8*c15-4279605957210423960.0_r8*c14-2915791867689908985.0_r8*c13&
                 +973459027265184775.0_r8*c12+3730399203798277107.0_r8*c11+3042903732882256788.0_r8*c10&
                 +396117195396233010.0_r8*c9-1328396918245677270.0_r8*c8-1187627541412476150.0_r8*c7&
                 -291199621343628524.0_r8*c6+211393950537241899.0_r8*c5+181498437006912315.0_r8*c4+&
                  19874016923212235.0_r8*c3-16503501533468100.0_r8*c2-2542411653780036.0_r8*c+&
                  155654367677916.0_r8)*zeta*c17/cdenom9*W20)/W4
        bk(2)=(1.0_r8/6720.0_r8*c*(9*c5+9*c4+9*c3+4*c+4)*(c-1.0_r8)*zeta*W4&
              +1.0_r8/302400.0_r8*(98*c7+196*c6+213*c5+83*c4-47*c3-96*c2-76*c-56)*c3*(c-1.0_r8)&
              *zeta/cdenom*W6+1.0_r8/38808000.0_r8*(4077*c9+12231.0_r8*c8+16916.0_r8*c7&
              +9978.0_r8*c6-3237.0_r8*c5-11410.0_r8*c4-10413.0_r8*c3-5592.0_r8*c2&
              -1828.0_r8*c+828)*c5*(c-1.0_r8)*zeta/cdenom2*W8&
              +1.0_r8/7567560000.0_r8*(271656.0_r8*c11+1086624.0_r8*c10+&
              1921545.0_r8*c9+1573920.0_r8*c8-107700.0_r8*c7-1684734.0_r8*c6&
              -1875571.0_r8*c5-975175.0_r8*c4-83470.0_r8*c3+273160.0_r8*c2+252628.0_r8*c&
              +22792.0_r8)*c7*(c-1.0_r8)*zeta/cdenom3*W10+1.0_r8/794593800000.0_r8*&
              (9926995.0_r8*c13+49634975.0_r8*c12+108971493.0_r8*c11+117996680.0_r8*c10&
              +23304685.0_r8*c9-112142745.0_r8*c8-164160315.0_r8*c7-99982062.0_r8*c6&
              -2998230.0_r8*c5+44249930.0_r8*c4+37346320.0_r8*c3+13501080.0_r8*c2&
              -2899976.0_r8*c-1250080.0_r8)*c9*(c-1.0_r8)*zeta/cdenom4*W12+1.0_r8/22513491000000.0_r8*&
              (98709912.0_r8*c15+592259472.0_r8*c14+1566631407.0_r8*c13+2144609910.0_r8*c12&
              +993484611.0_r8*c11-1648415217.0_r8*c10-3459671744.0_r8*c9-2675618409.0_r8*c8&
              -303627696.0_r8*c7+1332420910.0_r8*c6+1322416896.0_r8*c5+514721244.0_r8*c4&
              -73283846.0_r8*c3-186087792.0_r8*c2-37851504.0_r8*c+4317096.0_r8)*c11*(c-1.0_r8)*&
              zeta/cdenom5*W14+1.0_r8/1482175679985000000.0_r8*(2290875403815.0_r8*c17&
              +16036127826705.0_r8*c16+49847663782476.0_r8*c15+83192099634792.0_r8*c14+&
              60433694450568.0_r8*c13-42447122909352.0_r8*c12-150313744710496.0_r8*c11&
              -151930198189654.0_r8*c10-40647228995323.0_r8*c9+71949889801323.0_r8*c8+&
              94022483075388.0_r8*c7+44189908944534.0_r8*c6-5389203336702.0_r8*c5-&
              19139960492979.0_r8*c4-9728875988691.0_r8*c3+299051082664.0_r8*c2+&
              984902727172.0_r8*c+51375356260.0_r8)*c13*(c-1.0_r8)*zeta/cdenom6*W16+&
              1.0_r8/7410878399925000000.0_r8*(4048556587768.0_r8*c19+32388452702144.0_r8*c18&
              +115966656910053.0_r8*c17+229357243890444.0_r8*c16+226393507645548.0_r8*c15&
              -38068980041772.0_r8*c14-439693813180464.0_r8*c13-587811603844524.0_r8*c12&
              -270394076368940.0_r8*c11+226230760410240.0_r8*c10+439866482466809.0_r8*c9&
              +268458673161834.0_r8*c8-5466610953768.0_r8*c7-122882117855892.0_r8*c6&
              -81933726464274.0_r8*c5-10976641917615.0_r8*c4+15049928649936.0_r8*c3+&
              5555561805944.0_r8*c2-336943037932.0_r8*c-125602833664.0_r8)*c15*(c-1.0_r8)*&
              zeta/cdenom7*W18+1.0_r8/77554842455215125000000.0_r8*(15000147517319766.0_r8*c21&
              +135001327655877894.0_r8*c20+547902204206413110.0_r8*c19+1256083696114375320.0_r8*c18+&
              1564644811187692800.0_r8*c17+338396649780494412.0_r8*c16-2373267185032712572.0_r8*c15&
             -4279605957210423960.0_r8*c14-2915791867689908985.0_r8*c13+&
              973459027265184775.0_r8*c12+3730399203798277107.0_r8*c11+3042903732882256788.0_r8*c10&
              +396117195396233010.0_r8*c9-1328396918245677270.0_r8*c8-1187627541412476150.0_r8*c7&
              -291199621343628524.0_r8*c6+211393950537241899.0_r8*c5+181498437006912315.0_r8*c4&
              +19874016923212235.0_r8*c3-16503501533468100.0_r8*c2-2542411653780036.0_r8*c+&
              155654367677916.0_r8)*c17*(c-1.0_r8)*zeta/cdenom8*W20)/W4;
        bk(3)=(-1.0_r8/524160000.0_r8*c*zeta*(6059583.0_r8*c13+12119166.0_r8*c12+&
               18178749.0_r8*c11+9111312.0_r8*c10+95615.0_r8*c9-8920082.0_r8*c8-5698629.0_r8*c7&
               -2545036.0_r8*c6+608557.0_r8*c5+546510.0_r8*c4+1181699.0_r8*c3+1816888.0_r8*c2&
               +1229112.0_r8*c+614556.0_r8)/cdenom2*W10-1.0_r8/1816214400000.0_r8*(9727357705.0_r8*c15&
               +29182073115.0_r8*c14+54164855211.0_r8*c13+55512024253.0_r8*c12+33367820181.0_r8*c11&
               -1784732145.0_r8*c10-19542687620.0_r8*c9-20194526124.0_r8*c8-12274421967.0_r8*c7&
               -8011196494.0_r8*c6-7270887045.0_r8*c5-7771225740.0_r8*c4-7435063897.0_r8*c3&
               -6252124236.0_r8*c2-3350967852.0_r8*c-1258951720.0_r8)*zeta*c3/cdenom3*W12&
               -1.0_r8/51459408000000.0_r8*(144811913265.0_r8*c17+579247653060.0_r8*c16+1282754051665.0_r8&
               *c15+1728105185620.0_r8*c14+1451348623879.0_r8*c13+484314490976.0_r8*c12&
               -489620336905.0_r8*c11-945041189740.0_r8*c10-873657521970.0_r8*c9-631216670541.0_r8*c8&
               -445274929539.0_r8*c7-335069642665.0_r8*c6-245861617660.0_r8*c5-156878390545.0_r8*c4&
               -77977426823.0_r8*c3-17492885312.0_r8*c2+3279403240.0_r8*c+6170395620.0_r8)*zeta*c5/cdenom4*W14&
               -1.0_r8/20532303792000000.0_r8*(29048679511047.0_r8*c19+145243397555235.0_r8*c18+&
                377950239272970.0_r8*c17+611632330911882.0_r8*c16+633363321305990.0_r8*c15+&
                334989811415962.0_r8*c14-124793259618935.0_r8*c13-462893943016455.0_r8*c12&
                -532059257514943.0_r8*c11-409102997160979.0_r8*c10-246328251031409.0_r8*c9&
                -125738641603410.0_r8*c8-49813792932765.0_r8*c7-1352246462804.0_r8*c6+&
                26758233574976.0_r8*c5+35762086322945.0_r8*c4+30037846475850.0_r8*c3+&
                16475189210860.0_r8*c2+5900059068820.0_r8*c+687614243288.0_r8)*zeta*c7/cdenom5*W16)/(-W10);  
        IF (abs(zeta-beta)<1.e-10_r8) THEN
          A0=beta2*c+2.0_r8*beta*(3.0_r8+2.0_r8*beta2+2.0_r8/c2)*(zeta-beta)&
             /5.0_r8/be2p1
        ELSE 
          A0=xpowy((zeta2-beta2)/(1.0_r8+beta2*(1.0_r8-x2)),0.25_r8)
          A0=sqrt(sqrt((zeta-beta)*(zeta+beta)/(1.0_r8+beta2*(1.0_r8-x2))))
        ENDIF
      ELSE
        p=x/sqrt(beta*beta*(x*x-1.0_r8)-1.0_r8);        
        bk(1)=-1.0_r8/24.0_r8*(-5*beta4+3*W3*p*beta2+3*W2*beta2+5*W3*p3*beta2-5*beta2+&
              3*W2-3*W3*p)*zeta/(W4*be2p1);
        bk(2)=-1.0_r8/576.0_r8*zeta*beta2*(-5*beta4+3*W3*p*beta2+3*W2*beta2+5*W3*p3*beta2&
              -5*beta2+3*W2-3*W3*p)/(W4*be2p12);
        bk(3)=-1.0_r8/414720.0_r8*(-20790.0_r8*W5*p*beta8+34650.0_r8*W6*beta4*p4+30375.0_r8*W6&
              -32355.0_r8*W9*beta4*p-1276275.0_r8*beta10+2297295.0_r8*W2*beta8+765765.0_r8*W2&
              *beta10+194400.0_r8*beta2*W9*p-28875.0_r8*W6*beta6*p6+3240.0_r8*W8-20790.0_r8*W8*p4&
              *beta2+3645.0_r8*W7*p*beta4+6075.0_r8*W7*p3*beta2+12150.0_r8*W7*p3*beta4+&
               85725.0_r8*W6*beta2+30300.0_r8*W9*beta6*p3+369603.0_r8*W9*beta6*p5-765765.0_r8*W9*&
               beta4*p7+57750.0_r8*W3*p3*beta8+75.0_r8*W6*beta8-1108809.0_r8*W4*beta4-1108809.0_r8*&
               W4*beta6-20790.0_r8*W5*p*beta6-425425.0_r8*beta6+33075.0_r8*W6*beta4*p2+&
               28875.0_r8*W3*beta10*p3-6075.0_r8*W6*beta8*p2+35850.0_r8*W6*beta6+17325.0_r8*W3*beta10*p&
              -369603.0_r8*W4*beta2-32400.0_r8*W9*p-3645.0_r8*W7*p*beta2-3645.0_r8*W7*p+6075.0_r8*W7*p3&
              *beta6-564165.0_r8*W9*beta4*p3+91125.0_r8*W6*beta4+17325.0_r8*W3*beta8*p+17325.0_r8*W8*&
               beta6*p6-1300806.0_r8*W9*beta4*p5+33075.0_r8*W6*beta6*p2-425425.0_r8*beta12&
              -34650.0_r8*W5*p3*beta4-28875.0_r8*W6*beta8*p6-34650.0_r8*W6*beta8*p4+20790.0_r8*W8&
              *beta6*p4+17325.0_r8*W8*beta4*p6+3645.0_r8*W8*p2+425425.0_r8*W9*beta6*p9-6075.0_r8*W6&
              *beta2*p2-369603.0_r8*W4*beta8+765765.0_r8*W9*beta6*p7+28875.0_r8*beta6*W3*p3+&
               2297295.0_r8*W2*beta6-17325.0_r8*W3*beta6*p+765765.0_r8*W2*beta4-1276275.0_r8*beta8&
              -17325.0_r8*W3*beta4*p+45.0_r8*W8*beta6+3645.0_r8*W8*beta6*p2-30375.0_r8*W9*p3+&
               3645.0_r8*W7*p*beta6+20790.0_r8*W5*p*beta4+20790.0_r8*W5*p*beta2+&
               564165.0_r8*W9*beta2*p3+369603.0_r8*W9*beta2*p5-19845.0_r8*W8*beta2*p2&
              -19845.0_r8*W8*beta4*p2-3285.0_r8*W8*beta4-34650.0_r8*W5*p3*beta8-45.0_r8*W9*beta6*p&
              -69300.0_r8*W5*p3*beta6)*zeta/(-W10*(beta6+1.0_r8+3*beta2+3*beta4));
        ak(2)=1.0_r8/1152.0_r8*(594.0_r8*W2*beta6-W6*beta4+594.0_r8*W2*beta2-42.0_r8*W3*p*beta2&
              -270*W4*beta2-455.0_r8*beta8-72.0_r8*W6*beta2+81.0_r8*W6*p2+72.0_r8*W6-54.0_r8*W5*p*&
               beta4+54.0_r8*W5*p+70.0_r8*W3*beta6*p3+70.0_r8*W3*beta4*p3-910.0_r8*beta6-455.0_r8&
              *beta4-135.0_r8*W4+42.0_r8*W3*beta6*p+1188.0_r8*W2*beta4+385.0_r8*W6*&
               beta4*p6+462.0_r8*W6*beta4*p4-462.0_r8*W6*beta2*p4-522.0_r8*W6*beta2*p2&
              -90.0_r8*W5*p3*beta4-90*W5*p3*beta2-135.0_r8*W4*beta4+81.0_r8*W6*beta4*p2)/(-W6*be2p12);
        ak(3)=1.0_r8/414720.0_r8*(6930.0_r8*W6*beta4*p4+7128.0_r8*W6+1215.0_r8*W6*p2+1944.0_r8*W6*beta2&
              -2025.0_r8*W4*beta4-13650.0_r8*beta6-6825.0_r8*beta4+1215.0_r8*W6*beta4*p2&
              -4050.0_r8*W4*beta2-6930.0_r8*W6*beta2*p4+1003.0_r8*W6*beta4+810.0_r8*W5*p&
              -1350.0_r8*W5*p3*beta2-1350.0_r8*W5*p3*beta4-7830.0_r8*W6*beta2*p2+5775.0_r8*W6*beta4*p6&
              +1050.0_r8*beta6*W3*p3+8910.0_r8*W2*beta6+630.0_r8*W3*beta6*p+17820.0_r8*W2*beta4+&
               1050.0_r8*W3*p3*beta4-6825.0_r8*beta8-630.0_r8*W3*p*beta2+8910.0_r8*W2*beta2&
              -810.0_r8*W5*p*beta4-2025.0_r8*W4)*beta2/(-W6*(beta6+1.0_r8+3*beta2+3.0_r8*beta4));
        ak(4)=1.0_r8/39813120.0_r8*(4451265.0_r8*beta8*W12*p4-2430.0_r8*beta8*&
               W12*p2+748417428.0_r8*beta4*W12*p6+490305150.0_r8*beta4*W12*p4&
               -1022624460.0_r8*beta6*W12*p8-5740875.0_r8*W8-6107940.0_r8*beta6*W12*p2&
               +93484530.0_r8*beta4*W12*p2-178129440.0_r8*beta6*W12*p4-748417428.0_r8*&
                beta6*W12*p6+94110126.0_r8*beta8*W12*p6+349922430.0_r8*beta8*W12*p8&
               +349922430.0_r8*beta4*W12*p8+1632960.0_r8*W12+4914000.0_r8*W6*beta10*p2&
               -202076875.0_r8*beta8-11411400.0_r8*beta8*W3*p-10510500.0_r8*W6*beta10*p6&
               -808307500.0_r8*beta10-1871100.0_r8*beta8*W10*p4-328050.0_r8*beta8*W10*p2&
               -178143300.0_r8*W12*p4*beta2+111234708.0_r8*beta2*W6-1212461250.0_r8*beta12&
               -808307500.0_r8*beta14-202076875.0_r8*beta16+295650.0_r8*W10*beta4+299700.0_r8*W10&
               *beta6+11911900.0_r8*beta8*W9*p9-5255250.0_r8*W6*beta12*p6-6415200.0_r8*W8*beta8*p2&
               +28528500.0_r8*beta10*W3*p3-93486960.0_r8*beta2*W12*p2+6123600.0_r8*W12*p2&
               +1443420.0_r8*W8*beta10*p2-291600.0_r8*W10+113400.0_r8*beta2*W9*p&
               +1620.0_r8*W11*beta8*p-15315300.0_r8*W11*beta8*p9-27567540.0_r8*W11*beta8*p7&
               +493152660.0_r8*beta6*W2+443956032.0_r8*beta4*W6-396578754.0_r8*beta4*W4&
               -21680460.0_r8*beta2*W8-13305708.0_r8*W11*beta8*p5+6306300.0_r8*W6*beta8*p4&
               -1090800.0_r8*W11*beta8*p3+9509500.0_r8*beta8*W3*p3-5420844.0_r8*beta2*W7*p&
               -33162210.0_r8*beta4*W8-1559250.0_r8*W10*beta8*p6+28528500.0_r8*beta12*W3*p3&
               +4050.0_r8*beta8*W10-4027.0_r8*beta8*W12-9936.0_r8*beta6*W12+185910725.0_r8*beta8&
               *W12*p12+666425448.0_r8*beta6*W6-1586315016.0_r8*beta6*W4+1972610640.0_r8*beta8*W2&
               +446185740.0_r8*beta8*W12*p10-9486720.0_r8*beta2*W12-291600.0_r8*W10*beta2&
               -6306300.0_r8*W6*beta10*p4+4465125.0_r8*W12*p4-5255250.0_r8*W6*beta8*p6&
               +1458000.0_r8*W10*p2*beta2-5705700.0_r8*beta6*W3*p-19216440.0_r8*W11*beta2*p3&
               -1105650.0_r8*W6*beta12*p2+1606608.0_r8*beta4*W12-94121676.0_r8*beta2*W12*p6&
               +445935282.0_r8*W6*beta8+112244808.0_r8*W6*beta10+13650.0_r8*W6*beta12&
               -2379472524.0_r8*W4*beta8-1586315016.0_r8*W4*beta10-396578754.0_r8*W4*beta12&
               +2958915960.0_r8*beta10*W2+1972610640.0_r8*beta12*W2+493152660.0_r8*beta14*W2&
               +6306300.0_r8*beta6*W6*p4+4914000.0_r8*beta6*W6*p2-1105650.0_r8*beta4*W6*p2&
               +12039300.0_r8*W6*beta8*p2-6306300.0_r8*W6*beta12*p4+10602900.0_r8*W5*p*beta4&
               +21205800.0_r8*W5*beta6*p-21205800.0_r8*W5*p*beta10+5705700.0_r8*beta14*W3*p&
               -17671500.0_r8*beta6*W5*p3-53014500.0_r8*beta8*W5*p3-53014500.0_r8*beta10*W5*p3&
               +11411400.0_r8*beta12*W3*p+9509500.0_r8*beta14*W3*p3-17671500.0_r8*beta12*W5*p3&
               -10602900.0_r8*beta12*W5*p-13305708.0_r8*beta2*W11*p5-1559250.0_r8*W10*beta4*p6&
               +4536000.0_r8*beta4*W9*p+3516660.0_r8*beta6*W9*p-1417500.0_r8*beta8*W9*p&
               -1260.0_r8*beta10*W9*p+8232840.0_r8*W8*beta8*p4+8232840.0_r8*W8*beta10*p4&
               -24264360.0_r8*W8*beta6-7059555.0_r8*W8*beta8-17820.0_r8*W8*beta10-10841688.0_r8*beta4*W7*p&
               +10841688.0_r8*beta8*W7*p+5420844.0_r8*beta10*W7*p+11911900.0_r8*beta10*W9*p9&
               +21441420.0_r8*beta10*W9*p7+10348884.0_r8*beta10*W9*p5+848400.0_r8*beta10*W9*p3&
               -26073684.0_r8*beta8*W9*p5-15798720.0_r8*beta8*W9*p3-21441420.0_r8*beta6*W9*p7&
               +1871100.0_r8*W10*beta4*p4-26073684.0_r8*beta6*W9*p5-2551500.0_r8*beta6*W9*p3&
               -6415200.0_r8*W8*p2*beta4-15717240.0_r8*W8*p2*beta6+510300.0_r8*W9*p&
               +1443420.0_r8*W8*p2*beta2-1701000.0_r8*W9*beta2*p3+12394620.0_r8*W9*beta4*p3&
               +1871100.0_r8*W10*beta2*p4-5832000.0_r8*W11*p*beta2+1166400.0_r8*W11*beta6*p+&
                19219140.0_r8*W11*beta6*p3-5833620.0_r8*W11*p*beta4+9034740.0_r8*beta4*W7*p3+27104220.0_r8*&
                beta8*W7*p3+27104220.0_r8*beta6*W7*p3+9034740.0_r8*beta10*W7*p3+1166400.0_r8*W11*p+&
                1093500.0_r8*W11*p3-8232840.0_r8*beta6*W8*p4-8232840.0_r8*beta4*W8*p4+&
                10348884.0*beta4*W9*p5+6860700.0_r8*W8*beta10*p6+6860700.0_r8*W8*beta6*p6&
               +13721400.0_r8*W8*beta8*p6+27567540.0_r8*W11*beta4*p7-446185740.0_r8*beta6*W12*p10&
               +33523308.0_r8*W11*beta4*p5+33523308.0_r8*W11*beta6*p5-328050.0_r8*W10*p2&
               -1871100.0_r8*W10*beta6*p4-3118500.0_r8*W10*beta6*p6-15315300.0_r8*W11*beta6*p9&
               +3572100.0_r8*W10*p2*beta4+1458000.0_r8*W10*p2*beta6)/(W12*be2p14);              
         bk(4)=-1.0_r8/9953280.0_r8*zeta*(21384.0_r8*W8+30375.0_r8*W6-425425.0_r8*beta6&
               -1276275.0_r8*beta8+17325.0_r8*beta8*W3*p+20790.0_r8*W5*p*beta2-&
                1276275.0_r8*beta10+55485.0_r8*beta2*W6-369603.0_r8*beta2*W4-425425.0_r8*beta12&
               +765765.0_r8*beta4*W2+28875.0_r8*beta10*W3*p3+425425.0_r8*beta6*W9*p9&
               +28875.0_r8*beta6*W3*p3+203472.0_r8*beta2*W9*p+2297295.0_r8*beta6*W2&
               +45765.0_r8*beta4*W6-1108809.0_r8*beta4*W4+27216.0_r8*beta2*W8-34650.0_r8*W6*beta8*p4&
               +57750.0_r8*beta8*W3*p3-3645.0_r8*beta2*W7*p+8841.0_r8*beta4*W8-28875.0_r8*&
                beta6*W6*p6+17325.0_r8*beta10*W3*p+15640.0_r8*beta6*W6-1108809.0_r8*beta6*W4&
               +2297295.0_r8*beta8*W2-28875.0_r8*W6*beta8*p6-17325.0_r8*beta6*W3*p-&
                17325.0_r8*beta4*W3*p-5015.0_r8*W6*beta8-369603.0_r8*W4*beta8+765765.0_r8*beta10*W2&
               +33075.0_r8*beta6*W6*p2+34650.0_r8*beta4*W6*p4+33075.0_r8*beta4*W6*p2&
               -6075.0_r8*W6*beta8*p2+20790.0_r8*W5*p*beta4-20790.0_r8*W5*beta6*p&
               -20790.0_r8*W5*p*beta8-6075.0_r8*W6*p2*beta2-34650.0_r8*beta4*W5*p3&
               -69300.0_r8*beta6*W5*p3-34650.0_r8*beta8*W5*p3-26337.0_r8*beta4*W9*p+&
                3009.0_r8*beta6*W9*p+3009.0_r8*W8*beta6+3645.0_r8*beta4*W7*p+3645.0_r8*beta6*&
                W7*p+765765.0_r8*beta6*W9*p7+3645.0_r8*W8*p2-3645.0_r8*W7*p+369603.0_r8*beta6*W9*p5&
               +35390.0_r8*beta6*W9*p3-19845.0_r8*W8*p2*beta4+3645.0_r8*W8*p2*beta6&
               -50544.0_r8*W9*p-19845.0_r8*W8*p2*beta2+594405.0_r8*W9*beta2*p3&
               -549045.0_r8*W9*beta4*p3-30375.0_r8*W9*p3+12150.0_r8*beta4*W7*p3+6075.0_r8*beta2*W7*p3&
               +6075.0_r8*beta6*W7*p3-20790.0_r8*beta2*W8*p4+20790.0_r8*beta6*W8*p4&
               -1300806.0_r8*beta4*W9*p5+369603.0_r8*beta2*W9*p5-765765.0_r8*beta4*W9*p7&
               +17325.0_r8*W8*beta6*p6+17325.0_r8*beta4*W8*p6)*beta2/(W8*be2p14);
          A0=sqrt(sqrt((zeta-beta)*(zeta+beta)/(1.0_r8+beta2*(1.0_r8-x2))))
        ENDIF
      ENDIF
      Ac=A0*(ak(0)+ak(1)/mu+ak(2)/mu2+ak(3)/mu3+ak(4)/mu4)
      Bc=A0*(bk(0)+bk(1)/mu+bk(2)/mu2+bk(3)/mu3+bk(4)/mu4)
      CALL dkia(1,mu*zeta,tau,dki,dkid,ierro)  
      fac2=1.0_r8
      DO n=1,mu
        fac2=sqq*fac2
      ENDDO
      IF (ierro==0) THEN
        pmmu=(exp(-mu*lambda)*fac2)*2.0_r8*gammah(mu)/&
           (sqrt(2.0_r8*pi)*pi)*(cosh(pi*tau)*(Ac*dki-Bc*dkid))
      ELSE
        pmmu=0.0_r8
        ierr=1
      ENDIF
      END SUBROUTINE expakia

      SUBROUTINE tauex1(x,mu,tau,pmmu,ierr)
    ! --------------------------------------------------
    ! Calculation of P^(mu)_(-1/2+i*tau)(x) by using 
    ! an asymptotic expansion for tau large, x>1.
    ! This expansion is used for m=0,1
    ! --------------------------------------------------
    ! Inputs:
    !   x ,    argument of the functions
    !   mu,    function parameter (integer value)
    !   tau,   function parameter (real value)
    ! Outputs:
    !   pmmu , P^(mu)_(-1/2+i*tau)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! ---------------------------------------------------------
     USE Someconstants
     REAL(r8), INTENT(IN) :: x
     INTEGER, INTENT(IN) :: mu
     REAL(r8), INTENT(IN) :: tau
     REAL(r8), INTENT(OUT) :: pmmu
     INTEGER,  INTENT(OUT) :: ierr
     REAL(r8) :: beta,s,x2,x4,x6,x8,x10,x12,x14,x16,x18,&
                 mu2,mu3,mu4,mu5,mu6,mu7,&
                 s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,&
                 s12,s13,s14,s15,s16,s17,s18,&
                 bn(0:18),fac2,tau2,sumt,T, poch
     INTEGER :: n, NN
     !argcosh(x) = ln( x + sqrt(x2 - 1) )
     ierr=0
     beta=log(x+sqrt((x-1.0_r8)*(x+1.0_r8))); s= sinh(beta);
     x2=x*x
     x4=x2*x2
     x6=x4*x2
     x8=x6*x2
     x10=x8*x2
     x12=x10*x2
     x14=x12*x2
     x16=x14*x2
     x18=x16*x2
     mu2=mu*mu
     mu3=mu2*mu
     mu4=mu3*mu
     mu5=mu4*mu
     mu6=mu5*mu
     mu7=mu6*mu
     s2=s*s
     s3=s2*s
     s4=s3*s
     s5=s4*s
     s6=s5*s
     s7=s6*s
     s8=s7*s
     s9=s8*s
     s10=s9*s
     s11=s10*s
     s12=s11*s
     s13=s12*s
     s14=s13*s
     s15=s14*s
     s16=s15*s
     s17=s16*s
     s18=s17*s
     bn(0)= 1.0_r8;
     bn(1)= (1.0_r8/4.0_r8)*(2*mu-1.0_r8)*x/s;
     bn(2)= (1.0_r8/96.0_r8)*(2*mu-1.0_r8)*(-8.0_r8+(6*mu-1.0_r8)*x2)/s2;
     bn(3)= (1.0_r8/384.0_r8)*(2*mu-1.0_r8)*x*((-1.0_r8+4*mu2)*x2+16.0_r8-16*mu)/s3;
     bn(4)= (1.0_r8/92160.0_r8)*(2*mu-1.0_r8)*((120.0_r8*mu3+60.0_r8*mu2-110.0_r8*mu+9)*x4+(-1008.0_r8-960*mu2+&
            1600*mu)*x2-576.0_r8+640*mu)/s4;
     bn(5)= (1.0_r8/368640.0_r8)*(2*mu-1.0_r8)*x*((19.0_r8-16*mu+64*mu3+48*mu4-88*mu2)*x4+&
          (832.0_r8-640*mu3+1280*mu2-1376*mu)*x2+1984-3072*mu+1280*mu2)/s5;
     bn(6)= (1.0_r8/185794560.0_r8)*(2*mu-1.0_r8)*((4718.0_r8*mu-5432.0_r8*mu2-5040.0_r8*mu3+&
          5040.0_r8*mu4+2016.0_r8*mu5-237.0_r8)*x6+(130368.0_r8*mu-79992.0_r8+80640.0_r8*mu3&
          -40320.0_r8*mu4-106176.0_r8*mu2)*x4+(-480960.0_r8+857472.0_r8*mu-575232.0_r8*mu2&
          +161280.0_r8*mu3)*x2-93696.0_r8+157696.0_r8*mu-71680.0_r8*mu2)/s6;
     bn(7)= (1.0_r8/743178240.0_r8)*(2*mu-1.0_r8)*x*((-825+1168*mu-5248*mu3-1200*mu4+3564*mu2+576*mu6&
            +2304*mu5)*x6+(-16128.0_r8*mu5-34944.0_r8*mu3+26880.0_r8*mu4-100848.0_r8*mu+52560.0_r8&
            +77952.0_r8*mu2)*x4+(-473088.0_r8*mu3+688320.0_r8+107520.0_r8*mu4+1053696.0_r8*mu2&
            -1319424.0_r8*mu)*x2+476160.0_r8-143360.0_r8*mu3+602112.0_r8*mu2-904192.0_r8*mu)/s7;
     bn(8)= (1.0_r8/356725555200.0_r8)*(2*mu-1.0_r8)*((17280.0_r8*mu7+100800.0_r8*mu6+10080.0_r8*mu5&
            -371280.0_r8*mu4+127400.0_r8*mu3+316500.0_r8*mu2-198638.0_r8*mu+7017)*x8+&
           (-3203808.0_r8+645120.0_r8*mu5-483840.0_r8*mu4-7962240.0_r8*mu2+3834880.0_r8*mu3&
            +7405952.0_r8*mu-645120.0_r8*mu6)*x6+(196922112.0_r8*mu-96582528.0_r8&
            -171555840.0_r8*mu2-31610880.0_r8*mu4+89241600.0_r8*mu3+6451200.0_r8*mu5)*x4&
            +(-158373888.0_r8+323796992.0_r8*mu-259031040.0_r8*mu2-17203200.0_r8*mu4&
            +100925440.0_r8*mu3)*x2-15495168.0_r8+31301632.0_r8*mu-22364160.0_r8*mu2+&
            5734400.0_r8*mu3)/s8;
     IF (mu==0) THEN
        bn(9)=-(1.0_r8/475634073600.0_r8)*x*(11813.0_r8*x8+444928.0_r8*x6+31916928.0_r8*x4&
              +105054208.0_r8*x2+34869248.0_r8)/s9;
        bn(10)=-(1.0_r8/1993133260800.0_r8)*(-677.0_r8*x10-216648.0_r8*x8&
             -28856448.0_r8*x6-174496768.0_r8*x4-131616768.0_r8*x2-7766016.0_r8)/s10;
        bn(11)=-(1.0_r8/3720515420160.0_r8)*x*(-2117.0_r8*x10+50768.0_r8*x8+10389120.0_r8*x6&
              +108591104.0_r8*x4+156037120.0_r8*x2+30474240.0_r8)/s11;
        bn(12)=-(1.0_r8/48753634065776640.0_r8)*(308963.0_r8*x12-46284528.0_r8*x10&
               -23934223680.0_r8*x8-412895119360.0_r8*x6-1022292602880.0_r8*x4-441849151488.0_r8*x2&
               -17470062592.0_r8)/s12;
        bn(13)=-(1.0_r8/4875363406577664000.0_r8)*x*(64604977.0_r8*x12-64640832.0_r8*x10&
                +385965769920.0_r8*x8+10677616844800.0_r8*x6+42672499077120.0_r8*x4&
               +34108050505728.0_r8*x2+4390826278912.0_r8)/s13;
        bn(14)=-(1.0_r8/1053078495820775424000.0_r8)*(-131301607.0_r8*x14&
              -4830158024.0_r8*x12-12347747442624.0_r8*x10-540439580577280.0_r8*x8&
              -3331531451863040.0_r8*x6-4442551522787328.0_r8*x4-1243461370249216.0_r8*x2&
              -35275466080256.0_r8)/s14;
        bn(15)=-(1.0_r8/842462796656620339200.0_r8)*x*(-263101079.0_r8*x14+2290911920.0_r8*x12&
             +1357702779072.0_r8*x10+92884281238528.0_r8*x8+854940668407808.0_r8*x6&
              +1778949256445952.0_r8*x4+903564402884608.0_r8*x2+82471672610816.0_r8)/s15;
        bn(16)=-(1.0_r8/19839816510008721408000.0_r8)*(50531787.0_r8*x16-861894976.0_r8*x14&
             -4147032373504.0_r8*x12-434759541092352.0_r8*x10-5832986472079360.0_r8*x8&
             -18074981856182272.0_r8*x6-15014905154371584.0_r8*x4-2942638810464256.0_r8*x2&
             -62844313796608.0_r8)/s16;
        bn(17)=-(1.0_r8/2618855779321151225856000.0_r8)*x*(19449462373.0_r8*x16&
             -161684570624.0_r8*x14+67021045928704.0_r8*x12+10632633811566592.0_r8*x10&
             +204101191194664960.0_r8*x8+910188421801050112.0_r8*x6+1155710005185347584.0_r8*x4&
             +405671965785849856.0_r8*x2+27625990880493568.0_r8)/s17;
        bn(18)=-(1.0_r8/877735702997277044857896960000.0_r8)*(-46949081169401.0_r8*x18&
             +329626742029560.0_r8*x16-2545450046209384704.0_r8*x14-619275164999346444288.0_r8*x12&
             -16758930124749577347072.0_r8*x10-104795615204228060872704.0_r8*x8&
             -193967804415750297354240.0_r8*x6-109854268161822321278976.0_r8*x4&
             -15914370046613878996992.0_r8*x2-265167190484014071808.0_r8)/s18;
        NN=18
      ELSEIF (mu==1) THEN

        bn(9)=(1.0_r8/22649241600.0_r8)*x*(x8+2496.0_r8*x6+103296.0_r8*x4&
             +290816.0_r8*x2+86016.0_r8)/s9
        bn(10)=(1.0_r8/930128855040.0_r8)*(17.0_r8*x10-9864.0_r8*x8&
             -841088.0_r8*x6-4398080.0_r8*x4-3010560.0_r8*x2-163840.0_r8)/s10
        bn(11)=(1.0_r8/744103084032.0_r8)*x*(19.0_r8*x10+576.0_r8*x8&
             +119936.0_r8*x6+1089536.0_r8*x4+1437696.0_r8*x2+262144.0_r8)/s11
        bn(12)=(1.0_r8/1218840851644416000.0_r8)*(935917.0_r8*x12&
             -94321392.0_r8*x10-31752534720.0_r8*x8-482370160640.0_r8*x6&
             -1104539013120.0_r8*x4-449675722752.0_r8*x2-16881287168.0_r8)/s12
        bn(13)=(1.0_r8/43878270659198976000.0_r8)*x*(-20287103.0_r8*x12&
             +368883168.0_r8*x10+170216971200.0_r8*x8+4192944148480.0_r8*x6&
             +15584536473600.0_r8*x4+11804184281088.0_r8*x2+1452354568192.0_r8)/s13
        bn(14)=(1.0_r8/210615699164155084800.0_r8)*(-2452337.0_r8*x14&
             -62553640.0_r8*x12-113403269184.0_r8*x10-4414009697792.0_r8*x8&
             -25432626982912.0_r8*x6-32278932848640.0_r8*x4-8676144054272.0_r8*x2&
             -237500366848.0_r8)/s14
        bn(15)=(1.0_r8/3647025093751603200.0_r8)*x*(35371.0_r8*x14&
             -176736.0_r8*x12+254656320.0_r8*x10+15417456640.0_r8*x8&
             +133153689600.0_r8*x6+264630042624.0_r8*x4+129530331136.0_r8*x2&
             +11450449920.0_r8)/s15
        bn(16)=(1.0_r8/654713944830287806464000.0_r8)*(141886453.0_r8*x16&
             -2026320704.0_r8*x14-5525282752256.0_r8*x12-518298332622848.0_r8*x10&
             -6542031953469440.0_r8*x8-19418164349370368.0_r8*x6-15587804196110336.0_r8*x4&
             -2967073282064384.0_r8*x2-61736581332992.0_r8)/s16
        bn(17)=(1.0_r8/30553317425413430968320000.0_r8)*x*&
             (-6293972187.0_r8*x16+55946317696.0_r8*x14+29304606351104.0_r8*x12&
             +4234306656018432.0_r8*x10+76628197541109760.0_r8*x8&
             +328094037924380672.0_r8*x6+403479956542193664.0_r8*x4&
             +137860560984211456.0_r8*x2+9167505905942528.0_r8)/s17
        bn(18)=(1.0_r8/13503626199958108382429184000.0_r8)*(-55660769987.0_r8*x18&
             +414280971000.0_r8*x16-1411977483166464.0_r8*x14-307978396577273856.0_r8*x12&
             -7874362709012717568.0_r8*x10-47365221274598965248.0_r8*x8&
             -85066924445175644160.0_r8*x6-46981216234673012736.0_r8*x4&
             -6658073400730189824.0_r8*x2-108757184552108032.0_r8)/s18
       NN=18 
     ELSE
       NN=8
     ENDIF
     sumt=0.0_r8; 
     DO n=0,NN
        T= cos(pi*(n-mu-0.5_r8)*0.5_r8+beta*tau)*pochham(mu+0.5_r8,n)*bn(n)/xpowy(tau,(n+mu+0.5_r8))
       sumt= sumt+T 
     ENDDO
     pmmu=sqrt(2.0_r8/pi/sinh(beta))*sumt
     fac2=1.0_r8  
     IF (mu>0) THEN
       tau2=tau*tau
       DO n=1,mu
         fac2=((mu-n+0.5_r8)*(mu-n+0.5_r8)+tau2)*fac2
       ENDDO
     ENDIF
     pmmu=fac2*pmmu
     END SUBROUTINE tauex1
    

     SUBROUTINE pmuint(x,mu,tau,pmmu,ierr)
    ! --------------------------------------------------
    ! Calculation of P^(mu)_(-1/2+i*tau)(x) by using 
    ! numerical quadrature, -1<x<1. 
    ! --------------------------------------------------
    ! Inputs:
    !   x ,    argument of the functions
    !   mu,    function parameter (integer value)
    !   tau,   function parameter (real value)
    ! Outputs:
    !   pmmu , P^(mu)_(-1/2+i*tau)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! ---------------------------------------------------------
     USE Someconstants
     REAL(r8), INTENT(IN) :: x
     INTEGER, INTENT(IN) :: mu
     REAL(r8), INTENT(IN) :: tau
     REAL(r8), INTENT(OUT) :: pmmu
     INTEGER,  INTENT(OUT) :: ierr
     real(r8) :: x2, tau2, muplushalf, beta, beta2, w, p, v, s0, sum
     COMMON/param/x2,beta,beta2,p,muplushalf	
     ierr=0
     x2= x*x; muplushalf= mu+0.5_r8; 
     beta= tau/muplushalf;
     tau2=tau*tau
     beta2= beta*beta;
     w=sqrt(1.0_r8+beta2*(1.0_r8-x2)); p= x/w; 
     v=4.0_r8*(1.0_r8+beta2)/(1.0_r8+p); 
     v=(exp(40.0_r8*2.31_r8/muplushalf)-1.0_r8)/v;
     s0=2.0_r8*log(sqrt(v)+sqrt(v+1.0_r8));
     CALL trapz(-s0, s0, sum)
     s0=acos(w*(1.0_r8-p*beta2)/(1.0_r8+beta2));
     pmmu= sum*exp(-tau*s0-muplushalf*log(w*(p+1.0_r8)/(1.0_r8+beta2)))/sqrttwopi;
     pmmu=cosh(pi*tau)*gammah(mu)*xpowy(1.0_r8-x2,mu/2.0_r8)*pmmu/pi
     END SUBROUTINE pmuint

     SUBROUTINE trapz(a,b,ti)
!-------------------------------------------------
! Implementation of an adaptative trapezoidal rule
!   inputs:
!          a, b: integration limits
!   output:
!         ti,  contribution 
!---------------------------------------
     USE Someconstants
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: b 
     REAL(r8), INTENT(OUT) :: ti
     REAL(r8) :: aa, bb, h, eps, delta, sum, xac, tin
     INTEGER :: n, ifi, i
     eps=epss
     IF (eps < 1.e-14_r8) eps=1.e-14_r8
     n=0
     h=b-a
     aa=a
     bb=b
     ti=0.5_r8*h*(fint(aa)+fint(bb))
     delta=1.0_r8+eps
11   n=n+1
     h=0.5_r8*h
     IF (n == 1) THEN
       ifi=1
     ELSE
       ifi=2*ifi
     END IF
     sum=0.0_r8
     DO  i=1,ifi
       xac=a+(2*i-1)*h
       sum=sum+fint(xac)
     ENDDO
     tin=0.5D0*ti+h*sum
     IF ((tin /= 0).AND.(n > 4)) THEN
       delta=ABS(1.d0-ti/tin)
     ENDIF
     ti=tin
     IF ((delta > eps).AND.(n < 20)) GO TO 11
     RETURN
     END SUBROUTINE trapz

     FUNCTION fint(s)
!------------------------------------
!   input:                                     
!      s,     integration abcissa              
!   output:                                    
!      fint,  function                         
!---------------------------------------------
     REAL(r8), INTENT(IN OUT)         :: s
     REAL(r8) :: fint, u, v, s1, es1, es1i, sh, sh2, ch, re, im
     REAL(r8) :: x2, beta, beta2, p, muplushalf
     COMMON/param/x2,beta,beta2,p,muplushalf
     IF (s==0) THEN 
       fint= 1 
     ELSE
       s1=s*0.5_r8
       es1=exp(s1); es1i=1.0_r8/es1;
       sh=(es1-es1i)*0.5_r8;
       ch=(es1+es1i)*0.5_r8;
       sh2= sh*sh;
       u=4.0_r8*sh2*(1.0_r8+beta2)*(1.0_r8+(1.0_r8+p*p*beta2)*sh2/(1.0_r8+p))/(1.0_r8+p);
       v= 2.0_r8*beta*(1.0_r8+p)*sh*ch/(p+1.0_r8+2.0_r8*sh2*(1.0_r8-p*beta2));
       re= 0.5_r8*log(1.0_r8+u);
       im= atan(v)-beta*s;
       fint= exp(-muplushalf*re)*cos(muplushalf*im)
     ENDIF
     END FUNCTION fint
 END MODULE Conical
