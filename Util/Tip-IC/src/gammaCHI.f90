  MODULE GammaCHI
  IMPLICIT NONE
  INTEGER, PARAMETER  :: r8 = KIND(0.0d0)
  PRIVATE
  PUBLIC  :: cdfgamC, invcdfgamC, cdfgamNC, invcdfgamNC, &
             errorfunction, inverfc, gamma, loggam,  &
             gamstar, quotgamm, recipgam, startijbes, hypfun
  CONTAINS 
     SUBROUTINE cdfgamC(ichi,ai,xi,p,q,ierr)
    ! ---------------------------------------------------------------
    ! Calculation of the central gamma (or chi-square distribution)
    ! functions P(a,x) and Q(a,x). 
    !
    ! The parameters a and x are related to the input arguments
    ! ai and xi as follows:
    !
    !      a=ai,   x=xi    for central gamma distributions.
    !      a=ai/2, x=xi/2  for central chi-square distributions.
    !                      The parameter ai corresponds to the degrees 
    !                      of freedom of the distribution function.
    !
    ! The parameter xi should be zero or a positive real number.
    !
    ! The parameter ai should be a positive real number. In order
    ! to avoid underflow/overflow problems, the admissible values
    ! of ai should be greater than 1.e-15.
    !
    ! The aimed relative accuracy is close to 1.0e-13 in
    ! most cases.  
    ! ---------------------------------------------------------------
    ! Inputs:
    !   ichi,  flag for the choice of gamma or 
    !          chi-square distribution function:
    !          ichi=1, computes the gamma distribution function.
    !          ichi=2, computes the central chi-square
    !                  distribution function. 
    !   ai ,    argument of the functions
    !   xi ,    argument of the functions
    ! Outputs:
    !   p,     function P(a,x)
    !   q,     function Q(a,x)  
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow/underflow problems. The function values 
    !                  (P(a,x) and Q(a,x)) are set to zero.
    !          ierr=2, any of the arguments of the function is 
    !                  out of range. The function values  
    !                  (P(a,x) and Q(a,x)) are set to zero.
    ! ----------------------------------------------------------------------
    ! Authors:
    !  Amparo Gil    (U. Cantabria, Santander, Spain)
    !                 e-mail: amparo.gil@unican.es
    !  Javier Segura (U. Cantabria, Santander, Spain)
    !                 e-mail: javier.segura@unican.es
    !  Nico M. Temme (CWI, Amsterdam, The Netherlands)
    !                 e-mail: nico.temme@cwi.nl
    ! -------------------------------------------------------------
    !  References: 
    !  1. A. Gil, J. Segura and N.M. Temme, GammaCHI: a package 
    !     for the inversion and computation of the gamma and 
    !     chi-square distribution functions (central and noncentral).
    !     Computer Physics Commun.
    !  2. "Efficient and accurate algorithms for the computation 
    !     and inversion of the incomplete gamma function ratios",    
    !     A. Gil, J. Segura and N.M. Temme, SIAM J Sci Comput (2012)
    ! -------------------------------------------------------------------
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: ai
    REAL(r8), INTENT(IN) :: xi
    INTEGER, INTENT(IN) :: ichi
    REAL(r8), INTENT(OUT) :: p
    REAL(r8), INTENT(OUT) :: q   
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: a, x, asmall, lnx, dp   
    ierr=0
    asmall=1.0e-15_r8
    IF (ichi==2) THEN
    ! chi-square distribution variables
      a=ai*0.5_r8
      x=xi*0.5_r8
    ELSE
      a=ai
      x=xi
    ENDIF
    IF (x<0.0_r8) THEN
      ierr=2
      p=0.0_r8
      q=0.0_r8
    ENDIF
    IF (a<asmall) THEN
      ierr=2
      p=0.0_r8
      q=0.0_r8
    ENDIF
    IF (ierr==0) THEN
      IF (x<dwarf) THEN
        lnx=log(dwarf)
      ELSE
        lnx=log(x)
      ENDIF
      IF (a>alfa(x)) THEN  
        dp=dompart(a,x,.false.);
        IF (dp<0) THEN
          ierr=1
          p=0; q=0;
        ELSE
          IF ((x <0.34_r8*a).OR.(a<16.5_r8)) THEN
            p=ptaylor(a,x,dp)
          ELSE
            p=pqasymp(a,x,dp,.true.)
          ENDIF
          q=1.0_r8-p
        ENDIF
      ELSE
        IF (a<-dwarf/lnx) THEN
          q=0.0_r8
        ELSE
          IF (x<1.0_r8) THEN
            dp=dompart(a,x,.true.)
            IF (dp<0) THEN
              ierr=1
              q=0; p=0;
            ELSE
              q=qtaylor(a,x,dp)
              p=1.0_r8-q
            ENDIF
          ELSE
            dp=dompart(a,x,.false.);
            IF (dp<0) THEN
              ierr=1
              p=0; q=0;
            ELSE 
              IF ((x>1.5_r8*a).OR.(a<12.0_r8)) THEN
                q=qfraction(a,x,dp)
              ELSE
                q=pqasymp(a,x,dp,.false.)
                IF (dp==0.0_r8) THEN
                  q=0.0_r8
                ENDIF
              ENDIF
              p=1.0_r8-q
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    END SUBROUTINE cdfgamC

    SUBROUTINE invcdfgamC(ichi,ai,p,q,xr,ierr)
    ! --------------------------------------------------------------
    !     Inversion routine for the central gamma distribution
    !---------------------------------------------------------------
    ! This routine computes xr in the equations P(a,xr)=p 
    ! and Q(a,xr)=q with a as a given positive parameter;
    ! p and q satisfy p+q=1. 
    ! We invert the equation with min(p,q).
    !
    !  The parameter a is related to the input argument ai as:
    !         a=ai,     for central gamma distributions
    !         a=ai/2,   for central chi-square distributions.
    !                   The parameter ai corresponds to the degrees 
    !                   of freedom of the distribution function.
    ! 
    !  The parameter ai should be a positive real number. In order
    !  to avoid underflow/overflow problems, the admissible values
    !  of a and p (or q) should be greater than 1.e-15 and 1.e-250, 
    !  respectively.
    ! --------------------------------------------------------------
    ! Inputs:
    !   ichi,  flag for the choice of gamma or 
    !          chi-square distribution function:
    !          ichi=1, computes the inverse of the central
    !                  gamma distribution function.
    !          ichi=2, computes the inverse of the central 
    !                  chi-square distribution function. 
    !   ai,    argument of the functions
    !   p,     function value P(a,x)
    !   q,     function value Q(a,x)  
    ! Outputs:
    !   xr   , solution of the equations P(a,xr)=p and Q(a,xr)=q
    !          with a as a given positive parameter.
    !   ierr , error flag
    !          ierr=0,  computation succesful
    !          ierr=1,  overflow problems in one or more steps of the 
    !                   computation.
    !                   Either the value zero or the initial approximation 
    !                   to the root is given as output.
    !          ierr=2,  the number of iterations in the Newton method
    !                   reached the upper limit N=15. The last value
    !                   obtained for the root is given as output.
    !          ierr=3,  any of the arguments of the function is 
    !                   out of range. The value 0 is given as output. 
    ! ----------------------------------------------------------------------
    ! The claimed accuracy obtained using this inversion routine is near
    ! 1.e-12. 
    ! ----------------------------------------------------------------------
    !           METHODS OF COMPUTATION
    ! ----------------------------------------------------------------------
    ! The present code uses different methods of computation
    ! depending on the values of the input values: Taylor, asymptotic 
    ! expansions and high-order Newton methods.
    ! ---------------------------------------------------------------------- 
    ! Authors:
    ! 
    !  Amparo Gil    (U. Cantabria, Santander, Spain)
    !                 e-mail: amparo.gil@unican.es
    !  Javier Segura (U. Cantabria, Santander, Spain)
    !                 e-mail: javier.segura@unican.es
    !  Nico M. Temme (CWI, Amsterdam, The Netherlands)
    !                 e-mail: nico.temme@cwi.nl
    ! -------------------------------------------------------------------
    !  References: 
    !  1. A. Gil, J. Segura and N.M. Temme, GammaCHI: a package 
    !     for the inversion and computation of the gamma and 
    !     chi-square distribution functions (central and noncentral).
    !     Computer Physics Commun
    !  2. A. Gil, J. Segura and N.M. Temme. Efficient and accurate 
    !     algorithms for the computation and inversion 
    !     of the incomplete gamma function ratios. SIAM J Sci Comput.  
    !     (2012) 34(6), A2965-A2981
    ! ------------------------------------------------------------------
    USE Someconstants    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ichi
    REAL(r8), INTENT(IN) :: ai
    REAL(r8), INTENT(IN) :: p
    REAL(r8), INTENT(IN) :: q 
    REAL(r8), INTENT(OUT) :: xr
    REAL(r8) :: a, asmall, porq, s, dlnr, logr, r, a2, a3, a4, ap1,& 
                ap12, ap13, ap14, ap2, ap22, ap3, ainv, x0,&
                ck(1:5), b, eta, L, L2, L3, L4, b2, b3, x, &
                x2, t, px, qx, y, vmin, vdif, lgama, fp, invfp
    INTEGER,  INTENT(OUT) :: ierr
    INTEGER :: k, n, m, ierrf
    LOGICAL :: pcase
    ierr=0
    asmall=1.0e-15_r8
    IF (ichi==2) THEN
    ! chi-square distribution variables
      a=ai*0.5_r8      
    ELSE
      a=ai
    ENDIF 
    IF (a<asmall) THEN
      ierr=3
      xr=0.0_r8
    ENDIF
    IF ((p>1.0_r8).OR.(q>1.0_r8)) THEN
      ierr=3
      xr=0.0_r8
    ENDIF
    IF ((p<=1.e-250_r8).OR.(q<=1.e-250_r8)) THEN
      ierr=3
      xr=0.0_r8
    ENDIF
    IF (p+q>1.0_r8+1.e-8_r8) THEN
      ierr=3
      xr=0.0_r8
    ENDIF    
    IF (ierr==0) THEN
      IF (p<0.5_r8) THEN
        pcase=.true.
        porq=p
        s=-1 
      ELSE
        pcase=.false.
        porq=q 
        s=1 
      ENDIF
      logr=(1.0_r8/a)*(log(p)+loggam(a+1.0_r8))
      k=0
      IF (logr <log(0.2_r8*(1.0_r8+a))) THEN
        r=exp(logr)
        m=0
        a2=a*a
        a3=a2*a
        a4=a3*a
        ap1=a+1.0_r8
        ap12=ap1*ap1
        ap13=ap1*ap12
        ap14=ap12*ap12
        ap2=a+2.0_r8
        ap22=ap2*ap2
        ap3=a+3.0_r8
        ck(1)= 1.0_r8;
        ck(2)= 1.0_r8/ap1;
        ck(3)=0.5_r8*(3.0_r8*a+5.0_r8)/(ap12*ap2);
        ck(4)=0.333333333333333333333333333333_r8*(31.0_r8+8.0_r8*a2+33.0_r8*a)&
            /(ap13*ap2*ap3);
        ck(5)=0.0416666666666666666666666666667_r8*(2888.0_r8+1179.0_r8*a3+125.0_r8*a4+&
            3971.0_r8*a2+5661.0_r8*a)/(ap14*ap22*ap3*(a+4.0_r8)); 
        x0=r*(1.0_r8+r*(ck(2)+r*(ck(3)+r*(ck(4)+r*ck(5)))))
        lgama=loggam(a)
        k=1
      ENDIF
      IF ((a<10.0_r8).AND.(k==0)) THEN
        vmin=min(0.02_r8,exp(-1.5_r8*a)/gamma(a))
        vdif=abs(q-vmin)
        IF ((q<vmin).AND.(vdif>1.e-30_r8)) THEN
          m=0
          b=1.0_r8-a; 
          b2=b*b;
          b3=b2*b;
          eta=sqrt(-2.0_r8/a*log(q*gamstar(a)*sqrttwopi/sqrt(a)));
          x0=a*lambdaeta(eta);
          L=log(x0); 
          IF (x0>5.0_r8) THEN
            L2=L*L;
            L3=L2*L;
            L4=L3*L
            r=1.0_r8/x0;
            ck(1)=L-1.0_r8;
            ck(2)=(3*b-2*b*L+L2-2*L+2)*0.5_r8;
            ck(3)=(24*b*L-11*b2-24*b-6*L2+12*L-12-9*b*L2+6*b2*L+2*L3)&
                  *0.166666666666666666666666666667_r8;
            ck(4)=(-12*b3*L+84*b*L2-114*b2*L+72+36*L2+3*L4&
                   -72*L+162*b-168.0_r8*b*L-12*L3+25*b3&
                   -22*b*L3+36*b2*L2+120*b2)*0.0833333333333333333333333333333_r8;
            x0=x0-L+b*r*(ck(1)+r*(ck(2)+r*(ck(3)+r*ck(4))))
          ELSE
            r=1.0_r8/x0;
            L2=L*L;
            ck(1)=L-1.0_r8;
            IF ((L-b*r*ck(1))<x0) THEN
              x0=x0-L+b*r*ck(1);    
            ENDIF
          ENDIF
          lgama=loggam(a) 
          k=1
        ENDIF
      ENDIF
      IF ((abs(porq-0.5_r8)< 1.0e-5_r8).AND.(k==0)) THEN
        m=0
        ainv=1.0_r8/a;
        x0=a-0.333333333333333333333333333333_r8+&
           (0.0197530864197530864197530864198_r8 &
           +0.00721144424848128551832255535959_r8*ainv)*ainv;
        lgama=loggam(a)
        k=1
      ENDIF
      IF ((abs(a-1.0_r8)<1.0e-4_r8).AND.(k==0)) THEN
        m=0
        IF (pcase) THEN
          x0=-log(1.0_r8-p) 
        ELSE
          x0=-log(q) 
        ENDIF
        lgama=loggam(a)
        k=1
      ENDIF
      IF ((a<1.0_r8).AND.(k==0)) THEN
        m=0 
        IF (pcase) THEN 
          x0=exp((1.0_r8/a)*(log(porq)+loggam(a+1.0_r8)))
        ELSE
          x0=exp((1.0_r8/a)*(log(1.0_r8-porq)+loggam(a+1.0_r8)))
        ENDIF
        lgama=loggam(a)
        k=1
      ENDIF
      IF (k==0) THEN
        m=1  
        ainv=1.0_r8/a;
        r=inverfc(2*porq);
        eta=s*r/sqrt(a*0.5_r8);
        IF (r<giant) THEN
          eta=eta+(eps1(eta)+(eps2(eta)+eps3(eta)*ainv)*ainv)*ainv;
          x0=a*lambdaeta(eta)
          y=eta
          fp=-sqrt(a/twopi)*exp(-0.5_r8*a*y*y)/(gamstar(a));
          invfp=1.0_r8/fp;
        ELSE
          x0=0.0_r8
          xr=x0
          ierr=1
        ENDIF
      ENDIF
      x=x0;
      IF (ierr==0) THEN 
        t=1
        n=1;
        a2=a*a
        a3=a2*a
        ! Implementation of the high order Newton-like method
        DO WHILE ((t>1.0e-15_r8).AND.(n< 15))
          x=x0;
          x2=x*x
          IF (m==0) THEN
            dlnr=(1.0_r8-a)*log(x)+x+lgama;
            IF (dlnr>log(giant)) THEN
              n=20
              ierr=1
            ELSE
              r=exp(dlnr);
              IF (pcase) THEN 
                CALL cdfgamC(1,a,x,px,qx,ierrf)
                ck(1)=-r*(px-p);  
              ELSE
                CALL cdfgamC(1,a,x,px,qx,ierrf)
                ck(1)=r*(qx-q);
              ENDIF   
              r=ck(1);
              IF (a>0.1_r8) THEN
                ck(2)=0.5_r8*(x-a+1.0_r8)/x;
                ck(3)=0.166666666666666666666666666667_r8*&
                      (2*x2-4*x*a+4*x+2*a2-3*a+1)/x2;
                x0=x+r*(1.0_r8+r*(ck(2)+r*ck(3)));
              ELSE  
                IF (a>0.05_r8) THEN
                  ck(2)=0.5_r8*(x-a+1.0_r8)/x;
                  x0=x+r*(1.0_r8+r*(ck(2)));
                ELSE
                  x0=x+r;
                ENDIF
              ENDIF 
            ENDIF       
          ELSE
            r=-invfp*x
            IF (pcase) THEN 
              CALL cdfgamC(1,a,x,px,qx,ierrf)
              ck(1)=-r*(px-p);  
            ELSE
              CALL cdfgamC(1,a,x,px,qx,ierrf)
              ck(1)=r*(qx-q);
            ENDIF
            r=ck(1);
            IF (a>0.1_r8) THEN
              ck(2)=0.5_r8*(x-a+1.0_r8)/x;
              ck(3)=0.166666666666666666666666666667_r8*&
                   (2.0_r8*x2-4.0_r8*x*a+4.0_r8*x+2.0_r8*a2-3.0_r8*a+1.0_r8)/x2;
              x0=x+r*(1.0_r8+r*(ck(2)+r*ck(3)));
            ELSE  
              IF (a>0.05_r8) THEN
                ck(2)=0.5_r8*(x-a+1.0_r8)/x;
                x0=x+r*(1.0_r8+r*(ck(2)));
              ELSE
                x0=x+r;
              ENDIF
            ENDIF 
          ENDIF
          t=abs(x/x0-1.0_r8);
          n=n+1; 
          x=x0
        ENDDO
        IF (n==15) ierr=2
        xr=x
      ENDIF
    ENDIF
    IF (ichi==2) xr=xr*2.0_r8
    END SUBROUTINE invcdfgamC

    SUBROUTINE cdfgamNC(ichi,mui,xi,yi,p,q,ierr)
    ! ----------------------------------------------------------------
    ! Calculation of the noncentral gamma or noncentral
    !   chi-square distribution functions P_mu(x,y) 
    !   and Q_mu(x,y).
    ! 
    ! The parameters mu, x and y are related to the input 
    ! arguments mui, xi and yi as follows:
    !
    !   mu=mui,   x=xi,   y=yi    for noncentral gamma distributions.
    !
    !   mu=mui/2, x=xi/2, y=yi/2  for noncentral chi-square
    !                             distributions. The parameter
    !                             mui corresponds to the degrees
    !                             of freedom; xi is the noncentrality
    !                             parameter; yi is the limit of 
    !                             integration.
    !                                                      
    ! In order to avoid, overflow/underflow problems in IEEE double
    !  precision arithmetic, the admissible parameter ranges 
    !  for computation are:
    !
    !        0<=x<=10000,   0<=y<=10000,   0.5<=mu<=10000 
    !
    !  The aimed relative accuracy is close to 1.0e-11 in the 
    !  previous parameter domain.  
    ! -------------------------------------------------------------
    ! Inputs:
    !   ichi,  flag for the choice of gamma or 
    !          chi-square distribution function:
    !          ichi=1, computes the noncentral gamma 
    !                  distribution function.
    !          ichi=2, computes the noncentral chi-square
    !                  distribution function.  
    !   mui ,   argument of the functions
    !   xi ,    argument of the functions
    !   yi ,    argument of the functions
    ! Outputs:
    !   p,     function P_mu(x,y)
    !   q,     function Q_mu(x,y)  
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, Underflow problems. The function values 
    !                  are set to zero and one.
    !          ierr=2, any of the arguments of the function is 
    !                  out of range. The function values 
    !                  (P_mu(x,y) and Q_mu(x,y)) are set to zero.
    ! --------------------------------------------------------------------
    !           METHODS OF COMPUTATION
    ! --------------------------------------------------------------------
    ! The present code uses different methods of computation
    ! depending on the values of mu, x and y: series expansions in terms
    ! of incomplete gamma functions, integral representations, 
    ! asymptotic expansions, and use of three-term homogeneous 
    ! recurrence relations.
    ! -------------------------------------------------------------------- 
    ! Authors:
    ! 
    !  Amparo Gil    (U. Cantabria, Santander, Spain)
    !                 e-mail: amparo.gil@unican.es
    !  Javier Segura (U. Cantabria, Santander, Spain)
    !                 e-mail: javier.segura@unican.es
    !  Nico M. Temme (CWI, Amsterdam, The Netherlands)
    !                 e-mail: nico.temme@cwi.nl
    ! ---------------------------------------------------------------
    !  References: 
    !  1. A. Gil, J. Segura and N.M. Temme, GammaCHI: a package 
    !     for the inversion and computation of the gamma and 
    !     chi-square distribution functions (central and noncentral).
    !     Computer Physics Commun
    !  2. A. Gil, J. Segura and N.M. Temme, Computation of the Marcum 
    !      Q-function. ACM Trans Math Soft (2014)
    !  3. A. Gil, J. Segura and N.M. Temme. Efficient and accurate 
    !     algorithms for the computation and inversion 
    !     of the incomplete gamma function ratios. SIAM J Sci Comput.  
    !     (2012) 34(6), A2965-A2981
    ! ---------------------------------------------------------------
    USE Someconstants    
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: ichi   
    REAL(r8), INTENT(IN) :: mui
    REAL(r8), INTENT(IN) :: xi
    REAL(r8), INTENT(IN) :: yi      
    REAL(r8), INTENT(OUT) :: p
    REAL(r8), INTENT(OUT) :: q   
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: mu, x ,y
    REAL(r8) :: b,w,xii,y0,y1,mulim
    ierr=0
    IF (ichi==2) THEN
    ! chi-square distribution variables
      mu=mui*0.5_r8
      x=xi*0.5_r8
      y=yi*0.5_r8
    ELSE
      mu=mui
      x=xi
      y=yi
    ENDIF
    IF (((x>10000.0_r8).OR.(y>10000.0_r8)).OR.(mu>10000.0_r8)) THEN
      p=0.0_r8
      q=0.0_r8
      ierr=2
    ENDIF
    IF (((x<0.0_r8).OR.(y<0.0_r8)).OR.(mu<0.5_r8)) THEN
      p=0.0_r8
      q=0.0_r8
      ierr=2
    ENDIF
    IF (ierr==0) THEN
      mulim=135.0_r8
      b=1.0_r8
      w=b*sqrt(4*x+2*mu); xii=2*sqrt(x*y);
      y0=x+mu-w; y1=x+mu+w; 
      IF ((y >x+mu).AND.(x<30)) THEN
      ! Series for Q in terms of ratios of Gamma functions 
        CALL qser(mu,x,y,p,q,ierr)
      ELSEIF ((y <= x+mu).AND.(x < 30)) THEN
      ! Series for P in terms of ratios of Gamma functions 
        CALL pser(mu,x,y,p,q,ierr)
      ELSEIF  ((mu*mu<2*xii).AND.(xii>30)) THEN
      ! Asymptotic expansion for xy large 
        CALL pqasyxy(mu,x,y,p,q,ierr); 
      ELSEIF (((mu>=mulim).AND.(y0 <= y)).AND.(y <= y1)) THEN
      ! Asymptotic expansion for mu large
        CALL pqasymu(mu,x,y,p,q,ierr)
      ELSEIF (((y <= y1).AND.(y > x+mu)).AND.(mu < mulim)) THEN
      ! Recurrence relation for Q
        CALL qrec(mu,x,y,p,q,ierr)  
      ELSEIF (((y >= y0).AND.(y <= x+mu)).AND.(mu < mulim)) THEN
      ! Recurrence relation for P
        CALL prec(mu,x,y,p,q,ierr)
      ELSE
        IF (mu<=1) THEN
          ! One step of the Recurrence Relation
          CALL rec1step(mu,x,y,p,q,ierr)
        ELSE
          ! Integral representation
          CALL cdfPQtrap(mu,x,y,p,q,ierr)
        ENDIF
      ENDIF
    ENDIF
    IF (ierr==0) THEN
      IF (p<1.e-290_r8) THEN 
        p=0.0_r8; q=1.0_r8;
        ierr=1
      ENDIF
      IF (q<1.e-290_r8) THEN
        q=0.0_r8; p=1.0_r8;
        ierr=1         
      ENDIF  
    ENDIF
    END SUBROUTINE cdfgamNC 

    SUBROUTINE invcdfgamNC(ichi,icho,mui,p,q,yxi,xy,ierr)
    ! ----------------------------------------------------------------
    !     Inversion routine for the noncentral gamma/chi-square 
    !             cumulative distribution function                
    !-----------------------------------------------------------------
    ! This routine computes, depending on the value of
    ! icho:
    !    a)  x in the equations Q_mu(x,y)=q, P_mu(x,y)=p
    !    b)  y in the equations P_mu(x,y)=p, Q_mu(x,y)=q
    !         with (mu,y,p,q) or (mu,x,p,q) as given parameters,
    !         respectively.
    ! In a) and b), we invert the equation with min(p,q).
    !
    ! The parameters mu and x (or y) are related to the input 
    ! arguments mui and yxi as follows:
    !
    !   mu=mui,   x=yxi    (or y=yxi)   for noncentral gamma cumulative 
    !                                   distributions
    !
    !   mu=mui/2, x=yxi/2  (or y=yxi/2) for noncentral chi-square 
    !                                   cumulative distributions
    !                                   The parameter mui corresponds 
    !                                   to the degrees of freedom; yxi is 
    !                                   either the noncentrality
    !                                   parameter or the limit of 
    !                                   integration.
    !
    ! In statistics, the inversion of  Q_mu(x,y)=q with respect 
    ! to x corresponds to the problem of inverting the noncentral
    ! gamma (or chi-square) distribution function with respect 
    ! to the non-centrality parameter given the upper tail 
    ! probability.
    !
    ! On the other hand, the inversion of P_mu(x,y)=p with respect 
    ! to y with fixed x corresponds to the problem of computing the 
    ! p-quantiles of the distribution function.
    !
    !  The admissible parameter ranges for computation are:
    !  
    !  a) 1.e-20<p<1,  1.e-20<q<1,  p+q=1
    !
    !  b) Inversion with respect to x:
    !
    !      0<=y<=10000,    0.5<=mu<=10000,  
    !              
    !  c) Inversion with respect to y:
    !
    !      0<=x<=10000,    0.5<=mu<=10000.
    !              
    ! ----------------------------------------------------------------
    ! Inputs:
    !   ichi,  flag for the choice of gamma or 
    !          chi-square distribution function:
    !          ichi=1, computes the inverse of the noncentral 
    !                  gamma distribution function.
    !          ichi=2, computes the inverse of the noncentral 
    !                  chi-square distribution function.
    !   icho,  choice of the inversion process.
    !          If icho=1, x is computed in the
    !                     equations Q_mu(x,y)=q, P_mu(x,y)=p.
    !          If icho=2, y is computed in the
    !                     equations Q_mu(x,y)=q, P_mu(x,y)=p.
    !   mui ,  parameter of the function 
    !   p,     function value P_mu(x,y)
    !
    !   q,     function value Q_mu(x,y)
    !                
    !   yxi,   If icho=1, yxi is y in the
    !                     equations P_mu(x,y)=p, Q_mu(x,y)=q.
    !          If icho=2, yxi is x in the 
    !                     equation P_mu(x,y)=p, Q_mu(x,y)=q.
    ! Outputs:
    !   xy,    If icho=1, xy is x in the
    !                     equation  P_mu(x,y)=p, Q_mu(x,y)=q.
    !          If icho=2, xy is y in the 
    !                     equation P_mu(x,y)=p, Q_mu(x,y)=q. 
    !   ierr , error flag
    !          ierr=0, computation succesful.
    !          ierr=1, in the inversion process with respect to the 
    !                  non-centrality parameter (x-variable), the input 
    !                  values (yx,qp) do not make sense (i.e. the 
    !                  inequality qp<Q_mu(0,yx) holds). The value 0 
    !                  is given as output.
    !          ierr=2, at least one of the gamma distribution 
    !                  function values needed in the inversion process
    !                  cannot be correctly computed. Either 0 or the
    !                  last value obtained for the root are given
    !                  as output.   
    !          ierr=3, the number of iterations in the secant method
    !                  reached the upper limit N=25. The last value
    !                  obtained for the root is given as output.
    !          ierr=4, any of the arguments of the function is 
    !                  out of range. The value 0 is given as output. 
    !----------------------------------------------------------------------
    ! The claimed accuracy obtained using this inversion routine is close
    ! to the value of eps.
    ! ----------------------------------------------------------------------
    !           METHODS OF COMPUTATION
    ! ----------------------------------------------------------------------
    ! The present code uses different methods of computation
    ! depending on the values of the input values: asymptotic expansions 
    ! and/or secant method.
    ! ---------------------------------------------------------------------- 
    ! Authors:
    ! 
    !  Amparo Gil    (U. Cantabria, Santander, Spain)
    !                 e-mail: amparo.gil@unican.es
    !  Javier Segura (U. Cantabria, Santander, Spain)
    !                 e-mail: javier.segura@unican.es
    !  Nico M. Temme (CWI, Amsterdam, The Netherlands)
    !                 e-mail: nico.temme@cwi.nl
    ! -------------------------------------------------------------------
    !  References: 
    !  1. A. Gil, J. Segura and N.M. Temme, GammaCHI: a package 
    !     for the inversion and computation of the gamma and 
    !     chi-square distribution functions (central and noncentral).
    !     Computer Physics Commun
    !  2. A. Gil, J. Segura and N.M. Temme, The asymptotic and numerical 
    !     inversion of the Marcum Q-function. Studies in Applied
    !     Mathematics (2014)
    !  3. A. Gil, J. Segura and N.M. Temme, Computation of the Marcum 
    !      Q-function. ACM Trans Math Soft (2014)
    !  4. A. Gil, J. Segura and N.M. Temme. Efficient and accurate 
    !     algorithms for the computation and inversion 
    !     of the incomplete gamma function ratios. SIAM J Sci Comput.  
    !     (2012) 34(6), A2965-A2981
    ! ------------------------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ichi
    INTEGER, INTENT(IN) :: icho
    REAL(r8), INTENT(IN) :: mui
    REAL(r8), INTENT(IN) :: q
    REAL(r8), INTENT(IN) :: p
    REAL(r8), INTENT(IN) :: yxi
    REAL(r8), INTENT(OUT) :: xy  
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: eps,  mu, mu2, mu3, yx, y0, x0, p0, q0, porq, s, &
                qi, y1, y2, zeta0, zeta1, x1, x2, x, y, xs, xn, yn,&
                delta, q1, q2, p1, p2, zeta2, zeta3  
    INTEGER :: icon, ierr1, ierr2, ierr3
    LOGICAL :: pcase
    ierr=0
    eps=1.e-12_r8
    IF ((p<1.e-20_r8).OR.(q<1.e-20_r8)) THEN
      xy=0.0_r8
      ierr=4
    ENDIF
    IF (ichi==2) THEN
    ! chi-square distribution variables
      mu=mui*0.5_r8      
      yx=yxi*0.5_r8
    ELSE
      mu=mui
      yx=yxi
    ENDIF   
    IF ((yx>10000).OR.(yx<0.0_r8)) THEN
      xy=0.0_r8      
      ierr=4
    ENDIF    
    IF ((mu>10000).OR.(mu<0.5_r8)) THEN
      xy=0.0_r8      
      ierr=4
    ENDIF
    IF ((icho==1).AND.(ierr==0)) THEN
      IF (q>1.0_r8+1.e-30_r8) THEN
        xy=giant
      ELSE
        IF (p<0.5_r8) THEN
          pcase=.true.
          porq=p
          s=-1 
        ELSE
          pcase=.false.
          porq=q 
          s=1 
        ENDIF
        qi=q
        ! ------------------------------------------------------
        ! Inversion with respect to the x variable in the
        ! equation Q_mu(x,y)=q  
        ! Calculation of the noncentrality parameter.
        !-------------------------------------------------------
        y0=yx
        ! We have first to check that qi >=q_mu(0,y)
        CALL cdfgamC(1,mu,y0,p0,q0,ierr3)
        IF (ierr3==0) THEN
          IF (qi < q0) THEN
            xy=0.0_r8
            ierr=1 
          ELSE
            IF (mu<5.0_r8) THEN
              zeta0=s*inverfc(2.0_r8*porq)*sqrt(2.0_r8);   
              x=xetay(zeta0,mu,yx);
              IF (sqrt(x*yx)>10.0_r8) THEN
                mu2=mu*mu
                mu3=mu2*mu
                zeta0=zeta0/sqrt(mu);
                y=yx/mu;
                x=x/mu;
                CALL zeta3xy(icho,zeta0,x,y,zeta1,zeta2,zeta3)        
                zeta0=zeta0+zeta1/mu+zeta2/mu2+zeta3/mu3
                x=mu*xzetay(zeta0,y);  
              ENDIF
              IF (x>0.0_r8) THEN
                x1=x
                x2=x+0.01_r8
              ELSE
                x1=1.e-6_r8
                x2=1.e-1_r8
              ENDIF
              CALL cdfgamNC(1,mu,x1,y0,p1,q1,ierr1)
              CALL cdfgamNC(1,mu,x2,y0,p2,q2,ierr2)
              IF ((ierr1==0).AND.(ierr2==0)) THEN 
                delta=1.0_r8
                icon=0
                DO WHILE ((delta>eps).AND.(icon<25))
                  icon=icon+1
                  IF (pcase) THEN
                    xn=x2-(x2-x1)*(p2-p)/(p2-p1)
                  ELSE
                    xn=x2-(x2-x1)*(q2-q)/(q2-q1)
                  ENDIF 
                  IF (xn<0) xn=1.e-10_r8
                  x1=x2
                  q1=q2
                  p1=p2
                  x2=xn
                  CALL cdfgamNC(1,mu,x2,y0, p2,q2,ierr1)
                  IF (ierr1==0) THEN
                    IF (pcase) THEN
                      delta=abs(1.0_r8-p2/p)
                      IF(p1-p2==0.0_r8) delta=0.0_r8
                    ELSE
                      delta=abs(1.0_r8-q2/q)
                      IF(q1-q2==0.0_r8) delta=0.0_r8
                    ENDIF
                  ELSE
                    ierr=2
                    delta=0.0_r8 
                  ENDIF
                ENDDO
                xy=xn
                IF (icon==25) ierr=3
              ELSE
                xy=0.0_r8
                ierr=2
              ENDIF    
            ELSE
            !-----------------
            ! Asymptotics
            !-----------------
              mu2=mu*mu
              mu3=mu2*mu
              zeta0=s*inverfc(2.0_r8*porq)*sqrt(2.0_r8/mu); 
              y=y0/mu;
              x=xzetay(zeta0, y);
              CALL zeta3xy(icho,zeta0,x,y,zeta1,zeta2,zeta3)          
              zeta0=zeta0+zeta1/mu+zeta2/mu2+zeta3/mu3
              x=mu*xzetay(zeta0,y);
              xn=x
              IF (x>0.0_r8) THEN            
                x1=x 
                x2=x+0.01_r8
              ELSE
                x1=1.e-6_r8
                x2=1.e-1_r8
              ENDIF
              CALL cdfgamNC(1,mu,x1,y0, p1,q1,ierr1)
              CALL cdfgamNC(1,mu,x2,y0, p2,q2,ierr2)
              IF ((ierr1==0).AND.(ierr2==0)) THEN
                delta=1.0_r8
                icon=0
                DO WHILE ((delta>eps).AND.(icon<25))
                  icon=icon+1
                  IF (pcase) THEN
                    xn=x2-(x2-x1)*(p2-p)/(p2-p1)
                  ELSE
                    xn=x2-(x2-x1)*(q2-q)/(q2-q1)
                  ENDIF 
                  IF (xn<0) xn=1.e-10_r8
                  x1=x2
                  q1=q2
                  p1=p2
                  x2=xn
                  CALL cdfgamNC(1,mu,x2,y0, p2,q2,ierr1)
                  IF (ierr1==0) THEN
                    IF (pcase) THEN
                      delta=abs(1.0_r8-p2/p)
                      IF(p1-p2==0.0_r8) delta=0.0_r8
                    ELSE
                      delta=abs(1.0_r8-q2/qi)
                      IF(q1-q2==0.0_r8) delta=0.0_r8
                    ENDIF
                  ELSE
                    ierr=2
                    delta=0.0_r8 
                  ENDIF
                ENDDO
                xy=xn
                IF (icon==25) ierr=3
              ELSE
                xy=xn
                ierr=2
              ENDIF
            ENDIF
            xy=xn        
          ENDIF  
        ELSE
          xy=0.0_r8
          ierr=2
        ENDIF
      ENDIF       
    ELSEIF (ierr==0) THEN
      ! --------------------------------------------
      ! Inversion with respect to the y variable 
      ! calculation of the noncentrality parameter)
      !---------------------------------------------
      IF (abs(q)<dwarf) THEN
        xy=giant
      ELSE
        IF (p<0.5_r8) THEN
          pcase=.true.
          porq=p
          s=-1 
        ELSE
          pcase=.false.
          porq=q 
          s=1 
        ENDIF
        x0=yx
        IF (mu<5) THEN
          zeta0=s*inverfc(2.0_r8*porq)*sqrt(2.0_r8); 
          y=yetax(zeta0,mu,yx);            
          IF (sqrt(x0*y)>10.0_r8) THEN
            mu2=mu*mu
            mu3=mu2*mu      
            zeta0=zeta0/sqrt(mu);
            y=y/mu;
            x=yx/mu;
            CALL zeta3xy(icho,zeta0,x,y,zeta1,zeta2,zeta3)    
            zeta0=zeta0+zeta1/mu+zeta2/mu2;
            y=mu*yzetax(zeta0,x);  
          ENDIF
          IF (y>0.0_r8) THEN
            y1=y
            y2=y+0.01_r8
          ELSE
            y1=1.e-6_r8
            y2=1.e-1_r8
          ENDIF
          CALL cdfgamNC(1,mu,x0,y1, p1,q1,ierr1)
          CALL cdfgamNC(1,mu,x0,y2, p2,q2,ierr2)
          IF ((ierr1==0).AND.(ierr2==0)) THEN
            delta=1.0_r8
            icon=0
            DO WHILE((delta>eps).AND.(icon<25))
              icon=icon+1
              IF (pcase) THEN
                yn=y2-(y2-y1)/(p2-p1)*(p2-p)
              ELSE
                yn=y2-(y2-y1)/(q2-q1)*(q2-q)
              ENDIF
              IF(yn<0) yn=1.e-10_r8
              y1=y2
              p1=p2
              q1=q2
              y2=yn
              CALL cdfgamNC(1,mu,x0,y2,p2,q2,ierr1)
              IF (ierr1==0) THEN
                IF (pcase) THEN
                  delta=abs(1.0_r8-p2/p)
                  IF(p1-p2==0.0_r8) delta=0.0_r8
                ELSE
                  delta=abs(1.0_r8-q2/q)
                  IF(q1-q2==0.0_r8) delta=0.0_r8
                ENDIF
              ELSE
                ierr=2
                delta=0.0_r8 
              ENDIF
            ENDDO
            IF (icon==25) ierr=3 
            xy=yn
          ELSE
            ierr=2
            xy=0.0_r8
          ENDIF
        ELSE 
          x=yx
          xs=x/mu;
          zeta0=s*inverfc(2.0_r8*porq)*sqrt(2.0_r8/mu);    
          y=yzetax(zeta0,xs);
          CALL zeta3xy(icho,zeta0,xs,y,zeta1,zeta2,zeta3)
          zeta0=zeta0+zeta1/mu+zeta2/(mu*mu); 
          y= mu*yzetax(zeta0, xs);
          yn=y;   
          IF (y>0.0_r8) THEN
            y1=y
            y2=y+0.01_r8
          ELSE
            y1=1.e-6_r8
            y2=1.e-1_r8
          ENDIF
          CALL cdfgamNC(1,mu,x,y1, p1,q1,ierr1)
          CALL cdfgamNC(1,mu,x,y2, p2,q2,ierr2)
          IF ((ierr1==0).AND.(ierr2==0)) THEN
            delta=1.0_r8
            icon=0
            DO WHILE((delta>eps).AND.(icon<25))
              icon=icon+1
              IF (pcase) THEN
                yn=y2-(y2-y1)/(p2-p1)*(p2-p)
              ELSE
                yn=y2-(y2-y1)/(q2-q1)*(q2-q)
              ENDIF
              IF(yn<0) yn=1.e-10_r8
              y1=y2
              p1=p2
              q1=q2
              y2=yn
              CALL cdfgamNC(1,mu,x,y2,p2,q2,ierr1)
              IF (ierr1==0) THEN
                IF (pcase) THEN
                  delta=abs(1.0_r8-p2/p)
                  IF(p1-p2==0.0_r8) delta=0.0_r8
                ELSE
                  delta=abs(1.0_r8-q2/q)
                  IF(q1-q2==0.0_r8) delta=0.0_r8
                ENDIF
              ELSE
                ierr=2
                delta=0.0_r8 
              ENDIF
            ENDDO
            xy=yn
            IF (icon==25) ierr=3
          ELSE
            xy=0.0_r8
            ierr=2
          ENDIF 
        ENDIF
      ENDIF     
    ENDIF
    IF (ichi==2) xy=xy*2.0_r8
    END SUBROUTINE invcdfgamNC

    RECURSIVE FUNCTION errorfunction (x, erfcc, expo) RESULT(errfu)
    !-----------------------------------------------------------------
    ! Computation of: the error function, the complementary error 
    !                 function or the scaled complementary error
    !                 function. 
    !-----------------------------------------------------------------
    ! Inputs:
    !       x, real number.
    !       erfcc, logical variable.  
    !       expo,  logical variable.
    !
    ! The meaning of the logical variables erfcc and expo 
    !   is the following:
    !    When erfcc=.true. and expo=.false., the function computes 
    !         the complementary error function erfc(x). 
    !    When erfcc=.true., expo=.true. and x>0, the function computes
    !         the scaled complementary error function exp(x*x}erfc(x). 
    !    When erfcc=.false. and expo=.false., the function computes 
    !         the error function erf(x).
    !------------------------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, y, z, r(0:8), s(0:8), errfu
    LOGICAL erfcc, expo
    IF (erfcc) THEN
      IF (x < -6.5_r8) THEN
        y= 2.0_r8 
      ELSEIF (x < 0.0_r8) THEN
        y= 2.0_r8 - errorfunction(-x, .true., .false.) 
      ELSEIF (x == 0.0_r8) THEN
        y= 1.0_r8 
      ELSEIF (x < 0.5_r8) THEN
        IF (expo) THEN
          y=exp(x*x)
        ELSE
          y=1.0_r8
        ENDIF
        y=y*(1.0_r8-errorfunction(x, .false., .false.))
      ELSEIF (x < 4.0_r8) THEN
        IF (expo) THEN
          y= 1.0_r8 
        ELSE
          y= exp(-x*x)
        ENDIF
        r(0)= 1.230339354797997253e3_r8
        r(1)= 2.051078377826071465e3_r8
        r(2)= 1.712047612634070583e3_r8
        r(3)= 8.819522212417690904e2_r8
        r(4)= 2.986351381974001311e2_r8
        r(5)= 6.611919063714162948e1_r8
        r(6)= 8.883149794388375941_r8
        r(7)= 5.641884969886700892e-1_r8
        r(8)= 2.153115354744038463e-8_r8
        s(0)= 1.230339354803749420e3_r8
        s(1)= 3.439367674143721637e3_r8
        s(2)= 4.362619090143247158e3_r8
        s(3)= 3.290799235733459627e3_r8
        s(4)= 1.621389574566690189e3_r8
        s(5)= 5.371811018620098575e2_r8
        s(6)= 1.176939508913124993e2_r8
        s(7)= 1.574492611070983473e1_r8
        y=y*fractio(x,8,r,s)
      ELSE
        z=x*x
        IF (expo) THEN
          y=1.0_r8 
        ELSE
          y= exp(-z)
        ENDIF
        z=1.0_r8/z
        r(0)=6.587491615298378032e-4_r8
        r(1)=1.608378514874227663e-2_r8
        r(2)=1.257817261112292462e-1_r8
        r(3)=3.603448999498044394e-1_r8
        r(4)=3.053266349612323440e-1_r8
        r(5)=1.631538713730209785e-2_r8
        s(0)=2.335204976268691854e-3_r8
        s(1)=6.051834131244131912e-2_r8
        s(2)=5.279051029514284122e-1_r8
        s(3)=1.872952849923460472_r8
        s(4)=2.568520192289822421_r8
        y=y*(oneoversqrtpi-z*fractio(z,5,r,s))/x
      ENDIF
      errfu=y
    ELSE
      IF (x == 0.0_r8) THEN 
        y=0.0_r8
      ELSEIF (abs(x) > 6.5_r8) THEN 
        y=x/abs(x)
      ELSEIF (x > 0.5_r8) THEN
        y=1.0_r8 - errorfunction(x, .true., .false.) 
      ELSEIF (x < -0.5_r8) THEN
        y=errorfunction(-x, .true., .false.)-1.0_r8
      ELSE
        r(0)=3.209377589138469473e3_r8
        r(1)=3.774852376853020208e2_r8
        r(2)=1.138641541510501556e2_r8
        r(3)=3.161123743870565597e0_r8
        r(4)=1.857777061846031527e-1_r8
        s(0)=2.844236833439170622e3_r8
        s(1)=1.282616526077372276e3_r8
        s(2)=2.440246379344441733e2_r8
        s(3)=2.360129095234412093e1_r8
        z=x*x
        y=x*fractio(z,4,r,s)
      ENDIF  
      errfu= y
    ENDIF        
    END FUNCTION errorfunction


    RECURSIVE FUNCTION inverfc(x) RESULT(y)
    !---------------------------------------------
    ! Inverse of the complementary error function
    !---------------------------------------------
    ! Input:
    !       x, real positive number.
    !---------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) ::  x, y, y0, y02, h, r, f, fp, c1, c2, c3, c4, c5;
    IF (x==0.0_r8) THEN
      y=giant
    ELSE
      IF (x > 1) THEN
        y=-inverfc(2.0_r8-x)
      ELSE
        y0=0.70710678_r8*invq(x*0.5_r8);
        f=errorfunction(y0,.true.,.false.)-x;
        y02= y0*y0;
        fp=-2.0_r8/sqrtpi*exp(-y02);
        c1=-1.0_r8/fp;
        c2= y0;
        c3=(4.0_r8*y02+1.0_r8)*0.333333333333333333333333333333_r8;
        c4=y0*(12*y02+7.0_r8)*.166666666666666666666666666667_r8;
        c5=(8*y02+7)*(12*y02+1.0_r8)*0.333333333333333333333333333333_r8;
        r= f*c1;
        h=r*(1+r*(c2+r*(c3+r*(c4+r*c5))));
        y=y0+h
      ENDIF
    ENDIF
    END FUNCTION inverfc
   
    RECURSIVE FUNCTION gamma(x) RESULT(gam)
    !----------------------------------------
    ! Computation of the Euler gamma function 
    ! Gamma(x).
    !----------------------------------------
    !  Input:
    !       x, real number.
    !----------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, dw, gam, z
    INTEGER :: k, k1, n
    k=nint(x)
    k1=k-1
    IF (k==0) THEN
      dw=dwarf
    ELSE
      dw=machtol
    ENDIF
    IF ((k <= 0).AND.(abs(k - x)<= dw)) THEN
      IF (mod(k,2)>0) THEN
        ! k is odd
        gam=sign(1.0_r8,k-x)*giant
      ELSE
        ! k is even
        gam=sign(1.0_r8,x-k)*giant
      ENDIF
    ELSEIF (x<0.45_r8) THEN
      gam=pi/(sin(pi*x)*gamma(1.0_r8-x))
    ELSEIF ((abs(k-x)<dw).AND.(x<21.0_r8)) THEN
      gam=1;
      DO n=2,k1 
        gam=gam*n
      ENDDO
    ELSEIF ((abs(k-x-0.5_r8)<dw).AND.(x<21.0_r8)) THEN
      gam=sqrt(pi);
      DO n=1,k1 
        gam=gam*(n-0.5)
      ENDDO
    ELSEIF (x<3.0_r8) THEN
      IF (k>x) THEN
        k=k1
      ENDIF
      k1=3-k;
      z=k1+x;
      gam=gamma(z);
      DO n=1,k1 
        gam=gam/(z-n)
      ENDDO
    ELSE
      gam=sqrttwopi*exp(-x+(x-0.5_r8)*log(x)+stirling(x))
    ENDIF
    END FUNCTION gamma

    FUNCTION loggam(x)
    !------------------------------------
    ! Computation of the logarithm of
    ! the Gamma Function, log(gamma(x))
    !------------------------------------
    ! Input:
    !       x, real positive number.
    !------------------------------------
    USE Someconstants
    IMPLICIT NONE  
    REAL(r8) :: loggam, x
    IF (x>=3) THEN
      loggam=(x-0.5_r8)*log(x)- x+lnsqrttwopi+stirling(x)
    ELSEIF (x >= 2) THEN
      loggam=(x-2.0_r8)*(3.0_r8-x)*auxloggam(x-2.0_r8)+logoneplusx(x-2.0_r8)
    ELSEIF (x>=1) THEN
      loggam=(x-1.0_r8)*(2.0_r8-x)*auxloggam(x-1.0_r8)
    ELSEIF (x>0.5_r8) THEN
      loggam=x*(1.0_r8-x)*auxloggam(x)-logoneplusx(x-1.0_r8)
    ELSEIF (x>0.0_r8) THEN
      loggam=x*(1.0_r8-x)*auxloggam(x)-log(x)
    ELSE
      loggam=giant
    ENDIF
    END FUNCTION loggam
  
    FUNCTION gamstar(x)
    !-----------------------------------------------
    ! Computation of the regulated gamma function:
    !    gamstar(x)=exp(stirling(x)), x>0; or 
    !    gamma(x)/(exp(-x+(x-0.5)*ln(x))/sqrt(2*pi)
    !-----------------------------------------------
    ! Input: x, real positive number.
    !-----------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: gamstar, x
    IF (x>=3.0_r8) THEN
      gamstar=exp(stirling(x))
    ELSEIF (x>0.0_r8) THEN
      gamstar=gamma(x)/(exp(-x+(x-0.5_r8)*log(x))*sqrttwopi)
    ELSE
      gamstar=giant
    ENDIF
    END FUNCTION gamstar


    RECURSIVE FUNCTION quotgamm(x,y) RESULT(quotgam)
    !----------------------------------------------------
    ! Computation of the quotient of two gamma functions
    !  gamma(x)/gamma(y)
    !----------------------------------------------------
    USE Someconstants
    IMPLICIT NONE  
    REAL(r8) :: quotgam,x,y
    INTEGER :: n
    IF ((x <= 0.0_r8).OR.(y <= 0.0_r8)) THEN
	quotgam=qratio(x,y)
    ELSEIF (x>y) THEN
      quotgam=1.0_r8/quotgamm(y,x)
    ELSE
      n=int(y-x);
      IF (n==(y - x)) THEN
        quotgam=1.0_r8/shiftfact(x, n)
      ELSEIF (n>15) THEN
        quotgam=qratio(x,y)
      ELSEIF (n>0) THEN
        quotgam=quotgamm(x+n,y)/shiftfact(x,n)
      ELSEIF (x<26.0_r8) THEN
        quotgam=qratio(x,y)
      ELSE
        quotgam=qasym(x,y)
      ENDIF
    ENDIF
    END FUNCTION quotgamm

    FUNCTION qratio(x,y)
    USE Someconstants
    IMPLICIT NONE  			
    REAL(r8) :: qratio,x,y
    REAL(r8) :: b, c, q
    IF (((x<=0.0_r8).OR.(y<=0.0_r8)).OR.(y>1.5_r8*x)) THEN
      qratio=gamma(x)/gamma(y)
    ELSE
      b=y-x;
      c=b/x;
      q=lnec(c);
      q=exp(-b*log(x)-(b-0.5_r8)*log(1.0_r8+c)-x*q);
      qratio=q*gamstar(x)/gamstar(y)
    ENDIF
    END FUNCTION qratio

    FUNCTION qasym(x,y)
    USE Someconstants
    IMPLICIT NONE  			
    REAL(r8) :: qasym, x, y			
    REAL(r8) :: s, t, u, w, w2, r, r2, r3, r4, r5, r6, r7, cc(0:7)
    INTEGER :: k
    w=(x+y-1)/2.0_r8;
    w2=1.0_r8/(w*w);
    r=(x-y+1)/2.0_r8;
    r2=r*r;
    r3=r2*r;
    r4=r2*r2;
    r5=r4*r;
    r6=r5*r;
    r7=r6*r;
    cc(0)=1.0_r8;
    cc(1)=r/12.0_r8;
    cc(2)=r/1440.0_r8+r2/288.0_r8;
    cc(3)=r/90720.0_r8+r2/17280.0_r8+r3/10368.0_r8;
    cc(4)=r/4838400.0_r8+101*r2/87091200.0_r8+r3/414720.0_r8+r4/497664.0_r8;
    cc(5)=r/239500800.0_r8+13*r2/522547200.0_r8+61*r3/1045094400.0_r8+r4/14929920.0_r8+r5/29859840.0_r8;
    cc(6)=691*r/7846046208000.0_r8+7999*r2/14485008384000.0_r8+59*r3/41803776000.0_r8+143*r4/75246796800.0_r8&
          +r5/716636160.0_r8+r6/2149908480.0_r8;
    cc(7)=r/523069747200.0_r8 + 2357*r2/188305108992000.0_r8+5941.0_r8*r3/173820100608000.0_r8 &
          +11*r4/214990848000.0_r8+41*r5/902961561600.0_r8+r6/42998169600.0_r8+r7/180592312320.0_r8;
    s=1.0_r8;
    k=1;
    t=1.0_r8;
    u=1.0_r8
    DO WHILE ((abs(u)>machtol).AND.(k< 8))			
      t=-4.0_r8*w2*(k-r)*(k-r-0.5)*t;
      u=t*cc(k);
      s=s+u;
      k=k+1;
      qasym=s*exp((x-y)*log(w))
    ENDDO
    END FUNCTION qasym

    RECURSIVE FUNCTION shiftfact(x,n) RESULT(shift)
    USE Someconstants
    IMPLICIT NONE  			
    REAL(r8) :: shift, x, s
    INTEGER :: k,n 
    IF (n==0) THEN
      shift=1.0_r8
    ELSEIF (n<0) THEN
      shift=1.0_r8/shiftfact(x-n,n)
    ELSE
      s=1.0_r8;
      k=0;
      DO WHILE (k<n)
        s=s*(x+k);
        k=k+1
      ENDDO					
      shift=s
    ENDIF
    END FUNCTION shiftfact

    FUNCTION invq(x)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: invq, x, t
!  Abramowitx & Stegun 26.2.23; 
    t=sqrt(-2*log(x));
    t=t-(2.515517_r8+t*(0.802853_r8+t*0.010328_r8))/&
      (1.0_r8+t*(1.432788+t*(0.189269_r8+t*0.001308_r8)))
    invq=t
    END FUNCTION invq

    FUNCTION xetay(zeta,mu,y)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: zeta
    REAL(r8), INTENT(IN) :: y
    REAL(r8), INTENT(IN) :: mu
    REAL(r8) :: xetay, x, x0, mu2, z2, p, q, r, s
    INTEGER :: k
    mu2=mu*mu
    z2=0.5*zeta*zeta; 
    p=y*y; q=mu-y; 
    r=y-mu+mu*log(mu/y); 
    s=2*(z2-r);
    IF (y<=mu) THEN    
      p=q*q+p*s;
      IF (p>0) THEN
        x0=s/(q+sqrt(p)) 
      ELSE
        x0=y+z2+sqrt(mu2+4.0_r8*y*z2);
      ENDIF
    ELSE
      IF (zeta > 0) THEN 
        x0=0.0_r8 
      ELSE
        x0=y+z2+sqrt(mu2+4.0_r8*y*z2);
      ENDIF
    ENDIF 
    x=x0; 
    IF (x<0.0_r8) x=0.0_r8
    x0=fxeta(x,y,mu,z2);
    k=1;
    IF (x0<0) k=20
    DO WHILE ((abs(x/x0-1.0_r8)>1.0e-12_r8).AND.(k<20))
      x=x0; x0=fxeta(x,y,mu,z2); 
      IF (x0<0) k=20
      k= k+1;
    ENDDO 
    IF (k<20) x=x0
    xetay=x 
    END FUNCTION xetay

    FUNCTION xzetay(zeta,y)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: zeta
    REAL(r8), INTENT(IN) :: y
    REAL(r8) :: xzetay, tym1, x, x0, w, z, z2, p, q, r, s, t
    REAL(r8) :: yk(0:10), ck(0:10)
    INTEGER :: k
    tym1=2.0_r8*y-1.0_r8
    w=tym1*tym1; 
    z=zeta*sqrt(abs(2.0_r8*y-1.0_r8))/w;
    IF ((abs(z)<0.03_r8).AND.(y>=1.0_r8)) THEN
      yk(2)=y*y
      DO k=3,8
        yk(k)= y*yk(k-1)
      ENDDO
      ck(2)=0.3333333333333333333333333333_r8*(3.0_r8*y-1.0_r8);
      ck(3)=0.0277777777777777777777777777778_r8*(6.0_r8*y-1.0_r8);
      ck(4)=0.0037037037037037037037037037037_r8*(54.0_r8*yk(2)-9.0_r8*y+1.0_r8);
      ck(5)=0.000231481481481481481481481481481_r8*(1152.0_r8*yk(3)-36.0_r8*yk(2)+12.0_r8*y-1.0_r8);
      ck(6)=0.0000587889476778365667254556143445_r8*(6480.0_r8*yk(4)+1404.0_r8*yk(3)&
            -108.0_r8*yk(2)+15.0_r8*y-1.0_r8);
      ck(7)=1.83715461493239271017048794827e-7_r8*(3110400.0_r8*yk(5)+1741824.0_r8*yk(4)&
            -78408.0_r8*yk(3)+20196.0_r8*yk(2)-2502.0_r8*y+139.0_r8);
      ck(8)=0.00000489907897315304722712130119538_r8*(181440.0_r8*yk(6)+178848.0_r8*yk(5)+&
             7128.0_r8*yk(4)+1044.0_r8*yk(3)-198.0_r8*yk(2)+21.0_r8*y-1.0_r8);
      ck(9)=4.25267271975090905132057395432e-10_r8*(571.0_r8+144072.0_r8*yk(2)&
            -13704.0_r8*y+3344302080.0_r8*yk(7)+4961710080.0_r8*yk(6)&
            +805469184.0_r8*yk(5)-845856.0_r8*yk(3)+9283248.0_r8*yk(4));
      ck(10)=6.59808615912868313417010262003e-10_r8*(3527193600.0_r8*yk(8)&
             +7215350400.0_r8*yk(7)+2287963584.0_r8*yk(6)+75781008.0_r8*yk(5)&
             +3557520.0_r8*yk(4)-738720.0_r8*yk(3)+95580.0_r8*yk(2)&
             -7587.0_r8*y+281.0_r8);
      s=-z; t=1.0_r8; k=2; z2=z*z; 
      DO WHILE ((abs(t/s)>1.0e-15_r8).AND.(k < 11)) 
        t=ck(k)*z2; s=s+t; z2=z*z2; k=k+1
      ENDDO  
      x=y-1+w*s
    ELSE
      z2=0.5*zeta*zeta; 
      p=y*y; q=1.0_r8-y; 
      r=y-1-log(y); 
      s=2*(z2-r);
      IF (y<=1) THEN    
        p=q*q+p*s;
        IF (p>0) THEN
          x0=s/(q+sqrt(p)) 
        ELSE
          x0=y+z2+sqrt(1.0_r8+4.0_r8*y*z2);
        ENDIF
      ELSE
        IF (zeta > 0) THEN 
          x0=0 
        ELSE
          x0=y+z2+sqrt(1.0_r8+4.0_r8*y*z2);
        ENDIF
      ENDIF
      x=x0; 
      x0=fx(x,y,z2);
      k=1;
      DO WHILE ((abs(x/x0-1.0_r8)>1.0e-12_r8).AND.(k<20))
        x=x0; x0=fx(x,y,z2); k= k+1;
      ENDDO 
      x=x0
    ENDIF
    xzetay=x 
    END FUNCTION xzetay

    FUNCTION yzetax(zeta,x) 
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8) :: yzetax, zeta, x
    REAL(r8) :: txp1, y, y0, w, z, z2,  s, t, xk(0:100), ck(0:100)  
    INTEGER :: k
    txp1=2*x+1.0_r8 
    w=txp1*txp1; z=zeta*sqrt(txp1)/w;
    IF (abs(z)<0.05) THEN    
      xk(2)=x*x; 
      DO k=3,8 
        xk(k)=x*xk(k-1);
      ENDDO 
      ck(2)=0.3333333333333333333333333333_r8*(3.0_r8*x+1.0_r8);
      ck(3)=0.0277777777777777777777777777778_r8*(6.0_r8*x+1.0_r8);
      ck(4)=-0.0037037037037037037037037037037_r8*(54.0_r8*xk(2)+9.0_r8*x+1.0_r8);
      ck(5)=0.000231481481481481481481481481481_r8*(1152.0_r8*xk(3)+36.0_r8*xk(2)+12.0_r8*x+1.0_r8);
      ck(6)=-0.0000587889476778365667254556143445_r8*(6480.0_r8*xk(4)-1404.0_r8*xk(3)-108.0_r8*xk(2)-15*x-1);
      ck(7)=1.83715461493239271017048794827e-7_r8*(3110400.0_r8*xk(5)-1741824.0_r8*xk(4)&
            -78408.0_r8*xk(3)-20196.0_r8*xk(2)-2502.0_r8*x-139);
      ck(8)=-0.00000489907897315304722712130119538_r8*(181440.0_r8*xk(6)-178848.0_r8*xk(5)&
            +7128.0_r8*xk(4)-1044.0_r8*xk(3)-198.0_r8*xk(2)-21*x-1);
      ck(9)=4.25267271975090905132057395432e-10_r8*(3344302080.0_r8*xk(7)-4961710080.0_r8*xk(6)&
            +805469184.0_r8*xk(5)-9283248.0_r8*xk(4)-845856.0_r8*xk(3)&
            -144072.0_r8*xk(2)-13704.0_r8*x-571);
      ck(10)=-6.59808615912868313417010262003e-10_r8*(281.0_r8+7587.0_r8*x+95580.0_r8*xk(2)&
             +738720.0_r8*xk(3)+3557520.0_r8*xk(4)-75781008.0_r8*xk(5)+2287963584.0_r8*xk(6)&
             -7215350400.0_r8*xk(7)+3527193600.0_r8*xk(8));
      s=z; t=1.0_r8; k=2; z2=z*z; 
      DO WHILE ((abs(t/s)>1.0e-15_r8).AND.(k < 11)) 
        t=ck(k)*z2; s=s+t; z2=z*z2; k=k+1
      ENDDO
      y=x+1+w*s
    ELSE
      z2=0.5_r8*zeta*zeta;
      IF (zeta>0) THEN 
        y0=x+2 
      ELSE
        y0=exp(-z2);
      ENDIF
      y=0; k=1;
      DO WHILE ((abs(y/y0-1)>1.0e-8_r8).AND.(k<20))
        y=y0; y0=fy(x,y,z2); k= k+1;
      ENDDO 
      y= y0
    ENDIF
    yzetax=y
    END FUNCTION yzetax     

    FUNCTION yetax(zeta,mu,x)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8) :: yetax, zeta, mu, x
    REAL(r8) :: y, y0, z2 
    INTEGER :: k
    z2=0.5_r8*zeta*zeta;
    IF (zeta>0) THEN 
      y0=x+mu+1.0_r8 
    ELSE
      y0=exp(-z2);
    ENDIF
    y=0; k=1;
    DO WHILE ((abs(y/y0-1)>1.0e-8_r8).AND.(k<20))
      y=y0; y0=fyeta(x,y,mu,z2);
      IF (y0<0) k=20
      k= k+1;
    ENDDO
    IF (k<20) y= y0
    yetax=y
    END FUNCTION yetax 


    FUNCTION fxeta(x,y,mu,z2)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x,y,mu,z2
    REAL(r8) :: mu2, w, fxeta, f, fd
    mu2=mu*mu
    w=sqrt(mu2+4*x*y);
    fd=(w+mu-2.0_r8*y)/(mu+w);
    f=x+y-w+mu*log((mu+w)/(2.0_r8*y))-z2;
    fxeta=x-f/fd
    END FUNCTION fxeta

    FUNCTION fx(x,y,z2)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x,y,z2
    REAL(r8) :: w, fx, f, fd
    w=sqrt(1.0_r8+4*x*y);
    fd=(w+1.0_r8-2.0_r8*y)/(1.0_r8+w);
    f=x+y-w+log((1.0_r8+w)/(2.0_r8*y))-z2;
    fx=x-f/fd
    END FUNCTION fx

    FUNCTION fyeta(x,y,mu,z2)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x,y,mu,z2
    REAL(r8) :: mu2, w, fyeta, f, fd
    mu2=mu*mu
    w=sqrt(mu2+4*x*y);
    fd=((y-mu)*w+y*mu-2*x*y-mu2)/(y*(mu+w));
    f=x+y-w+mu*log((mu+w)/(2.0_r8*y))-z2;
    fyeta=y-f/fd
    END FUNCTION fyeta

    FUNCTION fy(x,y,z2)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x,y,z2
    REAL(r8) :: w, fy, f, fd
    w=sqrt(1.0_r8+4*x*y);
    fd=((y-1)*w+y-2*x*y-1)/(y*(1+w));
    f=x+y-w+log((1.0_r8+w)/(2*y))-z2;
    fy=y-f/fd
    END FUNCTION fy

    FUNCTION eta1xy(zeta,x,y,mu)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: zeta,x,y,mu
    REAL(r8) :: eta1xy, w, mu2, zeta1, f0
    mu2=mu*mu
    w=sqrt(mu2+4.0_r8*x*y);
    f0=zeta*(mu+w)/((-mu-w+2*y)*sqrt(w));
    zeta1=log(f0)/zeta
    eta1xy= zeta1
    END FUNCTION eta1xy

    FUNCTION zeta1xy(zeta,x,y)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: zeta,x,y
    REAL(r8) :: zeta1xy, txp1, w, z, zeta1, s, t, f0
    REAL(r8) :: zk(0:11), xk(0:11), ck(0:11)
    INTEGER k
    txp1=2.0_r8*x+1
    w=txp1*sqrt(txp1); z=zeta/w;
    IF (abs(z)<0.03_r8) THEN    
      xk(2)=x*x; zk(2)=z*z; 
      DO k=3,11  
        xk(k)=x*xk(k-1); 
        zk(k)=z*zk(k-1) 
      ENDDO
      ck(0)=-(1.0_r8/3.0_r8)*(3.0_r8*x+1);
      ck(1)=(1.0_r8/36.0_r8)*(36*xk(2)+6*x+1.0_r8)*z;
      ck(2)=-(1.0_r8/1620.0_r8)*(2160.0_r8*xk(3)-594.0_r8*xk(2)-9.0_r8*x-1)*zk(2);
      ck(3)=(1.0_r8/6480.0_r8)*(12960.0_r8*xk(4)-11664.0_r8*xk(3)+108.0_r8*xk(2)&
             -84.0_r8*x-7)*zk(3);
      ck(4)=-(1.0_r8/90720.0_r8)*(-2952.0_r8*xk(2)-375.0_r8*x+62676.0_r8*xk(3)+&
             290304.0_r8*xk(5)-482112.0_r8*xk(4)-25.0_r8) *zk(4);
      ck(5)=(1.0_r8/382725.0_r8)*(-1809.0_r8*xk(2)-198.0_r8*x+2041200.0_r8*xk(6)-11.0_r8-&
            72117.0_r8*xk(3)-5161320.0_r8*xk(5)+1588734.0_r8*xk(4))*zk(5);
      ck(6)=-(1.0_r8/73483200.0_r8)*(1147209.0_r8*x+54629.0_r8+10506132.0_r8*xk(2)+&
             53276346.0_r8*xk(3)+145651608.0_r8*xk(4)+26074872.0_r8*xk(5)-3096325440.0_r8*xk(6)&
             +17239858560.0_r8*xk(7))*zk(6); 
      ck(7)=(1.0_r8/9797760.0_r8)*(9864.0_r8*xk(2)+888.0_r8*x+562961664.0_r8*xk(6)+37.0_r8&
             +156764160.0_r8*xk(8)-707927040.0_r8*xk(7)+75456.0_r8*xk(3)-100476288.0_r8*xk(5)&
              +3282768.0_r8*xk(4))*zk(7); 
      ck(8)=-(1.0_r8/3491921664000.0_r8)*(1085205600.0_r8*xk(2)+85989897.0_r8*x-&
              173727743003136.0_r8*xk(6)+3184811.0_r8+99325771776000.0_r8*xk(9)&
              -557433645465600.0_r8*xk(8)+609525121720320.0_r8*xk(7)+&
              8540105400.0_r8*xk(3)+11923016301072.0_r8*xk(5)-8705251440.0_r8*xk(4))*zk(8);
      ck(9)=(1.0_r8/15913705500.0_r8)*(523341.0_r8*xk(2)+36930.0_r8*x+386373372840.0_r8*xk(6)&
             +1231.0_r8+814781721600.0_r8*xk(10)-5512856630400.0_r8*xk(9)+7843246148160.0_r8*xk(8)&
             -3236596377408.0_r8*xk(7)+4431105.0_r8*xk(3)-8476062588.0_r8*xk(5)+4677750.0_r8*xk(4))&
             *zk(9);
      ck(10)=-(1.0_r8/84737299046400.0_r8)*(-1405644714.0_r8*xk(2)-90601269.0_r8*x&
              -570910822987680.0_r8*xk(6)-2745493.0_r8-62881400385699840.0_r8*xk(10)&
              +111932154011713536.0_r8*xk(9)-62582386500956160.0_r8*xk(8)&
              +11503114507558656.0_r8*xk(7)-13895395416.0_r8*xk(3)+1634520610800.0_r8*xk(5)&
              -100596777120.0_r8*xk(4)+7888272202137600.0_r8*xk(11))*zk(10);
      s=ck(0)+ck(1); t= 1; k= 2;
      DO WHILE ((abs(t)>1.0e-15_r8).AND.(k < 11)) 
        t=ck(k); s=s+t; k= k+1 
      ENDDO
      zeta1=s/w
    ELSE 
      w=sqrt(1+4.0_r8*x*y);
      f0=zeta*(1.0_r8+2*x+w)/(2*(y-x-1.0_r8)*sqrt(w));
      zeta1=log(f0)/zeta
    ENDIF
    zeta1xy= zeta1
    END FUNCTION zeta1xy

    SUBROUTINE zeta3xy(iflag,zeta0,x,y,zeta1,zeta2,zeta3)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN)  :: zeta0, x, y
    REAL(r8), INTENT(OUT) :: zeta1, zeta2, zeta3
    REAL(r8) :: xx2, x3, x4, pp, p, pp2, p3, p4, p5, p6, p7, p8, p9, p10, y2, &
                zeta02, zeta03, zeta04, zeta05, zeta06, zeta12, zeta13, dd, dd4
    REAL(r8) :: s, t
    REAL(r8) :: a1, a12,a13,a14,a15,a16,a17,a18,a19,a110,a111,a112,a113,&
                a114,a115,a116,a117,a118,a119,a120,a121,a122,a123,a124,a125,&
                a126,a127,a128,a129,a130,a131,a132,a133,a134,a135,a136,a137,a138,&
                a139,a140,a141,a142,a143,a144,a145
    REAL(r8) :: zk(0:13), z1k(0:13), z2k(0:10), z3k(0:10)
    REAL(r8) :: denom
    INTEGER, INTENT(IN) :: iflag
    INTEGER k
    zeta1=zeta1xy(zeta0,x,y)
    pp=1.0_r8/sqrt(1.0_r8+4.0_r8*x*y)
    IF ((abs(zeta0)>=0.1_r8).AND.(abs(1.0_r8-pp)>1.e-4_r8)) THEN
      p=pp;           
      pp2=pp*pp
      p3=pp2*pp
      p4=p3*pp
      p5=p4*pp
      p6=p5*pp
      p7=p6*pp
      p8=p7*pp  
      p9=p8*pp
      p10=p9*pp
      zeta02=zeta0*zeta0
      zeta03=zeta02*zeta0
      zeta04=zeta03*zeta0
      zeta05=zeta04*zeta0
      zeta06=zeta05*zeta0
      zeta12=zeta1*zeta1
      zeta13=zeta12*zeta1   
      IF (iflag==1) THEN
        xx2=x*x     
        x3=xx2*x
        x4=x3*x
        dd=(2.0_r8*pp*x+pp-1.0_r8)*(2.0_r8*pp*x+pp-1.0_r8)
        dd4=dd*dd
        zeta2=-(1.0_r8/24.0_r8)*(-24.0_r8+48.0_r8*pp-24.0_r8*pp2-&
              96.0_r8*pp2*xx2+96.0_r8*pp*x+38.0_r8*zeta02*p3+&
              5.0_r8*zeta02*p5-22.0_r8*zeta02*p4+12.0_r8*zeta02*zeta12+&
              48.0_r8*zeta02*zeta12*x*pp2-48.0_r8*zeta02*zeta12*x*pp-&
              24.0_r8*zeta1*zeta03*x*p4+24.0_r8*zeta1*zeta03*x*pp2&
              -36.0_r8*zeta03*zeta1*pp2+12.0_r8*zeta02*zeta12*pp2+&
              12.0_r8*zeta03*zeta1*pp+36.0_r8*zeta02*x*pp2-12.0_r8*zeta02*x*p3&
              +48.0_r8*pp2*xx2*zeta02*zeta12-24.0_r8*pp*zeta02*zeta12+9.0_r8*zeta02*pp&
              -30.0_r8*zeta02*pp2-96.0_r8*x*pp2-12.0_r8*zeta03*zeta1*p4&
              +36.0_r8*zeta03*zeta1*p3-44.0_r8*zeta02*x*p4+&
              20.0_r8*zeta02*x*p5+20.0_r8*p5*xx2*zeta02-12.0_r8*p3*xx2*zeta02)&
              /(zeta03*((2.0_r8*pp*x+pp-1.0_r8)*(2.0_r8*pp*x+pp-1.0_r8)));
        IF (zeta0>0.1_r8) THEN
          zeta3=(1.0_r8/48.0_r8)*(-120.0_r8+480*p+2880*p3*x+960*pp*x+240*p4*zeta02&
                -24*zeta02*p-240*zeta02*p3-2880*pp2*xx2+&
                120*pp2*zeta02+5760*p3*xx2-120*zeta02*p5+24*p6*zeta02-2880*p4*xx2&
                -3840*p4*x3+3840*p3*x3&
                -960*p4*x-720*pp2+480*p3-120*p4-2880*pp2*x+576*zeta02*p4*x+48*zeta02*pp2*x&
                -288*zeta02*p3*x+&
                144*zeta02*p6*x-480*zeta02*p5*x+288*zeta02*p6*xx2+96*zeta02*p4*xx2+&
                96*zeta02*p3*xx2-480*zeta02*p5*xx2&
                -192*p5*x3*zeta03*zeta1+320*zeta03*p7*x3*zeta1+240*zeta03*p7*xx2*zeta1&
                -96*p5*x4*zeta03*zeta1+80*zeta1*zeta03*x*p7&
                +1536*zeta0*zeta1*x3*p3+1152*zeta0*zeta1*x*p3+2304*zeta0*zeta1*xx2*p3+&
                160*zeta03*x4*p7*zeta1-240*zeta1*zeta05*x3*p9&
                -48*zeta0*zeta1*p4-168*p5*x3*zeta04+174*zeta03*p5*zeta1-90*zeta05*zeta1*p5+&
                30*zeta05*zeta1*p3-6*zeta05*zeta1*pp2&
                -204*p5*xx2*zeta04-780*zeta04*x*p4+420*p4*xx2*zeta04+192*p6*x3*zeta02&
                -365*zeta04*p7+305*zeta04*p6&
                -95*zeta04*p9+253*zeta04*p8+15*zeta04*p10-1920*p4*x4-145*zeta04*p5+35*zeta04*p4&
                -48*zeta0*zeta1-256*zeta03*zeta1*p4&
                +214*zeta03*zeta1*p3-288*zeta0*zeta1*pp2-96*zeta03*zeta1*pp2+18*zeta03*zeta1*p+&
                192*zeta0*zeta1*p-120*pp2*zeta04*zeta12&
                -12*pp2*zeta06*zeta12-96*p4*zeta06*zeta12+60*p3*zeta06*zeta12+&
                240*p3*zeta04*zeta12-96*zeta03*pp*zeta13+144*zeta03*pp2*zeta13&
                +24*zeta03*zeta13+24*pp*zeta04*zeta12+120*p5*zeta04*zeta12+24*p8*zeta06*zeta12+&
                24*p5*zeta06*zeta12-84*p7*zeta06*zeta12+&
                84*p6*zeta06*zeta12-24*p6*zeta04*zeta12+24*zeta03*p4*zeta13&
                -96*p4*zeta04*zeta12*xx2+480*p5*zeta04*zeta12*xx2+&
                192*p4*zeta04*zeta12*x3+288*p3*zeta04*zeta12*x+480*p5*zeta06*zeta12*x+&
                72*p3*zeta06*zeta12*x-336*p4*zeta06*zeta12*x&
                +96*p4*zeta06*zeta12*xx2-192*p6*zeta06*zeta12*xx2+96*p8*zeta06*zeta12*xx2&
                -192*p6*zeta04*zeta12*x3&
                -576*zeta03*p3*zeta13*x-192*zeta03*pp*zeta13*x+576*zeta03*pp2*zeta13*xx2&
                +576*zeta03*p4*zeta13*xx2+&
                768*zeta03*p4*zeta13*x3-768*zeta03*p3*zeta13*x3+576*zeta03*pp2*zeta13*x&
                +384*zeta03*p4*zeta13*x4&
                -576*p4*zeta04*zeta12*x-48*pp2*zeta04*zeta12*x-96*p3*zeta04*zeta12*xx2&
                -144*p6*zeta04*zeta12*x+&
                96*p8*zeta06*zeta12*x-144*p6*zeta06*zeta12*x-168*p7*zeta06*zeta12*x&
                -96*zeta03*p3*zeta13+&
                192*zeta03*p4*zeta13*x+480*p5*zeta04*zeta12*x-288*p6*zeta04*zeta12*xx2&
                -240*p4*zeta04*zeta12&
                -1152*zeta03*xx2*p3*zeta13+600*zeta05*zeta1*p8*x-360*zeta05*zeta1*p9*xx2&
                -180*zeta05*zeta1*p9*x&
                +312*zeta05*zeta1*p4*xx2-744*zeta05*zeta1*p4*x+156*zeta05*zeta1*p3*x+&
                432*zeta05*zeta1*p7*xx2-912*zeta05*zeta1*p6*xx2&
                -444*zeta05*zeta1*p7*x-624*zeta05*zeta1*p6*x-48*zeta05*zeta1*p5*x3&
                -72*zeta05*zeta1*p5*xx2+1236*zeta05*zeta1*p5*x&
                +288*zeta1*zeta05*x3*p7-3*zeta04*p3+96*p6*x3*zeta04&
                -576*p8*x3*zeta04-288*p8*x4*zeta04+&
                162*zeta04*p3*x+10*zeta03*zeta1*p7+48*p6*x4*zeta04-64*zeta03*p6*zeta1&
                -294*zeta05*zeta1*p7+270*zeta05*zeta1*p6&
                -30*zeta05*zeta1*p9+150*zeta05*zeta1*p8-30*zeta05*zeta1*p4+&
                192*zeta0*zeta1*p3+144*zeta03*x*p3*zeta1+672*zeta03*p5*x*zeta1&
                +576*zeta03*xx2*p5*zeta1-384*zeta1*zeta0*x*p4-1536*p4*x3*zeta0*zeta1&
                -1152*p4*xx2*zeta0*zeta1-768*p4*x4*zeta0*zeta1&
                -1140*zeta04*xx2*p9+940*zeta04*x*p8-760*p9*x3*zeta04+&
                652*p8*xx2*zeta04+928*p7*x3*zeta04&
                -1400*p6*xx2*zeta04+360*zeta04*xx2*p10-570*zeta04*x*p9&
                -306*zeta04*x*p7+1312*zeta04*xx2*p7&
                -984*zeta04*x*p6+120*zeta04*p10*x+1418*zeta04*p5*x&
                -192*zeta02*p4*x3+480*p10*x3*zeta04+&
                240*p10*x4*zeta04+192*zeta03*xx2*p4*zeta1&
                -512*zeta03*x3*p6*zeta1+384*zeta03*x3*p4*zeta1+&
                600*zeta1*zeta05*xx2*p8-240*zeta03*xx2*p3*zeta1&
                -768*zeta03*p6*xx2*zeta1-384*zeta03*p6*zeta1*x&
                -512*zeta03*zeta1*p4*x-1152*zeta0*zeta1*pp2*xx2&
                -1152*zeta0*zeta1*pp2*x+384*zeta0*zeta1*pp*x)/(zeta05*dd4);
        ELSE
          zeta3=0.0_r8
        ENDIF
        IF (abs(zeta2)>abs(10*zeta1)) zeta2=0.0_r8
        IF (abs(zeta3)>abs(5*zeta2)) zeta3=0.0_r8
      ELSE
        y2=y*y
        dd=(2.0_r8*pp*y+pp-1.0_r8)*(2.0_r8*pp*y+pp-1.0_r8)
        dd4=dd*dd
        zeta2=-(1.0_r8/24.0_r8)*(-24.0_r8-24.0_r8*pp2-20.0_r8*zeta02*y*p5+&
             60.0_r8*zeta02*y*p3+4.0_r8*zeta02*y*p4+&
             36.0_r8*zeta02*y*pp2+5.0_r8*zeta02*p5+2.0_r8*zeta02*p3-3.0_r8*zeta02*pp+&
             10.0_r8*zeta02*p4-6.0_r8*zeta02*pp2-96.0_r8*y2*pp2+96.0_r8*y*pp2-48.0_r8*pp+&
             96.0_r8*y*pp+12.0_r8*zeta02*zeta12&
             +48.0_r8*zeta02*zeta12*y2*pp2-48.0_r8*zeta02*zeta12*y*pp2&
             -48.0_r8*zeta02*zeta12*y*pp&
             -24.0_r8*p4*zeta1*zeta03*y+36.0_r8*zeta02*y2*p3-48.0_r8*zeta02*y2*p4+&
              20.0_r8*zeta02*y2*p5+72.0_r8*pp2*zeta1*zeta03*y+&
             48.0_r8*p3*zeta1*zeta03*y+12.0_r8*p4*zeta03*zeta1+&
             12.0_r8*p3*zeta1*zeta03-12.0_r8*pp*zeta1*zeta03-12.0_r8*pp2*zeta1*zeta03&
             +12.0_r8*zeta02*zeta12*pp2+24.0_r8*zeta02*zeta12*pp)/(zeta03*&
             (-pp+2.0_r8*y*pp-1.0_r8)*(-pp+2.0_r8*y*pp-1.0_r8));
        zeta3=0.0_r8
      ENDIF    
    ELSE
     IF (abs(zeta0)<0.1_r8) THEN

        zk(1)=zeta0
        DO k=2,13  
          zk(k)=zeta0*zk(k-1) 
        ENDDO       
      IF (iflag==1) THEN 
        a1=-sqrt(2.0_r8*y-1.0_r8); 
      ELSE
        a1=sqrt(2.0_r8*x+1.0_r8);
      ENDIF
      a12=a1*a1;a13=a12*a1;a14=a13*a1;a15=a14*a1;a16=a15*a1;
      a17=a16*a1;a18=a17*a1;a19=a18*a1;a110=a19*a1;
      a111=a110*a1;a112=a111*a1;a113=a112*a1;a114=a113*a1;
      a115=a114*a1;a116=a115*a1;a117=a116*a1;a118=a117*a1;
      a119=a118*a1;
      a120=a119*a1;a121=a120*a1;a122=a121*a1;a123=a122*a1;
      a124=a123*a1;a125=a124*a1;a126=a125*a1;a127=a126*a1;
      a128=a127*a1;a129=a128*a1;
      a130=a129*a1;a131=a130*a1;a132=a131*a1;a133=a132*a1;
      a134=a133*a1;a135=a134*a1;a136=a135*a1;a137=a136*a1;
      a138=a137*a1;a139=a138*a1;
      a140=a139*a1;a141=a140*a1;a142=a141*a1;a143=a142*a1;
      a144=a143*a1;a145=a144*a1;
      IF (iflag==1) THEN          
        z1k(0)=(1.0_r8/6.0_r8)*(-1.0_r8+3.0_r8*a12)/a13;
        z1k(1)=-(1.0_r8/36.0_r8)*(-7.0_r8+9.0_r8*a14-3.0_r8*a12)/a16;
        z1k(2)=(1.0_r8/3240.0_r8)*(-27*a14-1305*a12-830+540*a16)/a19;
        z1k(3)=-(1.0_r8/6480.0_r8)*(378.0_r8*a16-3861.0_r8*a14-5790.0_r8*a12&
                 -2330.0_r8+810.0_r8*a18)/a112;
        z1k(4)=(1.0_r8/181440.0_r8)*(21384.0_r8*a18-316680.0_r8*a12&
                -362187.0_r8*a14+18144.0_r8*a110-137403.0_r8*a16-95228.0_r8)/a115;
        z1k(5)=-(1.0_r8/3061800.0_r8)*(-2687337.0_r8*a18-9958725.0_r8*a12&
                 -15661485.0_r8*a14+255150.0_r8*a112+524880.0_r8*a110-11068785.0_r8*a16&
               -2409050.0_r8)/a118;
        z1k(6)=(1.0_r8/32659200.0_r8)*(-188597241.0_r8*a18-193437300.0_r8*a12&
                -383271210.0_r8*a14+7186320.0_r8*a112-31075812.0_r8*a110&
                +2332800.0_r8*a114-381779055.0_r8*a16-39276200.0_r8)/a121;
        z1k(7)=-(1.0_r8/9797760.0_r8)*(-224575173.0_r8*a18-104286000.0_r8*a12&
                -248164308.0_r8*a14-9558324.0_r8*a112-83064852.0_r8*a110&
                +2591352.0_r8*a114-315202860.0_r8*a16+612360.0_r8*a116-18277000.0_r8)/a124;
        z1k(8)=(1.0_r8/10158317568000.0_r8)*(-109954260438210.0_r8*a18+&
                46249674715125.0_r8*a12+62781588191700.0_r8*a14&
                -54608288635548.0_r8*a112-119825060843610.0_r8*a110&
                -7563831922620.0_r8*a114-8702169317100.0_r8*a16&
                +453523384125.0_r8*a118+1757028090825.0_r8*a116+10636061704625.0_r8)/a127;
        z1k(9)=-(1.0_r8/127309644000.0_r8)*(-4279315309500.0_r8*a12-13522731279750.0_r8*a14&
                 -24318275837625.0_r8*a16-27193050218565.0_r8*a18-19306988291745.0_r8*a110&
                 -8416870188111.0_r8*a112-1966263434982.0_r8*a114-110547732420.0_r8*a116&
                 +43701836850.0_r8*a118+6365482200.0_r8*a120-590307086600.0_r8)/a130; 
        z1k(10)=(1.0_r8/169474598092800.0_r8)*(-10054160837120400.0_r8*a12&
                 -35647819563676320.0_r8*a14-73341753727362120.0_r8*a16-96411337398782406.0_r8*a18&
                 -83837472483930165.0_r8*a110-48030703419545847.0_r8*a112-17258239538608068.0_r8*a114&
                 -3328592071143456.0_r8*a116-124745787998208.0_r8*a118+64150951173120.0_r8*a120&
                 +7703390822400.0_r8*a122-1254367500495200.0_r8)/a133; 
        z1k(11)=-(1.0_r8/127673385840000.0_r8)*(-13328252707770000.0_r8*a12&
                 -52369797400420500.0_r8*a14&
                 -121214416407084000.0_r8*a16-182880107334960750.0_r8*a18-187699522401641250.0_r8*a110&
                 -132420813236216505.0_r8*a112-63077110684718130.0_r8*a114-19108750611117258.0_r8*a116&
                 -3096024909924960.0_r8*a118-70413390440550.0_r8*a120+52537017231000.0_r8*a122&
                 +5319724410000.0_r8*a124-1518227841130000.0_r8)/a136; 
        z1k(12)=(1.0_r8/61010855313408000.0_r8)*(-11182030123135872000.0_r8*a12&
                 -48211055728647411600.0_r8*a14-123923771254537254000.0_r8*a16&
                 -210854486895115763700.0_r8*a18-249105487192355593800.0_r8*a110&
                 -208187581439825642655.0_r8*a112-122734345760116505415.0_r8*a114&
                 -49700657941195125384.0_r8*a116&
                 -12918774774646679328.0_r8*a118-1782426573908893440.0_r8*a120&
                 -19225434558297600.0_r8*a122+26994120536678400.0_r8*a124&
                 +2346571358208000.0_r8*a126-1172073510107192000.0_r8)/a139; 
        denom=135.0_r8*a16-81.0_r8*a14+45.0_r8*a12-155.0_r8 
        z2k(0)=-(1.0_r8/3240.0_r8)*denom/a19;
        z2k(1)=(1.0_r8/2592.0_r8)*(243.0_r8*a18-108.0_r8*a16-81.0_r8*a14-420.0_r8*a12&
                -505.0_r8)/a112;
        z2k(2)=-(1.0_r8/816480.0_r8)*(115668.0_r8*a110-10125.0_r8*a18-168777.0_r8*a16&
                -422163.0_r8*a14-808185.0_r8*a12-450646.0_r8)/a115;
        z2k(3)=(1.0_r8/2099520.0_r8)*(386370.0_r8*a112+169128.0_r8*a110-1119501.0_r8*a18&
               -3223665.0_r8*a16-7046865.0_r8*a14-7557660.0_r8*a12-2817130.0_r8)/a118; 
        z2k(4)=-(1.0_r8/97977600.0_r8)*(21695040.0_r8*a114+24086160.0_r8*a112&
               -96673419.0_r8*a110-379667007.0_r8*a18-930893985.0_r8*a16&
               -1394843625.0_r8*a14-1038337440.0_r8*a12-293913460.0_r8)/a121; 
        z2k(5)=(1.0_r8/1763596800.0_r8)*(449166060.0_r8*a116+859525992.0_r8*a114&
               -2656196064.0_r8*a112-14784723612.0_r8*a110-41801254713.0_r8*a18&
               -77951791260.0_r8*a16-85566907944.0_r8*a14-48921809160.0_r8*a12&
               -11223491860.0_r8)/a124; 
        z2k(6)=-(1.0_r8/1745960832000.0_r8)*(496767427200.0_r8*a118+1408198977360.0_r8*a116&
              -3493441585260.0_r8*a114-27999657809937.0_r8*a112-92926054500345.0_r8*a110&
              -205740653936895.0_r8*a18-294351301812645.0_r8*a16-252775949608710.0_r8*a14&
              -117580863321900.0_r8*a12-22746727187800.0_r8)/a127; 
        z2k(7)=(1.0_r8/185177664000.0_r8)*(57686149080.0_r8*a120+222823107360.0_r8*a118&
              -436405833540.0_r8*a116-5162127807816.0_r8*a114-20227890756213.0_r8*a112&
              -52060153504905.0_r8*a110-91113392194965.0_r8*a18-103966409573460.0_r8*a16&
              -73162123900740.0_r8*a14-28729989825600.0_r8*a12-4813736065880.0_r8)/a130; 
        z2k(8)=-(1.0_r8/4575814148505600.0_r8)*(1538064514022400.0_r8*a122+7672351342360320.0_r8*a120&
              -11167805088902592.0_r8*a118-205652023257710016.0_r8*a116-950381517372586299.0_r8*a114&
              -2810929213567245771.0_r8*a112-5804053176780884025.0_r8*a110&
              -8234830657985188257.0_r8*a18-7765555911833079360.0_r8*a16&
              -4626732255591405960.0_r8*a14-1573659776235346800.0_r8*a12-232808295847072400.0_r8)/a133; 
        z2k(9)=(1.0_r8/9899597917440000.0_r8)*(3551202490374000.0_r8*a124+22031130267554400.0_r8*a122&
             -20837102148613440.0_r8*a120-675301583963071080.0_r8*a118-3664495983529292058.0_r8*a116&
             -12356951269780552050.0_r8*a114-29435189701611164625.0_r8*a112-49859198031571298370.0_r8*a110&
             -58946499123740602110.0_r8*a18-47245278542339908800.0_r8*a16&
             -24406920061093785600.0_r8*a14-7326842720724450000.0_r8*a12-971150252248330000.0_r8)/a136; 
        z2k(10)=-(1.0_r8/1372744244551680000.0_r8)*(521075978809344000.0_r8*a126+&
              3909422480008780800.0_r8*a124-1633235164439558400.0_r8*a122-135284080062035250720.0_r8*a120&
              -857295093769286465124.0_r8*a118-3274109931036233003247.0_r8*a116&
              -8863662891070952959995.0_r8*a114-17449997056857395695665.0_r8*a112&
              -24805449338445590346675.0_r8*a110-25028886072322240272450.0_r8*a18&
              -17423719285466831337000.0_r8*a16-7946267334428713399200.0_r8*a14&
              -2135976241246185702000.0_r8*a12-256579957380860372000.0_r8)/a139; 
        z3k(0)=(1.0_r8/816480.0_r8)*(36645.0_r8*a12-18873.0_r8*a18+51849.0_r8*a14-&
               70119.0_r8*a16+23814.0_r8*a110-72268.0_r8)/a115;
        z3k(1)=-(1.0_r8/20995200.0_r8)*(-8069175.0_r8*a12-9594126.0_r8*a18+12125430.0_r8*a14&
               -601695.0_r8*a16-725355.0_r8*a110+1840725.0_r8*a112-13017775.0_r8)/a118;
        z3k(2)=(1.0_r8/146966400.0_r8)*(-642384960.0_r8*a12+32359635.0_r8*a14-182419209.0_r8*a18&
               -197892153.0_r8*a110&
               +277231815.0_r8*a16+3090960.0_r8*a112-393014140.0_r8+25981560.0_r8*a114)/a121;
        z3k(3)=-(1.0_r8/5290790400.0_r8)*(-125006633430.0_r8*a12+10913802705.0_r8*a18&
               -77775185187.0_r8*a14&
               +24578935920.0_r8*a16+1163129706.0_r8*a114-30234584880.0_r8*a110&
               -15750047187.0_r8*a112+1570244130.0_r8*a116-48417507530.0_r8)/a124;
        z3k(4)=(1.0_r8/1047576499200.0_r8)*(-99583081889940.0_r8*a12+14149065233457.0_r8*a18&
               -115129891436076.0_r8*a14&
               -34747982243667.0_r8*a16-5763893468058.0_r8*a114-7234652444283.0_r8*a110&
               -17998039415961.0_r8*a112&
               +686547606456.0_r8*a116+468534132000.0_r8*a118-28570914236480.0_r8)/a127;
        z3k(5)=-(1.0_r8/246903552000.0_r8)*(-79763065718700.0_r8*a12&
                -16013022932655.0_r8*a18-130493283697140.0_r8*a14&
                -90942208981260.0_r8*a16-10200929842908.0_r8*a114+2512219926825.0_r8*a110&
                +154901667060.0_r8*a120&
               -11180851682679.0_r8*a112-2208929530140.0_r8*a116+354198515580.0_r8*a118&
               -18298525303660.0_r8)/a130;
        z3k(6)=(1.0_r8/5719767685632000.0_r8)*(-5618405733128588400.0_r8*a12-5882173102258447995.0_r8*a18&
               -11729493633609952200.0_r8*a14-12103105143601672200.0_r8*a16-882178514895765087.0_r8*a114&
               -911634009726294855.0_r8*a110+15298094789389440.0_r8*a120-419312125712563083.0_r8*a112&
               -491823047868672312.0_r8*a116-75180895923051360.0_r8*a118+4782670826169600.0_r8*a122&
               -1078217577276585200.0_r8)/a133;
        z3k(7)=-(1.0_r8/39598391669760000.0_r8)*(-109241601084281219400.0_r8*a12&
               -273755617744412105280.0_r8*a18&
               -275531080361356997100.0_r8*a14-369234543309244628400.0_r8*a16-16629926930547479010.0_r8*a114&
               -104248218659100772410.0_r8*a110-700142882381607420.0_r8*a120-21434621063709196035.0_r8*a112&
               -16119995582024586972.0_r8*a116-6389833770645023520.0_r8*a118+178205363298894600.0_r8*a122&
               +42471918559639800.0_r8*a124-18066934415886340600.0_r8)/a136;
        z3k(8)=(1.0_r8/4118232733655040000.0_r8)*(-30058294461554586318000.0_r8*a12&
               -140901213435715105015500.0_r8*a18&
               -88530627902844856856400.0_r8*a14-144616754379461084721000.0_r8*a16-7905053018941478835435.0_r8*a114&
               -81513990419502162468375.0_r8*a110-1152647132166760552560.0_r8*a120-27145640835037997266995.0_r8*a112&
               +5500015182539136000.0_r8*a126-5904155800054927152621.0_r8*a116-3802193002950503796882.0_r8*a118&
               -89958010436156592000.0_r8*a122+29001477679986355200.0_r8*a124-4375459586172885448000.0_r8)/a139;
        z3k(9)=-(1.0_r8/4633011825361920000.0_r8)*(-85391848625055033024600.0_r8*a12&
               -657624808266209998555770.0_r8*a18-286937627528945308701240.0_r8*a14&
               -550253819731820080247280.0_r8*a16-76885144694468451642015.0_r8*a114&
               -501113443972115359383570.0_r8*a110-8731660442937289348014.0_r8*a120&
               -240951516125845943628639.0_r8*a112&
               +48345076518821304480.0_r8*a126-28316641973836032929619.0_r8*a116-18073829022079601109456.0_r8*a118&
               -2112862299988788781488.0_r8*a122-113834507284119452040.0_r8*a124+7524133830529617600.0_r8*a128&
               -11114038570390896171200.0_r8)/a142;
        z3k(10)=(1.0_r8/1323188177323364352000000.0_r8)*(-59322584268580541162995440000.0_r8*a12&
                -691457486921045899947309280500.0_r8*a18-223595825592435993168270246000.0_r8*a14&
                -490908387679243133903442594000.0_r8*a16-171020750418261506307267347625.0_r8*a114&
                -647976462663409836197469228000.0_r8*a110-12206949275396521902047114736.0_r8*a120&
                -406244910341218086344121040125.0_r8*a112-32295304702451043090432000.0_r8*a126&
                -55282508132808557750323745025.0_r8*a116-22729370374617460059789308025.0_r8*a118&
                -4703198025098354603242138200.0_r8*a122-934171582904834946460056000.0_r8*a124&
                +19604322782640697627008000.0_r8*a128+2563130465125906185216000.0_r8*a130&
                -6988238080689917196734008000.0_r8)/a145;
       ELSE
         z1k(0)=-(1.0_r8/6.0_r8)*(3.0_r8*a12-1.0_r8)/a13;
         z1k(1)=(1.0_r8/36.0_r8)*(9*a14-15.0_r8*a12+7.0_r8)/a16;
         z1k(2)=-(1.0_r8/3240.0_r8)*(540.0_r8*a16-1917.0_r8*a14+2205.0_r8*a12-830.0_r8)/a19;
         z1k(3)=(1.0_r8/6480.0_r8)*(810.0_r8*a18-4698.0_r8*a16+9261.0_r8*a14-7710.0_r8*a12&
                +2330.0_r8)/a112;
         z1k(4)=-(1.0_r8/181440.0_r8)*(381360.0_r8*a12-591507.0_r8*a14+438165.0_r8*a16&
                -150984.0_r8*a18+18144.0_r8*a110-95228.0_r8)/a115;
         z1k(5)=(1.0_r8/3061800.0_r8)*(-11369925.0_r8*a12+21709485.0_r8*a14-21255885.0_r8*a16&
               +11073267.0_r8*a18-2821230.0_r8*a110+255150.0_r8*a112+2409050.0_r8)/a118;
         z1k(6)=-(1.0_r8/32659200.0_r8)*(213857700.0_r8*a12-485973810.0_r8*a14+593130195.0_r8*a16&
               -415312191.0_r8*a18+163829628.0_r8*a110-32587920.0_r8*a112+2332800.0_r8*a114&
               -39276200.0_r8)/a121;
         z1k(7)=(1.0_r8/9797760.0_r8)*(-112887600.0_r8*a12+297838548.0_r8*a14-436001580.0_r8*a16&
               +384287733.0_r8*a18-206353980.0_r8*a110+64657116.0_r8*a112-10429560.0_r8*a114+612360.0_r8*a116&
               +18277000.0_r8)/a124;
         z1k(8)=-(1.0_r8/6983843328000.0_r8)*(141307987590000.0_r8*a12-424805883301800.0_r8*a14&
               +725840751913200.0_r8*a16-772229389384035.0_r8*a18+526083694083585.0_r8*a110&
               -226625680966608.0_r8*a112+58331119524480.0_r8*a114-7846872019200.0_r8*a116&
               +387991296000.0_r8*a118-20443712366000.0_r8)/a127;
         z1k(9)=(1.0_r8/127309644000.0_r8)*(-4516095853500.0_r8*a12+15243943695750.0_r8*a14&
               -29792346039225.0_r8*a16+37162390768965.0_r8*a18-30723150555405.0_r8*a110&
               +16899523555761.0_r8*a112-6027938529228.0_r8*a114+1306793604780.0_r8*a116&
               -149793206850.0_r8*a118+6365482200.0_r8*a120+590307086600.0_r8)/a130;
         z1k(10)=-(1.0_r8/169474598092800.0_r8)*(10525154504264400.0_r8*a12-39423889481296320.0_r8*a14&
                +86765073240966120.0_r8*a16-124208571575574486.0_r8*a18+120861414343735065.0_r8*a110&
                -81043795131556257.0_r8*a112+37111541665177980.0_r8*a114-11221784750483424.0_r8*a116&
                +2089073823123456.0_r8*a118-207552534174720.0_r8*a120+7703390822400.0_r8*a122&
                -1254367500495200.0_r8)/a133;
         z1k(11)=(1.0_r8/127673385840000.0_r8)*(-13865529287610000.0_r8*a12+57079362420580500.0_r8*a14&
                -139720236603948000.0_r8*a16+225822732500400750.0_r8*a18-252922506605618850.0_r8*a110&
                +200288672908829505.0_r8*a112-112243862958457230.0_r8*a114+43771264478848998.0_r8*a116&
                -11426583533685840.0_r8*a118+1855110460918950.0_r8*a120-161971347951000.0_r8*a122&
                +5319724410000.0_r8*a124+1518227841130000.0_r8)/a136;
         z1k(12)=-(1.0_r8/61010855313408000.0_r8)*(11575128683790672000.0_r8*a12&
                -51951107748591651600.0_r8*a14+140021157313351314000.0_r8*a16&
                -252228200979738683700.0_r8*a18+319682393781675724800.0_r8*a110&
                -292147632621771581655.0_r8*a112+193923227055129156585.0_r8*a114&
                -92862978773912211384.0_r8*a116+31385770990121594592.0_r8*a118-7177250648648367360.0_r8*a120&
                +1029090082244774400.0_r8*a122-79870195141632000.0_r8*a124+2346571358208000.0_r8*a126&
                -1172073510107192000.0_r8)/a139;
         z2k(0)=-(1.0_r8/3240.0_r8)*(675.0_r8*a16-999.0_r8*a14+225.0_r8*a12+155)/a19;
         z2k(1)=(1.0_r8/2592.0_r8)*(891.0_r8*a18-2808.0_r8*a16+2403.0_r8*a14+12.0_r8*a12-505.0_r8)/a112;
         z2k(2)=-(1.0_r8/816480.0_r8)*(360612.0_r8*a110-1847367.0_r8*a18+2973159.0_r8*a16&
                -1409877.0_r8*a14-529305.0_r8*a12+450646.0_r8)/a115;
         z2k(3)=(1.0_r8/2099520.0_r8)*(1086210.0_r8*a112-7963596.0_r8*a110+19603053.0_r8*a18&
                -19160091.0_r8*a16+3104595.0_r8*a14+6145380.0_r8*a12-2817130.0_r8)/a118;
         z2k(4)=-(1.0_r8/97977600.0_r8)*(56687040.0_r8*a114-551882160.0_r8*a112+1875163689.0_r8*a110&
                -2848377393.0_r8*a18+1668615795.0_r8*a16+420389865.0_r8*a14&
                -914516400.0_r8*a12+293913460.0_r8)/a121;
         z2k(5)=(1.0_r8/1763596800.0_r8)*(1110514860.0_r8*a116-13665128328.0_r8*a114&
               +60227454744.0_r8*a112-126426509100.0_r8*a110+126738043767.0_r8*a18&
               -31045914900.0_r8*a16-50559078024.0_r8*a14+44844180360.0_r8*a12-11223491860.0_r8)/a124;
         z2k(6)=-(1.0_r8/1745960832000.0_r8)*(1175752195200.0_r8*a118-17657759442960.0_r8*a116&
               +96804669865860.0_r8*a114-262891435135923.0_r8*a112+376197022690335.0_r8*a110&
               -237036075843105.0_r8*a18-54081767113605.0_r8*a16+185058317296710.0_r8*a14&
               -110315414573100.0_r8*a12+22746727187800.0_r8)/a127;
         z2k(7)=(1.0_r8/185177664000.0_r8)*(131757214680.0_r8*a120-2352880625760.0_r8*a118&
               +15566372327340.0_r8*a116-52398360836136.0_r8*a114+98321547399957.0_r8*a112&
               -97554313727955.0_r8*a110+26026236868275.0_r8*a18+49064508958500.0_r8*a16&
               -59349145009140.0_r8*a14+27358014415200.0_r8*a12-4813736065880.0_r8)/a130;
         z2k(8)=-(1.0_r8/4575814148505600.0_r8)*(3409988483865600.0_r8*a122-70961952197364480.0_r8*a120&
               +553703220184315968.0_r8*a118-2241437046929907840.0_r8*a116+5244963547824369369.0_r8*a114&
               -7095472642735153677.0_r8*a112+4486300046287926867.0_r8*a110+1337585653909522017.0_r8*a18&
               -4915423617047143200.0_r8*a16+3978271574265348360.0_r8*a14-1513747067086254000.0_r8*a12&
               +232808295847072400.0_r8)/a133;
         z2k(9)=(1.0_r8/9899597917440000.0_r8)*(7676034955974000.0_r8*a124-183183385919954400.0_r8*a122&
               +1655721130272604560.0_r8*a120-7880259888799266600.0_r8*a118+22235847209559202842.0_r8*a116&
               -38230186388999393250.0_r8*a114+36747733201815133575.0_r8*a112-8636628178136012670.0_r8*a110&
               -24623457047307646110.0_r8*a18+34546349053518271200.0_r8*a16-21767456292393873600.0_r8*a14&
               +7098994800259362000.0_r8*a12-971150252248330000.0_r8)/a136;
         z2k(10)=-(1.0_r8/1372744244551680000.0_r8)*(1101852389965824000.0_r8*a126-29762759848147430400.0_r8*a124&
               +307134113531054227200.0_r8*a122-1688681289865388980320.0_r8*a120+5608250074097827426764.0_r8*a118&
               -11746912410175672315053.0_r8*a116+15025452822410836851045.0_r8*a114-9019215315411973606935.0_r8*a112&
               -4425894762564289432275.0_r8*a110+14384768308354215936450.0_r8*a18-13855276123879399209000.0_r8*a16&
               +7263110850633113767200.0_r8*a14-2080655317262245542000.0_r8*a12+256579957380860372000.0_r8)/a139;
       ENDIF
       s=z1k(0); t= 1; k= 1;
       DO WHILE ((abs(t)>1.0e-15_r8).AND.(k < 13)) 
         t=z1k(k)*zk(k); s=s+t; k= k+1 
       ENDDO
       zeta1=s;
       s=z2k(0); t= 1; k= 1; 
       DO WHILE ((abs(t)>1.0e-15_r8).AND.(k < 11)) 
         t=z2k(k)*zk(k); s=s+t; k= k+1 
       ENDDO
       zeta2=s;
       IF (iflag==1) THEN      
         s=z3k(0); t= 1; k= 1;
         DO WHILE ((abs(t)>1.0e-15_r8).AND.(k < 11)) 
           t=z3k(k)*zk(k); s=s+t; k= k+1  
         ENDDO
         zeta3=s;       
       ELSE
         zeta3=0.0_r8 
       ENDIF   
       IF (abs(zeta2)>abs(10*zeta1)) zeta2=0.0_r8
       IF (abs(zeta3)>abs(5*zeta2)) zeta3=0.0_r8
      ELSE
        zeta2=0.0_r8
        zeta3=0.0_r8
      ENDIF
     ENDIF
     IF (abs(zeta1)<epss) zeta1=0.0_r8
     IF (abs(zeta2)<epss) zeta2=0.0_r8
     IF (abs(zeta3)<epss) zeta3=0.0_r8
     END SUBROUTINE zeta3xy

    FUNCTION fc(pnu,z)
!-----------------------------------------------------------
!   Evaluation of the cf for the ratio Ipnu(z)/Ipnu-1(z) 
!   We use Lentz-Thompson algorithm.
!-----------------------------------------------------------
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: pnu
    REAL(r8), INTENT(IN) :: z
    REAL(r8) :: fc,b,a,c0,d0,delta
    INTEGER :: m   
    m=0
    b=2.0_r8*pnu/z
    a=1.0_r8
    fc=dwarf
    c0=fc
    d0=0.0_r8
    delta=0.0_r8
    DO WHILE((ABS(delta-1.0_r8)>epss).AND.(m<500)) 
      d0=b+a*d0
      IF(ABS(d0)<dwarf) d0=dwarf
      c0=b+a/c0
      IF(ABS(c0)<dwarf) c0=dwarf
      d0=1.0_r8/d0
      delta=c0*d0
      fc=fc*delta
      m=m+1
      a=1.0_r8  
      b=2.0_r8*(pnu+m)/z  
    ENDDO 
    END FUNCTION fc

    FUNCTION nmax(mu,x,y)
   !---------------------------------------------------------
   ! Computes a starting value for the backward summation of 
   ! the series in pmuxyser
   !---------------------------------------------------------
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(IN) :: y
    REAL(r8) :: lneps, c, n, n1
    INTEGER :: nmax
    lneps=-36.0_r8;
    c=-mu+y-mu*log(y)+lneps;
    n=10+2*(-mu+sqrt(mu*mu+4*x*y)); n1= 0;
    DO WHILE ((abs(n-n1)>1).AND.(n>0)) 
      n1=n;
      n=-(log(mu+n)*mu-2*n+c)/(log(n/(x*y))+log(mu+n));
    ENDDO
    IF (n<0) n=0
    nmax=int(n)+1
    END FUNCTION nmax

    FUNCTION factor(x,n)
    IMPLICIT NONE
    REAL(r8) :: x, facto, factor
    INTEGER :: n, i
    facto=1
    DO i=1,n
      facto=facto*(x/i)
    ENDDO
    factor=facto
    END FUNCTION factor

    FUNCTION pol(fjkm,d,v)
    IMPLICIT NONE
    REAL(r8) :: pol,v,s,fjkm(0:32)
    INTEGER :: d,m
    m=d; s= fjkm(d);
    DO WHILE (m > 0)
      m=m-1; s=s*v + fjkm(m)
    ENDDO
    pol=s
    END FUNCTION pol

    SUBROUTINE fjkproc16(u,fjk)
    IMPLICIT NONE
    REAL(r8) :: u
    REAL(r8) :: fjk(0:16,0:16), fjkm(0:32) 
    REAL(r8) :: un(0:64), v
    INTEGER :: j,k,d,n
    un(1)= u; v= u*u; un(2)= v; 
    DO n=2,64 
      un(n)= u*un(n-1);
    ENDDO	
    fjk(0,0)=1.0_r8;
    fjkm(0)=0.50000000000000000000_r8;
    fjkm(1)=0.16666666666666666667_r8;
    j= 1; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.12500000000000000000_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=0.20833333333333333333_r8;
    j= 2; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.62500000000000000000e-1_r8;
    fjkm(1)=-0.54166666666666666667e-1_r8;
    fjkm(2)=-0.31250000000000000000_r8;
    fjkm(3)=0.28935185185185185185_r8;
    j= 3; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.39062500000000000000e-1_r8;
    fjkm(1)=0.83333333333333333333e-1_r8;
    fjkm(2)=0.36631944444444444444_r8;
    fjkm(3)=-0.83333333333333333333_r8;
    fjkm(4)=0.42390046296296296296_r8;
    j= 4; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.27343750000000000000e-1_r8;
    fjkm(1)=-0.10145089285714285714_r8;
    fjkm(2)=-0.38281250000000000000_r8;
    fjkm(3)=1.6061921296296296296_r8;
    fjkm(4)=-1.7903645833333333333_r8;
    fjkm(5)=0.64144483024691358025_r8;
    j= 5; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.20507812500000000000e-1_r8;
    fjkm(1)=0.11354166666666666667_r8;
    fjkm(2)=0.36983072916666666667_r8;
    fjkm(3)=-2.5763888888888888889_r8;
    fjkm(4)=4.6821108217592592593_r8;
    fjkm(5)=-3.5607638888888888889_r8;
    fjkm(6)=0.99199861754115226337_r8;
    j= 6; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.16113281250000000000e-1_r8;
    fjkm(1)=-0.12196955605158730159_r8;
    fjkm(2)=-0.33297526041666666667_r8;
    fjkm(3)=3.7101836350859788360_r8;
    fjkm(4)=-9.7124626253858024691_r8;
    fjkm(5)=11.698143727494855967_r8;
    fjkm(6)=-6.8153513213734567901_r8;
    fjkm(7)=1.5583573120284636488_r8;
    j= 7; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.13092041015625000000e-1_r8;
    fjkm(1)=0.12801339285714285714_r8;
    fjkm(2)=0.27645252046130952381_r8;
    fjkm(3)=-4.9738777281746031746_r8;
    fjkm(4)=17.501935105096726190_r8;
    fjkm(5)=-29.549479166666666667_r8;
    fjkm(6)=26.907133829250257202_r8;
    fjkm(7)=-12.754267939814814815_r8;
    fjkm(8)=2.4771798425577632030_r8;
    j= 8; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.10910034179687500000e-1_r8;
    fjkm(1)=-0.13242874035415539322_r8;
    fjkm(2)=-0.20350690569196428571_r8;
    fjkm(3)=6.3349384739790013228_r8;
    fjkm(4)=-28.662114811111800044_r8;
    fjkm(5)=63.367483364421434083_r8;
    fjkm(6)=-79.925485618811085391_r8;
    fjkm(7)=58.757341382271304870_r8;
    fjkm(8)=-23.521455678429623200_r8;
    fjkm(9)=3.9743166454849898231_r8;
    j= 9; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.92735290527343750000e-2_r8;
    fjkm(1)=0.13569064670138888889_r8;
    fjkm(2)=0.11668911254144670021_r8;
    fjkm(3)=-7.7625075954861111111_r8;
    fjkm(4)=43.784562625335567773_r8;
    fjkm(5)=-121.31910738398368607_r8;
    fjkm(6)=198.20121981295421734_r8;
    fjkm(7)=-200.43673900016432327_r8;
    fjkm(8)=123.80342757950794259_r8;
    fjkm(9)=-42.937783937667895519_r8;
    fjkm(10)=6.4238224989853211488_r8;
    j= 10; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.80089569091796875000e-2_r8;
    fjkm(1)=-0.13811212730852318255_r8;
    fjkm(2)=-0.18036238655211433532e-1_r8;
    fjkm(3)=9.2275853445797140866_r8;
    fjkm(4)=-63.433189058657045718_r8;
    fjkm(5)=213.60596888977804302_r8;
    fjkm(6)=-432.96183396641609600_r8;
    fjkm(7)=563.58282810729226948_r8;
    fjkm(8)=-476.64858951490111802_r8;
    fjkm(9)=254.12602383553942414_r8;
    fjkm(10)=-77.797248335368675787_r8;
    fjkm(11)=10.446593930548512362_r8;
    j= 11; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.70078372955322265625e-2_r8;
    fjkm(1)=0.13990718736965074856_r8;
    fjkm(2)=-0.90802493534784075207e-1_r8;
    fjkm(3)=-10.703046719402920575_r8;
    fjkm(4)=88.139055705916160082_r8;
    fjkm(5)=-352.55365414896970073_r8;
    fjkm(6)=860.26747669490580229_r8;
    fjkm(7)=-1381.3884907075539460_r8;
    fjkm(8)=1497.5262381375579615_r8;
    fjkm(9)=-1089.5695395426785795_r8;
    fjkm(10)=511.32054028583482617_r8;
    fjkm(11)=-140.15612725058882506_r8;
    fjkm(12)=17.075450695147740963_r8;
    j= 12; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.61992406845092773438e-2_r8;
    fjkm(1)=-0.14122658948520402530_r8;
    fjkm(2)=0.20847570003254474385_r8;
    fjkm(3)=12.163573370672875144_r8;
    fjkm(4)=-118.39689039212288489_r8;
    fjkm(5)=552.67487991757989471_r8;
    fjkm(6)=-1587.5976792806460534_r8;
    fjkm(7)=3052.8623335041016490_r8;
    fjkm(8)=-4067.0706975337409188_r8;
    fjkm(9)=3781.4312193762993828_r8;
    fjkm(10)=-2415.5306966669781670_r8;
    fjkm(11)=1012.7298787459738084_r8;
    fjkm(12)=-251.37116645382357870_r8;
    fjkm(13)=28.031797071713952493_r8;
    j= 13; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.55350363254547119141e-2_r8;
    fjkm(1)=0.14217921713372686883_r8;
    fjkm(2)=-0.33386405994128191388_r8;
    fjkm(3)=-13.585546133738642981_r8;
    fjkm(4)=154.66282442015249412_r8;
    fjkm(5)=-830.71069930083407076_r8;
    fjkm(6)=2761.0291182562342601_r8;
    fjkm(7)=-6219.8351157050681259_r8;
    fjkm(8)=9888.1927799238643295_r8;
    fjkm(9)=-11266.694472611704499_r8;
    fjkm(10)=9175.5017581920039296_r8;
    fjkm(11)=-5225.7429703251833306_r8;
    fjkm(12)=1980.4053574007652015_r8;
    fjkm(13)=-449.21570290311749301_r8;
    fjkm(14)=46.189888661376921323_r8;
    j= 14; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.49815326929092407227e-2_r8;
    fjkm(1)=-0.14284537645361756992_r8;
    fjkm(2)=0.46603144339556328292_r8;
    fjkm(3)=14.946922390174830635_r8;
    fjkm(4)=-197.35300817536964730_r8;
    fjkm(5)=1205.6532423474201960_r8;
    fjkm(6)=-4572.9473467250314847_r8;
    fjkm(7)=11865.183572985041892_r8;
    fjkm(8)=-22026.784993357215819_r8;
    fjkm(9)=29873.206689727991728_r8;
    fjkm(10)=-29749.925047590507307_r8;
    fjkm(11)=21561.076414337110462_r8;
    fjkm(12)=-11081.438701085531999_r8;
    fjkm(13)=3832.1051284526998677_r8;
    fjkm(14)=-800.40791995840375064_r8;
    fjkm(15)=76.356879052900946470_r8;
    j= 15; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.45145140029489994049e-2_r8;
    fjkm(1)=0.14328537051348750262_r8;
    fjkm(2)=-0.60418804953366830816_r8;
    fjkm(3)=-16.227111168548122706_r8;
    fjkm(4)=246.84286168977111296_r8;
    fjkm(5)=-1698.7528990888950203_r8;
    fjkm(6)=7270.2387180078329370_r8;
    fjkm(7)=-21434.839860240815288_r8;
    fjkm(8)=45694.866035689911070_r8;
    fjkm(9)=-72195.530107556687632_r8;
    fjkm(10)=85409.022842807474925_r8;
    fjkm(11)=-75563.234444869051891_r8;
    fjkm(12)=49344.501227769590532_r8;
    fjkm(13)=-23110.149147008741710_r8;
    fjkm(14)=7349.7909384681412957_r8;
    fjkm(15)=-1422.6485707704091767_r8;
    fjkm(16)=126.58493346342458430_r8;
    j= 16; k=0; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.12500000000000000000_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-0.20833333333333333333_r8;
    j= 0; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.62500000000000000000e-1_r8;
    fjkm(1)=0.14583333333333333333_r8;
    fjkm(2)=0.52083333333333333333_r8;
    fjkm(3)=-0.65972222222222222222_r8;
    j= 1; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.46875000000000000000e-1_r8;
    fjkm(1)=-0.25000000000000000000_r8;
    fjkm(2)=-0.69791666666666666667_r8;
    fjkm(3)=2.5000000000000000000_r8;
    fjkm(4)=-1.6059027777777777778_r8;
    j= 2; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.39062500000000000000e-1_r8;
    fjkm(1)=0.34218750000000000000_r8;
    fjkm(2)=0.72916666666666666667_r8;
    fjkm(3)=-5.6712962962962962963_r8;
    fjkm(4)=8.1640625000000000000_r8;
    fjkm(5)=-3.5238233024691358025_r8;
    j= 3; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.34179687500000000000e-1_r8;
    fjkm(1)=-0.42708333333333333333_r8;
    fjkm(2)=-0.59798177083333333333_r8;
    fjkm(3)=10.208333333333333333_r8;
    fjkm(4)=-24.385308159722222222_r8;
    fjkm(5)=22.482638888888888889_r8;
    fjkm(6)=-7.3148750964506172840_r8;
    j= 4; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.30761718750000000000e-1_r8;
    fjkm(1)=0.50665457589285714286_r8;
    fjkm(2)=0.29326171875000000000_r8;
    fjkm(3)=-16.044663008432539683_r8;
    fjkm(4)=56.156774450231481481_r8;
    fjkm(5)=-82.372823832947530864_r8;
    fjkm(6)=56.160933883101851852_r8;
    fjkm(7)=-14.669405462319958848_r8;
    j= 5; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.28198242187500000000e-1_r8;
    fjkm(1)=-0.58203125000000000000_r8;
    fjkm(2)=0.19236328125000000000_r8;
    fjkm(3)=23.032335069444444444_r8;
    fjkm(4)=-110.33599717881944444_r8;
    fjkm(5)=227.74508101851851852_r8;
    fjkm(6)=-243.01300676761831276_r8;
    fjkm(7)=131.66775173611111111_r8;
    fjkm(8)=-28.734679254811814129_r8;
    j= 6; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.26184082031250000000e-1_r8;
    fjkm(1)=0.65396069723462301587_r8;
    fjkm(2)=-0.86386369977678571429_r8;
    fjkm(3)=-30.956497628348214286_r8;
    fjkm(4)=194.54890778287588183_r8;
    fjkm(5)=-527.74348743041776896_r8;
    fjkm(6)=780.79702721113040123_r8;
    fjkm(7)=-656.29672278886959877_r8;
    fjkm(8)=295.22178492918917181_r8;
    fjkm(9)=-55.334928257039108939_r8;
    j= 7; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.24547576904296875000e-1_r8;
    fjkm(1)=-0.72297712053571428571_r8;
    fjkm(2)=1.7246239871070498512_r8;
    fjkm(3)=39.546341145833333333_r8;
    fjkm(4)=-316.99617299397786458_r8;
    fjkm(5)=1081.5824590773809524_r8;
    fjkm(6)=-2074.4037171674994144_r8;
    fjkm(7)=2398.8177766525205761_r8;
    fjkm(8)=-1664.7533222350236974_r8;
    fjkm(9)=640.37285196437757202_r8;
    fjkm(10)=-105.19241070496638071_r8;
    j= 8; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.23183822631835937500e-1_r8;
    fjkm(1)=0.78948182770700165720_r8;
    fjkm(2)=-2.7769457196432446677_r8;
    fjkm(3)=-48.483511725054617511_r8;
    fjkm(4)=486.26944746794524016_r8;
    fjkm(5)=-2023.8687997794445650_r8;
    fjkm(6)=4819.5203340475451309_r8;
    fjkm(7)=-7173.8455521540386687_r8;
    fjkm(8)=6815.5497547693867415_r8;
    fjkm(9)=-4029.1488859965138299_r8;
    fjkm(10)=1353.9765257894256969_r8;
    fjkm(11)=-197.95866455017786514_r8;
    j= 9; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.22024631500244140625e-1_r8;
    fjkm(1)=-0.85378706190321180556_r8;
    fjkm(2)=4.0223697649378354858_r8;
    fjkm(3)=57.408728524667245370_r8;
    fjkm(4)=-711.17788174487525874_r8;
    fjkm(5)=3529.7027186963924024_r8;
    fjkm(6)=-10126.360073656459287_r8;
    fjkm(7)=18593.571833843032106_r8;
    fjkm(8)=-22636.974191769862737_r8;
    fjkm(9)=18256.758136740546277_r8;
    fjkm(10)=-9401.7390963482140877_r8;
    fjkm(11)=2805.1309521324368189_r8;
    fjkm(12)=-369.51173382133760772_r8;
    j= 10; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.21023511886596679688e-1_r8;
    fjkm(1)=0.91614229084450603921_r8;
    fjkm(2)=-5.4618897365663646792_r8;
    fjkm(3)=-65.927090730899790043_r8;
    fjkm(4)=1000.5846459432155138_r8;
    fjkm(5)=-5819.5104958244309663_r8;
    fjkm(6)=19669.512303383909626_r8;
    fjkm(7)=-43248.109984766956219_r8;
    fjkm(8)=64686.833219925562541_r8;
    fjkm(9)=-66644.721592787005810_r8;
    fjkm(10)=46700.224876576248105_r8;
    fjkm(11)=-21305.241783200054197_r8;
    fjkm(12)=5716.0564388863858560_r8;
    fjkm(13)=-685.13376643364457482_r8;
    j= 11; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.20147532224655151367e-1_r8;
    fjkm(1)=-0.97675104265089158888_r8;
    fjkm(2)=7.0960977732520869186_r8;
    fjkm(3)=73.612404446931816840_r8;
    fjkm(4)=-1363.2529599857205103_r8;
    fjkm(5)=9163.5749605651064785_r8;
    fjkm(6)=-35864.832027517679572_r8;
    fjkm(7)=92376.947946883135629_r8;
    fjkm(8)=-164834.83784046136147_r8;
    fjkm(9)=208040.74214864238663_r8;
    fjkm(10)=-185789.85543054601897_r8;
    fjkm(11)=115101.05980498683057_r8;
    fjkm(12)=-47134.497747406170101_r8;
    fjkm(13)=11488.455724263405622_r8;
    fjkm(14)=-1263.2564781342309572_r8;
    j= 12; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.19372627139091491699e-1_r8;
    fjkm(1)=1.0357822247248985273_r8;
    fjkm(2)=-8.9252867409545555669_r8;
    fjkm(3)=-80.010760848208515926_r8;
    fjkm(4)=1807.7010112763722305_r8;
    fjkm(5)=-13886.239779926939526_r8;
    fjkm(6)=62074.025784691064765_r8;
    fjkm(7)=-184145.36258906485362_r8;
    fjkm(8)=383582.03973738516554_r8;
    fjkm(9)=-576038.14802241650407_r8;
    fjkm(10)=628833.57177487054480_r8;
    fjkm(11)=-495573.34466362548232_r8;
    fjkm(12)=275136.26556037481059_r8;
    fjkm(13)=-102207.94135045741701_r8;
    fjkm(14)=22823.518594320784899_r8;
    fjkm(15)=-2318.1664194368241801_r8;
    j= 13; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.18680747598409652710e-1_r8;
    fjkm(1)=-1.0933780262877533843_r8;
    fjkm(2)=10.949523517716476548_r8;
    fjkm(3)=84.643531757174997082_r8;
    fjkm(4)=-2342.0651134541361650_r8;
    fjkm(5)=20369.770814493525849_r8;
    fjkm(6)=-102837.36670370684414_r8;
    fjkm(7)=346648.70123159565251_r8;
    fjkm(8)=-828961.75261733863230_r8;
    fjkm(9)=1449716.2069081751493_r8;
    fjkm(10)=-1879301.2063630471735_r8;
    fjkm(11)=1807331.1927210534022_r8;
    fjkm(12)=-1274496.0479553760313_r8;
    fjkm(13)=641015.88735632964269_r8;
    fjkm(14)=-217895.35698164739918_r8;
    fjkm(15)=44894.191761385746325_r8;
    fjkm(16)=-4236.6734164587391641_r8;
    j= 14; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.18058056011795997620e-1_r8;
    fjkm(1)=1.1496596058234620572_r8;
    fjkm(2)=-13.168702574894088581_r8;
    fjkm(3)=-87.009904621222723801_r8;
    fjkm(4)=2973.9704746271364310_r8;
    fjkm(5)=-29057.863221639334277_r8;
    fjkm(6)=164134.77797509230729_r8;
    fjkm(7)=-621775.71008157787326_r8;
    fjkm(8)=1684395.4811437516714_r8;
    fjkm(9)=-3374103.7733925392782_r8;
    fjkm(10)=5084254.0101088150007_r8;
    fjkm(11)=-5797111.7426650438469_r8;
    fjkm(12)=4981685.5893221501493_r8;
    fjkm(13)=-3178510.0651440248351_r8;
    fjkm(14)=1461172.7927229457558_r8;
    fjkm(15)=-457795.06653758949338_r8;
    fjkm(16)=87552.178162658627517_r8;
    fjkm(17)=-7715.5318619797584859_r8;
    j= 15; k=1; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.70312500000000000000e-1_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-0.40104166666666666667_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=0.33420138888888888889_r8;
    j= 0; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.10546875000000000000_r8;
    fjkm(1)=0.15234375000000000000_r8;
    fjkm(2)=1.4036458333333333333_r8;
    fjkm(3)=-1.6710069444444444444_r8;
    fjkm(4)=-1.8381076388888888889_r8;
    fjkm(5)=2.0609085648148148148_r8;
    j= 1; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.13183593750000000000_r8;
    fjkm(1)=-0.42187500000000000000_r8;
    fjkm(2)=-2.8623046875000000000_r8;
    fjkm(3)=8.0208333333333333333_r8;
    fjkm(4)=1.0777994791666666667_r8;
    fjkm(5)=-14.036458333333333333_r8;
    fjkm(6)=8.0904586226851851852_r8;
    j= 2; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.15380859375000000000_r8;
    fjkm(1)=0.79892578125000000000_r8;
    fjkm(2)=4.5903320312500000000_r8;
    fjkm(3)=-22.751985677083333333_r8;
    fjkm(4)=14.934624565972222222_r8;
    fjkm(5)=42.526662567515432099_r8;
    fjkm(6)=-65.691460503472222222_r8;
    fjkm(7)=25.746658387988683128_r8;
    j= 3; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.17303466796875000000_r8;
    fjkm(1)=-1.2773437500000000000_r8;
    fjkm(2)=-6.3615722656250000000_r8;
    fjkm(3)=50.011935763888888889_r8;
    fjkm(4)=-73.559339735243055556_r8;
    fjkm(5)=-70.026331018518518519_r8;
    fjkm(6)=271.34066056616512346_r8;
    fjkm(7)=-242.71375868055555556_r8;
    fjkm(8)=72.412718470695087449_r8;
    j= 4; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.19033813476562500000_r8;
    fjkm(1)=1.8524126325334821429_r8;
    fjkm(2)=7.9220947265625000000_r8;
    fjkm(3)=-94.174715169270833333_r8;
    fjkm(4)=221.09830050998263889_r8;
    fjkm(5)=13.578712293836805556_r8;
    fjkm(6)=-765.03722541714891975_r8;
    fjkm(7)=1204.5108913845486111_r8;
    fjkm(8)=-777.22725008740837191_r8;
    fjkm(9)=187.66711848589945559_r8;
    j= 5; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.20619964599609375000_r8;
    fjkm(1)=-2.5202636718750000000_r8;
    fjkm(2)=-8.9979495239257812500_r8;
    fjkm(3)=159.65856119791666667_r8;
    fjkm(4)=-527.02200527615017361_r8;
    fjkm(5)=337.21907552083333333_r8;
    fjkm(6)=1618.7873626708984375_r8;
    fjkm(7)=-4211.0382245852623457_r8;
    fjkm(8)=4434.1363497656886306_r8;
    fjkm(9)=-2259.2768162856867284_r8;
    fjkm(10)=458.84770992088928362_r8;
    j= 6; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.22092819213867187500_r8;
    fjkm(1)=3.2776113237653459821_r8;
    fjkm(2)=9.3003209795270647321_r8;
    fjkm(3)=-250.77115683984504175_r8;
    fjkm(4)=1087.8174260457356771_r8;
    fjkm(5)=-1404.3028911260910976_r8;
    fjkm(6)=-2563.9444452795962738_r8;
    fjkm(7)=11622.495969086321293_r8;
    fjkm(8)=-17934.344163614699543_r8;
    fjkm(9)=14479.313178892270800_r8;
    fjkm(10)=-6122.6397406024697386_r8;
    fjkm(11)=1074.0188194633057124_r8;
    j= 7; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.23473620414733886719_r8;
    fjkm(1)=-4.1216033935546875000_r8;
    fjkm(2)=-8.5290844099862234933_r8;
    fjkm(3)=371.57630452473958333_r8;
    fjkm(4)=-2030.4431650042155432_r8;
    fjkm(5)=3928.2498148600260417_r8;
    fjkm(6)=2472.6031768756442600_r8;
    fjkm(7)=-26784.706192883150077_r8;
    fjkm(8)=57707.467479758203765_r8;
    fjkm(9)=-65779.375284558624561_r8;
    fjkm(10)=43428.755429000357378_r8;
    fjkm(11)=-15731.921483001918296_r8;
    fjkm(12)=2430.2098720207426574_r8;
    j= 8; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.24777710437774658203_r8;
    fjkm(1)=5.0497246547178788619_r8;
    fjkm(2)=6.3753494648706345331_r8;
    fjkm(3)=-525.77763195628211612_r8;
    fjkm(4)=3515.3901490500364354_r8;
    fjkm(5)=-9098.9465854134383025_r8;
    fjkm(6)=1501.4499341175879961_r8;
    fjkm(7)=52968.402569427411743_r8;
    fjkm(8)=-156962.17039551999834_r8;
    fjkm(9)=237710.55444046526561_r8;
    fjkm(10)=-217889.76091367761240_r8;
    fjkm(11)=122183.02599420558225_r8;
    fjkm(12)=-38765.366765003690558_r8;
    fjkm(13)=5352.0219072834891857_r8;
    j= 9; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.26016595959663391113_r8;
    fjkm(1)=-6.0597300529479980469_r8;
    fjkm(2)=-2.5233336282830660035_r8;
    fjkm(3)=716.61609867398701017_r8;
    fjkm(4)=-5739.3531518108567233_r8;
    fjkm(5)=18710.056678357368214_r8;
    fjkm(6)=-15227.052778022872591_r8;
    fjkm(7)=-90410.693429278463984_r8;
    fjkm(8)=374173.78606362103489_r8;
    fjkm(9)=-726678.85908967426430_r8;
    fjkm(10)=867690.63613487545936_r8;
    fjkm(11)=-669330.97296360068003_r8;
    fjkm(12)=326923.43164094028666_r8;
    fjkm(13)=-92347.975136788220980_r8;
    fjkm(14)=11528.702830431737704_r8;
    j= 10; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.27199168503284454346_r8;
    fjkm(1)=7.1495960926884537810_r8;
    fjkm(2)=-3.3482232633271774688_r8;
    fjkm(3)=-946.77886875050071077_r8;
    fjkm(4)=8937.5223755513871158_r8;
    fjkm(5)=-35337.290730230231910_r8;
    fjkm(6)=49459.110954680523755_r8;
    fjkm(7)=130172.10642681313787_r8;
    fjkm(8)=-799759.43622037315810_r8;
    fjkm(9)=1951258.3781149473349_r8;
    fjkm(10)=-2917568.5809228727915_r8;
    fjkm(11)=2902364.7717490216619_r8;
    fjkm(12)=-1939009.9106363828308_r8;
    fjkm(13)=839993.29007310418975_r8;
    fjkm(14)=-213947.07574365748020_r8;
    fjkm(15)=24380.364047003816130_r8;
    j= 11; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.28332467190921306610_r8;
    fjkm(1)=-8.3174842023230218268_r8;
    fjkm(2)=11.564968830087189329_r8;
    fjkm(3)=1218.3176746274854345_r8;
    fjkm(4)=-13385.506466376132701_r8;
    fjkm(5)=62540.301444639907312_r8;
    fjkm(6)=-122382.91365675340023_r8;
    fjkm(7)=-143089.29056273054521_r8;
    fjkm(8)=1554707.5273250397863_r8;
    fjkm(9)=-4718379.1202435790866_r8;
    fjkm(10)=8605382.3882172381886_r8;
    fjkm(11)=-10599895.885891960945_r8;
    fjkm(12)=9077838.6492033673420_r8;
    fjkm(13)=-5358306.8036736424269_r8;
    fjkm(14)=2087192.8797093735160_r8;
    fjkm(15)=-484205.51887813298356_r8;
    fjkm(16)=50761.444989589644273_r8;
    j= 12; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.29422177467495203018_r8;
    fjkm(1)=9.5617123863640260863_r8;
    fjkm(2)=-22.456083202924761739_r8;
    fjkm(3)=-1532.5752010328979216_r8;
    fjkm(4)=19400.899387993296528_r8;
    fjkm(5)=-105087.87893355958488_r8;
    fjkm(6)=262967.92353732330762_r8;
    fjkm(7)=61613.709412718183861_r8;
    fjkm(8)=-2770349.5737553264845_r8;
    fjkm(9)=10454840.708341905254_r8;
    fjkm(10)=-22841191.197018135498_r8;
    fjkm(11)=33883152.513403664925_r8;
    fjkm(12)=-35702419.891311431439_r8;
    fjkm(13)=26908321.887260639130_r8;
    fjkm(14)=-14241918.449935481857_r8;
    fjkm(15)=5042187.1772326832703_r8;
    fjkm(16)=-1074259.7908211056482_r8;
    fjkm(17)=104287.72699173731366_r8;
    j= 13; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.30472969519905745983_r8;
    fjkm(1)=-10.880732823367957230_r8;
    fjkm(2)=36.353598566849979929_r8;
    fjkm(3)=1890.1183163422473995_r8;
    fjkm(4)=-27344.503966626558254_r8;
    fjkm(5)=169206.00701162634518_r8;
    fjkm(6)=-515265.57929420780237_r8;
    fjkm(7)=248217.18376755761876_r8;
    fjkm(8)=4532205.9690912962460_r8;
    fjkm(9)=-21492951.405026820809_r8;
    fjkm(10)=55557247.108151935296_r8;
    fjkm(11)=-97285892.160889723465_r8;
    fjkm(12)=122607925.76741209887_r8;
    fjkm(13)=-113234217.41453912066_r8;
    fjkm(14)=76312885.330872101297_r8;
    fjkm(15)=-36634431.269047918012_r8;
    fjkm(16)=11891540.519643040965_r8;
    fjkm(17)=-2342835.9225964451203_r8;
    fjkm(18)=211794.47349942484210_r8;
    j= 14; k=2; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.73242187500000000000e-1_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-0.89121093750000000000_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=1.8464626736111111111_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-1.0258125964506172840_r8;
    j= 0; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.18310546875000000000_r8;
    fjkm(1)=0.23193359375000000000_r8;
    fjkm(2)=4.0104492187500000000_r8;
    fjkm(3)=-4.6045898437500000000_r8;
    fjkm(4)=-12.002007378472222222_r8;
    fjkm(5)=13.232982494212962963_r8;
    fjkm(6)=8.7194070698302469136_r8;
    fjkm(7)=-9.4032821341306584362_r8;
    j= 1; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.32043457031250000000_r8;
    fjkm(1)=-0.87890625000000000000_r8;
    fjkm(2)=-10.464160156250000000_r8;
    fjkm(3)=26.736328125000000000_r8;
    fjkm(4)=29.225667317708333333_r8;
    fjkm(5)=-103.40190972222222222_r8;
    fjkm(6)=17.131070360725308642_r8;
    fjkm(7)=92.323133680555555556_r8;
    fjkm(8)=-50.991434481899434156_r8;
    j= 2; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.48065185546875000000_r8;
    fjkm(1)=2.1109008789062500000_r8;
    fjkm(2)=21.025415039062500000_r8;
    fjkm(3)=-89.876334092881944444_r8;
    fjkm(4)=-15.284450954861111111_r8;
    fjkm(5)=411.04911024305555556_r8;
    fjkm(6)=-389.32152566792052469_r8;
    fjkm(7)=-293.92095419801311728_r8;
    fjkm(8)=567.72315884813850309_r8;
    fjkm(9)=-213.02470796338270176_r8;
    j= 3; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.66089630126953125000_r8;
    fjkm(1)=-4.0893554687500000000_r8;
    fjkm(2)=-36.043276468912760417_r8;
    fjkm(3)=229.67792968750000000_r8;
    fjkm(4)=-150.64704827202690972_r8;
    fjkm(5)=-1115.0236545138888889_r8;
    fjkm(6)=2175.9328758333936150_r8;
    fjkm(7)=-176.78170412165637860_r8;
    fjkm(8)=-2817.5643744553721654_r8;
    fjkm(9)=2651.5545930587705761_r8;
    fjkm(10)=-757.67687847693870206_r8;
    j= 4; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.85916519165039062500_r8;
    fjkm(1)=6.9690023149762834821_r8;
    fjkm(2)=55.378133392333984375_r8;
    fjkm(3)=-495.03058466109018477_r8;
    fjkm(4)=733.55991770426432292_r8;
    fjkm(5)=2262.8678469622576678_r8;
    fjkm(6)=-7898.5588740407684703_r8;
    fjkm(7)=5695.1407000317985629_r8;
    fjkm(8)=7718.7923309734328784_r8;
    fjkm(9)=-16089.784052072184022_r8;
    fjkm(10)=10424.896645958040967_r8;
    fjkm(11)=-2413.3719004256171872_r8;
    j= 5; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=1.0739564895629882812_r8;
    fjkm(1)=-10.898971557617187500_r8;
    fjkm(2)=-78.354169082641601562_r8;
    fjkm(3)=948.65802978515625000_r8;
    fjkm(4)=-2219.4645106141832140_r8;
    fjkm(5)=-3394.0875061035156250_r8;
    fjkm(6)=22215.581371235788604_r8;
    fjkm(7)=-30531.682191548916538_r8;
    fjkm(8)=-6945.0954169277301051_r8;
    fjkm(9)=63170.236743335697713_r8;
    fjkm(10)=-72433.473359744190134_r8;
    fjkm(11)=36368.490166893057699_r8;
    fjkm(12)=-7090.9841426397698721_r8;
    j= 6; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-1.3040900230407714844_r8;
    fjkm(1)=16.023568312327067057_r8;
    fjkm(2)=103.72332413083031064_r8;
    fjkm(3)=-1667.3156506674630301_r8;
    fjkm(4)=5406.6561977448034539_r8;
    fjkm(5)=2872.9832533515445770_r8;
    fjkm(6)=-52104.157882492630570_r8;
    fjkm(7)=110439.44251210509668_r8;
    fjkm(8)=-46659.137282813036883_r8;
    fjkm(9)=-173333.53977304078369_r8;
    fjkm(10)=339551.86235951279889_r8;
    fjkm(11)=-281181.85781749484440_r8;
    fjkm(12)=116143.52270798282713_r8;
    fjkm(13)=-19586.901426503340492_r8;
    j= 7; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=1.5486069023609161377_r8;
    fjkm(1)=-22.482862472534179688_r8;
    fjkm(2)=-129.63729095714432853_r8;
    fjkm(3)=2741.6385669817243304_r8;
    fjkm(4)=-11510.430102675971531_r8;
    fjkm(5)=3157.6331450774177672_r8;
    fjkm(6)=105738.58606273177860_r8;
    fjkm(7)=-320687.78222286271460_r8;
    fjkm(8)=312110.04133755429291_r8;
    fjkm(9)=306706.99777237116391_r8;
    fjkm(10)=-1199602.2626751819183_r8;
    fjkm(11)=1489876.7581807900958_r8;
    fjkm(12)=-983301.03812460934208_r8;
    fjkm(13)=346445.22525468589947_r8;
    fjkm(14)=-51524.795648340968513_r8;
    j= 8; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-1.8067080527544021606_r8;
    fjkm(1)=30.413155585075869705_r8;
    fjkm(2)=153.62532423010894230_r8;
    fjkm(3)=-4275.6818557748761872_r8;
    fjkm(4)=22280.234407631843178_r8;
    fjkm(5)=-22294.360392132099974_r8;
    fjkm(6)=-188892.07658652729984_r8;
    fjkm(7)=800265.57729686724177_r8;
    fjkm(8)=-1211192.8548070343556_r8;
    fjkm(9)=-77428.408184713489979_r8;
    fjkm(10)=3343683.8379703140094_r8;
    fjkm(11)=-6075462.1554119391419_r8;
    fjkm(12)=5742939.8025630234344_r8;
    fjkm(13)=-3178262.5289756645302_r8;
    fjkm(14)=978732.98065558879521_r8;
    fjkm(15)=-130276.59845140693203_r8;
    j= 9; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=2.0777142606675624847_r8;
    fjkm(1)=-39.947360754013061523_r8;
    fjkm(2)=-172.57638879468381540_r8;
    fjkm(3)=6386.1869555628867376_r8;
    fjkm(4)=-40128.233950133856041_r8;
    fjkm(5)=67957.947814703914434_r8;
    fjkm(6)=297268.90885718166849_r8;
    fjkm(7)=-1779345.5277624794845_r8;
    fjkm(8)=3703482.9515239098067_r8;
    fjkm(9)=-1986101.2546910898185_r8;
    fjkm(10)=-7335848.9571808003709_r8;
    fjkm(11)=20236466.148311260729_r8;
    fjkm(12)=-26009069.048248407006_r8;
    fjkm(13)=20168378.697199375155_r8;
    fjkm(14)=-9655403.7938215681211_r8;
    fjkm(15)=2644939.5099481697170_r8;
    fjkm(16)=-318773.08892039496616_r8;
    j= 10; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-2.3610389325767755508_r8;
    fjkm(1)=51.215318048867833364_r8;
    fjkm(2)=182.72462206625770697_r8;
    fjkm(3)=-9201.6026100350086393_r8;
    fjkm(4)=68268.256270370096570_r8;
    fjkm(5)=-162141.74207057405274_r8;
    fjkm(6)=-402709.29424480088095_r8;
    fjkm(7)=3603004.6646962551842_r8;
    fjkm(8)=-9740822.7436915944519_r8;
    fjkm(9)=10107345.074822039354_r8;
    fjkm(10)=11498843.003383757375_r8;
    fjkm(11)=-56841651.345428293012_r8;
    fjkm(12)=97219891.251414595951_r8;
    fjkm(13)=-99519441.572391483730_r8;
    fjkm(14)=65942266.170395132237_r8;
    fjkm(15)=-27893470.527717198286_r8;
    fjkm(16)=6888375.1431181415309_r8;
    fjkm(17)=-758786.31484749532876_r8;
    j= 11; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=2.6561687991488724947_r8;
    fjkm(1)=-64.344060508074698510_r8;
    fjkm(2)=-179.63738093291836674_r8;
    fjkm(3)=12860.884276726402663_r8;
    fjkm(4)=-110864.10416322604021_r8;
    fjkm(5)=338921.66960863282529_r8;
    fjkm(6)=430704.46372470858983_r8;
    fjkm(7)=-6738355.6173354573243_r8;
    fjkm(8)=22959038.837731227338_r8;
    fjkm(9)=-34599901.818598601926_r8;
    fjkm(10)=-5491093.1482636375735_r8;
    fjkm(11)=136348100.13563343194_r8;
    fjkm(12)=-311391327.48073317359_r8;
    fjkm(13)=406747852.87490380879_r8;
    fjkm(14)=-350465400.65634558429_r8;
    fjkm(15)=203621539.64411279745_r8;
    fjkm(16)=-77285292.825523293342_r8;
    fjkm(17)=17387623.032021543609_r8;
    fjkm(18)=-1764164.5657772609975_r8;
    j= 12; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-2.9626498144352808595_r8;
    fjkm(1)=79.458041370067243966_r8;
    fjkm(2)=158.20533868869907341_r8;
    fjkm(3)=-17512.092404496780068_r8;
    fjkm(4)=173186.23779115767550_r8;
    fjkm(5)=-648825.90238130558651_r8;
    fjkm(6)=-226118.57798798103100_r8;
    fjkm(7)=11744385.639221317992_r8;
    fjkm(8)=-49655262.440949257658_r8;
    fjkm(9)=97992370.806143206674_r8;
    fjkm(10)=-45563229.252833811893_r8;
    fjkm(11)=-276725901.55879139753_r8;
    fjkm(12)=874903955.44068001049_r8;
    fjkm(13)=-1431639430.0141678479_r8;
    fjkm(14)=1544324559.7308983592_r8;
    fjkm(15)=-1156507831.0309397511_r8;
    fjkm(16)=599846625.33396072076_r8;
    fjkm(17)=-206711172.70469868149_r8;
    fjkm(18)=42729154.354849580701_r8;
    fjkm(19)=-4019188.6691200667599_r8;
    j= 13; k=3; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.11215209960937500000_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-2.3640869140625000000_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=8.7891235351562500000_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-11.207002616222993827_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=4.6695844234262474280_r8;
    j= 0; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-0.39253234863281250000_r8;
    fjkm(1)=0.46730041503906250000_r8;
    fjkm(2)=13.002478027343750000_r8;
    fjkm(3)=-14.578535970052083333_r8;
    fjkm(4)=-65.918426513671875000_r8;
    fjkm(5)=71.777842203776041667_r8;
    fjkm(6)=106.46652485411844136_r8;
    fjkm(7)=-113.93785993160043724_r8;
    fjkm(8)=-53.700220869401845422_r8;
    fjkm(9)=56.813277151686010374_r8;
    j= 1; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.88319778442382812500_r8;
    fjkm(1)=-2.2430419921875000000_r8;
    fjkm(2)=-40.888863372802734375_r8;
    fjkm(3)=99.291650390625000000_r8;
    fjkm(4)=222.92270863850911458_r8;
    fjkm(5)=-632.81689453125000000_r8;
    fjkm(6)=-205.55324667471426505_r8;
    fjkm(7)=1232.7702877845293210_r8;
    fjkm(8)=-339.12856875133121946_r8;
    fjkm(9)=-728.45517005449459877_r8;
    fjkm(10)=393.21792165601858550_r8;
    j= 2; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-1.6191959381103515625_r8;
    fjkm(1)=6.5174388885498046875_r8;
    fjkm(2)=97.292139053344726562_r8;
    fjkm(3)=-384.72013047112358941_r8;
    fjkm(4)=-422.46132278442382812_r8;
    fjkm(5)=2925.8162224946198640_r8;
    fjkm(6)=-1437.2810672241964458_r8;
    fjkm(7)=-5929.6163575631600839_r8;
    fjkm(8)=6678.9649706318545243_r8;
    fjkm(9)=1992.7516382310725152_r8;
    fjkm(10)=-5560.5995012212682653_r8;
    fjkm(11)=2034.9551693024277015_r8;
    j= 3; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=2.6311933994293212891_r8;
    fjkm(1)=-14.813423156738281250_r8;
    fjkm(2)=-194.76567316055297852_r8;
    fjkm(3)=1116.2957621256510417_r8;
    fjkm(4)=214.74175742997063531_r8;
    fjkm(5)=-9500.2007904052734375_r8;
    fjkm(6)=12733.852428636433166_r8;
    fjkm(7)=15619.117721871584041_r8;
    fjkm(8)=-43856.442195416477973_r8;
    fjkm(9)=16041.189890575016477_r8;
    fjkm(10)=30538.376827393703173_r8;
    fjkm(11)=-31457.433732481486840_r8;
    fjkm(12)=8757.4502329231489542_r8;
    j= 4; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-3.9467900991439819336_r8;
    fjkm(1)=28.973631262779235840_r8;
    fjkm(2)=345.96240515708923340_r8;
    fjkm(3)=-2698.2640026051657540_r8;
    fjkm(4)=1749.1663194396760729_r8;
    fjkm(5)=24230.291604531833104_r8;
    fjkm(6)=-57186.682706525590685_r8;
    fjkm(7)=-12268.269917924904529_r8;
    fjkm(8)=179642.68522044022878_r8;
    fjkm(9)=-184075.59647969633791_r8;
    fjkm(10)=-55836.464134952713487_r8;
    fjkm(11)=219854.46366368396092_r8;
    fjkm(12)=-146898.32628401899970_r8;
    fjkm(13)=33116.007471226346158_r8;
    j= 5; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=5.5912859737873077393_r8;
    fjkm(1)=-51.149418354034423828_r8;
    fjkm(2)=-562.32073248028755188_r8;
    fjkm(3)=5740.3790382385253906_r8;
    fjkm(4)=-8505.8885195685227712_r8;
    fjkm(5)=-51098.161945523156060_r8;
    fjkm(6)=189688.56368664133696_r8;
    fjkm(7)=-98986.676113505422333_r8;
    fjkm(8)=-524313.33320720157996_r8;
    fjkm(9)=1006412.2572891519230_r8;
    fjkm(10)=-338656.53584266578056_r8;
    fjkm(11)=-879242.30314162838225_r8;
    fjkm(12)=1184856.1356540792633_r8;
    fjkm(13)=-599009.59593194338847_r8;
    fjkm(14)=113723.03789882673771_r8;
    j= 6; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-7.5881738215684890747_r8;
    fjkm(1)=83.791322236259778341_r8;
    fjkm(2)=852.65149199664592743_r8;
    fjkm(3)=-11106.114063047180100_r8;
    fjkm(4)=25906.550896742895797_r8;
    fjkm(5)=90628.038560826472504_r8;
    fjkm(6)=-518627.00526554010007_r8;
    fjkm(7)=625235.13813439073187_r8;
    fjkm(8)=1093177.6471805288836_r8;
    fjkm(9)=-3890056.2699064931220_r8;
    fjkm(10)=3443257.5304835279133_r8;
    fjkm(11)=1688397.8063561002636_r8;
    fjkm(12)=-6072278.7538110165925_r8;
    fjkm(13)=5368522.3053911687123_r8;
    fjkm(14)=-2206353.9977704553128_r8;
    fjkm(15)=362368.26917284610367_r8;
    j= 7; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=9.9594781408086419106_r8;
    fjkm(1)=-129.64064240455627441_r8;
    fjkm(2)=-1221.6473410353064537_r8;
    fjkm(3)=19961.318612462793078_r8;
    fjkm(4)=-64135.377206358956439_r8;
    fjkm(5)=-132491.34865988838862_r8;
    fjkm(6)=1231383.3446542421759_r8;
    fjkm(7)=-2354589.9435372119185_r8;
    fjkm(8)=-1264272.3547066582332_r8;
    fjkm(9)=11829877.236705665039_r8;
    fjkm(10)=-17640228.849921599961_r8;
    fjkm(11)=3679036.5985685346058_r8;
    fjkm(12)=21328290.463957701877_r8;
    fjkm(13)=-31765842.068457922539_r8;
    fjkm(14)=21552604.862053478433_r8;
    fjkm(15)=-7505720.5013225646891_r8;
    fjkm(16)=1087467.9477654193166_r8;
    j= 8; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-12.725999846588820219_r8;
    fjkm(1)=191.72191139924424616_r8;
    fjkm(2)=1668.3300609576205413_r8;
    fjkm(3)=-33822.376037367149478_r8;
    fjkm(4)=139670.53397437636223_r8;
    fjkm(5)=142413.44221576977052_r8;
    fjkm(6)=-2615894.9488370680107_r8;
    fjkm(7)=7028008.7411020993735_r8;
    fjkm(8)=-1692308.6869349767646_r8;
    fjkm(9)=-29470781.812749969812_r8;
    fjkm(10)=66963160.307047821471_r8;
    fjkm(11)=-48089009.180540108686_r8;
    fjkm(12)=-47220345.652008289171_r8;
    fjkm(13)=138176622.72569342840_r8;
    fjkm(14)=-141477318.49446414033_r8;
    fjkm(15)=78991076.066288083382_r8;
    fjkm(16)=-23950277.790642797164_r8;
    fjkm(17)=3106959.7999206284827_r8;
    j= 9; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=15.907499808236025274_r8;
    fjkm(1)=-273.33609867841005325_r8;
    fjkm(2)=-2184.4474298407246048_r8;
    fjkm(3)=54603.014533953543693_r8;
    fjkm(4)=-277734.45331034886641_r8;
    fjkm(5)=-41109.662796960089573_r8;
    fjkm(6)=5064705.8054347585059_r8;
    fjkm(7)=-18090940.629099192262_r8;
    fjkm(8)=15917528.891696545807_r8;
    fjkm(9)=60220437.637860917599_r8;
    fjkm(10)=-208561974.26419501420_r8;
    fjkm(11)=250941409.40257928779_r8;
    fjkm(12)=12084648.536915454989_r8;
    fjkm(13)=-458976381.31741616556_r8;
    fjkm(14)=699975914.26074700776_r8;
    fjkm(15)=-563757375.34166870146_r8;
    fjkm(16)=269426949.85344351478_r8;
    fjkm(17)=-72497863.184361287773_r8;
    fjkm(18)=8519623.3256649401434_r8;
    j= 10; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-19.522840673744212836_r8;
    fjkm(1)=378.05442703043402114_r8;
    fjkm(2)=2752.8286152184838570_r8;
    fjkm(3)=-84659.009316768248795_r8;
    fjkm(4)=515277.57789757005885_r8;
    fjkm(5)=-327436.12858763123432_r8;
    fjkm(6)=-9039645.4021864355266_r8;
    fjkm(7)=41795541.145581633576_r8;
    fjkm(8)=-61523469.344794095368_r8;
    fjkm(9)=-95072522.336003533273_r8;
    fjkm(10)=557312492.76991978150_r8;
    fjkm(11)=-963898631.17215744429_r8;
    fjkm(12)=480963077.05278739701_r8;
    fjkm(13)=1131925258.2353560419_r8;
    fjkm(14)=-2762861092.1224474487_r8;
    fjkm(15)=3066633312.6192085228_r8;
    fjkm(16)=-2065683582.7903266052_r8;
    fjkm(17)=866733914.71761334168_r8;
    fjkm(18)=-209952808.47963646972_r8;
    fjkm(19)=22561861.306890567863_r8;
    j= 11; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=23.590099147440923844_r8;
    fjkm(1)=-509.71270865877158940_r8;
    fjkm(2)=-3345.7051051560481552_r8;
    fjkm(3)=126830.08496875773140_r8;
    fjkm(4)=-904536.17796320887184_r8;
    fjkm(5)=1241459.9239568200231_r8;
    fjkm(6)=14964746.535519588726_r8;
    fjkm(7)=-88697323.877097900818_r8;
    fjkm(8)=182496348.04247934870_r8;
    fjkm(9)=82548357.373675412652_r8;
    fjkm(10)=-1305701152.6740455738_r8;
    fjkm(11)=3075905875.9322221769_r8;
    fjkm(12)=-2915784132.1314453282_r8;
    fjkm(13)=-1631935529.3260648957_r8;
    fjkm(14)=8923557290.4467172403_r8;
    fjkm(15)=-13522339111.332256776_r8;
    fjkm(16)=12138202400.912393639_r8;
    fjkm(17)=-7081644626.2422883477_r8;
    fjkm(18)=2655510812.9195669082_r8;
    fjkm(19)=-585530475.83660861349_r8;
    fjkm(20)=57986597.253985419492_r8;
    j= 12; k=4; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.22710800170898437500_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-7.3687943594796316964_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=42.534998745388454861_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-91.818241543240017361_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=84.636217674600734632_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-28.212072558200244877_r8;
    j= 0; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-1.0219860076904296875_r8;
    fjkm(1)=1.1733913421630859375_r8;
    fjkm(2)=47.897163336617606027_r8;
    fjkm(3)=-52.809692909604027158_r8;
    fjkm(4)=-361.54748933580186632_r8;
    fjkm(5)=389.90415516606083623_r8;
    fjkm(6)=964.09153620402018229_r8;
    fjkm(7)=-1025.3036972328468605_r8;
    fjkm(8)=-1057.9527209325091829_r8;
    fjkm(9)=1114.3768660489096727_r8;
    fjkm(10)=409.07505209390355072_r8;
    fjkm(11)=-427.88310046603704731_r8;
    j= 1; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=2.8104615211486816406_r8;
    fjkm(1)=-6.8132400512695312500_r8;
    fjkm(2)=-175.59265831538609096_r8;
    fjkm(3)=412.65248413085937500_r8;
    fjkm(4)=1483.6983865298922100_r8;
    fjkm(5)=-3828.1498870849609375_r8;
    fjkm(6)=-3429.1824372044316045_r8;
    fjkm(7)=12120.007883707682292_r8;
    fjkm(8)=557.04779563126740632_r8;
    fjkm(9)=-15403.791616777333703_r8;
    fjkm(10)=5099.3321148946942616_r8;
    fjkm(11)=6770.8974139680587706_r8;
    fjkm(12)=-3602.9167662868229396_r8;
    j= 2; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-6.0893332958221435547_r8;
    fjkm(1)=23.218954324722290039_r8;
    fjkm(2)=480.30594648633684431_r8;
    fjkm(3)=-1808.5316684886387416_r8;
    fjkm(4)=-3878.5356589824434311_r8;
    fjkm(5)=19896.567837257637549_r8;
    fjkm(6)=-442.43697992960611979_r8;
    fjkm(7)=-68889.990792852959025_r8;
    fjkm(8)=51994.291933598341765_r8;
    fjkm(9)=82310.686911069807202_r8;
    fjkm(10)=-107791.27622674358562_r8;
    fjkm(11)=-11421.030640241631355_r8;
    fjkm(12)=61775.622629784098705_r8;
    fjkm(13)=-22242.802900370862046_r8;
    j= 3; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=11.417499929666519165_r8;
    fjkm(1)=-60.543208122253417969_r8;
    fjkm(2)=-1092.4312783437115805_r8;
    fjkm(3)=5864.8337360927036830_r8;
    fjkm(4)=6502.6598836863797808_r8;
    fjkm(5)=-73117.673620733634505_r8;
    fjkm(6)=62464.986490075687042_r8;
    fjkm(7)=248344.19895160816334_r8;
    fjkm(8)=-446788.55343178424816_r8;
    fjkm(9)=-141685.28980603760980_r8;
    fjkm(10)=805685.00855625677338_r8;
    fjkm(11)=-411181.55351158250234_r8;
    fjkm(12)=-353321.84981767407842_r8;
    fjkm(13)=410732.51135669781511_r8;
    fjkm(14)=-112357.72180097660343_r8;
    j= 4; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-19.409749880433082581_r8;
    fjkm(1)=133.60695303976535797_r8;
    fjkm(2)=2184.1347855359315872_r8;
    fjkm(3)=-15685.927751365060709_r8;
    fjkm(4)=-3330.0494749048683378_r8;
    fjkm(5)=213065.39140775687165_r8;
    fjkm(6)=-371035.73548135295431_r8;
    fjkm(7)=-595658.10351999312306_r8;
    fjkm(8)=2217706.3121620208025_r8;
    fjkm(9)=-928359.76150830112939_r8;
    fjkm(10)=-3462387.4565158783367_r8;
    fjkm(11)=4492508.5831562094441_r8;
    fjkm(12)=-105953.60990151918538_r8;
    fjkm(13)=-3174045.4228780972346_r8;
    fjkm(14)=2222890.1148558130258_r8;
    fjkm(15)=-492012.66653936007240_r8;
    j= 5; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=30.732103977352380753_r8;
    fjkm(1)=-262.68730902671813965_r8;
    fjkm(2)=-3966.7030987024307251_r8;
    fjkm(3)=36616.605837684018271_r8;
    fjkm(4)=-23288.020921949948583_r8;
    fjkm(5)=-522073.19210540329968_r8;
    fjkm(6)=1445873.6105443313563_r8;
    fjkm(7)=729826.91993359621660_r8;
    fjkm(8)=-8027322.7404775209228_r8;
    fjkm(9)=9022069.9722413070898_r8;
    fjkm(10)=8528377.7669713558429_r8;
    fjkm(11)=-26111911.326974580072_r8;
    fjkm(12)=15072848.600502055062_r8;
    fjkm(13)=11547035.062352154444_r8;
    fjkm(14)=-20141460.694713124158_r8;
    fjkm(15)=10381853.494410238157_r8;
    fjkm(16)=-1934247.3992962518385_r8;
    j= 6; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-46.098155966028571129_r8;
    fjkm(1)=474.28179555572569370_r8;
    fjkm(2)=6683.6986173737261977_r8;
    fjkm(3)=-77185.928108340199369_r8;
    fjkm(4)=113831.12007369411134_r8;
    fjkm(5)=1111273.3131467255535_r8;
    fjkm(6)=-4487654.8822124313622_r8;
    fjkm(7)=1363280.0193113290821_r8;
    fjkm(8)=22934079.022569534587_r8;
    fjkm(9)=-45086888.891430235790_r8;
    fjkm(10)=-1310272.8912292084566_r8;
    fjkm(11)=103829405.25391139096_r8;
    fjkm(12)=-124354846.65650933807_r8;
    fjkm(13)=7129538.3762123397968_r8;
    fjkm(14)=107358912.41929520843_r8;
    fjkm(15)=-104902766.60439072992_r8;
    fjkm(16)=43358616.238781106380_r8;
    fjkm(17)=-6986431.7916780392488_r8;
    j= 7; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=66.266099201166070998_r8;
    fjkm(1)=-801.85159036517143250_r8;
    fjkm(2)=-10599.673972529735017_r8;
    fjkm(3)=150238.26124378282197_r8;
    fjkm(4)=-349014.85897753891605_r8;
    fjkm(5)=-2089251.7495501184712_r8;
    fjkm(6)=11932847.810754978267_r8;
    fjkm(7)=-12233248.989355019522_r8;
    fjkm(8)=-52996346.810335384350_r8;
    fjkm(9)=167552806.49381000405_r8;
    fjkm(10)=-104151238.67453537869_r8;
    fjkm(11)=-295472139.38679802840_r8;
    fjkm(12)=638921130.05027917750_r8;
    fjkm(13)=-364575119.55069248400_r8;
    fjkm(14)=-321938848.50186568760_r8;
    fjkm(15)=670099675.87405621186_r8;
    fjkm(16)=-477068925.07477830810_r8;
    fjkm(17)=165792634.22539301473_r8;
    fjkm(18)=-23563863.859185525714_r8;
    j= 8; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-92.036248890508431941_r8;
    fjkm(1)=1286.5459820115674202_r8;
    fjkm(2)=15984.246642014751810_r8;
    fjkm(3)=-274251.94134559012498_r8;
    fjkm(4)=875242.48789234256927_r8;
    fjkm(5)=3479369.5294348938478_r8;
    fjkm(6)=-28243123.382389417092_r8;
    fjkm(7)=48978753.833933594038_r8;
    fjkm(8)=96403146.255130900766_r8;
    fjkm(9)=-511553591.82361285211_r8;
    fjkm(10)=629980523.20384634815_r8;
    fjkm(11)=530948403.41019250637_r8;
    fjkm(12)=-2455930387.0192782105_r8;
    fjkm(13)=2650615136.5125114389_r8;
    fjkm(14)=-13787083.107153494438_r8;
    fjkm(15)=-2934135354.5042859293_r8;
    fjkm(16)=3427343175.1586684272_r8;
    fjkm(17)=-1959752975.9949664819_r8;
    fjkm(18)=590135160.40330437780_r8;
    fjkm(19)=-75099321.778257988527_r8;
    j= 9; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=124.24893600218638312_r8;
    fjkm(1)=-1977.9098049255553633_r8;
    fjkm(2)=-23091.371137548782696_r8;
    fjkm(3)=474842.97883978903678_r8;
    fjkm(4)=-1940160.8026042800714_r8;
    fjkm(5)=-5051698.1764753913085_r8;
    fjkm(6)=60910728.961121054906_r8;
    fjkm(7)=-150451158.23316204445_r8;
    fjkm(8)=-115862597.99920025974_r8;
    fjkm(9)=1340568362.9608932867_r8;
    fjkm(10)=-2534626023.8559564861_r8;
    fjkm(11)=57223508.996001180928_r8;
    fjkm(12)=7462542511.1368897210_r8;
    fjkm(13)=-12803010436.906857334_r8;
    fjkm(14)=6529418551.9917203148_r8;
    fjkm(15)=8333347633.6309463361_r8;
    fjkm(16)=-17879701023.370781023_r8;
    fjkm(17)=15384544080.383156452_r8;
    fjkm(18)=-7429525037.6309918827_r8;
    fjkm(19)=1979364564.1715841600_r8;
    fjkm(20)=-228201703.20311712289_r8;
    j= 10; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-163.78268836651841411_r8;
    fjkm(1)=2934.5753597465244149_r8;
    fjkm(2)=32133.674298098814287_r8;
    fjkm(3)=-786448.50173515116761_r8;
    fjkm(4)=3940574.6676366506081_r8;
    fjkm(5)=6023482.0215492578926_r8;
    fjkm(6)=-121576300.80985078356_r8;
    fjkm(7)=396478314.88388992406_r8;
    fjkm(8)=-28867608.104982563453_r8;
    fjkm(9)=-3079357130.2226958330_r8;
    fjkm(10)=8249384124.8426910359_r8;
    fjkm(11)=-5275292901.0490022350_r8;
    fjkm(12)=-17864885603.458945318_r8;
    fjkm(13)=48480263067.875486448_r8;
    fjkm(14)=-46292296797.993831733_r8;
    fjkm(15)=-6808810264.3074322272_r8;
    fjkm(16)=70205428541.430035549_r8;
    fjkm(17)=-89866265315.798467190_r8;
    fjkm(18)=62728346454.981427111_r8;
    fjkm(19)=-26379975109.497266807_r8;
    fjkm(20)=6313975478.5070403854_r8;
    fjkm(21)=-665761463.93251599995_r8;
    j= 11; k=5; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=0.57250142097473144531_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-26.491430486951555525_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=218.19051174421159048_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-699.57962737613254123_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=1059.9904525279998779_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-765.25246814118164230_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=212.57013003921712286_r8;
    j= 0; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-3.1487578153610229492_r8;
    fjkm(1)=3.5304254293441772461_r8;
    fjkm(2)=198.68572865213666643_r8;
    fjkm(3)=-216.34668231010437012_r8;
    fjkm(4)=-2072.8098615700101096_r8;
    fjkm(5)=2218.2702027328178365_r8;
    fjkm(6)=8045.1657148255242242_r8;
    fjkm(7)=-8511.5521330762792517_r8;
    fjkm(8)=-14309.871109127998352_r8;
    fjkm(9)=15016.531410813331604_r8;
    fjkm(10)=11861.413256188315456_r8;
    fjkm(11)=-12371.581568282436550_r8;
    fjkm(12)=-3719.9772756862996501_r8;
    fjkm(13)=3861.6906957124443986_r8;
    j= 1; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=10.233462899923324585_r8;
    fjkm(1)=-24.045059680938720703_r8;
    fjkm(2)=-830.55504153881754194_r8;
    fjkm(3)=1907.3829950605119978_r8;
    fjkm(4)=9817.0755057463759468_r8;
    fjkm(5)=-24000.956291863274953_r8;
    fjkm(6)=-37145.398656393453558_r8;
    fjkm(7)=109134.42187067667643_r8;
    fjkm(8)=44836.131085879493643_r8;
    fjkm(9)=-222597.99503087997437_r8;
    fjkm(10)=21083.102663859050460_r8;
    fjkm(11)=208148.67133440140671_r8;
    fjkm(12)=-75945.993209761297570_r8;
    fjkm(13)=-72698.984473412256018_r8;
    fjkm(14)=38306.908850817252349_r8;
    j= 2; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-25.583657249808311462_r8;
    fjkm(1)=94.002347901463508606_r8;
    fjkm(2)=2561.4464542163269860_r8;
    fjkm(3)=-9323.5958246609994343_r8;
    fjkm(4)=-30925.007683879997995_r8;
    fjkm(5)=137806.75057795568307_r8;
    fjkm(6)=66832.114908046586804_r8;
    fjkm(7)=-695211.42408942898710_r8;
    fjkm(8)=297044.24306208789349_r8;
    fjkm(9)=1456689.0310313083630_r8;
    fjkm(10)=-1349408.4412036920976_r8;
    fjkm(11)=-1160751.9107670259288_r8;
    fjkm(12)=1778986.1326650806502_r8;
    fjkm(13)=2732.4118798791034334_r8;
    fjkm(14)=-771855.42780552482418_r8;
    fjkm(15)=274755.25810395356345_r8;
    j= 3; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=54.365271655842661858_r8;
    fjkm(1)=-276.51818633079528809_r8;
    fjkm(2)=-6505.3076391667127609_r8;
    fjkm(3)=33393.909001989024026_r8;
    fjkm(4)=70377.367684318314469_r8;
    fjkm(5)=-559956.42247611134141_r8;
    fjkm(6)=214189.08434628603635_r8;
    fjkm(7)=2932546.1434609688359_r8;
    fjkm(8)=-3873550.9334950489425_r8;
    fjkm(9)=-5169809.5626455059758_r8;
    fjkm(10)=12515387.720161636405_r8;
    fjkm(11)=-485287.10103841189331_r8;
    fjkm(12)=-14696506.049911874334_r8;
    fjkm(13)=8973612.6112505443054_r8;
    fjkm(14)=4358025.7113579717181_r8;
    fjkm(15)=-5899263.9630258568617_r8;
    fjkm(16)=1593568.9458830170786_r8;
    j= 4; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-103.29401614610105753_r8;
    fjkm(1)=679.49573871586471796_r8;
    fjkm(2)=14406.592158034443855_r8;
    fjkm(3)=-97833.485413427407644_r8;
    fjkm(4)=-109874.73447139263153_r8;
    fjkm(5)=1806898.4781330669971_r8;
    fjkm(6)=-2228617.5688649368428_r8;
    fjkm(7)=-9039218.0698707036945_r8;
    fjkm(8)=22913970.357558080752_r8;
    fjkm(9)=6014580.7747564135778_r8;
    fjkm(10)=-67365551.082731008652_r8;
    fjkm(11)=49347945.008100392566_r8;
    fjkm(12)=61291772.218955972273_r8;
    fjkm(13)=-101851641.61990357673_r8;
    fjkm(14)=18812228.438435360796_r8;
    fjkm(15)=48878028.306514326695_r8;
    fjkm(16)=-36319210.680616361669_r8;
    fjkm(17)=7931540.8655372143613_r8;
    j= 5; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=180.76452825567685068_r8;
    fjkm(1)=-1472.1516226977109909_r8;
    fjkm(2)=-28789.438034501904622_r8;
    fjkm(3)=248412.49281514968191_r8;
    fjkm(4)=57720.374439080618322_r8;
    fjkm(5)=-4919773.7285477433167_r8;
    fjkm(6)=10556839.574755387407_r8;
    fjkm(7)=20546468.524262875642_r8;
    fjkm(8)=-95782021.413901062575_r8;
    fjkm(9)=47859753.794019524423_r8;
    fjkm(10)=248947938.76422519634_r8;
    fjkm(11)=-397302112.77505450021_r8;
    fjkm(12)=-60373613.593196943601_r8;
    fjkm(13)=619895830.67946774788_r8;
    fjkm(14)=-460542977.57691540068_r8;
    fjkm(15)=-133881283.62288371220_r8;
    fjkm(16)=360816116.57296145253_r8;
    fjkm(17)=-191228273.50596204944_r8;
    fjkm(18)=35131056.264643928958_r8;
    j= 6; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-296.97029642004054040_r8;
    fjkm(1)=2903.8678610353963450_r8;
    fjkm(2)=53085.574017544759304_r8;
    fjkm(3)=-566162.21564224131681_r8;
    fjkm(4)=351303.05079310148650_r8;
    fjkm(5)=11717363.296337798053_r8;
    fjkm(6)=-37248885.401731669600_r8;
    fjkm(7)=-29556393.543845627395_r8;
    fjkm(8)=317963019.15514055453_r8;
    fjkm(9)=-391675101.18705789558_r8;
    fjkm(10)=-632571201.89067213266_r8;
    fjkm(11)=2001265322.6458411313_r8;
    fjkm(12)=-1008064787.2696644192_r8;
    fjkm(13)=-2363398838.2753953344_r8;
    fjkm(14)=3743079200.5912006268_r8;
    fjkm(15)=-1091651847.4629542050_r8;
    fjkm(16)=-1902208028.0358040356_r8;
    fjkm(17)=2133928476.4409109596_r8;
    fjkm(18)=-893289789.98112876745_r8;
    fjkm(19)=141870657.61208999896_r8;
    j= 7; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=464.01608815631334437_r8;
    fjkm(1)=-5325.3068480961956084_r8;
    fjkm(2)=-91725.505981153156193_r8;
    fjkm(3)=1185316.6466349283157_r8;
    fjkm(4)=-1734253.4162956622921_r8;
    fjkm(5)=-24965643.657687719930_r8;
    fjkm(6)=109865301.06964371257_r8;
    fjkm(7)=-7842536.5990090553082_r8;
    fjkm(8)=-882004248.16533042144_r8;
    fjkm(9)=1812569161.1229405097_r8;
    fjkm(10)=796190115.17482420785_r8;
    fjkm(11)=-7547640676.2891254543_r8;
    fjkm(12)=8606208381.3162846761_r8;
    fjkm(13)=4634326771.3840251843_r8;
    fjkm(14)=-19504767652.151161478_r8;
    fjkm(15)=15458811432.485826518_r8;
    fjkm(16)=3233232293.8508888375_r8;
    fjkm(17)=-14292761227.388542723_r8;
    fjkm(18)=10872211035.524945516_r8;
    fjkm(19)=-3794154076.5815443275_r8;
    fjkm(20)=531367092.46942384383_r8;
    j= 8; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-696.02413223447001656_r8;
    fjkm(1)=9211.7329484450783639_r8;
    fjkm(2)=150176.80785295595602_r8;
    fjkm(3)=-2316795.9016969799608_r8;
    fjkm(4)=5353301.7267099139240_r8;
    fjkm(5)=48240128.320645392607_r8;
    fjkm(6)=-285088568.66606943443_r8;
    fjkm(7)=236539658.57255855849_r8;
    fjkm(8)=2091063546.0424140229_r8;
    fjkm(9)=-6478408107.9479218188_r8;
    fjkm(10)=2114947617.4171696310_r8;
    fjkm(11)=22484222108.370436724_r8;
    fjkm(12)=-43504182331.312579142_r8;
    fjkm(13)=9073306866.7430378915_r8;
    fjkm(14)=72378013713.663446159_r8;
    fjkm(15)=-104755002543.75143509_r8;
    fjkm(16)=35288894490.708636111_r8;
    fjkm(17)=58426452630.062379587_r8;
    fjkm(18)=-83650899080.625595930_r8;
    fjkm(19)=49569134901.994243944_r8;
    fjkm(20)=-14909719423.420583328_r8;
    fjkm(21)=1869289195.4875346262_r8;
    j= 9; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=1009.2349917399815240_r8;
    fjkm(1)=-15188.489623096204014_r8;
    fjkm(2)=-234912.74468029359442_r8;
    fjkm(3)=4278033.8338821994157_r8;
    fjkm(4)=-13570280.995019535225_r8;
    fjkm(5)=-85095459.392767211454_r8;
    fjkm(6)=669951675.09023006543_r8;
    fjkm(7)=-1041016026.5703218480_r8;
    fjkm(8)=-4241350570.5261899700_r8;
    fjkm(9)=19536358548.670235798_r8;
    fjkm(10)=-19158789931.456781880_r8;
    fjkm(11)=-52588961930.566950783_r8;
    fjkm(12)=168223251386.27505861_r8;
    fjkm(13)=-131543443714.17931735_r8;
    fjkm(14)=-183278984759.32788630_r8;
    fjkm(15)=500731039137.68584765_r8;
    fjkm(16)=-395175402811.80868032_r8;
    fjkm(17)=-83866621459.668331313_r8;
    fjkm(18)=443918936157.05017610_r8;
    fjkm(19)=-420471377306.84165961_r8;
    fjkm(20)=207052924562.97132855_r8;
    fjkm(21)=-54907932888.507130529_r8;
    fjkm(22)=6236056730.2635893277_r8;
    j= 10; k=6; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=1.7277275025844573975_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-108.09091978839465550_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=1200.9029132163524628_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-5305.6469786134031084_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=11655.393336864533248_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-13586.550006434137439_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=8061.7221817373093845_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-1919.4576623184069963_r8;
    j= 0; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-11.230228766798973083_r8;
    fjkm(1)=12.382047101855278015_r8;
    fjkm(2)=918.77281820135457175_r8;
    fjkm(3)=-990.83343139361767542_r8;
    fjkm(4)=-12609.480588771700859_r8;
    fjkm(5)=13410.082530915935834_r8;
    fjkm(6)=66320.587232667538855_r8;
    fjkm(7)=-69857.685218409807594_r8;
    fjkm(8)=-169003.20338453573209_r8;
    fjkm(9)=176773.46560911208759_r8;
    fjkm(10)=224178.07510616326774_r8;
    fjkm(11)=-233235.77511045269270_r8;
    fjkm(12)=-149141.86036214022361_r8;
    fjkm(13)=154516.34181663176320_r8;
    fjkm(14)=39348.882077527343424_r8;
    fjkm(15)=-40628.520519072948089_r8;
    j= 1; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=42.113357875496149063_r8;
    fjkm(1)=-96.752740144729614258_r8;
    fjkm(2)=-4309.3875268953187125_r8;
    fjkm(3)=9728.1827809555189950_r8;
    fjkm(4)=67131.493914289162272_r8;
    fjkm(5)=-158519.18454455852509_r8;
    fjkm(6)=-361549.21741861661275_r8;
    fjkm(7)=965627.75010763936573_r8;
    fjkm(8)=791368.90269480066167_r8;
    fjkm(9)=-2797294.4008474879795_r8;
    fjkm(10)=-473067.29978352049251_r8;
    fjkm(11)=4157484.3019688460562_r8;
    fjkm(12)=-742925.21875958646140_r8;
    fjkm(13)=-3063454.4290601775661_r8;
    fjkm(14)=1186992.6183777028865_r8;
    fjkm(15)=886789.43999110403230_r8;
    fjkm(16)=-463948.91246287829107_r8;
    j= 2; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-119.32118064723908901_r8;
    fjkm(1)=426.72709654457867146_r8;
    fjkm(2)=14774.672950860112906_r8;
    fjkm(3)=-52455.440737222809167_r8;
    fjkm(4)=-242280.57068559899926_r8;
    fjkm(5)=994102.67395229967167_r8;
    fjkm(6)=1032233.4449197463691_r8;
    fjkm(7)=-6741567.1085603777512_r8;
    fjkm(8)=646495.83215296654790_r8;
    fjkm(9)=20680668.518424851831_r8;
    fjkm(10)=-13425393.794712164658_r8;
    fjkm(11)=-29950004.715897617246_r8;
    fjkm(12)=32133326.488920312047_r8;
    fjkm(13)=16877471.905929211282_r8;
    fjkm(14)=-30879035.210339582752_r8;
    fjkm(15)=1961435.1350279426026_r8;
    fjkm(16)=10741165.112229910651_r8;
    fjkm(17)=-3791244.3495002070744_r8;
    j= 3; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=283.38780403719283640_r8;
    fjkm(1)=-1397.7315495908260345_r8;
    fjkm(2)=-41380.460209263255820_r8;
    fjkm(3)=205570.88069305533455_r8;
    fjkm(4)=656958.69278146053058_r8;
    fjkm(5)=-4406726.3053560412498_r8;
    fjkm(6)=-460386.15389300137896_r8;
    fjkm(7)=31745963.553354091997_r8;
    fjkm(8)=-30185144.899501059008_r8;
    fjkm(9)=-92885893.761842946947_r8;
    fjkm(10)=159671874.67639934532_r8;
    fjkm(11)=89037317.391947147287_r8;
    fjkm(12)=-324112057.62064508289_r8;
    fjkm(13)=75417289.288447113978_r8;
    fjkm(14)=275169534.60077554540_r8;
    fjkm(15)=-191141129.57979682342_r8;
    fjkm(16)=-56632059.761422146328_r8;
    fjkm(17)=92789782.492575658213_r8;
    fjkm(18)=-24828398.690560814574_r8;
    j= 4; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-595.11438847810495645_r8;
    fjkm(1)=3784.5764889410929754_r8;
    fjkm(2)=100370.50080364884343_r8;
    fjkm(3)=-654419.09384311857985_r8;
    fjkm(4)=-1371381.6192330607878_r8;
    fjkm(5)=15498343.812751787009_r8;
    fjkm(6)=-11665144.592299357907_r8;
    fjkm(7)=-112601561.49426607138_r8;
    fjkm(8)=219262318.25667364074_r8;
    fjkm(9)=254104121.37964987144_r8;
    fjkm(10)=-1006396481.7103226659_r8;
    fjkm(11)=253036623.98729691889_r8;
    fjkm(12)=1831393810.4978639927_r8;
    fjkm(13)=-1774529065.4266491466_r8;
    fjkm(14)=-1003472637.9292914864_r8;
    fjkm(15)=2290280875.9786741918_r8;
    fjkm(16)=-651246794.10887221841_r8;
    fjkm(17)=-803784902.98109705890_r8;
    fjkm(18)=640483342.29369123263_r8;
    fjkm(19)=-138440607.21363135318_r8;
    j= 5; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=1140.6359112497011665_r8;
    fjkm(1)=-8957.2682584379799664_r8;
    fjkm(2)=-218407.33629125532461_r8;
    fjkm(3)=1794837.6930232931017_r8;
    fjkm(4)=2046161.0563529025012_r8;
    fjkm(5)=-45986041.828079905965_r8;
    fjkm(6)=74544641.814482732603_r8;
    fjkm(7)=314850611.45725349592_r8;
    fjkm(8)=-1044136171.4521382041_r8;
    fjkm(9)=-187106978.33355739111_r8;
    fjkm(10)=4438953626.0956817754_r8;
    fjkm(11)=-4408582954.0483059788_r8;
    fjkm(12)=-6264272163.4397192983_r8;
    fjkm(13)=14149702487.577580987_r8;
    fjkm(14)=-2713149740.6519136094_r8;
    fjkm(15)=-14233558941.890554329_r8;
    fjkm(16)=12753481214.719287040_r8;
    fjkm(17)=934842459.40433410209_r8;
    fjkm(18)=-6845048171.9341116529_r8;
    fjkm(19)=3754053882.0127951637_r8;
    fjkm(20)=-682202534.28377278514_r8;
    j= 6; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-2036.8498415173235117_r8;
    fjkm(1)=19163.284363303755526_r8;
    fjkm(2)=436370.36938456366105_r8;
    fjkm(3)=-4395730.5804739226086_r8;
    fjkm(4)=-1112960.9840974107984_r8;
    fjkm(5)=119491789.40923461317_r8;
    fjkm(6)=-304912988.84344882540_r8;
    fjkm(7)=-691963837.01336538937_r8;
    fjkm(8)=3885386169.7360802683_r8;
    fjkm(9)=-2313847981.3260245039_r8;
    fjkm(10)=-14816925759.772210986_r8;
    fjkm(11)=28337891715.905708458_r8;
    fjkm(12)=7872387353.1924133326_r8;
    fjkm(13)=-72435569735.091307437_r8;
    fjkm(14)=59410896148.189366465_r8;
    fjkm(15)=46141874867.831978016_r8;
    fjkm(16)=-105411711029.14681739_r8;
    fjkm(17)=45764643115.283153298_r8;
    fjkm(18)=33619493138.706044340_r8;
    fjkm(19)=-45524891856.627263052_r8;
    fjkm(20)=19398990085.810093364_r8;
    fjkm(21)=-3046176001.4829695650_r8;
    j= 7; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=3437.1841075604834259_r8;
    fjkm(1)=-37884.316448339552153_r8;
    fjkm(2)=-813572.82790964112831_r8;
    fjkm(3)=9844700.1196742646463_r8;
    fjkm(4)=-5932665.5737596173789_r8;
    fjkm(5)=-278602941.80981757775_r8;
    fjkm(6)=999702331.75599910068_r8;
    fjkm(7)=1082145758.5676864833_r8;
    fjkm(8)=-12097872233.934345656_r8;
    fjkm(9)=16269568192.896168430_r8;
    fjkm(10)=37442906272.629884505_r8;
    fjkm(11)=-127678482632.39934914_r8;
    fjkm(12)=56533169829.379619075_r8;
    fjkm(13)=263705962704.22363017_r8;
    fjkm(14)=-430880513772.40015440_r8;
    fjkm(15)=28685028305.645455607_r8;
    fjkm(16)=544657720548.16430202_r8;
    fjkm(17)=-543385080747.85240194_r8;
    fjkm(18)=33319140131.863860742_r8;
    fjkm(19)=311025095335.01861076_r8;
    fjkm(20)=-257492440682.78870598_r8;
    fjkm(21)=90635479554.844098597_r8;
    fjkm(22)=-12545989968.390205031_r8;
    j= 8; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-5537.6855066252232973_r8;
    fjkm(1)=70275.288364441469184_r8;
    fjkm(2)=1432194.8403133915934_r8;
    fjkm(3)=-20502591.084346145880_r8;
    fjkm(4)=29691158.739889488095_r8;
    fjkm(5)=592508323.24904278254_r8;
    fjkm(6)=-2831998971.7799303680_r8;
    fjkm(7)=-487772628.07968000257_r8;
    fjkm(8)=32638786024.059947790_r8;
    fjkm(9)=-70496897875.469508762_r8;
    fjkm(10)=-62531875035.153643030_r8;
    fjkm(11)=456706039700.77564846_r8;
    fjkm(12)=-503779673552.18388727_r8;
    fjkm(13)=-650638245764.84731771_r8;
    fjkm(14)=2119375307389.9958522_r8;
    fjkm(15)=-1373635599107.5234068_r8;
    fjkm(16)=-1758554601545.7817139_r8;
    fjkm(17)=3640007756944.3988399_r8;
    fjkm(18)=-1912987782878.0613666_r8;
    fjkm(19)=-1044942776057.9478291_r8;
    fjkm(20)=2082243082925.6114167_r8;
    fjkm(21)=-1292209352199.4454475_r8;
    fjkm(22)=389815335189.77376153_r8;
    fjkm(23)=-48292926381.689492854_r8;
    j= 9; k=7; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=6.0740420012734830379_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-493.91530477308801242_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=7109.5143024893637214_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-41192.654968897551298_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=122200.46498301745979_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-203400.17728041553428_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=192547.00123253153236_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-96980.598388637513489_r8;
    fjkm(15)=0.0_r8;
    fjkm(16)=20204.291330966148643_r8;
    j= 0; k=8; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-45.555315009551122785_r8;
    fjkm(1)=49.604676343733444810_r8;
    fjkm(2)=4692.1953953443361180_r8;
    fjkm(3)=-5021.4722651930614596_r8;
    fjkm(4)=-81759.414478627682797_r8;
    fjkm(5)=86499.090680287258611_r8;
    fjkm(6)=556100.84208011694252_r8;
    fjkm(7)=-583562.61205938197672_r8;
    fjkm(8)=-1894107.2072367706267_r8;
    fjkm(9)=1975574.1838921155999_r8;
    fjkm(10)=3559503.1024072718499_r8;
    fjkm(11)=-3695103.2205942155394_r8;
    fjkm(12)=-3754666.5240343648810_r8;
    fjkm(13)=3883031.1915227192359_r8;
    fjkm(14)=2085082.8653557065400_r8;
    fjkm(15)=-2149736.5976147982157_r8;
    fjkm(16)=-474800.84627770449312_r8;
    fjkm(17)=488270.37383168192555_r8;
    j= 1; k=8; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=193.61008879059227183_r8;
    fjkm(1)=-437.33102409169077873_r8;
    fjkm(2)=-24389.798720089893322_r8;
    fjkm(3)=54330.683525039681367_r8;
    fjkm(4)=481258.52318321001006_r8;
    fjkm(5)=-1109084.2311883407405_r8;
    fjkm(6)=-3433050.7548587226633_r8;
    fjkm(7)=8650457.5434684857726_r8;
    fjkm(8)=11004225.300068311602_r8;
    fjkm(9)=-33238526.475380749062_r8;
    fjkm(10)=-15303078.309507955098_r8;
    fjkm(11)=69562860.629902112723_r8;
    fjkm(12)=1830924.9239440239572_r8;
    fjkm(13)=-80869740.517663243591_r8;
    fjkm(14)=18943271.994495349280_r8;
    fjkm(15)=49072182.784650581825_r8;
    fjkm(16)=-19806771.899029389669_r8;
    fjkm(17)=-12122574.798579689186_r8;
    fjkm(18)=6307948.1226220563244_r8;
    j= 2; k=8; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-613.09861450354219414_r8;
    fjkm(1)=2147.8571771753195208_r8;
    fjkm(2)=91956.399098661058815_r8;
    fjkm(3)=-320284.80771632295052_r8;
    fjkm(4)=-1938565.2131506522862_r8;
    fjkm(5)=7533108.6348155997448_r8;
    fjkm(6)=12364512.209251265036_r8;
    fjkm(7)=-65358277.938386196419_r8;
    fjkm(8)=-16530840.331176116137_r8;
    fjkm(9)=269292234.00676317160_r8;
    fjkm(10)=-105780946.98529899276_r8;
    fjkm(11)=-575008738.51905591292_r8;
    fjkm(12)=462747211.13824444405_r8;
    fjkm(13)=619754381.87898314676_r8;
    fjkm(14)=-755460241.94008607635_r8;
    fjkm(15)=-252322946.06096464110_r8;
    fjkm(16)=569416784.91969405599_r8;
    fjkm(17)=-61357261.820865861243_r8;
    fjkm(18)=-164974352.55837953060_r8;
    fjkm(19)=57850732.229668052938_r8;
    j= 3; k=8; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=1609.3838630717982596_r8;
    fjkm(1)=-7751.9961041252827272_r8;
    fjkm(2)=-281313.57426896920515_r8;
    fjkm(3)=1362914.7891508336068_r8;
    fjkm(4)=5966457.1081110733990_r8;
    fjkm(5)=-36097318.956064416329_r8;
    fjkm(6)=-22309783.293551422257_r8;
    fjkm(7)=336158425.14225148967_r8;
    fjkm(8)=-205028522.03506721003_r8;
    fjkm(9)=-1388339192.4100137759_r8;
    fjkm(10)=1836112835.6997822732_r8;
    fjkm(11)=2556726458.1368042032_r8;
    fjkm(12)=-5652778880.4580993178_r8;
    fjkm(13)=-1072688790.0156425156_r8;
    fjkm(14)=8223828086.6764744334_r8;
    fjkm(15)=-2962125430.6175380614_r8;
    fjkm(16)=-5363095697.4586512056_r8;
    fjkm(17)=4151384158.2283357977_r8;
    fjkm(18)=758617226.29066830002_r8;
    fjkm(19)=-1589602926.9007581937_r8;
    fjkm(20)=422197436.26031767718_r8;
    j= 4; k=8; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-3701.5828850651359971_r8;
    fjkm(1)=22929.423816498228916_r8;
    fjkm(2)=740942.45462626213339_r8;
    fjkm(3)=-4683368.4794967535629_r8;
    fjkm(4)=-14688082.295896812171_r8;
    fjkm(5)=136985890.98054216279_r8;
    fjkm(6)=-37063877.251336542111_r8;
    fjkm(7)=-1319571104.7402945920_r8;
    fjkm(8)=2005453574.1427636671_r8;
    fjkm(9)=4917316117.1324589056_r8;
    fjkm(10)=-13428166320.891447862_r8;
    fjkm(11)=-3828355991.3598264145_r8;
    fjkm(12)=37778220592.815991623_r8;
    fjkm(13)=-20635506272.066951931_r8;
    fjkm(14)=-47082248514.737817855_r8;
    fjkm(15)=56649997707.122923550_r8;
    fjkm(16)=14169411429.388481941_r8;
    fjkm(17)=-52510965464.207838905_r8;
    fjkm(18)=18673363121.781645026_r8;
    fjkm(19)=14082226269.939926879_r8;
    fjkm(20)=-12159500780.523353877_r8;
    fjkm(21)=2607014902.9539700770_r8;
    j= 5; k=8; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=7711.6310105523666607_r8;
    fjkm(1)=-58858.321154496479721_r8;
    fjkm(2)=-1741907.2159207465870_r8;
    fjkm(3)=13794087.929804007718_r8;
    fjkm(4)=28916480.065464689455_r8;
    fjkm(5)=-437983696.24595980730_r8;
    fjkm(6)=494945307.55363422412_r8;
    fjkm(7)=4180226129.2961517207_r8;
    fjkm(8)=-10954467146.581302332_r8;
    fjkm(9)=-11060935369.012122216_r8;
    fjkm(10)=67274495072.654389414_r8;
    fjkm(11)=-33385597946.547532147_r8;
    fjkm(12)=-168738543918.69430001_r8;
    fjkm(13)=236353638119.12099768_r8;
    fjkm(14)=123426482948.35286349_r8;
    fjkm(15)=-460569493399.15116382_r8;
    fjkm(16)=183685386812.60549147_r8;
    fjkm(17)=323426896220.20397235_r8;
    fjkm(18)=-345059927788.03296106_r8;
    fjkm(19)=18066712637.598243570_r8;
    fjkm(20)=137638821647.93765226_r8;
    fjkm(21)=-78528723144.419087956_r8;
    fjkm(22)=14147149999.271829167_r8;
    j= 6; k=8; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-14872.431234636707131_r8;
    fjkm(1)=135739.23591067179473_r8;
    fjkm(2)=3743563.4066693378186_r8;
    fjkm(3)=-36117123.586437547311_r8;
    fjkm(4)=-41857416.064054280688_r8;
    fjkm(5)=1225475415.8400151136_r8;
    fjkm(6)=-2467842162.9074118554_r8;
    fjkm(7)=-10944143126.183310329_r8;
    fjkm(8)=45219896350.687109257_r8;
    fjkm(9)=3706147531.3941538260_r8;
    fjkm(10)=-259638349679.49799671_r8;
    fjkm(11)=331511238218.36097150_r8;
    fjkm(12)=503674442974.85219625_r8;
    fjkm(13)=-1481164218255.8010841_r8;
    fjkm(14)=369277357070.34485189_r8;
    fjkm(15)=2339646237464.7011180_r8;
    fjkm(16)=-2569885298111.3265694_r8;
    fjkm(17)=-667794049596.58740652_r8;
    fjkm(18)=2906153064933.9250617_r8;
    fjkm(19)=-1575977334606.5289263_r8;
    fjkm(20)=-578377418663.53101806_r8;
    fjkm(21)=1021516684285.5146634_r8;
    fjkm(22)=-444821917816.52114438_r8;
    fjkm(23)=69214137882.703873029_r8;
    j= 7; k=8; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=26956.281612779031676_r8;
    fjkm(1)=-287752.74431752909550_r8;
    fjkm(2)=-7479028.3293311244821_r8;
    fjkm(3)=86134990.009811768863_r8;
    fjkm(4)=23679004.644932133834_r8;
    fjkm(5)=-3077161905.3885988089_r8;
    fjkm(6)=9103660424.3905427045_r8;
    fjkm(7)=23518972784.280263949_r8;
    fjkm(8)=-154703115855.86457830_r8;
    fjkm(9)=108277078595.12757644_r8;
    fjkm(10)=805587125662.51652400_r8;
    fjkm(11)=-1805165807562.2091925_r8;
    fjkm(12)=-672890155245.49490410_r8;
    fjkm(13)=6669143665397.1378050_r8;
    fjkm(14)=-6132703262892.6427237_r8;
    fjkm(15)=-7384886877294.0818151_r8;
    fjkm(16)=17989877009001.983669_r8;
    fjkm(17)=-6680625953666.2171081_r8;
    fjkm(18)=-14315112877022.322823_r8;
    fjkm(19)=17853687377733.207721_r8;
    fjkm(20)=-3800827233430.3173356_r8;
    fjkm(21)=-6918410846679.1440993_r8;
    fjkm(22)=6365772360949.3812603_r8;
    fjkm(23)=-2267586042740.4274751_r8;
    fjkm(24)=310920009576.22258316_r8;
    j= 8; k=8; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=24.380529699556063861_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-2499.8304818112096241_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=45218.768981362726273_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-331645.17248456357783_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=1268365.2733216247816_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-2813563.2265865341107_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=3763271.2976564039964_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-2998015.9185381067501_r8;
    fjkm(15)=0.0_r8;
    fjkm(16)=1311763.6146629772007_r8;
    fjkm(17)=0.0_r8;
    fjkm(18)=-242919.18790055133346_r8;
    j= 0; k=9; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-207.23450244622654282_r8;
    fjkm(1)=223.48818891259725206_r8;
    fjkm(2)=26248.220059017701053_r8;
    fjkm(3)=-27914.773713558507469_r8;
    fjkm(4)=-565234.61226703407842_r8;
    fjkm(5)=595380.45825460922926_r8;
    fjkm(6)=4808855.0010261718786_r8;
    fjkm(7)=-5029951.7826825475971_r8;
    fjkm(8)=-20928027.009806808897_r8;
    fjkm(9)=21773603.858687892085_r8;
    fjkm(10)=52050919.691850881048_r8;
    fjkm(11)=-53926628.509575237122_r8;
    fjkm(12)=-77147061.601956281926_r8;
    fjkm(13)=79655909.133727217924_r8;
    fjkm(14)=67455358.167107401877_r8;
    fjkm(15)=-69454035.446132806377_r8;
    fjkm(16)=-32138208.559242941417_r8;
    fjkm(17)=33012717.635684926217_r8;
    fjkm(18)=6437358.4793646103367_r8;
    fjkm(19)=-6599304.6046316445590_r8;
    j= 1; k=9; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=984.36388661957607837_r8;
    fjkm(1)=-2194.2476729600457475_r8;
    fjkm(2)=-149715.34984220301505_r8;
    fjkm(3)=329977.62359907967038_r8;
    fjkm(4)=3636074.9553359345392_r8;
    fjkm(5)=-8229815.9546080161817_r8;
    fjkm(6)=-32850375.705398849013_r8;
    fjkm(7)=79594841.396295258680_r8;
    fjkm(8)=140766384.09976010426_r8;
    fjkm(9)=-388119773.63641718318_r8;
    fjkm(10)=-302391232.58882834949_r8;
    fjkm(11)=1069154026.1028829621_r8;
    fjkm(12)=267438889.51147761435_r8;
    fjkm(13)=-1738631339.5172586463_r8;
    fjkm(14)=117013574.77418801057_r8;
    fjkm(15)=1654904787.0330349261_r8;
    fjkm(16)=-452792004.09905362650_r8;
    fjkm(17)=-852646349.53093518044_r8;
    fjkm(18)=354479824.94387953335_r8;
    fjkm(19)=183646906.05281680809_r8;
    fjkm(20)=-95153470.227211795243_r8;
    j= 2; k=9; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-3445.2736031685162743_r8;
    fjkm(1)=11875.044917870854988_r8;
    fjkm(2)=615370.50617503436198_r8;
    fjkm(3)=-2111012.2880743031723_r8;
    fjkm(4)=-16085470.193130487741_r8;
    fjkm(5)=60142798.465327082670_r8;
    fjkm(6)=138071565.42237886418_r8;
    fjkm(7)=-645362876.57211229968_r8;
    fjkm(8)=-403042168.77936625406_r8;
    fjkm(9)=3392353592.8461412643_r8;
    fjkm(10)=-459521271.54126358493_r8;
    fjkm(11)=-9731832243.4128846845_r8;
    fjkm(12)=5664098916.1399176504_r8;
    fjkm(13)=15628587453.068137370_r8;
    fjkm(14)=-14586127233.110452133_r8;
    fjkm(15)=-13112286301.064475188_r8;
    fjkm(16)=18076138583.500999115_r8;
    fjkm(17)=3815036344.9012045022_r8;
    fjkm(18)=-11188129037.135692765_r8;
    fjkm(19)=1563031125.3210441483_r8;
    fjkm(20)=2774182673.1720275815_r8;
    fjkm(21)=-967769239.01720317824_r8;
    j= 3; k=9; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=9905.1616091094842886_r8;
    fjkm(1)=-46820.775577189124306_r8;
    fjkm(2)=-2040487.1579820312439_r8;
    fjkm(3)=9691735.5725892393225_r8;
    fjkm(4)=54815060.692073363805_r8;
    fjkm(5)=-309390686.57090219229_r8;
    fjkm(6)=-360017948.08027628199_r8;
    fjkm(7)=3579817310.0615042649_r8;
    fjkm(8)=-975311608.49901937973_r8;
    fjkm(9)=-19324823653.635730322_r8;
    fjkm(10)=19615966368.548239543_r8;
    fjkm(11)=52313543472.729543241_r8;
    fjkm(12)=-87319509031.594965641_r8;
    fjkm(13)=-62880701946.001201560_r8;
    fjkm(14)=187577712375.82047424_r8;
    fjkm(15)=-5014883987.8038807144_r8;
    fjkm(16)=-210187382319.75343210_r8;
    fjkm(17)=95721572764.807326435_r8;
    fjkm(18)=109424352395.85297165_r8;
    fjkm(19)=-93578868051.240421500_r8;
    fjkm(20)=-10054064393.166929203_r8;
    fjkm(21)=29497575770.435656525_r8;
    fjkm(22)=-7788016225.4016704591_r8;
    j= 4; k=9; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-24762.904022773710722_r8;
    fjkm(1)=150202.12250011534820_r8;
    fjkm(2)=5795836.3854749992775_r8;
    fjkm(3)=-35748503.915727151813_r8;
    fjkm(4)=-151894819.53780369009_r8;
    fjkm(5)=1257834576.5897312294_r8;
    fjkm(6)=283767191.88351472079_r8;
    fjkm(7)=-15254386392.322461649_r8;
    fjkm(8)=17551211231.567450237_r8;
    fjkm(9)=79224659176.160664514_r8;
    fjkm(10)=-169127446206.19948897_r8;
    fjkm(11)=-159965663833.05011647_r8;
    fjkm(12)=667618088176.46828695_r8;
    fjkm(13)=-96638419442.062469711_r8;
    fjkm(14)=-1307530021252.0141810_r8;
    fjkm(15)=998989571129.40670193_r8;
    fjkm(16)=1171444346499.5230642_r8;
    fjkm(17)=-1737201461305.1621935_r8;
    fjkm(18)=-114746520425.41058268_r8;
    fjkm(19)=1243665816084.6793073_r8;
    fjkm(20)=-512525557301.93515073_r8;
    fjkm(21)=-261796114655.33666678_r8;
    fjkm(22)=247688439610.96543843_r8;
    fjkm(23)=-52756420815.901269809_r8;
    j= 5; k=9; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=55716.534051240849124_r8;
    fjkm(1)=-415613.72155461745024_r8;
    fjkm(2)=-14629090.289521179515_r8;
    fjkm(3)=112517365.15692088569_r8;
    fjkm(4)=349569807.61470724572_r8;
    fjkm(5)=-4300749690.4008170962_r8;
    fjkm(6)=2780201489.2284702260_r8;
    fjkm(7)=53008111096.391616511_r8;
    fjkm(8)=-112833892526.88952297_r8;
    fjkm(9)=-236800084652.33408957_r8;
    fjkm(10)=948058960807.18952080_r8;
    fjkm(11)=15893204800.580178374_r8;
    fjkm(12)=-3467524908336.7843972_r8;
    fjkm(13)=3164309467171.3659456_r8;
    fjkm(14)=5582435675279.1817146_r8;
    fjkm(15)=-10541249137249.730392_r8;
    fjkm(16)=-1168774181536.2029753_r8;
    fjkm(17)=14413943369875.125772_r8;
    fjkm(18)=-7960463497916.1060492_r8;
    fjkm(19)=-7323822211977.2888688_r8;
    fjkm(20)=9405979314966.6882022_r8;
    fjkm(21)=-1275663773372.9196006_r8;
    fjkm(22)=-2930439830152.3544240_r8;
    fjkm(23)=1747630840980.1348510_r8;
    fjkm(24)=-312613977240.16973767_r8;
    j= 6; k=9; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-115412.82053471318747_r8;
    fjkm(1)=1027788.8260207103833_r8;
    fjkm(2)=33623754.973679064933_r8;
    fjkm(3)=-313584768.22511895708_r8;
    fjkm(4)=-659891574.79280520771_r8;
    fjkm(5)=12850764070.127691360_r8;
    fjkm(6)=-19440345035.583477731_r8;
    fjkm(7)=-155138555406.99154257_r8;
    fjkm(8)=516911193990.92260102_r8;
    fjkm(9)=451041007491.51124256_r8;
    fjkm(10)=-4074544650661.9635765_r8;
    fjkm(11)=3101195098601.7687528_r8;
    fjkm(12)=13186035652463.635190_r8;
    fjkm(13)=-24971415775329.260965_r8;
    fjkm(14)=-10847664957177.134063_r8;
    fjkm(15)=66205224018265.472705_r8;
    fjkm(16)=-37500185472771.584184_r8;
    fjkm(17)=-69777795804669.385671_r8;
    fjkm(18)=99311193873506.569397_r8;
    fjkm(19)=-492168373279.20159968_r8;
    fjkm(20)=-80134019127845.984412_r8;
    fjkm(21)=51160427044311.739411_r8;
    fjkm(22)=9045641917569.2467301_r8;
    fjkm(23)=-24123027505367.165512_r8;
    fjkm(24)=10768904399045.846700_r8;
    fjkm(25)=-1663085461560.5466581_r8;
    j= 7; k=9; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=110.01714026924673817_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-13886.089753717040532_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=308186.40461266239848_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-2785618.1280864546890_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=13288767.166421818329_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-37567176.660763351308_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=66344512.274729026665_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-74105148.211532657748_r8;
    fjkm(15)=0.0_r8;
    fjkm(16)=50952602.492664642206_r8;
    fjkm(17)=0.0_r8;
    fjkm(18)=-19706819.118432226927_r8;
    fjkm(19)=0.0_r8;
    fjkm(20)=3284469.8530720378211_r8;
    j= 0; k=10; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-1045.1628325578440126_r8;
    fjkm(1)=1118.5075927373418381_r8;
    fjkm(2)=159690.03216774596612_r8;
    fjkm(3)=-168947.42533689065981_r8;
    fjkm(4)=-4160516.4622709423795_r8;
    fjkm(5)=4365974.0653460506451_r8;
    fjkm(6)=43177080.985340047679_r8;
    fjkm(7)=-45034159.737397684138_r8;
    fjkm(8)=-232553425.41238182077_r8;
    fjkm(9)=241412603.52332969965_r8;
    fjkm(10)=732559944.88488535051_r8;
    fjkm(11)=-757604729.32539425138_r8;
    fjkm(12)=-1426407013.9066740733_r8;
    fjkm(13)=1470636688.7564934244_r8;
    fjkm(14)=1741470982.9710174571_r8;
    fjkm(15)=-1790874415.1120392289_r8;
    fjkm(16)=-1299291363.5629483763_r8;
    fjkm(17)=1333259765.2247248044_r8;
    fjkm(18)=541937525.75688624049_r8;
    fjkm(19)=-555075405.16917439177_r8;
    fjkm(20)=-96891860.665625115724_r8;
    fjkm(21)=99081507.234339807604_r8;
    j= 1; k=10; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=5487.1048709286810663_r8;
    fjkm(1)=-12101.885429617141199_r8;
    fjkm(2)=-991438.75239470139087_r8;
    fjkm(3)=2166230.0015798583230_r8;
    fjkm(4)=28994419.876786743130_r8;
    fjkm(5)=-64719144.968659103681_r8;
    fjkm(6)=-321629835.31147623339_r8;
    fjkm(7)=757688130.83951567540_r8;
    fjkm(8)=1749409837.5100643555_r8;
    fjkm(9)=-4544758370.9162618687_r8;
    fjkm(10)=-5113992851.9544763313_r8;
    fjkm(11)=15778214197.520607549_r8;
    fjkm(12)=7774473545.9444870052_r8;
    fjkm(13)=-33570323211.012887492_r8;
    fjkm(14)=-3804246527.4759322626_r8;
    fjkm(15)=44463088926.919594649_r8;
    fjkm(16)=-5920634247.3331925357_r8;
    fjkm(17)=-35768726949.850578829_r8;
    fjkm(18)=10834752690.813605970_r8;
    fjkm(19)=16001937124.166968265_r8;
    fjkm(20)=-6803368741.9070923418_r8;
    fjkm(21)=-3054556963.3569951737_r8;
    fjkm(22)=1577229794.0273014954_r8;
    j= 2; k=10; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-21033.902005226610754_r8;
    fjkm(1)=71551.022388357981754_r8;
    fjkm(2)=4410889.4214898588524_r8;
    fjkm(3)=-14945521.653227529678_r8;
    fjkm(4)=-139310283.22771754330_r8;
    fjkm(5)=506112018.78208167499_r8;
    fjkm(6)=1519598088.8133635683_r8;
    fjkm(7)=-6552067582.0564506220_r8;
    fjkm(8)=-6692370095.4282068536_r8;
    fjkm(9)=42444451633.414257582_r8;
    fjkm(10)=5560288488.7766793217_r8;
    fjkm(11)=-155035909620.11651130_r8;
    fjkm(12)=57543525805.054081844_r8;
    fjkm(13)=335136386281.40047184_r8;
    fjkm(14)=-242028870673.91835183_r8;
    fjkm(15)=-425158049601.57076545_r8;
    fjkm(16)=447281929473.23048294_r8;
    fjkm(17)=285539009675.26705395_r8;
    fjkm(18)=-446532174700.75534828_r8;
    fjkm(19)=-56114815650.416797381_r8;
    fjkm(20)=234199738711.39910784_r8;
    fjkm(21)=-38367666879.807869654_r8;
    fjkm(22)=-50717346515.577689017_r8;
    fjkm(23)=17618025541.849480837_r8;
    j= 3; k=10; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=65730.943766333158607_r8;
    fjkm(1)=-305976.00327882005331_r8;
    fjkm(2)=-15752484.913133094466_r8;
    fjkm(3)=73626048.766219891723_r8;
    fjkm(4)=517727813.38636845015_r8;
    fjkm(5)=-2779778621.0311819654_r8;
    fjkm(6)=-4848484539.3006943808_r8;
    fjkm(7)=38867748248.486701209_r8;
    fjkm(8)=2816113040.3833130041_r8;
    fjkm(9)=-261989764937.19714608_r8;
    fjkm(10)=193832179511.35043238_r8;
    fjkm(11)=941830771354.52418439_r8;
    fjkm(12)=-1252062600825.4805385_r8;
    fjkm(13)=-1788586580037.9230755_r8;
    fjkm(14)=3715406759129.6400037_r8;
    fjkm(15)=1338049108553.0784173_r8;
    fjkm(16)=-6078517409833.2670243_r8;
    fjkm(17)=1046877894952.0462047_r8;
    fjkm(18)=5489626110236.1661278_r8;
    fjkm(19)=-2928719118680.3861544_r8;
    fjkm(20)=-2340856193318.6421567_r8;
    fjkm(21)=2206204676067.3123808_r8;
    fjkm(22)=119173943641.45927347_r8;
    fjkm(23)=-589883942966.21075926_r8;
    fjkm(24)=154983207892.81174977_r8;
    j= 4; k=10; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-177473.54816909952824_r8;
    fjkm(1)=1058104.4704696212045_r8;
    fjkm(2)=47977824.897952560186_r8;
    fjkm(3)=-290115715.66147843083_r8;
    fjkm(4)=-1577194622.2444813090_r8;
    fjkm(5)=12040320836.774941282_r8;
    fjkm(6)=8938496432.7761692719_r8;
    fjkm(7)=-177729149538.54760085_r8;
    fjkm(8)=142764258530.36399735_r8;
    fjkm(9)=1190892990607.9998843_r8;
    fjkm(10)=-2059313111949.9121859_r8;
    fjkm(11)=-3727699663610.1911429_r8;
    fjkm(12)=10889256225145.522967_r8;
    fjkm(13)=3246312869985.9515053_r8;
    fjkm(14)=-29702334678008.056483_r8;
    fjkm(15)=12354361209407.777556_r8;
    fjkm(16)=43335047349400.917011_r8;
    fjkm(17)=-41504037524268.002306_r8;
    fjkm(18)=-28336804589555.309295_r8;
    fjkm(19)=52945932725332.519399_r8;
    fjkm(20)=-2984104941157.1313565_r8;
    fjkm(21)=-30618403352283.411013_r8;
    fjkm(22)=14099792752244.009722_r8;
    fjkm(23)=5138575314225.4923364_r8;
    fjkm(24)=-5394419195595.0379138_r8;
    fjkm(25)=1142750145697.5795143_r8;
    j= 5; k=10; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=428894.40807532385991_r8;
    fjkm(1)=-3139491.0587579525708_r8;
    fjkm(2)=-129346362.91769878378_r8;
    fjkm(3)=971658116.23636815026_r8;
    fjkm(4)=4056885717.6650654269_r8;
    fjkm(5)=-43776518263.068506188_r8;
    fjkm(6)=6886509445.7336568338_r8;
    fjkm(7)=666192636622.63242699_r8;
    fjkm(8)=-1146605162893.2335873_r8;
    fjkm(9)=-4150515189210.3819757_r8;
    fjkm(10)=12899393468582.902633_r8;
    fjkm(11)=7742402532790.6606612_r8;
    fjkm(12)=-63390570866544.420475_r8;
    fjkm(13)=31574846906545.697293_r8;
    fjkm(14)=155509830706127.21122_r8;
    fjkm(15)=-197328514020298.70121_r8;
    fjkm(16)=-161495194791368.91656_r8;
    fjkm(17)=431593958485279.11169_r8;
    fjkm(18)=-57955664355845.420187_r8;
    fjkm(19)=-444813800230214.83356_r8;
    fjkm(20)=304547176185415.82363_r8;
    fjkm(21)=164529366533754.02760_r8;
    fjkm(22)=-262168848376474.26253_r8;
    fjkm(23)=51430539333046.601717_r8;
    fjkm(24)=65933562346693.241986_r8;
    fjkm(25)=-41287526582645.050139_r8;
    fjkm(26)=7341963962580.3111652_r8;
    j= 6; k=10; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=551.33589612202058561_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-84005.433603024085289_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=2243768.1779224494292_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-24474062.725738728468_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=142062907.79753309519_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-495889784.27503030925_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=1106842816.8230144683_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-1621080552.1083370752_r8;
    fjkm(15)=0.0_r8;
    fjkm(16)=1553596899.5705800562_r8;
    fjkm(17)=0.0_r8;
    fjkm(18)=-939462359.68157840255_r8;
    fjkm(19)=0.0_r8;
    fjkm(20)=325573074.18576574902_r8;
    fjkm(21)=0.0_r8;
    fjkm(22)=-49329253.664509961973_r8;
    j= 0; k=11; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-5789.0269092812161489_r8;
    fjkm(1)=6156.5841733625632060_r8;
    fjkm(2)=1050067.9200378010661_r8;
    fjkm(3)=-1106071.5424398171230_r8;
    fjkm(4)=-32534638.579875516724_r8;
    fjkm(5)=34030484.031823816343_r8;
    fjkm(6)=403822034.97468901972_r8;
    fjkm(7)=-420138076.79184817203_r8;
    fjkm(8)=-2628163794.2543622609_r8;
    fjkm(9)=2722872399.4527176577_r8;
    fjkm(10)=10165740577.638121340_r8;
    fjkm(11)=-10496333767.154808213_r8;
    fjkm(12)=-24903963378.517825536_r8;
    fjkm(13)=25641858589.733168515_r8;
    fjkm(14)=39716473526.654258344_r8;
    fjkm(15)=-40797193894.726483060_r8;
    fjkm(16)=-41170317838.620371488_r8;
    fjkm(17)=42206049105.000758192_r8;
    fjkm(18)=26774677250.924984473_r8;
    fjkm(19)=-27400985490.712703408_r8;
    fjkm(20)=-9929978762.6658553451_r8;
    fjkm(21)=10147027478.789699178_r8;
    fjkm(22)=1603200744.0965737641_r8;
    fjkm(23)=-1636086913.2062470721_r8;
    j= 1; k=11; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=33286.904728366992856_r8;
    fjkm(1)=-72776.338288106717300_r8;
    fjkm(2)=-7048423.0820374073034_r8;
    fjkm(3)=15288988.915750383523_r8;
    fjkm(4)=243935418.08573977628_r8;
    fjkm(5)=-538504362.70138786302_r8;
    fjkm(6)=-3246894911.6396827767_r8;
    fjkm(7)=7489063194.0760509112_r8;
    fjkm(8)=21666937100.705365161_r8;
    fjkm(9)=-53983904963.062576171_r8;
    fjkm(10)=-80910564664.877465851_r8;
    fjkm(11)=229101080335.06400288_r8;
    fjkm(12)=172760876423.44066571_r8;
    fjkm(13)=-610977234886.30398648_r8;
    fjkm(14)=-187937135374.72033958_r8;
    fjkm(15)=1053702358870.4190989_r8;
    fjkm(16)=18639458829.443774842_r8;
    fjkm(17)=-1174519256075.3585225_r8;
    fjkm(18)=213630362751.48244186_r8;
    fjkm(19)=817332252922.97321022_r8;
    fjkm(20)=-266086886489.81593243_r8;
    fjkm(21)=-322968489592.27962303_r8;
    fjkm(22)=139744842706.19027127_r8;
    fjkm(23)=55347422611.580177333_r8;
    fjkm(24)=-28497920919.101275948_r8;
    j= 2; k=11; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-138695.43636819580357_r8;
    fjkm(1)=466698.94436858890046_r8;
    fjkm(2)=33739004.101605793511_r8;
    fjkm(3)=-113151831.10727233791_r8;
    fjkm(4)=-1262494179.8855194512_r8;
    fjkm(5)=4485679097.1221702841_r8;
    fjkm(6)=16876412838.414698204_r8;
    fjkm(7)=-68734331237.685231727_r8;
    fjkm(8)=-99319007190.867271605_r8;
    fjkm(9)=535143925100.64569760_r8;
    fjkm(10)=219559829042.68087919_r8;
    fjkm(11)=-2402169771537.2297740_r8;
    fjkm(12)=385225420494.61274582_r8;
    fjkm(13)=6605716163805.3629860_r8;
    fjkm(14)=-3536807746028.4626688_r8;
    fjkm(15)=-11320974870399.200872_r8;
    fjkm(16)=9487604757721.1762395_r8;
    fjkm(17)=11727696199159.893852_r8;
    fjkm(18)=-13727281916428.067929_r8;
    fjkm(19)=-6409955842808.2983124_r8;
    fjkm(20)=11468445778721.253331_r8;
    fjkm(21)=723160235150.79794197_r8;
    fjkm(22)=-5214971540434.5399686_r8;
    fjkm(23)=952560905703.62661138_r8;
    fjkm(24)=1001898723474.6755508_r8;
    fjkm(25)=-346817425242.52751961_r8;
    j= 3; k=11; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=468097.09774266083704_r8;
    fjkm(1)=-2151404.5559841446618_r8;
    fjkm(2)=-129076990.59865050915_r8;
    fjkm(3)=595316208.51428773472_r8;
    fjkm(4)=5064273602.5906671292_r8;
    fjkm(5)=-26185510869.073114681_r8;
    fjkm(6)=-61837037264.905576505_r8;
    fjkm(7)=433386998010.23271090_r8;
    fjkm(8)=186443883288.37312122_r8;
    fjkm(9)=-3537220973754.0018335_r8;
    fjkm(10)=1681730347233.8501195_r8;
    fjkm(11)=15988837394678.875714_r8;
    fjkm(12)=-16993822946683.547492_r8;
    fjkm(13)=-41340811227869.288548_r8;
    fjkm(14)=67591520066179.679649_r8;
    fjkm(15)=56621717984129.588648_r8;
    fjkm(16)=-149792388677379.47570_r8;
    fjkm(17)=-20183328379251.847125_r8;
    fjkm(18)=196651762467137.74506_r8;
    fjkm(19)=-54964390301576.370525_r8;
    fjkm(20)=-147653284026647.88789_r8;
    fjkm(21)=88925712863728.440410_r8;
    fjkm(22)=52472318964377.745533_r8;
    fjkm(23)=-54571161032830.893735_r8;
    fjkm(24)=-776744894101.81078918_r8;
    fjkm(25)=12653076888080.966521_r8;
    fjkm(26)=-3310861678129.4432166_r8;
    j= 4; k=11; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-1357481.5834537164274_r8;
    fjkm(1)=7977851.1535147665703_r8;
    fjkm(2)=419517728.92817665197_r8;
    fjkm(3)=-2495517054.7161424403_r8;
    fjkm(4)=-16720780562.050493234_r8;
    fjkm(5)=120296564778.19436333_r8;
    fjkm(6)=154404724838.89607546_r8;
    fjkm(7)=-2110123804399.9206379_r8;
    fjkm(8)=974123695127.88869692_r8;
    fjkm(9)=17443881430511.946630_r8;
    fjkm(10)=-24448434409173.656268_r8;
    fjkm(11)=-73513952038406.365892_r8;
    fjkm(12)=169708094306922.00526_r8;
    fjkm(13)=139392884698962.95542_r8;
    fjkm(14)=-607701552209729.47437_r8;
    fjkm(15)=42904363912061.406899_r8;
    fjkm(16)=1241046770899817.3524_r8;
    fjkm(17)=-768125556097572.10187_r8;
    fjkm(18)=-1403171934960438.0242_r8;
    fjkm(19)=1622188475923667.7103_r8;
    fjkm(20)=657623934547413.88722_r8;
    fjkm(21)=-1631882122398223.2875_r8;
    fjkm(22)=237083343503533.13124_r8;
    fjkm(23)=785835171244720.50574_r8;
    fjkm(24)=-396417689594574.31797_r8;
    fjkm(25)=-105868095076065.58561_r8;
    fjkm(26)=125179414423474.77661_r8;
    fjkm(27)=-26396909127729.654202_r8;
    j= 5; k=11; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=3038.0905109223842686_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-549842.32757228868713_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=17395107.553978164538_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-225105661.88941527780_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=1559279864.8792575133_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-6563293792.6192843320_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=17954213731.155600080_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-33026599749.800723140_r8;
    fjkm(15)=0.0_r8;
    fjkm(16)=41280185579.753973955_r8;
    fjkm(17)=0.0_r8;
    fjkm(18)=-34632043388.158777923_r8;
    fjkm(19)=0.0_r8;
    fjkm(20)=18688207509.295824922_r8;
    fjkm(21)=0.0_r8;
    fjkm(22)=-5866481492.0518472276_r8;
    fjkm(23)=0.0_r8;
    fjkm(24)=814789096.11831211495_r8;
    j= 0; k=12; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-34938.040875607419089_r8;
    fjkm(1)=36963.434549555675268_r8;
    fjkm(2)=7422871.4222258972763_r8;
    fjkm(3)=-7789432.9739407564011_r8;
    fjkm(4)=-269624167.08666155034_r8;
    fjkm(5)=281220905.45598032670_r8;
    fjkm(6)=3939349083.0647673616_r8;
    fjkm(7)=-4089419524.3243775468_r8;
    fjkm(8)=-30405957365.145521510_r8;
    fjkm(9)=31445477275.065026519_r8;
    fjkm(10)=141110816541.31461314_r8;
    fjkm(11)=-145486345736.39413603_r8;
    fjkm(12)=-421924022682.15660188_r8;
    fjkm(13)=433893498502.92700194_r8;
    fjkm(14)=842178293619.91844007_r8;
    fjkm(15)=-864196026786.45225550_r8;
    fjkm(16)=-1135205103443.2342838_r8;
    fjkm(17)=1162725227163.0702664_r8;
    fjkm(18)=1021645279950.6839487_r8;
    fjkm(19)=-1044733308876.1231340_r8;
    fjkm(20)=-588678536542.81848505_r8;
    fjkm(21)=601137341549.01570167_r8;
    fjkm(22)=196527129983.73688212_r8;
    fjkm(23)=-200438117645.10478028_r8;
    fjkm(24)=-28925012912.200080081_r8;
    fjkm(25)=29468205642.945621491_r8;
    j= 1; k=12; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=218362.75547254636931_r8;
    fjkm(1)=-473942.11970389194590_r8;
    fjkm(2)=-53559985.272697166145_r8;
    fjkm(3)=115466888.79018062430_r8;
    fjkm(4)=2162702487.2919505639_r8;
    fjkm(5)=-4731469254.6820607544_r8;
    fjkm(6)=-33930459549.835830283_r8;
    fjkm(7)=76986136366.180025009_r8;
    fjkm(8)=271095146839.75321729_r8;
    fjkm(9)=-654897543249.28815561_r8;
    fjkm(10)=-1244130265844.5028996_r8;
    fjkm(11)=3321026659065.3578720_r8;
    fjkm(12)=3434492363731.4649585_r8;
    fjkm(13)=-10772528238693.360048_r8;
    fjkm(14)=-5553407245149.3813559_r8;
    fjkm(15)=23184673024360.107644_r8;
    fjkm(16)=4148109873524.0835034_r8;
    fjkm(17)=-33519510690760.226852_r8;
    fjkm(18)=1766187462911.1875877_r8;
    fjkm(19)=32207800350987.663468_r8;
    fjkm(20)=-7064569616534.6127663_r8;
    fjkm(21)=-19734747129816.391118_r8;
    fjkm(22)=6780185269401.9041713_r8;
    fjkm(23)=6981112975541.6982009_r8;
    fjkm(24)=-3063627371132.2565100_r8;
    fjkm(25)=-1085299076029.5917371_r8;
    fjkm(26)=557485489473.28346831_r8;
    j= 2; k=12; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-982632.39962645866188_r8;
    fjkm(1)=3276416.0527937831379_r8;
    fjkm(2)=274430595.86767397376_r8;
    fjkm(3)=-912438377.61048048104_r8;
    fjkm(4)=-11979589299.502912885_r8;
    fjkm(5)=41816043571.570359444_r8;
    fjkm(6)=191330729548.58232223_r8;
    fjkm(7)=-746943185635.96708583_r8;
    fjkm(8)=-1416198224023.6464721_r8;
    fjkm(9)=6856956393274.5884258_r8;
    fjkm(10)=4829515891796.9969935_r8;
    fjkm(11)=-36877687475484.603376_r8;
    fjkm(12)=-2050957082915.0942624_r8;
    fjkm(13)=124378250253033.78820_r8;
    fjkm(14)=-44312193894090.533579_r8;
    fjkm(15)=-271207689501172.62955_r8;
    fjkm(16)=179586202876957.77361_r8;
    fjkm(17)=381618869449173.21877_r8;
    fjkm(18)=-359490061117882.96881_r8;
    fjkm(19)=-330393612742248.46582_r8;
    fjkm(20)=427952236370799.92495_r8;
    fjkm(21)=148052565489017.04595_r8;
    fjkm(22)=-307266140839707.15754_r8;
    fjkm(23)=-4687264091266.0864760_r8;
    fjkm(24)=123260326898726.61625_r8;
    fjkm(25)=-24376238900981.871642_r8;
    fjkm(26)=-21272360948501.370513_r8;
    fjkm(27)=7341892911307.8818418_r8;
    j= 3; k=12; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=3562042.4486459126493_r8;
    fjkm(1)=-16196313.687936474068_r8;
    fjkm(2)=-1119546446.8249808309_r8;
    fjkm(3)=5105896697.6570577310_r8;
    fjkm(4)=51476060928.630380287_r8;
    fjkm(5)=-258447407741.23544355_r8;
    fjkm(6)=-779649080899.64056288_r8;
    fjkm(7)=4982179421789.6942804_r8;
    fjkm(8)=4004696303571.8542110_r8;
    fjkm(9)=-48152095258426.146318_r8;
    fjkm(10)=10206156605278.263985_r8;
    fjkm(11)=264344253278752.35082_r8;
    fjkm(12)=-218798857956960.84499_r8;
    fjkm(13)=-868622560677829.56823_r8;
    fjkm(14)=1161601525559524.2658_r8;
    fjkm(15)=1687979001191903.2190_r8;
    fjkm(16)=-3349241418127553.8805_r8;
    fjkm(17)=-1649878672532120.5740_r8;
    fjkm(18)=5897451215285124.6365_r8;
    fjkm(19)=-125952827188585.92633_r8;
    fjkm(20)=-6433283969334751.7469_r8;
    fjkm(21)=2343935079955608.5977_r8;
    fjkm(22)=4106725553395393.7642_r8;
    fjkm(23)=-2735980005540132.5807_r8;
    fjkm(24)=-1230142327980830.8915_r8;
    fjkm(25)=1417490532431040.0017_r8;
    fjkm(26)=-23384761374805.900562_r8;
    fjkm(27)=-289892454526107.40352_r8;
    fjkm(28)=75592403781851.468191_r8;
    j= 4; k=12; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=18257.755474293174691_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-3871833.4425726126206_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=143157876.71888898129_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-2167164983.2237950935_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=17634730606.834969383_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-87867072178.023265677_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=287900649906.15058872_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-645364869245.37650328_r8;
    fjkm(15)=0.0_r8;
    fjkm(16)=1008158106865.3820948_r8;
    fjkm(17)=0.0_r8;
    fjkm(18)=-1098375156081.2233068_r8;
    fjkm(19)=0.0_r8;
    fjkm(20)=819218669548.57732864_r8;
    fjkm(21)=0.0_r8;
    fjkm(22)=-399096175224.46649796_r8;
    fjkm(23)=0.0_r8;
    fjkm(24)=114498237732.02580995_r8;
    fjkm(25)=0.0_r8;
    fjkm(26)=-14679261247.695616661_r8;
    j= 0; k=13; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-228221.94342866468364_r8;
    fjkm(1)=240393.78041152680010_r8;
    fjkm(2)=56141584.917302882999_r8;
    fjkm(3)=-58722807.212351291413_r8;
    fjkm(4)=-2362104965.8616681913_r8;
    fjkm(5)=2457543550.3409275122_r8;
    fjkm(6)=40092552189.640209230_r8;
    fjkm(7)=-41537328845.122739292_r8;
    fjkm(8)=-361511977440.11687235_r8;
    fjkm(9)=373268464511.34018528_r8;
    fjkm(10)=1977009124005.5234777_r8;
    fjkm(11)=-2035587172124.2056548_r8;
    fjkm(12)=-7053565922700.6894237_r8;
    fjkm(13)=7245499689304.7898162_r8;
    fjkm(14)=17102169035002.477337_r8;
    fjkm(15)=-17532412281166.061672_r8;
    fjkm(16)=-28732506045663.389701_r8;
    fjkm(17)=29404611450240.311097_r8;
    fjkm(18)=33500442260477.310858_r8;
    fjkm(19)=-34232692364531.459729_r8;
    fjkm(20)=-26624606760328.763181_r8;
    fjkm(21)=27170752540027.814733_r8;
    fjkm(22)=13768818045244.094179_r8;
    fjkm(23)=-14034882162060.405178_r8;
    fjkm(24)=-4179185677218.9420633_r8;
    fjkm(25)=4255517835706.9592699_r8;
    fjkm(26)=565151558036.28124143_r8;
    fjkm(27)=-574937732201.41165254_r8;
    j= 1; k=13; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=1540498.1181434866146_r8;
    fjkm(1)=-3322911.4963213577938_r8;
    fjkm(2)=-433313348.25129661430_r8;
    fjkm(3)=929240026.21742702895_r8;
    fjkm(4)=20173953055.394385937_r8;
    fjkm(5)=-43806310275.980028275_r8;
    fjkm(6)=-367752562201.24170098_r8;
    fjkm(7)=823522693625.04213554_r8;
    fjkm(8)=3453452850623.2709660_r8;
    fjkm(9)=-8147245540357.7558550_r8;
    fjkm(10)=-18967395863304.498472_r8;
    fjkm(11)=48502623842268.842654_r8;
    fjkm(12)=64652876623215.013090_r8;
    fjkm(13)=-187135422438997.88267_r8;
    fjkm(14)=-137928375585894.45832_r8;
    fjkm(15)=487895841149504.63648_r8;
    fjkm(16)=171009666849543.97695_r8;
    fjkm(17)=-877097552972882.42245_r8;
    fjkm(18)=-74254863627598.106482_r8;
    fjkm(19)=1089588154832573.5204_r8;
    fjkm(20)=-116085557257555.85969_r8;
    fjkm(21)=-919163347233503.76274_r8;
    fjkm(22)=228872931917376.68922_r8;
    fjkm(23)=502861180782827.78742_r8;
    fjkm(24)=-180138187046491.99093_r8;
    fjkm(25)=-160984522251228.28879_r8;
    fjkm(26)=71472589051967.572740_r8;
    fjkm(27)=22899647546405.161991_r8;
    fjkm(28)=-11739127546959.248774_r8;
    j= 2; k=13; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-7445740.9043601853037_r8;
    fjkm(1)=24634048.351746637287_r8;
    fjkm(2)=2366020806.0262900637_r8;
    fjkm(3)=-7808648260.6425043412_r8;
    fjkm(4)=-118977141825.00262846_r8;
    fjkm(5)=409347654830.37623610_r8;
    fjkm(6)=2227886871742.0341399_r8;
    fjkm(7)=-8418196051278.7614582_r8;
    fjkm(8)=-19993114336865.499162_r8;
    fjkm(9)=89749606285710.604304_r8;
    fjkm(10)=91027650219177.795567_r8;
    fjkm(11)=-567320456048460.00800_r8;
    fjkm(12)=-158541319151866.40491_r8;
    fjkm(13)=2286994298876364.5224_r8;
    fjkm(14)=-406108119093234.18714_r8;
    fjkm(15)=-6109006264284954.9274_r8;
    fjkm(16)=3061509093597577.3617_r8;
    fjkm(17)=10949569443954370.550_r8;
    fjkm(18)=-8409227264715248.2782_r8;
    fjkm(19)=-12972671111828071.122_r8;
    fjkm(20)=13506857028733485.383_r8;
    fjkm(21)=9542099548494263.3202_r8;
    fjkm(22)=-13665198589678693.875_r8;
    fjkm(23)=-3500771616752033.2515_r8;
    fjkm(24)=8599505926079281.8598_r8;
    fjkm(25)=-192584364085979.71888_r8;
    fjkm(26)=-3085105840164255.6689_r8;
    fjkm(27)=648294322883069.58159_r8;
    fjkm(28)=483163296698761.31750_r8;
    fjkm(29)=-166336791576720.83219_r8;
    j= 3; k=13; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=118838.42625678325312_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-29188388.122220813403_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=1247009293.5127103248_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-21822927757.529223729_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=205914503232.41001569_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-1196552880196.1815990_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=4612725780849.1319668_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-12320491305598.287160_r8;
    fjkm(15)=0.0_r8;
    fjkm(16)=23348364044581.840938_r8;
    fjkm(17)=0.0_r8;
    fjkm(18)=-31667088584785.158403_r8;
    fjkm(19)=0.0_r8;
    fjkm(20)=30565125519935.320612_r8;
    fjkm(21)=0.0_r8;
    fjkm(22)=-20516899410934.437391_r8;
    fjkm(23)=0.0_r8;
    fjkm(24)=9109341185239.8989559_r8;
    fjkm(25)=0.0_r8;
    fjkm(26)=-2406297900028.5039611_r8;
    fjkm(27)=0.0_r8;
    fjkm(28)=286464035717.67904299_r8;
    j= 0; k=14; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-1604318.7544665739172_r8;
    fjkm(1)=1683544.3719710960859_r8;
    fjkm(2)=452420015.89442260775_r8;
    fjkm(3)=-471878941.30923648336_r8;
    fjkm(4)=-21822662636.472430684_r8;
    fjkm(5)=22654002165.480904234_r8;
    fjkm(6)=425547091271.81986272_r8;
    fjkm(7)=-440095709776.83934521_r8;
    fjkm(8)=-4427161819496.8153373_r8;
    fjkm(9)=4564438154985.0886811_r8;
    fjkm(10)=28118992684610.267576_r8;
    fjkm(11)=-28916694604741.055309_r8;
    fjkm(12)=-117624507411652.86515_r8;
    fjkm(13)=120699657932218.95313_r8;
    fjkm(14)=338813510903952.89689_r8;
    fjkm(15)=-347027171774351.75500_r8;
    fjkm(16)=-688776739315164.30766_r8;
    fjkm(17)=704342315344885.53495_r8;
    fjkm(18)=997513290420732.48968_r8;
    fjkm(19)=-1018624682810589.2619_r8;
    fjkm(20)=-1023931704917833.2405_r8;
    fjkm(21)=1044308455264456.7876_r8;
    fjkm(22)=728349929088172.52737_r8;
    fjkm(23)=-742027862028795.48563_r8;
    fjkm(24)=-341600294446496.21085_r8;
    fjkm(25)=347673188569989.47682_r8;
    fjkm(26)=95048767051125.906463_r8;
    fjkm(27)=-96652965651144.909104_r8;
    fjkm(28)=-11888257482283.680284_r8;
    fjkm(29)=12079233506095.466313_r8;
    j= 1; k=14; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=11631310.969882660899_r8;
    fjkm(1)=-24956069.513924483156_r8;
    fjkm(2)=-3719130469.3827566264_r8;
    fjkm(3)=7939241569.2440612457_r8;
    fjkm(4)=197650420583.57805736_r8;
    fjkm(5)=-426477178381.34693109_r8;
    fjkm(6)=-4137136219101.0505865_r8;
    fjkm(7)=9165629658162.2739663_r8;
    fjkm(8)=44999979919399.924736_r8;
    fjkm(9)=-104192738635599.46794_r8;
    fjkm(10)=-290053332678279.44824_r8;
    fjkm(11)=717931728117708.95938_r8;
    fjkm(12)=1184950942733150.9332_r8;
    fjkm(13)=-3238133498156090.6407_r8;
    fjkm(14)=-3148099361614567.8423_r8;
    fjkm(15)=10004238940145809.174_r8;
    fjkm(16)=5326672157182975.4416_r8;
    fjkm(17)=-21713978561461112.072_r8;
    fjkm(18)=-4997511985428331.4237_r8;
    fjkm(19)=33440445545533127.273_r8;
    fjkm(20)=429328409587666.98618_r8;
    fjkm(21)=-36372499368723031.528_r8;
    fjkm(22)=5419838346824587.4483_r8;
    fjkm(23)=27328510015364670.604_r8;
    fjkm(24)=-7462027883028047.7909_r8;
    fjkm(25)=-13500043636525530.253_r8;
    fjkm(26)=5000259547410615.2462_r8;
    fjkm(27)=3946328556046746.4962_r8;
    fjkm(28)=-1769166076587921.0596_r8;
    fjkm(29)=-517354048506128.35163_r8;
    fjkm(30)=264752449010576.61885_r8;
    j= 2; k=14; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=832859.30401628929898_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-234557963.52225152478_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=11465754899.448237157_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-229619372968.24646817_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=2485000928034.0853236_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-16634824724892.480519_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=74373122908679.144941_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-232604831188939.92523_r8;
    fjkm(15)=0.0_r8;
    fjkm(16)=523054882578444.65558_r8;
    fjkm(17)=0.0_r8;
    fjkm(18)=-857461032982895.05140_r8;
    fjkm(19)=0.0_r8;
    fjkm(20)=1026955196082762.4888_r8;
    fjkm(21)=0.0_r8;
    fjkm(22)=-889496939881026.44181_r8;
    fjkm(23)=0.0_r8;
    fjkm(24)=542739664987659.72270_r8;
    fjkm(25)=0.0_r8;
    fjkm(26)=-221349638702525.19597_r8;
    fjkm(27)=0.0_r8;
    fjkm(28)=54177510755106.049005_r8;
    fjkm(29)=0.0_r8;
    fjkm(30)=-6019723417234.0054450_r8;
    j= 0; k=15; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=-12076459.908236194835_r8;
    fjkm(1)=12631699.444247054368_r8;
    fjkm(2)=3870206398.1171501588_r8;
    fjkm(3)=-4026578373.7986511753_r8;
    fjkm(4)=-212116465639.79238740_r8;
    fjkm(5)=219760302239.42454551_r8;
    fjkm(6)=4707197145849.0525974_r8;
    fjkm(7)=-4860276727827.8835762_r8;
    fjkm(8)=-55912520880766.919782_r8;
    fjkm(9)=57569188166122.976664_r8;
    fjkm(10)=407553205759865.77271_r8;
    fjkm(11)=-418643088909794.09305_r8;
    fjkm(12)=-1970887757079997.3409_r8;
    fjkm(13)=2020469839019116.7709_r8;
    fjkm(14)=6629237688884787.8691_r8;
    fjkm(15)=-6784307576344081.1526_r8;
    fjkm(16)=-15953173918642561.995_r8;
    fjkm(17)=16301877173694858.432_r8;
    fjkm(18)=27867483571944089.170_r8;
    fjkm(19)=-28439124260599352.538_r8;
    fjkm(20)=-35429954264855305.864_r8;
    fjkm(21)=36114591062243814.190_r8;
    fjkm(22)=32466638305657465.126_r8;
    fjkm(23)=-33059636265578149.421_r8;
    fjkm(24)=-20895477102024899.324_r8;
    fjkm(25)=21257303545350005.806_r8;
    fjkm(26)=8964660367452270.4366_r8;
    fjkm(27)=-9112226793253953.9006_r8;
    fjkm(28)=-2302544207092007.0827_r8;
    fjkm(29)=2338662547595411.1154_r8;
    fjkm(30)=267877692066913.24230_r8;
    fjkm(31)=-271890841011735.91260_r8;
    j= 1; k=15; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    fjkm(0)=6252951.4934347970025_r8;
    fjkm(1)=0.0_r8;
    fjkm(2)=-2001646928.1917763315_r8;
    fjkm(3)=0.0_r8;
    fjkm(4)=110997405139.17901279_r8;
    fjkm(5)=0.0_r8;
    fjkm(6)=-2521558474912.8546213_r8;
    fjkm(7)=0.0_r8;
    fjkm(8)=31007436472896.461417_r8;
    fjkm(9)=0.0_r8;
    fjkm(10)=-236652530451649.25168_r8;
    fjkm(11)=0.0_r8;
    fjkm(12)=1212675804250347.4165_r8;
    fjkm(13)=0.0_r8;
    fjkm(14)=-4379325838364015.4378_r8;
    fjkm(15)=0.0_r8;
    fjkm(16)=11486706978449752.110_r8;
    fjkm(17)=0.0_r8;
    fjkm(18)=-22268225133911142.562_r8;
    fjkm(19)=0.0_r8;
    fjkm(20)=32138275268586241.200_r8;
    fjkm(21)=0.0_r8;
    fjkm(22)=-34447226006485144.698_r8;
    fjkm(23)=0.0_r8;
    fjkm(24)=27054711306197081.241_r8;
    fjkm(25)=0.0_r8;
    fjkm(26)=-15129826322457681.181_r8;
    fjkm(27)=0.0_r8;
    fjkm(28)=5705782159023670.8096_r8;
    fjkm(29)=0.0_r8;
    fjkm(30)=-1301012723549699.4268_r8;
    fjkm(31)=0.0_r8;	
    fjkm(32)=135522158703093.69029_r8;
    j= 0; k=16; d=j+2*k; fjk(j,k)= un(d)*pol(fjkm,d,v);
    END SUBROUTINE fjkproc16

    FUNCTION startingpser(mu,x,y)
    !Computes a starting value for the backward summation of the series in pser}
    USE Someconstants
    REAL(r8) :: mu, x, y
    REAL(r8) :: n, n1, lnx, lny, mulnmu
    INTEGER ::  startingpser, a, b
    mulnmu= mu*log(mu); lnx=log(x); lny=log(y);
    IF (x<2.0_r8) THEN
      n= x+5 
    ELSE
      n=1.5*x;
    ENDIF 
    n1=0; a=0; b= 0;
    DO WHILE (abs(n-n1)>1)
      n1=n;
      n=ps(mu,mulnmu,lnx,y,lny,n,a,b);
    ENDDO
    n=n+1;
    IF (mu+n>y) THEN
      IF (y>mu) THEN
        a=1 
      ELSE 
        b=1;
      ENDIF 
      n1=0;
      DO WHILE (abs(n-n1)>1)   
        n1=n;
        n=ps(mu,mulnmu,lnx,y,lny,n,a,b);
      ENDDO
    ENDIF  
    startingpser=int(n)+1
    END FUNCTION startingpser

    FUNCTION ps(mu,mulnmu,lnx,y,lny,n, a, b)
    USE Someconstants
    REAL(r8) :: mu,mulnmu,lneps,lnx,y,lny,ps, n, f
    INTEGER :: a, b 
    lneps=log(epss) 
    IF ((a==0).AND.(b==0)) THEN
      f=(n-lneps)/(log(n)-lnx) 
    ELSEIF ((a==0).AND.(b==1)) THEN
      f=(2*n-lneps+mulnmu-mu*log(mu+n))/(log(n)-lnx-lny+log(mu+n)) 
    ELSEIF ((a==1).AND.(b==0)) THEN
      f=(2*n-lneps-y+mu*lny-mu*log(mu+n)+mu)/(log(n)-lnx-lny+log(mu+n));
    ENDIF
    ps= f
    END FUNCTION ps
   
    SUBROUTINE hypfun(x,sinh,cosh);
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x,sinh,cosh,ss
    REAL(r8) :: ax, y, f, f2
    ax=abs(x);
    IF (ax<0.21_r8) THEN
      IF (ax<0.07_r8) THEN
	y=x*x
      ELSE
	y=x*x/9.0_r8;
      ENDIF
      f=2.0_r8+y*(y*28+2520.0_r8)/(y*(y+420)+15120.0_r8);
      f2=f*f;
      sinh=2*x*f/(f2-y);
      cosh=(f2+ y)/(f2-y);
      IF (ax>=0.07_r8) THEN
        ss=2.0_r8*sinh/3.0_r8
	f=ss*ss
	sinh=sinh*(1.0_r8+f/3.0_r8);
	cosh=cosh*(1.0_r8+f)
      ENDIF
    ELSE
      y=exp(x);
      f=1.0_r8/y;
      cosh=(y+f)/2.0_r8;
      sinh=(y-f)/2.0_r8
    END IF
    END SUBROUTINE hypfun        

    FUNCTION ignega(n, x, eps)
    !-----------------------------------------------------
    ! Computes the Upper Incomplete Gama(1/2-n,x), x >= 0, 
    ! n=0,1,2, ...
    !-----------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, eps, ignega    
    REAL(r8) :: a, delta, g, p, q, r, s, t, tau, ro
    INTEGER :: n,k
    a=0.5-n; delta=epss/100.0_r8;
    IF (x>1.5_r8) THEN 
      p= 0.0_r8;
      q=(x-1-a)*(x+1-a);
      r=4*(x+1-a);
      s=1-a;
      ro=0.0_r8;
      t=1.0_r8; k=0;
      g=1.0_r8;
      DO WHILE (((t/g)>delta).AND.(k<1000))
        p=p+s;
        q=q+r;
        r=r+8;
        s=s+2;
        tau=p*(1+ro);
        ro=tau/(q-tau);
        t=ro*t;
        g=g+t;
        k=k+1;
      ENDDO
      g=g*exp(a*log(x))/(x+1-a)
    ELSE
      t=1; s=1.0_r8/a; k=1;
      DO WHILE ((abs(t/s)>delta).AND.(k<1000))
        t=-x*t/k; s=s+t/(k+a); k= k+1
      ENDDO
      g=sqrtpi;
      DO k=1,n 
        g=g/(0.5_r8-k); 
      ENDDO
      IF (x>0.0_r8) THEN
        g=exp(x)*(g-exp(a*log(x))*s)
      ENDIF 
    ENDIF
    ignega= g
    END FUNCTION ignega

    FUNCTION startkbes(x,eps)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x,eps
    REAL(r8) :: y, p, q, r, s, c, del, r2
    INTEGER :: startkbes
    IF (eps<machtol) THEN
      del=-log(machtol/2.0_r8)
    ELSE
      del=-log(eps/2.0_r8);
    ENDIF
    p=x+del;
    q=p/x;
    IF (q<2.0_r8) THEN
      r=log((q+1.0_r8)/(q-1.0_r8))/2.0_r8;
      y=r*(1.0_r8+2.0_r8/(1.0_r8+r*(q+1.0_r8/q)))
    ELSE
      r=2.0_r8*x/p;
      r2=r*r
      y=r*(1.0_r8+r2*r2/45.0_r8)
    ENDIF
    CALL hypfun(y,s,c);
    startkbes=1+int(x/(2.0_r8*s*s))
    END FUNCTION  startkbes

    FUNCTION startijbes(x,n,t,eps)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x,eps
    REAL(r8) :: p, q, r, y, del
    INTEGER :: startijbes, n,t, s
    IF (x<=0.0_r8) THEN
      startijbes=0
    ELSE
      s=2*t-1;
    ENDIF
    IF (eps<machtol) THEN
      del=-log(machtol/2.0_r8)
    ELSE
      del=-log(eps/2.0_r8);
    ENDIF
    p=del/x-t;
    r=n/x;
    IF ((r>1.0_r8).OR.(t==1)) THEN
      q=sqrt(r*r+s);
      r=r*log(q+r)-q
    ELSE
      r=0;
    ENDIF
    q=del/(2.0_r8*x)+r;
    IF (p>q) THEN
      r=p
    ELSE
      r=q;
    ENDIF
    y=alfinv(t,r,p,q)
    CALL hypfun(y,p,q)
    IF (t==0) THEN
      s=int(x*q)+1
    ELSE
      s=int(x*p)+1;
    ENDIF
    IF (MOD(s,2)>0) THEN
      s=s+1;
    ENDIF
    startijbes=s
    END FUNCTION startijbes
	
    FUNCTION alfinv(t,r,p,q)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: alfinv, r, p, q, a, b, a2, lna
    INTEGER :: t
    IF ((t+r)<2.7_r8) THEN
      IF (t==0) THEN
	a=exp(log(3.0_r8*r)/3.0_r8);
	a2=a*a;
	b=a*(1.0_r8+a2*(-1.0_r8/30.0_r8+0.004312_r8* a2))
      ELSE
	a=sqrt(2.0_r8*(1.0_r8+r));
	a2=a*a;
	b=a/(1.0_r8+a2/8.0_r8)
      ENDIF
    ELSE
      a=log(0.7357589_r8*(r+t));
      lna=log(a)/a
      b=1.0_r8+a+log(a)*(1.0_r8/a-1.0_r8)+0.5_r8*lna*lna
    ENDIF
    DO WHILE (abs(a/b-1.0_r8)>1.0e-2_r8)			
      a=b;
      b=fi(a,r,t,q)
    ENDDO
    alfinv=b
    END FUNCTION alfinv

    FUNCTION falfa(al,r,t,df)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: falfa, al, r, df, sh, ch
    INTEGER :: t
    CALL hypfun(al,sh,ch)
    IF (t==1) THEN
      falfa=al*sh/ch-1.0_r8-r/ch;
      df=(sh+(al+r*sh)/ch)/ch
    ELSE
      falfa=al-(sh+r)/ch;
      df=sh*(r+sh)/(ch*ch)
    ENDIF
    END FUNCTION falfa

    FUNCTION fi(al,r,t,q)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: fi, al, r, q, p
    INTEGER :: t
    p=falfa(al,r,t,q);
    fi=al-p/q
    END FUNCTION fi

    FUNCTION recipgam(x,q,r)	
    USE Someconstants
    IMPLICIT NONE		
    !recipgam(x)=1/gamma(x+1)=1 + x * (q + x * r); -0.5<=x<=0.5}
    REAL(r8) :: recipgam,x,q,r
    REAL(r8) :: t, tx, c(0:8)
    IF (x==0) THEN
      q=0.5772156649015328606e-0_r8;
      r=-0.6558780715202538811e-0_r8
    ELSE
      c(0)=+1.142022680371167841_r8;
      c(1)=-6.5165112670736881e-3_r8;
      c(2)=-3.087090173085368e-4_r8;
      c(3)=+3.4706269649043e-6_r8;
      c(4)=-6.9437664487e-9_r8;
      c(5)=-3.67795399e-11_r8;
      c(6)=+1.356395e-13_r8;
      c(7)=+3.68e-17_r8;
      c(8)=-5.5e-19_r8;
      tx=2.0_r8*x
      t=2*tx*tx-1;
      q=chepolsum(8,t,c);
      c(0)=-1.270583625778727532_r8;
      c(1)=+2.05083241859700357e-2_r8;
      c(2)=-7.84761097993185e-5_r8;
      c(3)=-5.377798984020e-7_r8;
      c(4)=+3.8823289907e-9_r8;
      c(5)=-2.6758703e-12_r8;
      c(6)=-2.39860e-14_r8;
      c(7)=+3.80e-17_r8;
      c(8)=+4e-20_r8;
      r=chepolsum(8,t,c)
    END IF			
    recipgam=1+x*(q+x*r)
    END FUNCTION recipgam   
    
    FUNCTION xpowy(x,y)
    IMPLICIT NONE
    REAL(r8) :: x,y,xpowy
    xpowy=x**y
    END FUNCTION xpowy

    FUNCTION fractio(x,n,r,s)
    IMPLICIT NONE
    INTEGER n,k
    REAL(r8) :: x, fractio, r(0:8), s(0:8), a, b
    a=r(n); b=1
    DO k=n-1,0,-1 
      a=a*x+r(k); b=b*x+s(k) 
    ENDDO
    fractio=a/b
    END FUNCTION fractio
                        
    FUNCTION zetaxy(x,y)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, y, zetaxy
    REAL(r8) :: ck(0:10), w, z, z2, zeta, S, t
    REAL(r8) :: x2, x3, x4, x5, x6, x7, x8, x9, x10, x2p1
    INTEGER :: k 
    z= (y-x-1.0_r8);
    x2=x*x
    x3=x2*x
    x4=x3*x
    x5=x4*x
    x6=x5*x
    x7=x6*x
    x8=x7*x
    x9=x8*x
    x10=x9*x
    x2p1=2*x+1.0_r8
    IF (abs(z)<0.05_r8) THEN    
      ck(0)=1.0_r8;
      ck(1)=-(1.0_r8/3.0_r8)*(3*x+1);
      ck(2)=(1.0_r8/36.0_r8)*(72*x2+42*x+7);
      ck(3)=-(1.0_r8/540.0_r8)*(2700*x3+2142*x2+657*x+73);
      ck(4)=(1.0_r8/12960.0_r8)*(1331+15972*x+76356*x2+177552*x3+181440.0_r8*x4);
      ck(5)=-(1.0_r8/272160.0_r8)*(22409+336135.0_r8*x+2115000.0_r8*x2+7097868.0_r8*x3+13105152.0_r8*x4&
             +11430720.0_r8*x5);
      ck(6)=(1.0_r8/5443200.0_r8)*(6706278.0_r8*x+52305684.0_r8*x2+228784392.0_r8*x3&
            +602453376.0_r8*x4+935038080.0_r8*x5+718502400.0_r8*x6+372571.0_r8);
      ck(7)=-(1.0_r8/16329600.0_r8)*(953677.0_r8+20027217.0_r8*x+186346566.0_r8*x2+1003641768.0_r8*x3&
             +3418065864.0_r8*x4+7496168976.0_r8*x5+10129665600.0_r8*x6+7005398400.0_r8*x7);
      ck(8)=(1.0_r8/783820800.0_r8)*(39833047.0_r8+955993128.0_r8*x+1120863744000.0_r8*x8&
             +10332818424.0_r8*x2+66071604672.0_r8*x3+275568952176.0_r8*x4&
             +776715910272.0_r8*x5+1472016602880.0_r8*x6+1773434373120.0_r8*x7);
      ck(9)=-(1.0_r8/387991296000.0_r8)*(17422499659.0_r8+470407490793.0_r8*x+3228423729868800.0_r8*x8+&
             1886413681152000.0_r8*x9+5791365522720.0_r8*x2+42859969263000.0_r8*x3+211370902874640.0_r8*x4+&
             726288467241168.0_r8*x5+1759764571151616.0_r8*x6+2954947944510720.0_r8*x7);
      ck(10)=(1.0_r8/6518253772800.0_r8)*(261834237251.0_r8+7855027117530.0_r8*x+&
              200149640441008128.0_r8*x8+200855460151664640.0_r8*x9+109480590367948800.0_r8*x10+&
              108506889674064.0_r8*x2+912062714644368.0_r8*x3+5189556987668592.0_r8*x4+&
              21011917557260448.0_r8*x5+61823384007654528.0_r8*x6+132131617757148672.0_r8*x7);
      z2=z/(x2p1*x2p1);  
      S=1; t=1; k=1;
      DO WHILE ((abs(t)> 1.0e-15_r8).AND.(k < 11)) 
        t=ck(k)*z2**k
        S= S+t; k= k+1 
      ENDDO
      zeta=-z/sqrt(x2p1)*S 
    ELSE
      w=sqrt(1.0_r8+4*x*y);
      zeta=sqrt(2.0_r8*(x+y-w-log(2.0_r8*y/(1.0_r8+w))));
      IF (x+1<y) THEN
        zeta=-zeta 
      ENDIF
    ENDIF
    zetaxy=zeta
    END FUNCTION zetaxy
   
    FUNCTION chepolsum(n,t,ak)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: t, ak(0:8)
    INTEGER, INTENT(IN) :: n
    REAL(r8) :: chepolsum, u0, u1, u2, s, tt;
    INTEGER :: k
    u0=0; u1=0; k=n; tt=t+t;
    DO WHILE (k>=0)
      u2=u1; 
      u1=u0; 
      u0=tt*u1-u2+ ak(k); 
      k= k-1 
    ENDDO
    s=(u0-u2)/2.0_r8
    chepolsum=s
    END FUNCTION chepolsum

    FUNCTION oddchepolsum(n, x, ak)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x, ak(0:100)
    INTEGER, INTENT(IN) :: n
    REAL(r8) :: oddchepolsum, h, r, s, y;
    INTEGER :: k
    IF (n==0) THEN 
      s=ak(0)*x 
    ELSEIF  (n == 1) THEN
      s=x*(ak(0)+ak(1)*(4*x*x - 3)) 
    ELSE
      y=2*(2*x*x - 1); 
      r= ak(n); h= ak(n-1)+ r*y;  
      k=n-2;
      DO WHILE (k >= 0)
        s=r; r= h; h= ak(k)+r*y-s; 
        k= k-1 
      ENDDO
      s=x*(h-r)
    ENDIF
    oddchepolsum=s
    END FUNCTION oddchepolsum

    FUNCTION logoneplusx(t)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: t
    REAL(r8) ::  logoneplusx, ck(0:100), c, p, p2, pj, x, y
    INTEGER :: j
    IF ((-0.2928_r8<t).AND.(t<0.4142_r8)) THEN
      p= twoexp14
      p=(p-1.0_r8)/(p+1.0_r8);
      pj=p; ck(0)= pj; p2= p*p; j=1; c=1;
      DO WHILE ((abs(c)> 1.0e-20).AND.(j<1000))   
        pj=pj*p2; c=pj/(2.0_r8*j+1.0_r8); 
        ck(j)=c; j= j+1 
      ENDDO
      x=t/(2.0_r8+t)*(1.0_r8+p2)/(2.0_r8*p);
      y=4*oddchepolsum(j-1, x, ck)
    ELSE 
      y=log(1.0_r8+t) 
    ENDIF
    logoneplusx=y
    END FUNCTION logoneplusx

    FUNCTION xminsinx(x) 
    !  {(x-sin(x))/(x^3/6)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x
    REAL(r8) :: xminsinx, f, fk(0:8), t
    IF (abs(x)>1.0_r8) THEN
      f=6*(x-sin(x))/(x*x*x) 
    ELSE
      fk(0)=1.95088260487819821294e-0_r8;
      fk(1)=-0.244124470324439564863e-1_r8;
      fk(2)=0.14574198156365500e-3_r8;
      fk(3)=-0.5073893903402518e-6_r8;
      fk(4)=0.11556455068443e-8_r8;
      fk(5)=-0.185522118416e-11_r8;
      fk(6)=0.22117315e-14_r8;
      fk(7)=-0.2035e-17_r8;
      fk(8)=0.15e-20_r8;
      t=2*x*x-1.0_r8;
      f=chepolsum(8,t,fk)
    ENDIF
    xminsinx=f
    END FUNCTION xminsinx
    	
    FUNCTION trapsum(a, b, h, d, xis2, mu,wxis,ys)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a, b, h, d, xis2, mu,wxis, ys
    REAL(r8) :: trapsum, b0, inte, s, aa, bb
    s=0;
    b0=b
    IF (d==0.0_r8) THEN
      CALL integrand(a,b0,inte,xis2,mu,wxis,ys)
      s=inte/2.0_r8;
      aa=a+h; 
      bb=b-h/2.0_r8
    ELSE
      aa=a+d; 
      bb=b 
    ENDIF
    DO WHILE ((aa<bb).AND.(aa< b0))
      CALL integrand(aa,b0,inte,xis2,mu,wxis,ys)
      s=s+inte; 
      aa=aa+h 
    ENDDO
    trapsum=s*h
    END FUNCTION trapsum

    FUNCTION trap(a, b, e, xis2, mu, wxis, ys)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a, b, e, xis2, mu, wxis, ys
    REAL(r8) :: trap, h, p, q, v, nc
    h=(b-a)/8.0_r8;
    p=trapsum(a,b,h,0.0_r8,xis2,mu,wxis,ys);
    nc=0; v=1;
    DO WHILE(((v>e).AND.(nc<10)).OR.(nc<=2)) 
      nc=nc+1;
      q=trapsum(a,b,h,h/2.0_r8,xis2,mu,wxis,ys);
      IF (abs(q)>0.0_r8) THEN  
        v=abs(p/q-1) 
      ELSE
        v=0 
      ENDIF
      h=h/2.0_r8;
      p=(p+q)/2.0_r8
    ENDDO
    trap=p 
   END FUNCTION trap
 
    SUBROUTINE integrand(theta,b0,inte,xis2,mu,wxis,ys)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: theta,xis2,mu,wxis,ys
    REAL(r8), INTENT(INOUT) :: b0
    REAL(r8), INTENT(OUT) :: inte
    REAL(r8) :: eps, lneps, psitheta, p, ft, f, theta2, sintheta, costheta,&
                term1, term2, sinth, dr, s2, wx, ts, rtheta, xminsinxtheta
    eps=1.0e-16_r8; lneps=log(eps);
    IF (theta>b0) THEN
      f=0 
    ELSE IF (abs(theta)<1.0e-10_r8) THEN
      rtheta=(1.0_r8+wxis)/(2*ys);
      theta2=theta*theta;
      psitheta=-wxis*theta2*0.5_r8 
      f=rtheta/(1.0_r8-rtheta)*exp(mu*psitheta) 
    ELSE
      theta2=theta*theta; 
      sintheta=sin(theta);  
      costheta=cos(theta);
      ts=theta/sintheta; 
      s2=sintheta*sintheta; 
      wx=sqrt(ts*ts+xis2);
      xminsinxtheta=xminsinx(theta);
      p=xminsinxtheta*theta2*ts/6.0_r8; 
      term1=((p*(ts+1)-theta2-s2*xis2)/(costheta*wx+wxis));
      p=(p*(1.0_r8+(ts+1)/(wx+wxis))/(1+wxis));
      term2=-logoneplusx(p);
      p= term1+term2
      psitheta=p
      f=mu*psitheta;
      IF (f>lneps) THEN
        f=exp(f) 
      ELSE 
 	f=0; 
        b0=min(theta,b0) 
      ENDIF
      rtheta=(ts+wx)/(2*ys); 
      sinth=sin(theta/2.0_r8)
      p=(2*theta*sinth*sinth-xminsinxtheta*theta2*theta/6)/(2*ys*s2);
      dr=p*(1+ts/wx);
      p=((dr*sintheta+(costheta-rtheta)*rtheta)/(rtheta*(rtheta-2*costheta)+1));
      ft=p
      f= f*ft	
    ENDIF
    inte=f
    END SUBROUTINE integrand

    SUBROUTINE qser(mu,x,y,p,q,ierro)
    !----------------------------------------------------
    ! Computes the series expansion for Q.
    ! For computing the incomplete gamma functions we
    ! use the routine incgam included in the module
    ! IncgamFI. Reference: A. Gil, J. Segura and 
    ! NM Temme, Efficient and accurate algorithms for 
    ! the computation and inversion of the incomplete 
    ! gamma function ratios. SIAM J Sci Comput.  
    !----------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: mu, x, y, p, q
    REAL(r8) :: delta, lh0, h0, q0, xy
    REAL(r8) :: x1, p1, q1, t, k, S, a
    INTEGER :: ierr, n, m, ierro
    ierro=0
    ierr=0
    CALL cdfgamC(1,mu,y,p,q,ierr)
    q0=q; 
    lh0=mu*log(y)-y-loggam(mu+1.0_r8)
    IF ((lh0>log(dwarf)).AND.(x<100.0_r8)) THEN
      h0=exp(lh0);
      n=0; xy=x*y;  delta= epss/100.0_r8;
      DO WHILE ((q0/q>delta).AND.(n<1000))
        q0=x*(q0+h0)/(n+1.0_r8); 
        h0=xy*h0/((n+1.0_r8)*(mu+n+1)); 
        q=q+q0;
        n=n+1
      ENDDO
      q=exp(-x)*q; p=1.0_r8-q;
    ELSE
      ! Computing Q forward
      x1=y
      S=0; t=1; k=0;
      m=0
      DO WHILE ((k<10000).AND.(m==0)) 
        a=mu+k
        CALL  cdfgamC(1,a,x1,p1,q1,ierr)
        t=dompart(k,x,.false.)*q1
        S= S+t; k= k+1;     
        IF ((S==0.0_r8).AND.(k>150)) m=1
        IF (S>0) THEN
          IF (((t/S)<1.e-16_r8).AND.(k>10)) m=1
        ENDIF
      ENDDO
      IF (ierr==0) THEN
        q=S
        p=1-q 
      ELSE
        q=0.0_r8
        p=1.0_r8
        ierro=1
      ENDIF
    ENDIF
    END SUBROUTINE qser


    SUBROUTINE pser(mu,x,y,p,q,ierro)
    !----------------------------------------------
    ! Computes backward the series expansion for P
    ! For computing the incomplete gamma functions we
    ! use the routine incgam included in the module
    ! IncgamFI. Reference: A. Gil, J. Segura and 
    ! NM Temme, Efficient and accurate algorithms for 
    ! the computation and inversion of the incomplete 
    ! gamma function ratios. SIAM J Sci Comput. 
    !----------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: mu, x, y, p, q
    REAL(r8) :: lh0, h0, p1, s, xy
    REAL(r8) :: a, x1, q1, expo, facto, t
    INTEGER :: n, k, nnmax, ierr, ierro
    ierro=0
    ierr=0
    xy=x*y; 
    nnmax=startingpser(mu,x,y);
    n=1+nnmax;
    lh0=-x-y+n*log(x)+(n+mu)*log(y)-loggam(mu+n+1.0_r8)-loggam(n+1.0_r8)
    IF (lh0<log(dwarf)) THEN
      x1=y
      expo=exp(-x)
      facto=1.0_r8 
      S=0; t=1; k=0;
      k=startingpser(mu,x,y)+1;
      DO WHILE ((k>0).AND.(ierr==0)) 
        a=mu+k
        facto=factor(x,k)
        CALL  cdfgamC(1,a,x1,p1,q1,ierr)
        t=facto*p1
        S= S+t;
        k= k-1;  
      ENDDO
      IF (ierr==0) THEN
        CALL  cdfgamC(1,mu,x1,p1,q1,ierr)
        S=S+p1;
        p=S*expo 
        q=1.0_r8-p
      ELSE
        ierro=1
        p=0.0_r8; q=1.0_r8
      ENDIF
    ELSE
      h0=exp(lh0);
      CALL cdfgamC(1,mu+n,y,p,q,ierr)
      IF (ierr==0) THEN
        p1=p*exp(-x+n*log(x)-loggam(n+1.0_r8));
        p=0;
        DO WHILE (n>0)
          h0=h0*n*(mu+n)/xy;         
          p1=n*p1/x+h0;
          p= p+p1;
          n= n-1;
        ENDDO
        q=1.0_r8-p;
      ELSE
        ierro=1
        p=0.0_r8; q=1.0_r8
      ENDIF
    ENDIF
    END SUBROUTINE pser

    SUBROUTINE prec(mu,x,y,p,q,ierro)
    !----------------------------------------
    ! Computes the backward recursion for P
    !----------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: mu, x, y, p, q
    REAL(r8) :: b, nu, mur, cmu, xi, p0, p1
    INTEGER :: ierr, ierro, n, n1, n2, n3
    ierro=0
    ierr=0 
    b=1.0_r8; nu=y-x+b*b+b*sqrt(2*(x+y)+b*b);
    n1=int(mu); n2=int(nu)+2; n3= n2-n1;
    mur= mu+n3;
    xi=2.0_r8*sqrt(x*y);
    cmu=sqrt(y/x)*fc(mur,xi)
    ! Numerical quadrature 
    CALL cdfPQtrap(mur+1,x,y,p1,q,ierr)      
    CALL cdfPQtrap(mur+0,x,y,p0,q,ierr);
    IF (ierr==0) THEN
      DO n=0,n3-1
        p=((1+cmu)*p0-p1)/cmu;
        p1= p0; p0= p;
        cmu=y/(mur-n-1+x*cmu);
      ENDDO
      q=1-p;
    ELSE
      p=0.0_r8
      q=1.0_r8
      ierro=1
    ENDIF
    END SUBROUTINE prec

    SUBROUTINE qrec(mu,x,y,p,q,ierro)
    !------------------------------------------
    ! Computes the forward recursion for Q
    !------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: mu, x, y, p, q
    REAL(r8) :: b, c, nu, mur, xi, q0, q1, cmu(0:300)
    INTEGER :: ierr, ierro, n, n1, n2, n3
    ierro=0  
    ierr=0
    b=1.0_r8; nu= y-x+b*(b-sqrt(2*(x+y)+b*b));
    IF (nu<5) THEN
      IF (x<200.0_r8) THEN
        CALL qser(mu,x,y,p,q,ierr)
      ELSE
        CALL prec(mu,x,y,p,q,ierr)
      ENDIF
    ELSE
      n1=int(mu); n2=int(nu)-1; n3=n1-n2;
      mur= mu-n3;
      xi=2*sqrt(x*y);
      cmu(0)=sqrt(y/x)*fc(mu,xi)
      DO n=1,n3 
        cmu(n)=y/(mu-n+x*cmu(n-1));
      ENDDO
      ! Numerical quadrature 
      CALL cdfPQtrap(mur-1,x,y,p,q0,ierr)      
      CALL cdfPQtrap(mur+0,x,y,p,q1,ierr);
      IF (ierr==0) THEN
        DO n=1,n3
          c=cmu(n3+1-n); q=(1.0_r8+c)*q1-c*q0;
          q0=q1; q1= q
        ENDDO
        p=1.0_r8-q;
      ELSE
        q=0.0_r8
        p=1.0_r8
        ierro=1
      ENDIF
    ENDIF
    END SUBROUTINE qrec
   
    SUBROUTINE pqasyxy(mu,x,y,p,q,ierro)
! ----------------------------------------------------------
! Computes P_{\mu}(x,y) and Q_{\mu}(x,y)
! using an asymptotic expansion for large xi.
! ----------------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: mu, x, y, p, q
    REAL(r8) :: delta, c, psi, psi0, an, xi, sqxi, pq,&
                rho, mulrho, rhom, rhop, sigmaxi, er, nu, mu2,&
                rhomu, s, tnm1, phin(0:100), bn(0:100)
    INTEGER :: n, n0, nrec, ierro
    ierro=0
    IF (y>=x) THEN
      s=1.0_r8 
    ELSE
      s=-1.0_r8;
    ENDIF
    delta=epss/100.0_r8;
    xi=2*sqrt(x*y); sqxi=sqrt(xi); rho=sqrt(y/x); 
    sigmaxi=((y-x)*(y-x))/(x+y+xi);
    mulrho=mu*log(rho)
    IF ((mulrho<log(dwarf)).OR.(mulrho>log(giant))) THEN
      IF (s==1.0_r8) THEN 
        q=0.0_r8; p=1.0_r8
      ELSE
        p=0.0_r8; q=1.0_r8 
      ENDIF
      ierro=1
    ELSE 
      rhomu=exp(mulrho);
      er=errorfunction(sqrt(sigmaxi),.true.,.true.);
      psi0=0.5*rhomu*er/sqrt(rho);
      nu=2*mu-1; rhom=nu*(rho-1); rhop=2*(rho+1); mu2=4*mu*mu;
      c=s*rhomu/sqrt(8.0_r8*pi);
      an=sqxi; n=0; n0=100;
      bn(0)= 1;
      DO WHILE ((abs(bn(n))>delta).AND.(n<n0))
        n=n+1; 
        tnm1=2*n-1;
        an=(mu2-tnm1*tnm1)*an/(8*n*xi);
        bn(n)=an*(rhom-n*rhop)/(rho*(nu+2*n));
      ENDDO
      n0=n;
      nrec=int(sigmaxi)+1;
      IF (nrec>n0) nrec=n0;
      phin(nrec)=exp((nrec-0.5)*log(sigmaxi))*ignega(nrec,sigmaxi,epss);
      DO n=nrec+1,n0 
        phin(n)=(-sigmaxi*phin(n-1)+1)/(n-0.5);
      ENDDO
      DO n=nrec-1,1,-1 
        phin(n)=(1-(n+0.5)*phin(n+1))/sigmaxi;
      ENDDO
      pq=psi0;
      DO n=1,n0
        c=-c;
        psi=c*bn(n)*phin(n);
        pq=pq+psi;
      ENDDO
      pq=pq*exp(-sigmaxi);
      IF (s==1.0_r8) THEN 
        q= pq; p=1.0_r8-q
      ELSE
        p=pq; q=1.0_r8-p 
      ENDIF
    ENDIF
    END SUBROUTINE pqasyxy

    SUBROUTINE pqasymu(mu0,x0,y0,p,q,ierro)
    !--------------------------------------------------
    ! Computes P_(mu)(x, y) and Q_(mu)(x, y) by using 
    ! the large mu asymptotic expansion
    !--------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: mu0, x0, y0, p, q
    REAL(r8) :: mu, muk(0:16), x, y, r, s, u, zeta, zetaj, bk
    REAL(r8) :: lexpor, psik(-1:16), fjk(0:16,0:16) 
    INTEGER :: ierro, a, b, j, k, t
    ierro=0
    mu= mu0-1.0_r8; x=x0/mu; y= y0/mu;
    zeta= zetaxy(x,y);
    IF (zeta < 0.0_r8) THEN 
      a=1 
    ELSE
      a=-1;
    ENDIF
    u=1.0_r8/sqrt(2.0_r8*x+1.0_r8);
    CALL fjkproc16(u,fjk)
    zeta=a*zeta;
    r=zeta*sqrt(mu/2.0_r8);
    psik(0)=sqrt(pi/(2.0_r8*mu))*errorfunction(-r,.true.,.false.);  
    s=psik(0); 
    lexpor=-mu*0.5_r8*zeta*zeta
    IF ((lexpor<log(dwarf)).OR.(lexpor>log(giant))) THEN    
      IF (a==1) THEN 
        q=0.0_r8; 
        p=1.0_r8; 
      ELSE 
        p=0.0_r8; 
        q=1.0_r8;  
      ENDIF
      ierro=1
    ELSE
      r=exp(lexpor);
      psik(-1)=0.0_r8; muk(0)=1.0_r8;
      bk=s; k=1; zetaj= 1.0_r8;
      DO WHILE ((abs(bk/s)> 1.e-30_r8).AND.(k <= 16)) 
        muk(k)=mu*muk(k-1);
        psik(k)=((k-1)*psik(k-2)+r*zetaj)/mu;
        bk=0; b=1; zetaj=-zeta*zetaj;
        DO j=0,k  
          IF ((a==-1).AND.(b==-1)) THEN
            t=-1 
          ELSE
            t=1; 
          ENDIF
          b=-b;
          bk=bk+t*fjk(j,k-j)*psik(j)/muk(k-j); 
        ENDDO
        s=s+bk; 
        k=k+1;
      ENDDO
      r=sqrt(mu/(2.0_r8*pi))*s;
      IF (a==1) THEN 
        q=r; 
        p=1.0_r8-q 
      ELSE 
        p=r; 
        q=1.0_r8-p 
      ENDIF
    ENDIF
    END SUBROUTINE pqasymu


    SUBROUTINE cdfPQtrap(mu,x,y,p,q,ierr)
    !---------------------------------------------------------------
    ! Computes the P and Q using an integral representation which
    ! is approximated using the trapezoidal rule
    !---------------------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: mu
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(IN) :: y
    REAL(r8), INTENT(OUT) :: p
    REAL(r8), INTENT(OUT) :: q
    INTEGER, INTENT(INOUT) :: ierr 
    REAL(r8) ::  pq, a, b, zeta, epstrap, xs, ys, xis2, wxis
    xs=x/mu; ys=y/mu;
    xis2=4*xs*ys;
    wxis=sqrt(1.0_r8+xis2);
    a=0.0_r8; b=3.0_r8;   
    epstrap=1.0e-13_r8;
    pq=trap(a,b,epstrap,xis2,mu,wxis,ys); 
    zeta=zetaxy(xs,ys);
    IF ((-mu*0.5_r8*zeta*zeta)<log(dwarf)) THEN
      IF (y >x+mu) THEN
        p=1.0_r8
        q=0.0_r8
      ELSE
        p=0.0_r8
        q=1.0_r8 
      ENDIF
      ierr=1
    ELSE
      pq=pq*exp(-mu*0.5_r8*zeta*zeta)/pi;
      IF (zeta<0.0_r8) THEN
        q=pq; p=1-q 
      ELSE
        p=-pq; q= 1-p 
      ENDIF
    ENDIF
    END SUBROUTINE cdfPQtrap

    SUBROUTINE rec1step(mu,x,y,p,q,ierro)
    !-----------------------------------------
    ! 0<mu<1:
    ! Computes 1 step backward of the TTRR
    !-----------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: mu, x, y, p, q
    REAL(r8) :: mu1, mu2, cmu, xi, p0, p1, q0, q1
    INTEGER :: ierro, ierr1, ierr2 
    ierro=0
    mu2=mu+2
    mu1=mu+1
    xi=2.0_r8*sqrt(x*y);
    cmu=sqrt(y/x)*fc(mu1,xi)
    ! function values for the TTRR
    IF (y<x+mu) THEN
      CALL cdfgamNC(1,mu2,x,y,p1,q,ierr1)
      CALL cdfgamNC(1,mu1,x,y,p0,q,ierr2)
      IF ((ierr1==0).AND.(ierr2==0)) THEN
        p=((1+cmu)*p0-p1)/cmu;
        q=1-p;
      ELSE
        p=0.0_r8
        q=1.0_r8
        ierro=1
      ENDIF
    ELSE
      CALL cdfgamNC(1,mu2,x,y,p,q1,ierr1)
      CALL cdfgamNC(1,mu1,x,y,p,q0,ierr2)
      IF ((ierr1==0).AND.(ierr2==0)) THEN
        q=((1+cmu)*q0-q1)/cmu;
        p=1-q;
      ELSE
        p=1.0_r8
        q=0.0_r8
        ierro=1
      ENDIF
    ENDIF
    END SUBROUTINE rec1step

    RECURSIVE FUNCTION sinh(x,eps) RESULT(sinhh)
    USE Someconstants  
    IMPLICIT NONE
    !to compute hyperbolic function sinh (x)}
    REAL(r8) :: sinhh, x, eps
    REAL(r8) :: ax, e, t, x2, y
    INTEGER  :: u, k
    ax=abs(x);
    IF (x==0.0_r8) THEN
      y=0.0_r8
    ELSEIF (ax<0.12) THEN
      e=eps*0.1_r8;
      x2=x*x;
      y=1;
      t=1;
      u=0;
      k=1;
      DO WHILE(t>e)
        u=u+8*k-2;
        k=k+1;
        t=t*x2/u;
        y=y+t
      END DO
      y=x*y
    ELSEIF (ax<0.36_r8) THEN
      t=sinh(x*0.333333333333333333333333333333_r8,eps);
      y=t*(3.0_r8+4.0_r8*t*t);
    ELSE
      t=exp(x);
      y=(t-1.0_r8/t)*0.5_r8
    ENDIF
    sinhh=y
    END FUNCTION sinh

    FUNCTION exmin1(x,eps)
    USE Someconstants  
    IMPLICIT NONE
    !computes (exp(x)-1)/x 
    REAL(r8) :: exmin1, x, eps
    REAL(r8) :: t, y
    IF (x==0) THEN
      y=1.0_r8
    ELSEIF ((x<-0.69_r8).OR.(x>0.4_r8)) THEN
      y=(exp(x)-1.0_r8)/x
    ELSE
      t=x*0.5_r8;
      y=exp(t)*sinh(t,eps)/t
    ENDIF
    exmin1=y
    END FUNCTION exmin1

    FUNCTION exmin1minx(x,eps)
    USE Someconstants    
    IMPLICIT NONE
    !computes (exp(x)-1-x)/(0.5*x*x) 
    REAL(r8) :: exmin1minx, x, eps
    REAL(r8) :: t, t2, y
    IF (x==0) THEN
      y=1.0_r8
    ELSEIF (abs(x)>0.9_r8) THEN
      y=(exp(x)-1.0_r8-x)/(x*x*0.5_r8)
    ELSE
      t=sinh(x*0.5_r8,eps);
      t2=t*t;
      y=(2*t2+(2.0_r8*t*sqrt(1.0_r8+t2)-x))/(x*x*0.5_r8)
    ENDIF
    exmin1minx=y
    END FUNCTION exmin1minx

    FUNCTION lnec(x)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: lnec, ln1, x,  y0, z, e2, r, s
    !x>-1; lnec:=ln1:=ln(1+x)-x
    z=logoneplusx(x);
    y0=z-x;
    e2=exmin1minx(z,machtol);
    s=e2*z*z/2;
    r=(s+y0)/(s+1+z);
    ln1=y0-r*(6.0_r8-r)/(6.0_r8-4.0_r8*r);
    lnec=ln1
    END FUNCTION lnec

    FUNCTION alfa(x)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, alfa, lnx
    lnx=log(x)
    IF (x>0.25_r8) THEN
      alfa=x+0.25_r8
    ELSEIF (x>=dwarf) THEN
      alfa=-0.6931_r8/lnx
    ELSE
      alfa=-0.6931_r8/log(dwarf)
    ENDIF
    END FUNCTION alfa

    FUNCTION  dompart(a,x,qt)
    !dompart is approx. of  x^a*exp(-x)/gamma(a+1)   
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: dompart, a, x
    REAL(r8) :: lnx, c, dp, la, lambda, mu, r 
    LOGICAL :: qt
    lnx=log(x)
    IF (a<=1.0_r8) THEN                     
      r=-x+a*lnx
    ELSE
      IF (x==a) THEN
        r=0
      ELSE
        la=x/a
        r=a*(1.0_r8-la+log(la))
      ENDIF
      r=r-0.5_r8*log(6.2832_r8*a)
    ENDIF
    IF (r<explow) THEN
      dp=0.0_r8
    ELSE
      dp=exp(r)
    ENDIF
    IF (qt) THEN
      dompart=dp
    ELSE
      IF (a<8.0_r8) THEN
        dompart=exp(a*lnx-x)/gamma(a+1.0_r8)
      ELSE
        dompart=1.0_r8/(sqrt(a*twopi)*gamstar(a))
        lambda=x/a
        IF ((lambda>0.3_r8).AND.(lambda<2.36_r8)) THEN
          mu=lambda-1.0_r8
          c=lnec(mu);
          dompart=dompart*exp(a*c)          
        ELSE
          dompart=dompart*exp(a*log(lambda)+a-x)
        ENDIF
      ENDIF
    ENDIF
    END FUNCTION dompart
    
    FUNCTION  pqasymp (a,x,dp,p)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: pqasymp, a, x, dp
    REAL(r8) :: y, mu, eta, u, v 
    INTEGER :: s
    LOGICAL :: p
    IF (dp==0.0_r8) THEN
      IF (p) THEN
        pqasymp=0.0_r8
      ELSE
        pqasymp=1.0_r8
      ENDIF
    ELSE
      IF (p) THEN
        s=-1
      ELSE
        s=1
      ENDIF
      mu=(x-a)/a;
      y=-lnec(mu)
      IF (y<0) THEN
        eta=0.0_r8
      ELSE
        eta=sqrt(2.0_r8*y)
      ENDIF
      y=y*a;
      v=sqrt(abs(y));
      IF (mu<0.0_r8) THEN     
        eta=-eta
        v=-v
      ENDIF  
      u=0.5_r8*errorfunction(s*v,.true.,.false.);
      v=s*exp(-y)*saeta(a,eta)/sqrt(2.0_r8*pi*a);
      pqasymp=u+v
    ENDIF
    END FUNCTION pqasymp
                    
    FUNCTION saeta(a,eta)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: a, eta, saeta, y, s, t
    REAL(r8) :: eps, fm(0:26), bm(0:26)
    INTEGER  ::  m
    eps=epss
    fm(0)=1.0_r8;
    fm(1)=-1.0_r8/3.0_r8;
    fm(2)=1.0_r8/12.0_r8;
    fm(3)=-2.0_r8/135.0_r8;
    fm(4)=1.0_r8/864.0_r8;
    fm(5)=1.0_r8/ 2835.0_r8;
    fm(6)=-139.0_r8/777600.0_r8;
    fm(7)=1.0_r8/25515.0_r8;
    fm(8)=-571.0_r8/261273600.0_r8;
    fm(9)=-281.0_r8/151559100.0_r8;
    fm(10)=8.29671134095308601e-7_r8;
    fm(11)=-1.76659527368260793e-7_r8;
    fm(12)=6.70785354340149857e-9_r8;
    fm(13)=1.02618097842403080e-8_r8;
    fm(14)=-4.38203601845335319e-9_r8;
    fm(15)=9.14769958223679023e-10_r8;
    fm(16)=-2.55141939949462497e-11_r8;
    fm(17)=-5.83077213255042507e-11_r8;
    fm(18)=2.43619480206674162e-11_r8;
    fm(19)=-5.02766928011417559e-12_r8;
    fm(20)=1.10043920319561347e-13_r8;
    fm(21)=3.37176326240098538e-13_r8;
    fm(22)=-1.39238872241816207e-13_r8;
    fm(23)=2.85348938070474432e-14_r8;
    fm(24)=-5.13911183424257258e-16_r8;
    fm(25)=-1.97522882943494428e-15_r8;
    fm(26)= 8.09952115670456133e-16_r8;
    bm(25)=fm(26);
    bm(24)=fm(25);
    DO m=24,1,-1 
      bm(m-1)=fm(m)+(m+1)*bm(m+1)/a;
    ENDDO
    s=bm(0);
    t=s;
    y=eta;
    m=1;
    DO WHILE ((abs(t/s)>eps).AND.(m<25))
      t=bm(m)*y;
      s=s+t;
      m=m+1;
      y=y*eta
    ENDDO 
    saeta=s/(1.0_r8+bm(1)/a);
    END FUNCTION saeta

    FUNCTION qfraction(a,x,dp)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: qfraction, a, x, dp
    REAL(r8) :: eps, g, p, q, r, s, t, tau, ro
    eps=epss
    IF (dp==0) THEN
      q=0.0_r8
    ELSE
      p=0;
      q=(x-1.0_r8-a)*(x+1.0_r8-a);
      r=4*(x+1.0_r8-a);
      s=1.0_r8-a;
      ro=0.0_r8;
      t=1.0_r8;
      g=1.0_r8;
      DO WHILE(abs(t/g)>=eps)
        p=p+s;
        q=q+r;
        r=r+8;
        s=s+2;
        tau=p*(1.0_r8+ro);
        ro=tau/(q-tau);
        t=ro*t;
        g=g+t
      ENDDO
      q=(a/(x+1.0_r8-a))*g*dp; 
    ENDIF
    qfraction= q
    END FUNCTION qfraction

    FUNCTION qtaylor(a,x,dp)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: qtaylor, a, x, dp
    REAL(r8) :: eps, lnx, p, q, r, s, t, u, v
    eps=epss
    lnx=log(x)
    IF (dp==0) THEN
      q=0.0_r8
    ELSE
      r=a*lnx;
      q=r*exmin1(r,eps);  
      s=a*(1.0_r8-a)*auxgam(a); 
      q=(1-s)*q;
      u=s-q;            
      p=a*x;
      q=a+1;
      r=a+3;
      t=1.0_r8;
      v=1.0_r8;
      DO WHILE (abs(t / v) > eps)
        p=p+x;
        q=q+r;
        r=r+2;
        t=-p*t/q;
        v=v+t
      ENDDO
      v=a*(1.0_r8-s)*exp((a+1.0_r8)*lnx)*v/(a+1.0_r8);
      q=u+v
    ENDIF
    qtaylor=q
    END FUNCTION qtaylor

    FUNCTION ptaylor(a,x,dp)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: ptaylor,a,x,dp
    REAL(r8) :: eps,p,c,r
    eps=epss
    IF (dp==0) THEN
      p=0.0_r8
    ELSE
      p=1.0_r8
      c=1.0_r8
      r=a
      DO WHILE ((c/p)>eps)
        r=r+1
        c=x*c/r
        p=p+c
      ENDDO
      p=p*dp
     ENDIF
    ptaylor=p
    END FUNCTION ptaylor

    FUNCTION eps1(eta)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: eps1, eta, la, ak(0:4), bk(0:4) 
    IF (abs(eta)<1.0) THEN
      ak(0)=-3.333333333438e-1_r8;  bk(0)= 1.000000000000e+0_r8;     
      ak(1)=-2.070740359969e-1_r8;  bk(1)= 7.045554412463e-1_r8;     
      ak(2)=-5.041806657154e-2_r8;  bk(2)= 2.118190062224e-1_r8;     
      ak(3)=-4.923635739372e-3_r8;  bk(3)= 3.048648397436e-2_r8;     
      ak(4)=-4.293658292782e-5_r8;  bk(4)= 1.605037988091e-3_r8;     
      eps1=ratfun(eta,ak,bk)
    ELSE
      la=lambdaeta(eta);
      eps1=log(eta/(la-1.0_r8))/eta
    ENDIF
    END FUNCTION eps1
     FUNCTION eps2(eta)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: eta, eps2, x, lnmeta, ak(0:4), bk(0:4)
    IF (eta < -5.0) THEN
      x=eta*eta;
      lnmeta=log(-eta)
      eps2=(12.0_r8-x-6.0*(lnmeta*lnmeta))/(12.0*x*eta)
    ELSEIF (eta<-2.0) THEN
      ak(0)=-1.72847633523e-2_r8;  bk(0)=1.00000000000e+0_r8;     
      ak(1)= -1.59372646475e-2_r8;  bk(1)= 7.64050615669e-1_r8;     
      ak(2)= -4.64910887221e-3_r8;  bk(2)= 2.97143406325e-1_r8;     
      ak(3)= -6.06834887760e-4_r8;  bk(3)= 5.79490176079e-2_r8;     
      ak(4)= -6.14830384279e-6_r8;  bk(4)= 5.74558524851e-3_r8;     
      eps2= ratfun(eta,ak,bk)
    ELSEIF (eta < 2.0) THEN
      ak(0)=-1.72839517431e-2_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)=-1.46362417966e-2_r8;  bk(1)= 6.90560400696e-1_r8;     
      ak(2)=-3.57406772616e-3_r8;  bk(2)= 2.49962384741e-1_r8;     
      ak(3)=-3.91032032692e-4_r8;  bk(3)= 4.43843438769e-2_r8;     
      ak(4)=2.49634036069e-6_r8;   bk(4)= 4.24073217211e-3_r8;     
      eps2= ratfun(eta,ak,bk)
   ELSEIF (eta < 1000.0) THEN
      ak(0)= 9.99944669480e-1_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= 1.04649839762e+2_r8;  bk(1)= 1.04526456943e+2_r8;     
      ak(2)= 8.57204033806e+2_r8;  bk(2)= 8.23313447808e+2_r8;     
      ak(3)= 7.31901559577e+2_r8;  bk(3)= 3.11993802124e+3_r8;     
      ak(4)= 4.55174411671e+1_r8;  bk(4)= 3.97003311219e+3_r8;     
      x=1.0_r8/eta;
      eps2=ratfun(x,ak,bk)/(-12.0*eta)
    ELSE
      eps2=-1.0/(12.0*eta)
    ENDIF
    END FUNCTION

    FUNCTION eps3(eta)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: eta, eps3, eta3, x, y, ak(0:4), bk(0:4)
    IF (eta <-8.0) THEN
      x=eta*eta
      y=log(-eta)/eta;
      eps3=(-30.0+eta*y*(6.0_r8*x*y*y-12.0+x))/(12.0*eta*x*x)
    ELSEIF (eta <-4.0) THEN
      ak(0)= 4.95346498136e-2_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= 2.99521337141e-2_r8;  bk(1)= 7.59803615283e-1_r8;     
      ak(2)= 6.88296911516e-3_r8;  bk(2)= 2.61547111595e-1_r8;     
      ak(3)= 5.12634846317e-4_r8;  bk(3)= 4.64854522477e-2_r8;     
      ak(4)= -2.01411722031e-5_r8; bk(4)= 4.03751193496e-3_r8;     
      eps3=ratfun(eta,ak,bk)/(eta*eta)
    ELSEIF (eta <-2.0) THEN
      ak(0)=4.52313583942e-3_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)=1.20744920113e-3_r8;  bk(1)= 9.12203410349e-1_r8;     
      ak(2)=-7.89724156582e-5_r8; bk(2)= 4.05368773071e-1_r8;     
      ak(3)=-5.04476066942e-5_r8; bk(3)= 9.01638932349e-2_r8;     
      ak(4)=-5.35770949796e-6_r8; bk(4)= 9.48935714996e-3_r8;     
      eps3=ratfun(eta,ak,bk)
    ELSEIF  (eta < 2.0) THEN
      ak(0)= 4.39937562904e-3_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= 4.87225670639e-4_r8;  bk(1)= 7.94435257415e-1_r8;     
      ak(2)= -1.28470657374e-4_r8; bk(2)= 3.33094721709e-1_r8;     
      ak(3)= 5.29110969589e-6_r8;  bk(3)= 7.03527806143e-2_r8;     
      ak(4)= 1.57166771750e-7_r8;  bk(4)= 8.06110846078e-3_r8;     
      eps3= ratfun(eta,ak,bk)
    ELSEIF (eta < 10.0) THEN
      ak(0)= -1.14811912320e-3_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= -1.12850923276e-1_r8;  bk(1)= 1.42482206905e+1_r8;     
      ak(2)= 1.51623048511e+0_r8;   bk(2)= 6.97360396285e+1_r8;     
      ak(3)= -2.18472031183e-1_r8;  bk(3)= 2.18938950816e+2_r8;     
      ak(4)= 7.30002451555e-2_r8;   bk(4)= 2.77067027185e+2_r8;     
      x= 1.0_r8/eta;
      eps3= ratfun(x,ak,bk)/(eta*eta)
    ELSEIF (eta < 100.0) THEN
      ak(0)= -1.45727889667e-4_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= -2.90806748131e-1_r8;  bk(1)= 1.39612587808e+2_r8;     
      ak(2)= -1.33085045450e+1_r8;  bk(2)= 2.18901116348e+3_r8;     
      ak(3)= 1.99722374056e+2_r8;   bk(3)= 7.11524019009e+3_r8;     
      ak(4)= -1.14311378756e+1_r8;  bk(4)= 4.55746081453e+4_r8;     
      x= 1.0_r8/eta;
      eps3= ratfun(x,ak,bk)/(eta*eta)
    ELSE
     eta3=eta*eta*eta
     eps3=-log(eta)/(12.0*eta3)
    ENDIF
    END FUNCTION eps3

    FUNCTION lambdaeta(eta)
! lambdaeta is the positive number satisfying 
! eta^2/2=lambda-1-ln(lambda)
! with sign(lambda-1)=sign(eta); 
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: eta, lambdaeta, ak(6), q, r, s, L, la
    REAL(r8) :: L2, L3, L4, L5
    s=eta*eta*0.5_r8
    IF (eta==0) THEN
      la=1 
    ELSEIF (eta < -1) THEN 
      r=exp(-1-s);
      ak(1)=1.0_r8;
      ak(2)=1.0_r8;
      ak(3)=1.5_r8;
      ak(4)=2.66666666666666666666666666667_r8;
      ak(5)=5.20833333333333333333333333333_r8;
      ak(6)=10.8_r8;
      la=r*(ak(1)+r*(ak(2)+r*(ak(3)+r*(ak(4)+r*(ak(5)+r*ak(6))))))
    ELSEIF (eta<1) THEN 
      ak(1)=1.0_r8;
      ak(2)=0.333333333333333333333333333333_r8;
      ak(3)=0.0277777777777777777777777777778_r8;
      ak(4)=-0.00370370370370370370370370370370_r8;
      ak(5)=0.000231481481481481481481481481481_r8;
      ak(6)=0.0000587889476778365667254556143445_r8;
      r=eta;
      la=1.0_r8+r*(ak(1)+r*(ak(2)+r*(ak(3)+r*(ak(4)+r*(ak(5)+r*ak(6))))))
    ELSE
      r=11+s; L=log(r); la=r+L; r=1.0_r8/r;
      L2=L*L
      L3=L2*L
      L4=L3*L 
      L5=L4*L
      ak(1)= 1;
      ak(2)=(2-L)*0.5_r8;
      ak(3)=(-9*L+6+2*L2)*0.166666666666666666666666666667_r8;
      ak(4)= -(3*L3+36*L-22*L2-12)*0.0833333333333333333333333333333_r8;
      ak(5)=(60+350*L2-300*L-125*L3+12*L4)*0.0166666666666666666666666666667_r8;
      ak(6)=-(-120-274*L4+900*L-1700*L2+1125*L3+20*L5)*&
            0.00833333333333333333333333333333_r8;
      la=la+L*r*(ak(1)+r*(ak(2)+r*(ak(3)+r*(ak(4)+r*(ak(5)+r*ak(6))))))
    ENDIF
    r= 1;
    IF (((eta>-3.5_r8).AND.(eta<-0.03_r8)).OR.((eta>0.03_r8).AND.(eta<40.0_r8))) THEN
      r=1.0_r8; 
      q=la;
      DO WHILE (r>1.0e-8_r8)
        la=q*(s+log(q))/(q-1.0_r8);
        r=abs(q/la-1.0_r8); 
        q=la
      ENDDO
    ENDIF
    lambdaeta=la
    END FUNCTION lambdaeta

      RECURSIVE FUNCTION auxloggam(x) RESULT(auxloggamm)
    !function g in ln(Gamma(1+x))=x*(1-x)*g(x), 0<=x<=1}
    USE Someconstants
    IMPLICIT NONE    
    REAL(r8) :: auxloggamm, x
    REAL(r8) :: ak(0:25)
    REAL(r8) :: g, t
    IF (x<-1.0_r8) THEN 
      g=giant
    ELSEIF (abs(x)<=dwarf) THEN 
      g=-eulmasc
    ELSEIF (abs(x - 1.0_r8)<=machtol) THEN
      g=eulmasc-1.0_r8
    ELSEIF (x<0) THEN
      g=-(x*(1+x)*auxloggam(x+1.0_r8)+logoneplusx(x))/(x*(1.0_r8-x))
    ELSEIF (x<1) THEN
      ak(0)=-0.98283078605877425496_r8;
      ak(1)=0.7611416167043584304e-1_r8;
      ak(2)=-0.843232496593277796e-2_r8;
      ak(3)=0.107949372632860815e-2_r8;
      ak(4)=-0.14900748003692965e-3_r8;
      ak(5)=0.2151239988855679e-4_r8;
      ak(6)=-0.319793298608622e-5_r8;
      ak(7)=0.48516930121399e-6_r8;
      ak(8)=-0.7471487821163e-7_r8;
      ak(9)=0.1163829670017e-7_r8;
      ak(10)=-0.182940043712e-8_r8;
      ak(11)= 0.28969180607e-9_r8;
      ak(12)=-0.4615701406e-10_r8;
      ak(13)= 0.739281023e-11_r8;
      ak(14)= -0.118942800e-11_r8;
      ak(15)= 0.19212069e-12_r8;
      ak(16)= -0.3113976e-13_r8;
      ak(17)= 0.506284e-14_r8;
      ak(18)= -0.82542e-15_r8;
      ak(19)= 0.13491e-15_r8;
      ak(20)= -0.2210e-16_r8;
      ak(21)= 0.363e-17_r8;
      ak(22)= -0.60e-18_r8;
      ak(23)= 0.98e-19_r8;
      ak(24)= -0.2e-19_r8;
      ak(25)= 0.3e-20_r8;
      t=2*x-1;
      g=chepolsum(25, t, ak)
    ELSEIF (x<1.5) THEN
      g=(logoneplusx(x-1.0_r8)+(x-1.0_r8)*(2.0_r8-x)*auxloggam(x-1.0_r8))/(x*(1.0_r8-x));
    ELSE
      g=(log(x)+(x-1.0_r8)*(2.0_r8-x)*auxloggam(x-1.0_r8))/(x*(1.0_r8-x));
    ENDIF
    auxloggamm= g
    END FUNCTION auxloggam

    RECURSIVE FUNCTION auxgam(x) RESULT(auxgamm)
    USE Someconstants
    IMPLICIT NONE    
    !function g in 1/gamma(x+1)=1+x*(x-1)*g(x), -1<=x<=1}
    REAL(r8) :: auxgamm, x
    REAL(r8) :: t, dr(0:17)
    IF (x<0.0_r8) THEN
      auxgamm=-(1.0_r8+(1.0_r8+x)*(1.0_r8+x)*auxgam(1.0_r8+x))/(1.0_r8-x)
    ELSE
      dr(0)= -1.013609258009865776949_r8;
      dr(1)= 0.784903531024782283535e-1_r8;
      dr(2)= 0.67588668743258315530e-2_r8;
      dr(3)= -0.12790434869623468120e-2_r8;
      dr(4)= 0.462939838642739585e-4_r8;
      dr(5)= 0.43381681744740352e-5_r8;
      dr(6)= -0.5326872422618006e-6_r8;
      dr(7)= 0.172233457410539e-7_r8;
      dr(8)= 0.8300542107118e-9_r8;
      dr(9)= -0.10553994239968e-9_r8;
      dr(10)= 0.39415842851e-11_r8;
      dr(11)= 0.362068537e-13_r8;
      dr(12)= -0.107440229e-13_r8;
      dr(13)= 0.5000413e-15_r8;
      dr(14)= -0.62452e-17_r8;
      dr(15)= -0.5185e-18_r8;
      dr(16)= 0.347e-19_r8;
      dr(17)= -0.9e-21_r8;
      t=2*x-1.0_r8;
      auxgamm=chepolsum(17,t,dr);
     ENDIF
    END FUNCTION auxgam

    FUNCTION lngam1(x)
    !ln(gamma(1+x)), -1<=x<=1}
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, lngam1   
    lngam1=-logoneplusx(x*(x-1.0_r8)*auxgam(x))
    END FUNCTION lngam1

    FUNCTION stirling(x)            
    !Stirling series, function corresponding with}
    !asymptotic series for log(gamma(x))}
    !that is:  1/(12x)-1/(360x**3)...; x>= 3}
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: stirling, x, a(0:17), c(0:6), z
    IF (x<dwarf) THEN
      stirling=giant
    ELSEIF (x<1.0_r8) THEN
      stirling= lngam1(x)-(x+0.5_r8)*log(x)+x-lnsqrttwopi
    ELSEIF (x<2.0_r8) THEN
      stirling=lngam1(x-1.0_r8)-(x-0.5_r8)*log(x)+x-lnsqrttwopi
    ELSEIF (x<3.0_r8) THEN
      stirling=lngam1(x-2.0_r8)-(x-0.5_r8)*log(x)+x&
               -lnsqrttwopi+log(x-1)
    ELSEIF (x<12.0_r8) THEN
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
      stirling=chepolsum(17,z,a)/(12.0_r8*x);
    ELSE
      z=1.0_r8/(x*x);
      IF (x<1000.0_r8) THEN
        c(0)=0.25721014990011306473e-1_r8;
        c(1)=0.82475966166999631057e-1_r8;
        c(2)=-0.25328157302663562668e-2_r8;
        c(3)=0.60992926669463371e-3_r8;
        c(4)=-0.33543297638406e-3_r8;
        c(5)=0.250505279903e-3_r8;
        c(6)=0.30865217988013567769_r8;
        stirling=((((((c(5)*z+c(4))*z+c(3))*z+c(2))*z+c(1))*z+c(0))/(c(6)+z)/x)
      ELSE
        stirling=(((-z*0.000595238095238095238095238095238_r8+&
                 0.000793650793650793650793650793651_r8)*z&
                -0.00277777777777777777777777777778_r8)*z+&
                 0.0833333333333333333333333333333_r8)/x
      ENDIF
    ENDIF 
    END FUNCTION stirling

    FUNCTION ratfun(x,ak,bk)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: ratfun, x, ak(0:4), bk(0:4), p, q
    p= ak(0)+x*(ak(1)+x*(ak(2)+x*(ak(3)+x*ak(4))));
    q= bk(0)+x*(bk(1)+x*(bk(2)+x*(bk(3)+x*bk(4))));
    ratfun=p/q
    END FUNCTION ratfun
 END MODULE GammaCHI


 
