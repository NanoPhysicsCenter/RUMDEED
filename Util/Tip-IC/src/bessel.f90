  MODULE BesselJY
  ! Computation of the Bessel functions Ja(x), Ya(x) and 
  ! their derivatives
  IMPLICIT NONE
  INTEGER, PARAMETER  :: r8 = KIND(0.0d0)
  PRIVATE
  PUBLIC  :: bessJY, bessJYN, JYBesAiry, expax, deb2, deb1
  CONTAINS
     SUBROUTINE bessJY(a, x, jax, jaxp, yax, yaxp, ierr)
     !--------------------------------------------------
     ! Computation of the Bessel functions Ja(x), Ya(x)
     ! and their derivatives for real positive orders (a) 
     ! and arguments (x).
     !--------------------------------------------------
     ! Inputs:
     !   a ,    order of the Bessel functions
     !   x,     argument of the Bessel functions.
     ! Outputs:
     !   jax,     function Ja(x)
     !   jaxp,    first order derivative of Ja(x)
     !   yax,     function Ya(x)
     !   yaxp,    first order derivative of Ya(x)
     !   ierr ,   error flag
     !            ierr=0, computation succesful
     !            ierr=1, overflow/underflow problems. The function values 
     !                    are set to zero.
     !            ierr=2, any of the arguments of the function is 
     !                    out of range. The function values  
     !                    are set to zero.
     ! ----------------------------------------------------------------------
     ! Authors:
     !  Amparo Gil    (U. Cantabria, Santander, Spain)
     !                 e-mail: amparo.gil@unican.es
     !  Javier Segura (U. Cantabria, Santander, Spain)
     !                 e-mail: javier.segura@unican.es
     !  Nico M. Temme (CWI, Amsterdam, The Netherlands)
     !                 e-mail: nico.temme@cwi.nl
     ! -------------------------------------------------------------
     USE Someconstants  
     IMPLICIT NONE
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: x
     REAL(r8), INTENT(OUT) :: jax, jaxp
     REAL(r8), INTENT(OUT) :: yax, yaxp
     INTEGER,  INTENT(OUT) :: ierr
     REAL(r8) :: xlim
     ierr=0
     IF ((x<0.0_r8).OR.(a<0.0_r8)) THEN
       ! Argument out of range
       ierr=2
       jax=0.0_r8; yax=0.0_r8
       jaxp=0.0_r8; yaxp=0.0_r8
    ENDIF
    IF (a/=0) THEN
      xlim=(2.0_r8*a)/exp(1.0_r8)*xpowy(dwarf,1.0_r8/a)
      IF (x<xlim) THEN
        ! Underflow/overflow problems
        ierr=1
        jax=0.0_r8; yax=0.0_r8
        jaxp=0.0_r8; yaxp=0.0_r8
      ENDIF
   ENDIF
     IF (ierr==0) THEN
     !  CALL JYBesAiry(a, x, jax, jaxp, yax, yaxp, ierr)
       IF (a>15.0_r8) THEN
         IF (a<1500.0_r8) THEN
           IF ((x>0.55_r8*a).AND.(x<2.0_r8*a)) THEN
              ! Airy-type asymptotic expansion
             CALL JYBesAiry(a, x, jax, jaxp, yax, yaxp, ierr)
           ELSE   
             CALL bessjax(a, x, jax, jaxp)
             CALL bessyax(a, x, yax, yaxp) 
           ENDIF
         ELSE
           IF ((x>0.55_r8*a).AND.(x<2.3_r8*a-36.41)) THEN
              ! Airy-type asymptotic expansion
             CALL JYBesAiry(a, x, jax, jaxp, yax, yaxp, ierr)
           ELSE   
             CALL bessjax(a, x, jax, jaxp)
             CALL bessyax(a, x, yax, yaxp) 
           ENDIF 
         ENDIF
       ELSE
         CALL bessjax(a, x, jax, jaxp)
         CALL bessyax(a, x, yax, yaxp)
       ENDIF
     ENDIF
     END SUBROUTINE bessJY 


      SUBROUTINE bessJYN(a, x, jax, jaxp, yax, yaxp, ierr)
     !--------------------------------------------------
     ! Computation of the Bessel functions Ja(x), Ya(x)
     ! and their derivatives for real positive orders (a) 
     ! and arguments (x).
     !--------------------------------------------------
     ! Inputs:
     !   a ,    order of the Bessel functions
     !   x,     argument of the Bessel functions.
     ! Outputs:
     !   jax,     function Ja(x)
     !   jaxp,    first order derivative of Ja(x)
     !   yax,     function Ya(x)
     !   yaxp,    first order derivative of Ya(x)
     !   ierr ,   error flag
     !            ierr=0, computation succesful
     !            ierr=1, overflow/underflow problems. The function values 
     !                    are set to zero.
     !            ierr=2, any of the arguments of the function is 
     !                    out of range. The function values  
     !                    are set to zero.
     ! ----------------------------------------------------------------------
     ! Authors:
     !  Amparo Gil    (U. Cantabria, Santander, Spain)
     !                 e-mail: amparo.gil@unican.es
     !  Javier Segura (U. Cantabria, Santander, Spain)
     !                 e-mail: javier.segura@unican.es
     !  Nico M. Temme (CWI, Amsterdam, The Netherlands)
     !                 e-mail: nico.temme@cwi.nl
     ! -------------------------------------------------------------
     USE Someconstants  
     USE gammaCHI  
     IMPLICIT NONE
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: x
     REAL(r8), INTENT(OUT) :: jax, jaxp
     REAL(r8), INTENT(OUT) :: yax, yaxp
     INTEGER,  INTENT(OUT) :: ierr
     REAL(r8) :: xlim
     REAL(r8) :: eps, alim, aa, c, r, s, jap, yap, ia, iad, lneps
     INTEGER :: n
     eps=1.e-12_r8
     aa=abs(a);
     lneps=-log(eps);
     ierr=0
     IF ((x<0.0_r8).OR.(a<0.0_r8)) THEN
       ! Argument out of range
       ierr=2
       jax=0.0_r8; yax=0.0_r8
       jaxp=0.0_r8; yaxp=0.0_r8
     ENDIF
     IF (a/=0) THEN
       xlim=(2.0_r8*a)/exp(1.0_r8)*xpowy(1.e-200_r8,1.0_r8/a)   
       IF (x<xlim) THEN
         ! Underflow/overflow problems
         ierr=1
         jax=0.0_r8; yax=0.0_r8
         jaxp=0.0_r8; yaxp=0.0_r8
       ENDIF
     ENDIF
     IF (ierr==0) THEN
       IF (a>300.0_r8) THEN
         alim=2.3_r8*a
       ELSE
         alim=2.0_r8*a
       ENDIF
       IF (a>15.0_r8) THEN          
         IF ((x>0.55_r8*a).AND.(x<alim)) THEN
           ! Airy-type asymptotic expansion
           CALL JYBesAiry(a, x, jax, jaxp, yax, yaxp, ierr)
         ELSE   
           IF ((x*x<4.0_r8*(1.0_r8 + aa)).AND.(x<11.0_r8)) THEN
              ! Power series
             ia=exp(a*(log(x/(2.0_r8*a))+1.0_r8))&
                  /(sqrt(2.0_r8*pi*a)*gamstar(a)); 
             r=x*x/4.0_r8;
             c=1.0_r8;
             s=a;
             n=1;
             jax=0.0_r8;
             jap=0.0_r8;
             eps=1.e-15_r8
             DO WHILE (abs(s)>eps)
               jax=jax+c;
               jap=jap+s;
               c=-c*r/(n*(a+n));
               s=(a+2.0_r8*n)*c;
               n=n+1
             ENDDO		
             jax=jax*ia;
             jap=jap*ia/x
             jaxp=jap
             CALL bessyax(a, x, yax, yaxp) 
           ELSEIF (x>(a/5.3_r8)**2+20.0_r8) THEN 
             CALL expax(a, x, jax, jaxp, yax, yaxp,ierr)
           ELSEIF ((a>5.0_r8+2.0_r8*lneps*(0.4343_r8+0.025_r8*lneps)&
                 /(2.0_r8-3.6_r8*x/a)).AND.(x<a/2.0_r8)) THEN
             CALL deb1(a, x, jax, jap, yax, yaxp)
           ELSEIF ((x>5.0_r8+lneps*(0.4343_r8+0.01_r8*lneps)/&
                 (2.0_r8-3.4_r8*a/x)).AND.(alim<x)) THEN
             CALL deb2(a, x, jax, jaxp, yax, yaxp);
           ELSE
             IF ((x<a).AND.(a>120.0_r8)) THEN
               CALL deb1(a, x, jax, jap, yax, yaxp)
             ELSE
               CALL jrec(a, x, jax, jaxp);
               CALL bessyax(a, x, yax, yaxp)
             ENDIF 
           ENDIF    
         ENDIF
       ELSE     
         IF ((x*x<4.0_r8*(1.0_r8 + aa)).AND.(x<11.0_r8)) THEN
           IF (a==0.0_r8) THEN
             ia=1.0_r8
           ELSE    
             ia=exp(a*(log(x/(2.0_r8*a))+1.0_r8))&
                /(sqrt(2.0_r8*pi*a)*gamstar(a));
           ENDIF
           r=x*x/4.0_r8;
           c=1.0_r8;
           s=a;
           n=1;
           jax=0.0_r8;
           jap=0.0_r8;
           eps=1.e-15_r8
           DO WHILE ((abs(s)>eps).OR.(n==1))
             jax=jax+c;
             jap=jap+s;
             c=-c*r/(n*(a+n));
             s=(a+2.0_r8*n)*c;
             n=n+1
           ENDDO				
           jax=jax*ia;
           jap=jap*ia/x
           jaxp=jap
           CALL bessyax(a, x, yax, yaxp)
         ELSEIF (x>(a/5.3_r8)**2+20.0_r8) THEN 
           CALL expax(a, x, jax, jaxp, yax, yaxp,ierr)
         ELSEIF ((a>5.0_r8+2.0_r8*lneps*(0.4343_r8+0.025_r8*lneps)&
                 /(2.0_r8-3.6_r8*x/a)).AND.(x<a/2.0_r8)) THEN
           CALL deb1(a, x, jax, jap, yax, yaxp)
         ELSEIF ((x>5.0_r8+lneps*(0.4343_r8+0.01_r8*lneps)/&
                 (2.0_r8-3.4_r8*a/x)).AND.(alim<x)) THEN
           CALL deb2(a, x, jax, jaxp, yax, yaxp);
         ELSE
           IF ((x<a).AND.(a>120.0_r8)) THEN
             CALL deb1(a, x, jax, jap, yax, yaxp)
           ELSE
             CALL jrec(a, x, jax, jaxp);
             CALL bessyax(a, x, yax, yaxp)
           ENDIF 
         ENDIF    
       ENDIF
     ENDIF
     END SUBROUTINE bessJYN 



     RECURSIVE SUBROUTINE  bessjax(a, x, jax, jap)
     ! Output: Ja(x) and derivative
     USE gammaCHI
     USE Someconstants
     IMPLICIT NONE
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: x
     REAL(r8), INTENT(OUT) :: jax
     REAL(r8), INTENT(OUT) :: jap
     REAL(r8) :: eps, aa, c, r, s, yax, yap, ia, iad, lneps
     INTEGER :: n
     eps=1.e-12_r8
     aa=abs(a);
     lneps=-log(eps);
     IF (a<-0.5_r8) THEN
       CALL bessjax(-a, x, ia, iad);
       CALL bessyax(-a, x, yax, yap);
       s=sin(a*pi);
       c=cos(a*pi);
       jax=c*ia+s*yax;
       jap=c*iad+s*yap
     ELSEIF (aa==0.5_r8) THEN
       s=sin(x);
       c=cos(x);
       r=sqrt(2.0_r8/(pi*x))
       IF (a<0) THEN
         jax=c*r;
         jap=a*jax/x-r*s
       ELSE
         jax=r*s;
         jap=r*c-a*jax/x
       ENDIF  
     ELSEIF ((x*x<4.0_r8*(1.0_r8 + aa))) THEN
       ia=exp(a*(log(x/(2.0_r8*a))+1.0_r8))&
          /(sqrt(2.0_r8*pi*a)*gamstar(a)); 
       r=x*x/4.0_r8;
       c=1.0_r8;
       s=a;
       n=1;
       jax=0.0_r8;
       jap=0.0_r8;
       eps=1.e-15_r8
       DO WHILE (abs(s)>eps)
         jax=jax+c;
         jap=jap+s;
         c=-c*r/(n*(a+n));
         s=(a+2.0_r8*n)*c;
         n=n+1
       ENDDO				
       jax=jax*ia;
       jap=jap*ia/x    
     ELSEIF ((a>5.0_r8+2.0_r8*lneps*(0.4343_r8+0.025_r8*lneps)&
         /(2.0_r8-3.6_r8*x/a)).AND.(x<a/2.0_r8)) THEN
       CALL deb1(a, x, jax, jap, ia, iad)
     ELSEIF ((x>5.0_r8+lneps*(0.4343_r8+0.01_r8*lneps)/&
            (2.0_r8-3.4_r8*a/x)).AND.(a<x/2.0_r8)) THEN
       CALL deb2(a, x, jax, jap, ia, iad);
     ELSE
       CALL jrec(a, x, jax, jap)
     ENDIF
     END SUBROUTINE  bessjax



     RECURSIVE SUBROUTINE  bessyax(a, x, yax, yap)
     USE Someconstants
     USE gammaCHI  
     IMPLICIT NONE
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: x
     REAL(r8), INTENT(OUT) :: yax
     REAL(r8), INTENT(OUT) :: yap
     REAL(r8) :: eps, ya1, b, c, d, e, f, g, h, p, q, r, s, t
     INTEGER :: n, na
     LOGICAL :: rec, rev
     eps=1.0e-15_r8
     IF (a <-0.5_r8) THEN	  
       CALL bessjax(-a, x, f, g);
       CALL bessyax(-a, x, p, q);
       c=cos(a*pi);
       s=sin(a*pi);
       yax=-s*f+c*p;
       yap=-s*g+c*q
     ELSEIF (abs(a)==0.5_r8) THEN
       s=sin(x);
       c=cos(x);
       t=sqrt(2.0_r8/(pi*x));
       IF (a < 0) THEN
         yax=s*t;
         yap=a*yax/x+t*c		
       ELSE
         yax=-c*t;
         yap=s*t-a*yax/x
       ENDIF
     ELSE
       na=int(a+0.5_r8);
       r=a-na;
       IF (x<11.5_r8) THEN
         b=x/2.0_r8;
         d=-log(b);
         e=r*d;
         c=r*pi;
         IF (abs(c)<1.0e-5_r8) THEN
           c=1.0_r8 + c*c/6.0_r8
         ELSE
           c=c/sin(c);
         ENDIF   
         c=c/pi;
         IF (abs(e)<1.0e-5_r8) THEN
	   q=e*e;
           s=1.0_r8+q/6.0_r8;
           q=1.0_r8+q/2.0_r8
         ELSE
	   CALL hypfun(e, p, q);
           s=p/e
         ENDIF
         e=exp(e);
         g=e*recipgam(-r, p, t);
         e=r*r;
         f=2.0_r8*c*(-p*q+(1.0_r8+e*t)*s*d);
         p=g*c;
         q=1.0_r8/(g*pi);
         t=r*pi/2.0_r8;
         IF (abs(t)<1.0e-5_r8) THEN
           s=1.0_r8-t*t/6.0_r8
         ELSE
           s=sin(t)/t;
         ENDIF   
         s=pi*t*s*s;
         c=1.0_r8;
         d=-b*b;
         yax=f+s*q;
         ya1=p;
         n=1;
         h=ya1
         g=yax
         DO WHILE ((abs(h/ya1)+ abs(g/yax))> eps)
           f=(f*n+p+q)/(n*n-e);
           c=c*d/n;
           p=p/(n-r);
           q=q/(n+r);
           g=c*(f+s*q);
           h=c*p-n* g;
           yax=yax+g;
           ya1=ya1+h;
           n=n+1
         ENDDO   
         f=-yax;
         g=-ya1/b
       ELSE
         b=x-pi*(r+0.5_r8)/2.0_r8;
         c=cos(b);
         s=sin(b);
         d=sqrt(2.0_r8/(pi*x));
         CALL besspqax(r, x, p, q, b, h);
         f=d*(p*s+q*c);
         g=d*(b*c-h*s);
         g=r*f/x-g
       ENDIF
       b=2.0_r8/x;
       DO n=1,na				
         h=b*(r+n)*g-f;
         f=g;
         g=h
       ENDDO
       yax=f;
       yap=a*f/x-g
     ENDIF 
     END SUBROUTINE bessyax

     SUBROUTINE besspqax(a, x, pa, qa, ra, sa) 
     !computes P, Q, R, S, of formulas 9.25- 9.2.16 of A&S}
     USE Someconstants
     USE gammaCHI  
     IMPLICIT NONE
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: x
     REAL(r8), INTENT(OUT) :: pa
     REAL(r8), INTENT(OUT) :: qa
     REAL(r8), INTENT(OUT) :: ra
     REAL(r8), INTENT(OUT) :: sa
     REAL(r8) :: eps, b, c, d, e, f, g, h, p, q, r, s, t, pa1, qa1
     INTEGER :: n, na
     eps=1.0e-15_r8
     r=abs(a);
     na=int(a);
     r=r-na;
     IF (r==0.5_r8) THEN
       pa=1.0_r8;
       pa1=1.0_r8;
       qa=0.0_r8;
       qa1=1.0_r8/x
     ELSEIF (x >= 5) THEN
       c=0.25_r8-r*r;
       b=x+x;
       f=1.0_r8;
       g=1.0_r8;
       p=1.0_r8;
       q=0.0_r8;
       n=startpqbes(x, eps);
       DO WHILE (n > 0)
         t=(n+1.0_r8)*(2.0_r8-p)-2.0_r8;
         s=b+(n+1)*q;
         d=(n-1.0_r8+c/n)/(s*s+t*t);
         p=d*t;
         q=d*s;
         e=f;
         f=p*(e+1.0_r8)-g*q;
         g=q*(e+1.0_r8)+g*p;
         n=n-1.0_r8;
       ENDDO   
       f=1.0_r8+f;
       d=f*f+g*g;
       pa=f/d;
       qa=-g/d;
       d=r+0.5_r8-p;
       q=q+x;
       pa1=(pa*q-qa*d)/x;
       qa1=(qa*q+pa*d)/x
     ELSE 
       e=sqrt(pi*x/2.0_r8);
       t=x-pi*(r/2.0_r8+0.25_r8);
       c=cos(t);
       s=sin(t);
       CALL bessyax(r, x, p, q);
       CALL bessjax(r, x, f, g);
       d=r/x;
       q=d*p-q;
       g=d*f-g;
       pa=e*(s*p+c*f);
       qa=e*(c*p-s*f);
       pa1=e*(s*g-c*q);
       qa1=e*(c*g+s*q);
     ENDIF    
     t=2.0_r8/x;
     b=(r+1.0_r8)*t;
     DO n=1,na 
       c=pa-qa1*b;
       s=qa+pa1*b;
       pa=pa1;
       pa1=c;
       qa=qa1;
       qa1=s;
       b=b+t
     ENDDO
     ra=abs(a)*qa/x+pa1;
     sa=-abs(a)*pa/x+qa1;
     END SUBROUTINE besspqax 

     SUBROUTINE deb1(a,x,jax,jap,yax,yap)
     ! Debye-type expansion a > x; 9.3.7 and 9.3.8 of A&S
     ! See Luke (1975) for this modification of 
     !   computing coefficients
     ! Out: Ja(x), Ya(x) and their derivatives}
     USE Someconstants
     IMPLICIT NONE
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: x
     REAL(r8), INTENT(OUT) :: jax
     REAL(r8), INTENT(OUT) :: jap
     REAL(r8), INTENT(OUT) :: yax
     REAL(r8), INTENT(OUT) :: yap
     REAL(r8) :: eps, y, eta, t, u, v, w, z, zz, wr, t2, ck, c1, c2, c3
     INTEGER :: jj, jjj, k, m, n, sig
     eps=1.e-14_r8
     y=x/a;
     t2=1.0_r8/(1.0_r8-y*y);
     t=sqrt(t2);
     eta=log((1.0_r8+t)/(t*y))-1.0_r8/t;
     c1=1.0_r8;
     c2=0.0_r8;
     c3=0.0_r8;
     jj=0;
     jjj=0;
     u=t/(8.0_r8*a);
     w=4.0_r8*u;
     z=sqrt(t/(2.0_r8*pi*a));
     zz=z/(y*t);
     k=1;
     m=1;
     n=-3;
     jax=1.0_r8;
     jap=1.0_r8;
     sig=1;
     yax=1.0_r8;
     yap=1.0_r8;
     jj=0
     jjj=0
     DO WHILE ((jj/=30).AND.(jjj/=400))
       v=m*c1-n*c3*t2;
       ck=m*u*v/k;
       IF (abs(c3)<abs(ck)) THEN
         jj=jj+10
       ELSE
         jj=0;
       ENDIF
       IF (abs(ck)<eps) THEN
         jjj=jjj+100
       ELSE
         jjj=0;
         c3=c2;
         c2=c1;
         c1=ck;
         k=k+1;
         m=m+2;
         n=n+2;
         sig=-sig;
         v=ck-w*v;
         jax=jax+ck;
         jap=jap+v;
         yax=yax+sig*ck;
         yap=yap+sig*v;
         wr=abs((jax*yap+jap*yax)/2.0_r8-1.0_r8)
       ENDIF
     ENDDO
     v=exp(a*eta);
     jax=z*jax/v;
     yax=-2.0_r8*z*yax*v;
     jap=zz*jap/v;
     yap=2.0_r8*zz*yap*v;
     END SUBROUTINE deb1

     SUBROUTINE deb2(a,x,jax,jap,yax,yap)
     ! Debye-type expansion x > a; 9.3.7 and 9.3.16 of A&S
     ! See Luke (1975) for this modification of computing coefficients}
     ! Out: Ja(x), Ya(x) and their derivatives}
     USE Someconstants
     IMPLICIT NONE
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: x
     REAL(r8), INTENT(OUT) :: jax
     REAL(r8), INTENT(OUT) :: jap
     REAL(r8), INTENT(OUT) :: yax
     REAL(r8), INTENT(OUT) :: yap
     REAL(r8) :: eps, a2, x2, w2, ww, chi, c, k1, k2, p, q, &
                 r, s, u, v, w, wr, ek, e1, e2, e3, fk
     INTEGER :: jj, jjj, k, m, n, sig, t
     eps=1.e-14_r8
     a2=a*a;
     x2=x*x;
     w2= a2/(x2-a2);
     ww=sqrt(x2-a2);
     chi=ww+a*(atan(a/ww)-pi/2.0_r8)-pi/4.0_r8;
     c=sqrt(2.0_r8/pi);
     k2=sqrt(ww);
     k1=c/k2;
     k2=c*k2/x;
     e1= 1.0_r8;
     e2=0.0_r8;
     e3=0.0_r8;
     jj=0;
     jjj=0;
     u=1.0_r8/(8.0_r8*ww);
     w=4.0_r8*u;
     k=1;
     m=1;
     n=-3;
     p=1.0_r8;
     q=0.0_r8;
     r=1;
     s=0;
     sig=1;
     t=1.0_r8;
     jj=0
     jjj=0
     DO WHILE ((jj/=30).AND.(jjj/=400))
       v=n*e3*w2-m*e1;
       ek=m*u*v/k;
       IF (abs(e3)< abs(ek)) THEN
         jj=jj+10
       ELSE
         jj=0;
       ENDIF   
       IF (abs(ek)<eps) THEN
         jjj=jjj+100
       ELSE
         jjj=0;
       ENDIF   
       e3=e2;
       e2=e1;
       e1=ek;
       fk=ek-w*v;
       IF (sig==1) THEN
         q=q+t*ek;
         s=s+t*fk;
         t=-t
       ELSE
         p=p+t*ek;
         r=r+t*fk
       ENDIF
       k=k+1;
       m=m+2;
       n=n+2;
       sig=-sig;
       wr=abs(p*r+q*s-1.0_r8)
     ENDDO   
     c=cos(chi);
     u=sin(chi);
     jax=k1*(p*c-q*u);
     yax=k1*(p*u+q*c);
     jap=-k2*(r*u+s*c);
     yap=k2*(r*c-s*u);
     END SUBROUTINE deb2  	

     SUBROUTINE expax(a,x,jax,jap,yax,yap,ierr)
     ! Asymptotic expansions, x large  
     USE Someconstants
     IMPLICIT NONE
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: x
     REAL(r8), INTENT(OUT):: jax, jap, yax, yap
     REAL(r8) :: eps, w, cosw, sinw, fact, fnus, x2, err1, err1p, &
                 y1, y1p, a2, x2k, x2kp, acof, acofd
     INTEGER, INTENT(OUT):: ierr
     INTEGER :: k, l, m
     m=0
     eps=epss
     w=x-a*pihalf-pikwart
     x2=x*x
     cosw=cos(w)
     sinw=sin(w)
     fact=sqrt(1.0_r8/(pihalf*x))
     fnus=4.0_r8*a*a
     err1=1.0_r8; 
     err1p=1.0_r8;
     y1=0.0_r8; 
     y1p=0.0_r8;
     k=0
     a2=1.0_r8
     x2k=1.0_r8
     x2kp=x
     DO WHILE (((err1 > eps).OR.(err1p >eps)).AND.(m<500))
       acof=a2/x2k;
       k=k+1
       l=2*k
       a2=(fnus-(l-1.0_r8)*(l-1.0_r8))*a2/(8.0_r8*k)   
       acofd=a2/x2kp;
       k=k+1;
       l=2*k
       a2=(fnus-(l-1.0_r8)*(l-1.0_r8))*a2/(8.0_r8*k)
       x2k=-x2*x2k
       x2kp=-x2*x2kp
       y1=y1+acof; 
       y1p=y1p+acofd    
       err1=abs(acof/y1); 
       err1p=abs(acofd/y1p);
       IF (abs(a2)>giant) m=500
     ENDDO     
     IF (m==500) THEN
       ierr=1
       jax=0.0_r8; yax=0.0_r8;
       jap=0.0_r8; yap=0.0_r8
     ELSE
       jax=fact*(cosw*y1-sinw*y1p)
       yax=fact*(sinw*y1+cosw*y1p)
     ENDIF
     IF (ierr==0) THEN
       err1=1.0_r8; 
       err1p=1.0_r8;
       y1=1.0_r8 
       a2=(fnus+3.0_r8)/8.0_r8
       x2k=1.0_r8
       x2kp=x
       y1p=a2/x2kp
       a2=(fnus-1.0_r8)*(fnus+15.0_r8)/128.0_r8
       x2k=-x2*x2k
       x2kp=-x2*x2kp
       k=2
       DO WHILE (((err1 > eps).OR.(err1p >eps)).AND.(m<500))
         acof=a2/x2k;
         k=k+1
         l=2*k
         a2=(fnus-(l-3.0_r8)*(l-3.0_r8))*a2*(fnus+4*k*k-1.0_r8)&
              /((8.0_r8*k)*(fnus+4*(k-1)*(k-1)-1.0_r8))
         acofd=a2/x2kp;
         k=k+1
         l=2*k
         a2=(fnus-(l-3.0_r8)*(l-3.0_r8))*a2*(fnus+4*k*k-1.0_r8)&
              /((8.0_r8*k)*(fnus+4*(k-1)*(k-1)-1.0_r8))
         x2k=-x2*x2k
         x2kp=-x2*x2kp
         y1=y1+acof; 
         y1p=y1p+acofd    
         err1=abs(acof/y1); 
         err1p=abs(acofd/y1p);
         IF (abs(a2)>giant) m=500
       ENDDO
       IF (m==500) THEN
         ierr=1
         jax=0.0_r8; yax=0.0_r8;
         jap=0.0_r8; yap=0.0_r8
       ELSE
         jap=-fact*(sinw*y1+cosw*y1p)
         yap=fact*(cosw*y1-sinw*y1p)
       ENDIF
     ENDIF   
     END SUBROUTINE expax


     SUBROUTINE jrec(a, x, jax, jap)
     !Miller, with the standard recurrence relation}
     !Out: jax =Ja(x), jap its derivative}    
     USE Someconstants
     USE GammaCHI  
     IMPLICIT NONE
     REAL(r8), INTENT(IN) :: a
     REAL(r8), INTENT(IN) :: x
     REAL(r8), INTENT(OUT) :: jax
     REAL(r8), INTENT(OUT) :: jap
     REAL(r8) :: eps, j0, j1, j2, am, ap, x2, r, s, &
                 t, l, la, g1, g2, g3
     INTEGER :: k, m, nu, nu2, n, nn, na
     LOGICAL :: notsimple
     eps=1.0e-12_r8
     IF (x <= 0) THEN
       IF (a == 0.0_r8) THEN
         jax=1.0_r8					
       ELSE
         jax=0.0_r8
       ENDIF					
       jap=0.0_r8	
     ELSE
     ! Function startijbes
       na=int(a) 
       nu=startijbes(x,na+1,0,eps);
       nu2=nu+nu;
       x2=2.0_r8/x;
       k=1;
       ap=a+1.0_r8;
       am=a-1.0_r8;
       ! Logical variable
       notsimple=((a/=0.0_r8).AND.(a/=1.0_r8));
       IF (a==0.0_r8) THEN
         la=2
       ELSEIF (a==1.0_r8) THEN
         la=1
       ELSE
         g1=gamstar(nu+a);
         g2=gamstar(nu*1.0_r8);
         g3=gamstar(a);
         la=g1*exp(a*log(1.0_r8+nu/a)+nu*log(1.0_r8+a/nu))&
           /(g2*g3*sqrt(2*pi*(nu+a)*nu*a));
       ENDIF
       j2=0.0_r8;
       j1=1.0_r8;
       n= nu2;
       m= nu;
       s=0.0_r8;
       DO WHILE (n>0)
         j0=-j2+x2*(ap+n)*j1;
         IF (k==1) THEN
	   IF (a==0.0_r8) THEN
             s=s+2*j0
           ELSE
             s=s+la*(a+n)*j0;
           ENDIF
	   IF (m>1) THEN
	     IF (notsimple) THEN
	       la=la*m/(m+am)
             END IF;
           ENDIF
           m=m-1;
         ENDIF
         j2=j1;
         j1=j0;
         k=-k;
         n=n-1
       ENDDO
       j0=-j2+x2*ap*j1;
       s=s+j0;
       IF (a==0.0_r8) THEN
         s=1.0_r8/s
       ELSEIF (a==1.0_r8) THEN
         s=1.0_r8/(s*x2)
       ELSE   
         s=exp(a*(1.0_r8-log(a*x2)))/(s*g3*sqrt(2.0_r8*pi*a));
       ENDIF
       jax=s*j0;
       jap=s*(a/x*j0-j1)
     ENDIF		
     END SUBROUTINE jrec	

     FUNCTION startpqbes(x,eps)
     USE Someconstants
     USE GammaCHI  
     IMPLICIT NONE
     REAL(r8) :: x,eps
     REAL(r8) :: startpqbes, del, a, b, c, p, q, r, s, t
     IF (eps<mactol) THEN
       del=-log(mactol/2.0_r8)
     ELSE
       del=-log(eps/2.0_r8);
     ENDIF
     t=del/(2.0_r8*x);
     IF (t<0.5_r8) THEN
       s=1.0_r8/(2.0_r8*t);
       c=sqrt(1.0_r8+s*s);
       a=log(s+c)+s*c/(1.0_r8+c*c);
       a=a*(1.0_r8+7.16_r8*t)/(1.0_r8+5.33_r8*t)
     ELSE
       q=1.0_r8/t*t;
       a=(1.0_r8+q*q*(-7.0_r8/360.0_r8+4.1e-3_r8*q))/t
     ENDIF
     b=0.0_r8
     DO WHILE (abs(b/a-1.0_r8)>1.0e-2_r8)
       b=a;
       CALL hypfun(a,s,c);
       a=a+(a*c*s+s*s*(1.0_r8-2.0_r8*t*s))/(a*(1.0_r8+c*c));
     ENDDO  
     CALL hypfun(a, s, c);
     startpqbes=1.0_r8+int(x*c/s*s)
     END FUNCTION startpqbes

     SUBROUTINE phizeta(zeta, z, s, kmax)
     ! Defined in our SIAM book (8.67)
     ! Inputs: zeta, z
     ! Outputs: s, kmax
     USE Someconstants
     IMPLICIT NONE
     REAL(r8) :: zeta, z, s
     REAL(r8) :: d, phik(0:24), t, zetak;
     INTEGER :: k, kmax
     phik(0)=1.2599210498948731648_r8;
     phik(1)=.20000000000000000000_r8;
     phik(2)=0.20409442096733993248e-1_r8;
     phik(3)=-0.35597769346236098942e-2_r8;
     phik(4)=-0.19683982683982683983e-2_r8;
     phik(5)=-0.20673367061132873279e-3_r8;
     phik(6)=0.11681173009598189462e-3_r8;
     phik(7)=0.53095696007340665204e-4_r8;
     phik(8)=0.42673899863261117719e-5_r8;
     phik(9)=-0.41528636693801060913e-5_r8;
     phik(10)=-0.17568372831164108478e-5_r8;
     phik(11)=-0.11402786896104923682e-6_r8;
     phik(12)=0.15533576025831417179e-6_r8;
     phik(13)=0.63395269099583076477e-7_r8;
     phik(14)= 0.34894701838555374020e-8_r8;
     phik(15)=-0.60044875533414815786e-8_r8;
     phik(16)=-0.23993332052814288615e-8_r8;
     phik(17)=-0.11582676745605076908e-9_r8;
     phik(18)=0.23754898857892364962e-9_r8;
     phik(19)=0.93607463563942867781e-10_r8;
     phik(20)=0.40542001500194274074e-11_r8;
     phik(21)=-0.95614572050793794738e-11_r8;
     phik(22)= -0.37301511168150182683e-11_r8;
     phik(23)= -0.14727062811805787322e-12_r8;
     phik(24)= 0.39000756737548022749e-12_r8;
     s= phik(0); 
     t=s;  
     zetak= 1.0_r8; 
     k= 1;
     DO WHILE ((abs(t) > 1.0e-16_r8).AND.(k < 25))
       zetak= zeta*zetak;
       t=phik(k)*zetak;
       s=s+t;
       k= k+1
     ENDDO
     kmax= k+1;
     END SUBROUTINE phizeta

     FUNCTION chizeta(zeta)
     ! Defined in our SIAM book (8.73)
     USE Someconstants
     IMPLICIT NONE
     REAL(r8) :: zeta, chizeta
     REAL(r8) :: d, chik(0:24), s, t, zetak;
     INTEGER :: k
     chik(0)=.15874010519681994748_r8;
     chik(1)= 0.71995488565421323703e-2_r8;
     chik(2)= -0.12190476190476190476e-1_r8;
     chik(3)= -0.39822790851437230115e-2_r8;
     chik(4)= 0.27754184922651308731e-3_r8;
     chik(5)= 0.57958491395226089104e-3_r8;
     chik(6)= 0.15466235789565560764e-3_r8;
     chik(7)= -0.21638396018033400537e-4_r8;
     chik(8)= -0.27028766531753854913e-4_r8;
     chik(9)= -0.65331520523895438118e-5_r8;
     chik(10)= 0.11834212951681698844e-5_r8;
     chik(11)= 0.12449403861697797972e-5_r8;
     chik(12)= 0.28434185951562648920e-6_r8;
     chik(13)= -0.58918145279531299099e-7_r8;
     chik(14)= -0.56913494700904364692e-7_r8;
     chik(15)= -0.12536139854637631961e-7_r8;
     chik(16)= 0.28216007891494136557e-8_r8;
     chik(17)= 0.25901270610604348767e-8_r8;
     chik(18)= 0.55628394372865850242e-9_r8;
     chik(19)= -0.13250218032852282737e-9_r8;
     chik(20)= -0.11754702052918243107e-9_r8;
     chik(21)= -0.24776231260359780515e-10_r8;
     chik(22)= 0.61525731123646300928e-11_r8;
     chik(23)= 0.53249498993041565381e-11_r8;
     chik(24)= 0.11060987152463791205e-11_r8;
     s= chik(0); 
     t= s; 
     d= 1; 
     zetak= 1.0_r8; 
     k= 1;
     DO WHILE ((abs(t) > 1.0e-16_r8).AND.(k < 25))
       zetak= zeta*zetak; 
       t= chik(k)*zetak;  
       s= s + t;  
       k= k+1
     ENDDO
     chizeta=s
     END FUNCTION chizeta

     FUNCTION zetaz(z)
     ! Defined in (8.68)
     USE Someconstants
     IMPLICIT NONE
     REAL(r8) :: z  
     REAL(r8) :: zetaz, zeta, zetak(1:24), w, wk, s, t, d;
     INTEGER :: k
     zetak(1)=-1.0_r8;
     zetak(2)=0.300000000000000000000000000000_r8
     zetak(3)=-0.182857142857142857142857142857_r8
     zetak(4)=0.131682539682539682539682539683_r8
     zetak(5)=-0.102636487322201607915893630179_r8
     zetak(6)=0.0840417027417027417027417027417_r8
     zetak(7)=-0.0709607816028768409720790673172_r8
     zetak(8)=0.0611325306917572223694672674265_r8
     zetak(9)=-0.0533701347889438816199130714547_r8
     zetak(10)=0.0469919788555835658740312665375_r8
     zetak(11)=-0.0415806453450275544520787454077_r8
     zetak(12)=0.0368672753779680636097084350267_r8
     zetak(13)=-0.0326717513334618341955351074614_r8
     zetak(14)=0.0288696438349622209737935735201_r8
     zetak(15)=-0.0253729385600982588211253000563_r8
     zetak(16)=0.0221182862894696069441279175389_r8
     zetak(17)=-0.0190595626009876389331767879833_r8
     zetak(18)=0.0161629994691626861897781420917_r8
     zetak(19)=-0.0134039065323106692770387708064_r8
     zetak(20)=0.0107644051661864058879958722483_r8
     zetak(21)=-0.00823182508269234831534366561927_r8
     zetak(22)=0.00579754441493412939502240278547_r8
     zetak(23)=-0.00345613270412384967447374617551_r8
     zetak(24)=0.00120470443781096787386723520378_r8
     w=z-1.0_r8;
     IF (z==1.0_r8) THEN
       zeta=0.0_r8 
     ELSEIF (abs(w)<0.01_r8) THEN
       s=0.0_r8; 
       t=1.0_r8; 
       d=1.0_r8;  
       wk=1.0_r8; 
       k= 1;
       DO WHILE ((abs(t)>1.0e-16_r8).AND.(k < 25)) 
         wk=w*wk;
         t=zetak(k)*wk;
         s=s+t;
         k=k+1
       ENDDO
       zeta=twoexp13*s
     ELSEIF (z<1.0_r8) THEN 
       w=sqrt((1.0_r8-z)*(1.0_r8+z)); 
       zeta=xpowy(1.5_r8*(log((1.0_r8+w)/z)-w),twothird)
     ELSE
       w=sqrt((z-1.0_r8)*(z+1.0_r8)); 
       zeta=-xpowy(1.5_r8*(w-acos(1.0_r8/z)),twothird)
     ENDIF
     zetaz=zeta
     END FUNCTION zetaz

     SUBROUTINE JYBesAiry(a,x,jax, jaxp, yax, yaxp,ierr)
     ! Airy-type expansion for Bessel
     ! Output variables jax, jaxp, yax, yaxp are 
     ! the Besssel functions and derivatives
     USE Someconstants  
     USE AiryFunction
     IMPLICIT NONE
     REAL(r8) :: a,x
     REAL(r8) :: jax, jaxp, yax, yaxp
     REAL(r8) :: z, zeta, phiz, chi
     REAL(r8) :: w, wlim, sf, sg, sfp, sgp
     REAL(r8) :: a13, a23, a43, az, bz, cz, dz, &
                 ai, bi, aip, bip, phih
     INTEGER :: kmax, ierr
     ierr=0
     z=x/a;
     zeta=zetaz(z);
     a13=xpowy(a,1.0_r8/3.0_r8)
     a23= a13*a13; 
     a43=a*a13;
     w=a23*zeta;
     wlim=xpowy((-1.5_r8*log(dwarf)),2.0_r8/3.0_r8)
     IF (w>wlim) THEN
       ! Under/overflow problems
       ierr=1
       jax=0.0_r8
       jaxp=0.0_r8
       yax=0.0_r8 
       yaxp=0.0_r8 
     ELSE
       CALL phizeta(zeta, z, phiz, kmax);
       chi=chizeta(zeta);
       IF (kmax<5) THEN
         kmax= 5 
       ENDIF
       CALL recur(a,zeta,kmax,sf,sg,sfp,sgp);
       w= 1.0_r8/(a*a);
       az=sf; bz=sg;
       cz=chi*sf+sfp+zeta*sg;
       dz=sf+w*(chi*sg+sgp);
       w=a23*zeta;
       CALL aibi(w, ai, bi)
       CALL aibip(w, aip, bip)
       !Now compute the output Bessel functions J, Y and derivatives
       jax=phiz/a13*(ai*az+aip*bz/a43);
       yax=-phiz/a13*(bi*az+bip*bz/a43);
       phih= 2.0_r8/(z*phiz);
       jaxp=-phih*(ai*cz/a43+aip*dz/a23);
       yaxp=phih*(bi*cz/a43+bip*dz/a23);
     ENDIF
     END SUBROUTINE JYBesAiry

     SUBROUTINE recur(nu,zeta,kmax,sf,sg,sfp,sgp);
     USE Someconstants  
     IMPLICIT NONE
     REAL(r8) :: nu,zeta
     REAL(r8) :: sf,sg,sfp,sgp
     REAL(r8) :: c, f0(0:30), f1(0:30), g0(0:30), g1(0:30), &
                 s, mu, psik(0:24), zetk, p, q
     INTEGER :: j, k, iter, kmax
     psik(0)= 0.17998872141355330925e-1_r8;
     psik(1)= 0.26666666666666666667e-1_r8;
     psik(2)= 0.81284358134178674405e-2_r8;
     psik(3)= -0.25499940565393672258e-2_r8;
     psik(4)= -0.27185440409930205849e-2_r8;
     psik(5)= -0.64287929800235579817e-3_r8;
     psik(6)= 0.21800832196221116516e-3_r8;
     psik(7)= 0.19524603774687808301e-3_r8;
     psik(8)= 0.41595772741004484377e-4_r8;
     psik(9)= -0.14680083795856540042e-4_r8;
     psik(10)= -0.12159606341232918219e-4_r8;
     psik(11)= -0.24579927092988893928e-5_r8;
     psik(12)= 0.88115019306698039283e-6_r8;
     psik(13)= 0.70074891930340538162e-6_r8;
     psik(14)= 0.13746831574149545326e-6_r8;
     psik(15)= -0.49561583422955524129e-7_r8;
     psik(16)= -0.38479724126131225767e-7_r8;
     psik(17)= -0.74072962899962262038e-8_r8;
     psik(18)= 0.26753958521879461353e-8_r8;
     psik(19)= 0.20450384057905768059e-8_r8;
     psik(20)= 0.38862478591128364861e-9_r8;
     psik(21)= -0.14040517140485884664e-9_r8;
     psik(22)= -0.10616024263813946569e-9_r8;
     psik(23)= -0.19986725826057819121e-10_r8;
     psik(24)= 0.72186876857324009442e-11_r8; 
     c=twoexp13;
     DO k=0,kmax+2  
       f1(k)= 0; 
       g1(k)=psik(k)/(2.0_r8*k+1.0_r8) 
     ENDDO
     mu=2.0_r8*nu*nu;
     s=2.0_r8/mu;
     f1(0)=1.0_r8+s*(-1.0_r8/225.0_r8+s*(6.93735541e-4_r8&
          -s*3.5421197e-4_r8));
     g1(0)=c*(1.0_r8/70.0_r8+s*(-1213.0_r8/1023750.0_r8&
          +s*(4.3829180944898811e-4_r8-s*(3.7670439477105e-4_r8+&
          s*5.84533e-4_r8))));
     DO iter=1,3 
       DO k=0,kmax 
         p=0; q=0;
         DO j=0,k 
           p=p + f1(j)*psik(k-j); 
           q=q + g1(j)*psik(k-j)
         ENDDO
         f0(k)= p; g0(k)= q
       ENDDO
       DO k=kmax-1,1,-1 
         f1(k)=(g0(k-1)/k-(k+1)*g1(k+1))/mu;
         g1(k-1)=(f0(k-1)-k*(k+1)*f1(k+1))/&
                 (2.0_r8*k-1.0_r8)
       ENDDO
     ENDDO 
     sf=f1(0);sg=g1(0);      
     sfp=0.0_r8;sgp=0.0_r8;  
     zetk=1.0_r8; p=1.0_r8; q=1.0_r8; k= 1;
     DO WHILE ((k<kmax).AND.(abs(p)+abs(q)>1.0e-16_r8)) 
       p=f1(k)*zetk;
       q=g1(k)*zetk;
       sfp=sfp+k*p;
       sgp=sgp+k*q;
       zetk=zeta*zetk;
       sf=sf+p*zeta; 
       sg=sg+q*zeta;  
       k= k+1
     ENDDO
     END SUBROUTINE recur

     FUNCTION xpowy(x,y)
     IMPLICIT NONE
     REAL(r8) :: x,y,xpowy
     xpowy=x**y
     END FUNCTION xpowy
   END MODULE BesselJY




