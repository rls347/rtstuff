      SUBROUTINE PLEG(ANGLEO,GMU,M,NEXP,NZEN,YPLEG)
C----------------------------------------------------------------------C
C      Computation of renormalized associated Legendre polynomials     C  
C      of the form:                                                    C
C                                                                      C
C                m                     1/2 m                           C
C               Y (u) = [(l-m)!/(l+m)!]   P (u)                        C
C                l                         l                           C
C                                                                      C
C         where,                                                       C
C                 m                                                    C
C                P (u) = Associated Legendre polynomial of             C
C                 l      order l and degree m                          C
C                                                                      C
C                    u = Cosine of zenith angle                        C
C                                                                      C
C     Reference:                                                       C
C                                                                      C
C             Dave, J. V.,and Armstrong, B. H., 1970: Computations     C
C                  of High-order Associated Legendre Polynomials,      C
C                  JQSRT, 10, 557-562, 1970                            C
C----------------------------------------------------------------------C
C                 I N P U T    V A R I A B L E S :                     C
C----------------------------------------------------------------------C
C         ANGLEO   :    Cosine of solar zenith angle                   C
C         GMU      :    Gaussian quadrature points                     C
C         M        :    Index for degree of Legendre polynomial        C
C         NEXP     :    Tot. no. of polynomial expansion terms         C
C         NZEN     :    Tot. no. of quadrature points                  C
C----------------------------------------------------------------------C
C                O U T P U T    V A R I A B L E S :                    C
C----------------------------------------------------------------------C
C         YPLEG    :    Renormalized associated Legendre polynomials   C
C----------------------------------------------------------------------C
      REAL  GMU(NZEN),LM1,LM2,MU,YPLEG(NEXP,NZEN+2)
      SAVE  CM1, DM1
C
      IF (M .EQ. 0) THEN
        CM1 = 1.
        DM1 = 1.
      ELSE
        CM2 = -SQRT((2.*M-1.)/(2.*M))*CM1
        DM2 = -SQRT((2.*M+1.)/(2.*M))*DM1
      END IF
      DO 20,I=1,NZEN+2
        IF (I .LE. NZEN) MU = GMU(I)
        IF (I .EQ. NZEN+1) MU = ANGLEO
        IF (I .EQ. NZEN+2) MU = 0.
C----------------------------------------------------------------------C
C                  Compute initial values                              C
C----------------------------------------------------------------------C
        IF (M .EQ. 0) THEN
          YPLEG(1,I) = 1.
          YPLEG(2,I) = MU       
        ELSE  
          YPLEG(M+1,I) = CM2*(1-MU**2)**(M/2.)
          YPLEG(M+2,I) = DM2*MU*(1-MU**2)**(M/2.)
        END IF
        DO 10,IL=1,NEXP-2 
          L = IL-1
          IF (IL .LT. M+1) THEN
            YPLEG(IL,I) = 0.
            GO TO 10
          ELSE
C----------------------------------------------------------------------C
C          Use varying-degree recurrence formula to compute            C
C          successive values                                           C
C----------------------------------------------------------------------C
            LM2 = (L-M+2)*(L+M+2)
            LM1 = (L+M+1)*(L-M+1)
            YPLEG(IL+2,I) = ((2*L+3)*MU*YPLEG(IL+1,I)-SQRT(LM1)*        
     $            YPLEG(IL,I))/SQRT(LM2)
          END IF
  10    CONTINUE
  20  CONTINUE
      IF (M .GT. 0) THEN 
        CM1 = CM2
        DM1 = DM2
      END IF 
      RETURN        
      END

