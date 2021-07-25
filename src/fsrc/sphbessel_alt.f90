!
!   sphbessel_alt.f90
!   pdb2saxs
!
!   Created by Robert Rambo on 6/16/12.
!   Copyright 2012 Lawrence Berkeley national lab. All rights reserved.
!
SUBROUTINE SBESJY  (X,LMAX, J,Y,JP,YP, IFAIL )
    !---------------------------------------------------------------------
    !   REAL SPHERICAL BESSEL FUNCTIONS AND X DERIVATIVES
    !            j , y , j', y'                    FROM L=0 TO L=LMAX
    !        FOR REAL X > SQRT(ACCUR) (E.G. 1D-7)    AND INTEGER LMAX
    !
    !  J (L)  =      j/L/(X) STORES   REGULAR SPHERICAL BESSEL FUNCTION:
    !  JP(L)  = D/DX j/L/(X)            j(0) =  SIN(X)/X
    !  Y (L)  =      y/L/(X) STORES IRREGULAR SPHERICAL BESSEL FUNCTION:
    !  YP(L)  = D/DX y/L/(X)            y(0) = -COS(X)/X
    !
    !    IFAIL = -1 FOR ARGUMENTS OUT OF RANGE
    !          =  0 FOR ALL RESULTS SATISFACTORY
    !
    !   USING LENTZ-THOMPSON EVALUATION OF CONTINUED FRACTION CF1,
    !   AND TRIGONOMETRIC FORMS FOR L = 0 SOLUTIONS.
    !   LMAX IS LARGEST L NEEDED AND MUST BE <= MAXL, THE ARRAY INDEX.
    !   MAXL CAN BE DELETED AND ALL THE ARRAYS DIMENSIONED (0:*)
    !   SMALL IS MACHINE DEPENDENT, ABOUT SQRT(MINIMUM REAL NUMBER),
    !         SO 1D-150 FOR DOUBLE PRECISION ON VAX, PCS ETC.
    !   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
    !   IN OSCILLATING REGION X .GE.  [ SQRT{LMAX*(LMAX+1)} ]
    !   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8
    !   IS THE SMALLEST NUMBER WITH 1+ACC8.NE.1 FOR OUR WORKING PRECISION
    !   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
    !   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
    !   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE:
    !   THE VARIABLE ERR IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
    !
    !   NOTE: FOR X=1 AND L=100  J = 7.4 E-190     Y = -6.7+E186    1.4.94
    !---------------------------------------------------------------------
    !   AUTHOR :   A.R.BARNETT       MANCHESTER    12 MARCH 1990/95
    !                                AUCKLAND      12 MARCH 1991
    !---------------------------------------------------------------------
    IMPLICIT    NONE
    include "mkl_vml.fi"
    INTEGER     LIMIT, MAXL, LMAX, IFAIL, NFP, L
    PARAMETER ( LIMIT = 20000, MAXL = 1001 )
    DOUBLE PRECISION  J(0:MAXL), Y(0:MAXL), JP(0:MAXL), YP(0:MAXL)
    DOUBLE PRECISION  ZERO,ONE,TWO,THREE,SMALL, ACCUR, TK,SL, ERR
    DOUBLE PRECISION  XINV, CF1,DCF1, DEN, C,D, OMEGA, TWOXI
    DOUBLE PRECISION, INTENT(IN):: X
    real(kind=4), DIMENSION(1)  :: fA, vsin, vcos
    INTEGER, PARAMETER          :: vec_len = 1
    PARAMETER ( ZERO  = 0.0D0  , ONE   = 1.0D0 , TWO = 2.0D0 )
    PARAMETER ( SMALL = 1.D-150, THREE = 3.0D0 )
    COMMON /STEDE/    ERR,NFP       ! not required in code
    !-------

    ACCUR = 1.D-14                  ! suitable for double precision
    IFAIL = -1                      ! user to check on exit
    IF (X .LT. DSQRT(ACCUR) )       GOTO 50
    !-------TAKES CARE OF NEGATIVE X ... USE REFLECTION FORMULA
    !-------BEGIN CALCULATION OF CF1 UNLESS LMAX = 0, WHEN SOLUTIONS BELOW
    XINV  = ONE / X


    IF (LMAX .GT. 0) THEN
        TWOXI =     XINV + XINV
        SL  =  REAL(LMAX)* XINV     ! used also in do loop 3
        TK  =  TWO * SL  + XINV * THREE
        CF1 =  SL                   ! initial value of CF1
        DEN =  ONE                  ! unnormalised j(Lmax,x)
        IF ( ABS(CF1) .LT. SMALL ) CF1 = SMALL
        C   = CF1                   ! inverse ratio of A convergents
        D   = ZERO                  ! direct  ratio of B convergents
        DO L = 1,LIMIT
            C   = TK - ONE / C
            D   = TK - D
            IF ( ABS(C) .LT. SMALL ) C = SMALL
            IF ( ABS(D) .LT. SMALL ) D = SMALL
            D   = ONE / D
            DCF1= D   * C
            CF1 = CF1 * DCF1
            IF ( D .LT. ZERO ) DEN = - DEN
            IF ( ABS(DCF1 - ONE) .LE. ACCUR ) GOTO 20
            TK   = TK + TWOXI
            NFP  = L                 ! ie number in loop
        END DO
        GOTO 50             ! error exit, no convergence
        20       CONTINUE
        ERR = ACCUR * DSQRT(DBLE(NFP))    ! error estimate
        J (LMAX) = DEN      ! lower-case j's  really
        JP(LMAX) = CF1 * DEN
        !------ DOWNWARD RECURSION TO L=0  AS SPHERICAL BESSEL FUNCTIONS
        DO L =  LMAX, 1, -1
            J (L-1)  = (SL + XINV) * J(L)   + JP(L)
            SL  =  SL - XINV
            JP(L-1)  =  SL * J(L-1)          - J(L)
        END DO
        DEN = J(0)
    END IF                           ! end loop for Lmax GT 0

    !------ CALCULATE THE L=0 SPHERICAL BESSEL FUNCTIONS DIRECTLY
    fA(1) = X
    vsin(1) = 0.0
    vcos(1) = 0.0
    !WRITE (*,*) fA(1), vsin(1), vcos(1)
    !CALL vssincos(vec_len,fA,vsin,vcos)

    J(0)   =  XINV * DSIN(X)
    Y(0)   = -XINV * DCOS(X)
    !J (0)   =  XINV * vsin(1)
    !Y (0)   = -XINV * vcos(1)
    JP(0)   = -Y(0) - XINV * J(0)
    YP(0)   =  J(0) - XINV * Y(0)
    IF (LMAX .GT. 0) THEN
        OMEGA  =  J(0) / DEN
        SL  = ZERO
        DO L = 1 , LMAX
            J (L) = OMEGA * J (L)
            JP(L) = OMEGA * JP(L)
            Y (L) = SL * Y(L-1)   -   YP(L-1)
            SL  = SL + XINV
            YP(L) = Y(L-1)  -  (SL + XINV) * Y(L)
        END DO
    END IF
    IFAIL = 0                       ! calculations successful
    RETURN
    !---------------------------------------------------------------------
    !       ERROR TRAPS
    !---------------------------------------------------------------------
    50  IF (X .LT. ZERO) THEN
        WRITE(6,1000) X
    ELSEIF (X .EQ. ZERO) THEN
        IFAIL = 0
        J(0) = ONE
        DO L = 1, LMAX
            J(L) = ZERO     ! remaining arrays untouched
        END DO
    END IF
    1000   FORMAT(' X NEGATIVE !',1PE15.5,'    USE REFLECTION FORMULA'/)
    RETURN
END SUBROUTINE SBESJY
!---------------------------------------------------------------------
!       END OF SUBROUTINE SBESJY
!---------------------------------------------------------------------