!*DECK CFFTI1
      SUBROUTINE ZFFTI1 (N,WA,IFAC)
!***BEGIN PROLOGUE  CFFTI1
!***PURPOSE  Initialize a real and an integer work array for CFFTF1 and
!            CFFTB1.
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A2
!***TYPE      COMPLEX (RFFTI1-S, CFFTI1-C)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine CFFTI1 initializes the work arrays WA and IFAC which are
!  used in both CFFTF1 and CFFTB1.  The prime factorization of N and a
!  tabulation of the trigonometric functions are computed and stored in
!  IFAC and WA, respectively.
!
!  Input Parameter
!
!  N       the length of the sequence to be transformed
!
!  Output Parameters
!
!  WA      a real work array which must be dimensioned at least 2*N.
!
!  IFAC    an integer work array which must be dimensioned at least 15.
!
!          The same work arrays can be used for both CFFTF1 and CFFTB1
!          as long as N remains unchanged.  Different WA and IFAC arrays
!          are required for different values of N.  The contents of
!          WA and IFAC must not be changed between calls of CFFTF1 or
!          CFFTB1.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing references to intrinsic function FLOAT
!               to REAL, and
!           (c) changing definition of variable TPI by using
!               FORTRAN intrinsic function ATAN instead of a DATA
!               statement.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900131  Routine changed from subsidiary to user-callable.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CFFTI1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
!***FIRST EXECUTABLE STATEMENT  CFFTI1
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J.le.4) GO TO 102
      GO TO 103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRYH(4)+2*(J-4)
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR.eq.0) GO TO 105
      GO TO 101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI =  6.28318530717958647692D0
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.0D0
            WA(I) = 0.0D0
            LD = LD+L1
            FI = 0.0D0
            ARGLD = FLOAT(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.D0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
