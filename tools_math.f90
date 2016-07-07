! Copyright (c) 2013-2016 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>, Felix Kannemann
! <felix.kannemann@dal.ca>, Erin R. Johnson <erin.johnson@dal.ca>,
! Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider
! <hs7@post.queensu.ca>, and Axel D. Becke <axel.becke@dal.ca>
!
! postg is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
module tools_math
  
  ! spline and linpack
  public :: spline
  ! lebedev for postg
  public :: LD0006, LD0014, LD0026, LD0038, LD0050, LD0074, LD0086,&
            LD0110, LD0146, LD0170, LD0194, LD0230, LD0266, LD0302,&
            LD0350, LD0434, LD0590, LD0770, LD0974, LD1202, LD1454,&
            LD1730, LD2030, LD2354, LD2702, LD3074, LD3470, LD3890,&
            LD4334, LD4802, LD5294, LD5810

contains

  SUBROUTINE SPLINE (H,Y,A,B,C,N,BCR)
    !
    !     NATURAL CUBIC SPLINE THROUGH DATA POINTS ON THE UNIFORM R-MESH.
    !     RETURNS COEFFICIENTS A,B,C DEFINING THE CUBIC POLYNOMIAL TO THE
    !     RIGHT OF EACH KNOT (INCLUDING A(0),B(0),C(0)).
    !
    !     INPUT:
    !         H    - INTERVAL BETWEEN KNOTS
    !         Y    - ARRAY OF ORDINATE VALUES
    !         N    - THE NUMBER OF DATA POINTS
    !         BCR  - VALUE OF FUNCTION AT R=INFINITY
    !
    !     OUTPUT:
    !         A    - ARRAY OF COEFFICIENTS OF X
    !         B    - ARRAY OF COEFFICIENTS OF X**2
    !         C    - ARRAY OF COEFFICIENTS OF X**3
    !
    IMPLICIT REAL*8(A-H,O-Z)
    DIMENSION Y(0:N),A(0:N),B(0:N),C(0:N)
    DIMENSION D(N),E(N),BB(N)
    intent(in) :: H, Y, N, BCR
    intent(out) :: A, B, C

    THRD=1.D0/3.D0
    THRD2=2.D0*THRD
    HINV=1.D0/H
    H2INV=HINV*HINV

    DO I=1,N-1
       D(I)=4.D0
       E(I)=1.D0
       BB(I)=3.D0*H2INV*(Y(I-1)-2.D0*Y(I)+Y(I+1))
    enddo
    D(N)=4.D0
    E(N)=0.D0
    BB(N)=3.D0*H2INV*(Y(N-1)-2.D0*Y(N)+BCR)

    CALL DPTSL(N,D,E,BB)

    A(0)=HINV*(Y(1)-Y(0))-THRD*H*BB(1)
    B(0)=0.D0
    C(0)=THRD*HINV*BB(1)
    DO I=1,N-1
       A(I)=HINV*(Y(I+1)-Y(I))-THRD2*H*BB(I)-THRD*H*BB(I+1)
       B(I)=BB(I)
       C(I)=THRD*HINV*(BB(I+1)-BB(I))
    enddo
    A(N)=HINV*(BCR-Y(N))-THRD2*H*BB(N)
    B(N)=BB(N)
    C(N)=-THRD*HINV*BB(N)

  END SUBROUTINE SPLINE

  SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(1)
      REAL*8 A(LDA,1),Z(1)
      REAL*8 RCOND
!
!     DGECO FACTORS A REAL*8 MATRIX BY GAUSSIAN ELIMINATION
!     AND ESTIMATES THE CONDITION OF THE MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.
!
!     ON ENTRY
!
!        A       REAL*8(LDA, N)
!                THE MATRIX TO BE FACTORED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        RCOND   REAL*8
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!
!        Z       REAL*8(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DGEFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,DSIGN
!
!     INTERNAL VARIABLES
!
      REAL*8 EK,T,WK,WKM
      REAL*8 ANORM,S,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
!
!
!     COMPUTE 1-NORM OF A
!
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
!
!     FACTOR
!
      CALL DGEFA(A,LDA,N,IPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.
!
!     SOLVE TRANS(U)*W = E
!
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!
!     SOLVE TRANS(L)*Y = W
!
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!
      YNORM = 1.0D0
!
!     SOLVE L*V = Y
!
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
!     SOLVE  U*Z = V
!
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
  END subroutine

  SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO
      REAL*8 A(LDA,1)
!
!     DGEFA FACTORS A REAL*8 MATRIX BY GAUSSIAN ELIMINATION.
!
!     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
!
!     ON ENTRY
!
!        A       REAL*8(LDA, N)
!                THE MATRIX TO BE FACTORED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
!                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL,IDAMAX
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      INTEGER J,K,KP1,L,NM1
!
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
!
!        FIND L = PIVOT INDEX
!
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
!
!           INTERCHANGE IF NECESSARY
!
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
!
!           COMPUTE MULTIPLIERS
!
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
  END subroutine

  SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      REAL*8 A(LDA,1),B(1)
!
!     DGESL SOLVES THE REAL*8 SYSTEM
!     A * X = B  OR  TRANS(A) * X = B
!     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
!
!     ON ENTRY
!
!        A       REAL*8(LDA, N)
!                THE OUTPUT FROM DGECO OR DGEFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DGECO OR DGEFA.
!
!        B       REAL*8(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
!                            TRANS(A)  IS THE TRANSPOSE.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0
!        OR DGEFA HAS SET INFO .EQ. 0 .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DDOT
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      INTEGER K,KB,L,NM1
!
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B
!
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
!
!        NOW SOLVE  U*X = Y
!
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
!
!        NOW SOLVE TRANS(L)*X = Y
!
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      REAL*8 A(LDA,1),DET(2),WORK(1)
!
!     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
!     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
!
!     ON ENTRY
!
!        A       REAL*8(LDA, N)
!                THE OUTPUT FROM DGECO OR DGEFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DGECO OR DGEFA.
!
!        WORK    REAL*8(N)
!                WORK VECTOR.  CONTENTS DESTROYED.
!
!        JOB     INTEGER
!                = 11   BOTH DETERMINANT AND INVERSE.
!                = 01   INVERSE ONLY.
!                = 10   DETERMINANT ONLY.
!
!     ON RETURN
!
!        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE UNCHANGED.
!
!        DET     REAL*8(2)
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
!        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
!        INFO .EQ. 0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL,DSWAP
!     FORTRAN DABS,MOD
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      REAL*8 TEN
      INTEGER I,J,K,KB,KP1,L,NM1
!
!
!     COMPUTE DETERMINANT
!
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
!        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     COMPUTE INVERSE(U)
!
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM INVERSE(U)*INVERSE(L)
!
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
      INTEGER LDA,N,ML,MU,IPVT(1)
      REAL*8 ABD(LDA,1),Z(1)
      REAL*8 RCOND
!
!     DGBCO FACTORS A REAL*8 BAND MATRIX BY GAUSSIAN
!     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DGBFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DGBCO BY DGBSL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DGBCO BY DGBSL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DGBCO BY DGBDI.
!
!     ON ENTRY
!
!        ABD     REAL*8(LDA, N)
!                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
!                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
!                SEE THE COMMENTS BELOW FOR DETAILS.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!                LDA MUST BE .GE. 2*ML + MU + 1 .
!
!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.
!
!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!                0 .LE. ML .LT. N .
!
!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!                0 .LE. MU .LT. N .
!                MORE EFFICIENT IF  ML .LE. MU .
!
!     ON RETURN
!
!        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        RCOND   REAL*8
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!
!        Z       REAL*8(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     BAND STORAGE
!
!           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
!           WILL SET UP THE INPUT.
!
!                   ML = (BAND WIDTH BELOW THE DIAGONAL)
!                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-MU)
!                      I2 = MIN0(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
!           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
!           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
!           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
!           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
!           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
!
!     EXAMPLE..  IF THE ORIGINAL MATRIX IS
!
!           11 12 13  0  0  0
!           21 22 23 24  0  0
!            0 32 33 34 35  0
!            0  0 43 44 45 46
!            0  0  0 54 55 56
!            0  0  0  0 65 66
!
!      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN
!
!            *  *  *  +  +  +  , * = NOT USED
!            *  * 13 24 35 46  , + = USED FOR PIVOTING
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!           21 32 43 54 65  *
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DGBFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,MAX0,MIN0,DSIGN
!
!     INTERNAL VARIABLES
!
      REAL*8 EK,T,WK,WKM
      REAL*8 ANORM,S,SM,YNORM
      INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM
!
!
!     COMPUTE 1-NORM OF A
!
      ANORM = 0.0D0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(L,ABD(IS,J),1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
!
!     FACTOR
!
      CALL DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.
!
!     SOLVE TRANS(U)*W = E
!
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(ABD(M,K))) GO TO 30
            S = DABS(ABD(M,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (ABD(M,K) .EQ. 0.0D0) GO TO 40
            WK = WK/ABD(M,K)
            WKM = WKM/ABD(M,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
         MM = M
         IF (KP1 .GT. JU) GO TO 90
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + DABS(Z(J)+WKM*ABD(MM,J))
               Z(J) = Z(J) + WK*ABD(MM,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!
!     SOLVE TRANS(L)*Y = W
!
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(ML,N-K)
         IF (K .LT. N) Z(K) = Z(K) + DDOT(LM,ABD(M+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!
      YNORM = 1.0D0
!
!     SOLVE L*V = Y
!
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN0(ML,N-K)
         IF (K .LT. N) CALL DAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
!     SOLVE  U*Z = W
!
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(ABD(M,K))) GO TO 150
            S = DABS(ABD(M,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (ABD(M,K) .NE. 0.0D0) Z(K) = Z(K)/ABD(M,K)
         IF (ABD(M,K) .EQ. 0.0D0) Z(K) = 1.0D0
         LM = MIN0(K,M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL DAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)
  160 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
  END subroutine

  SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
      INTEGER LDA,N,ML,MU,IPVT(1),INFO
      REAL*8 ABD(LDA,1)
!
!     DGBFA FACTORS A REAL*8 BAND MATRIX BY ELIMINATION.
!
!     DGBFA IS USUALLY CALLED BY DGBCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!
!     ON ENTRY
!
!        ABD     REAL*8(LDA, N)
!                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
!                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
!                SEE THE COMMENTS BELOW FOR DETAILS.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!                LDA MUST BE .GE. 2*ML + MU + 1 .
!
!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.
!
!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!                0 .LE. ML .LT. N .
!
!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!                0 .LE. MU .LT. N .
!                MORE EFFICIENT IF  ML .LE. MU .
!     ON RETURN
!
!        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF
!                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.
!
!     BAND STORAGE
!
!           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
!           WILL SET UP THE INPUT.
!
!                   ML = (BAND WIDTH BELOW THE DIAGONAL)
!                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-MU)
!                      I2 = MIN0(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
!           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
!           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
!           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
!           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
!           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL,IDAMAX
!     FORTRAN MAX0,MIN0
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      INTEGER I,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
!
!
      M = ML + MU + 1
      INFO = 0
!
!     ZERO INITIAL FILL-IN COLUMNS
!
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
!
!        ZERO NEXT FILL-IN COLUMN
!
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE
!
!        FIND L = PIVOT INDEX
!
         LM = MIN0(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
!
!           INTERCHANGE IF NECESSARY
!
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
!
!           COMPUTE MULTIPLIERS
!
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
  END subroutine

  SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
      INTEGER LDA,N,ML,MU,IPVT(1),JOB
      REAL*8 ABD(LDA,1),B(1)
!
!     DGBSL SOLVES THE REAL*8 BAND SYSTEM
!     A * X = B  OR  TRANS(A) * X = B
!     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.
!
!     ON ENTRY
!
!        ABD     REAL*8(LDA, N)
!                THE OUTPUT FROM DGBCO OR DGBFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!
!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.
!
!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!
!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DGBCO OR DGBFA.
!
!        B       REAL*8(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
!                            TRANS(A)  IS THE TRANSPOSE.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0
!        OR DGBFA HAS SET INFO .EQ. 0 .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DDOT
!     FORTRAN MIN0
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
!
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE L*Y = B
!
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
!
!        NOW SOLVE  U*X = Y
!
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
!
!        NOW SOLVE TRANS(L)*X = Y
!
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DGBDI(ABD,LDA,N,ML,MU,IPVT,DET)
      INTEGER LDA,N,ML,MU,IPVT(1)
      REAL*8 ABD(LDA,1),DET(2)
!
!     DGBDI COMPUTES THE DETERMINANT OF A BAND MATRIX
!     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.
!     IF THE INVERSE IS NEEDED, USE DGBSL  N  TIMES.
!
!     ON ENTRY
!
!        ABD     REAL*8(LDA, N)
!                THE OUTPUT FROM DGBCO OR DGBFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!
!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.
!
!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!
!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DGBCO OR DGBFA.
!
!     ON RETURN
!
!        DET     REAL*8(2)
!                DETERMINANT OF ORIGINAL MATRIX.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
!                OR  DET(1) = 0.0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     FORTRAN DABS
!
!     INTERNAL VARIABLES
!
      REAL*8 TEN
      INTEGER I,M
!
!
      M = ML + MU + 1
      DET(1) = 1.0D0
      DET(2) = 0.0D0
      TEN = 10.0D0
      DO 50 I = 1, N
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)
         DET(1) = ABD(M,I)*DET(1)
!     ...EXIT
         IF (DET(1) .EQ. 0.0D0) GO TO 60
   10    IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
            DET(1) = TEN*DET(1)
            DET(2) = DET(2) - 1.0D0
         GO TO 10
   20    CONTINUE
   30    IF (DABS(DET(1)) .LT. TEN) GO TO 40
            DET(1) = DET(1)/TEN
            DET(2) = DET(2) + 1.0D0
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPOCO(A,LDA,N,RCOND,Z,INFO)
      INTEGER LDA,N,INFO
      REAL*8 A(LDA,1),Z(1)
      REAL*8 RCOND
!
!     DPOCO FACTORS A REAL*8 SYMMETRIC POSITIVE DEFINITE
!     MATRIX AND ESTIMATES THE CONDITION OF THE MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DPOFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DPOCO BY DPOSL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DPOCO BY DPOSL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DPOCO BY DPODI.
!     TO COMPUTE  INVERSE(A) , FOLLOW DPOCO BY DPODI.
!
!     ON ENTRY
!
!        A       REAL*8(LDA, N)
!                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
!                DIAGONAL AND UPPER TRIANGLE ARE USED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
!                WHERE  TRANS(R)  IS THE TRANSPOSE.
!                THE STRICT LOWER TRIANGLE IS UNALTERED.
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
!
!        RCOND   REAL*8
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.
!
!        Z       REAL*8(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!                IF  INFO .NE. 0 , Z  IS UNCHANGED.
!
!        INFO    INTEGER
!                = 0  FOR NORMAL RETURN.
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DPOFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,DREAL,DSIGN
!
!     INTERNAL VARIABLES
!
      REAL*8 EK,T,WK,WKM
      REAL*8 ANORM,S,SM,YNORM
      INTEGER I,J,JM1,K,KB,KP1
!
!
!     FIND NORM OF A USING ONLY UPPER HALF
!
      DO 30 J = 1, N
         Z(J) = DASUM(J,A(1,J),1)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + DABS(A(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
!
!     FACTOR
!
      CALL DPOFA(A,LDA,N,INFO)
      IF (INFO .NE. 0) GO TO 180
!
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!        SOLVE TRANS(R)*W = E
!
         EK = 1.0D0
         DO 50 J = 1, N
            Z(J) = 0.0D0
   50    CONTINUE
         DO 110 K = 1, N
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
            IF (DABS(EK-Z(K)) .LE. A(K,K)) GO TO 60
               S = A(K,K)/DABS(EK-Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = DABS(WK)
            SM = DABS(WKM)
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
            KP1 = K + 1
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + DABS(Z(J)+WKM*A(K,J))
                  Z(J) = Z(J) + WK*A(K,J)
                  S = S + DABS(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*A(K,J)
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
!
!        SOLVE R*Y = W
!
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 120
               S = A(K,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  130    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
!
         YNORM = 1.0D0
!
!        SOLVE TRANS(R)*V = Y
!
         DO 150 K = 1, N
            Z(K) = Z(K) - DDOT(K-1,A(1,K),1,Z(1),1)
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 140
               S = A(K,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/A(K,K)
  150    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
!
!        SOLVE R*Z = V
!
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 160
               S = A(K,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  170    CONTINUE
!        MAKE ZNORM = 1.0
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
!
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
  180 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPOFA(A,LDA,N,INFO)
      INTEGER LDA,N,INFO
      REAL*8 A(LDA,1)
!
!     DPOFA FACTORS A REAL*8 SYMMETRIC POSITIVE DEFINITE
!     MATRIX.
!
!     DPOFA IS USUALLY CALLED BY DPOCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR DPOCO) = (1 + 18/N)*(TIME FOR DPOFA) .
!
!     ON ENTRY
!
!        A       REAL*8(LDA, N)
!                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
!                DIAGONAL AND UPPER TRIANGLE ARE USED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
!                WHERE  TRANS(R)  IS THE TRANSPOSE.
!                THE STRICT LOWER TRIANGLE IS UNALTERED.
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
!
!        INFO    INTEGER
!                = 0  FOR NORMAL RETURN.
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DDOT
!     FORTRAN DSQRT
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      REAL*8 S
      INTEGER J,JM1,K
!     BEGIN BLOCK WITH ...EXITS TO 40
!
!
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            S = A(J,J) - S
!     ......EXIT
            IF (S .LE. 0.0D0) GO TO 40
            A(J,J) = DSQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPOSL(A,LDA,N,B)
      INTEGER LDA,N
      REAL*8 A(LDA,1),B(1)
!
!     DPOSL SOLVES THE REAL*8 SYMMETRIC POSITIVE DEFINITE
!     SYSTEM A * X = B
!     USING THE FACTORS COMPUTED BY DPOCO OR DPOFA.
!
!     ON ENTRY
!
!        A       REAL*8(LDA, N)
!                THE OUTPUT FROM DPOCO OR DPOFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        B       REAL*8(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
!        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
!        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
!        CORRECTLY AND  INFO .EQ. 0 .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL DPOCO(A,LDA,N,RCOND,Z,INFO)
!           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
!           DO 10 J = 1, P
!              CALL DPOSL(A,LDA,N,C(1,J))
!        10 CONTINUE
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DDOT
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      INTEGER K,KB
!
!     SOLVE TRANS(R)*Y = B
!
      DO 10 K = 1, N
         T = DDOT(K-1,A(1,K),1,B(1),1)
         B(K) = (B(K) - T)/A(K,K)
   10 CONTINUE
!
!     SOLVE R*X = Y
!
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   20 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPODI(A,LDA,N,DET,JOB)
      INTEGER LDA,N,JOB
      REAL*8 A(LDA,1)
      REAL*8 DET(2)
!
!     DPODI COMPUTES THE DETERMINANT AND INVERSE OF A CERTAIN
!     REAL*8 SYMMETRIC POSITIVE DEFINITE MATRIX (SEE BELOW)
!     USING THE FACTORS COMPUTED BY DPOCO, DPOFA OR DQRDC.
!
!     ON ENTRY
!
!        A       REAL*8(LDA, N)
!                THE OUTPUT  A  FROM DPOCO OR DPOFA
!                OR THE OUTPUT  X  FROM DQRDC.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        JOB     INTEGER
!                = 11   BOTH DETERMINANT AND INVERSE.
!                = 01   INVERSE ONLY.
!                = 10   DETERMINANT ONLY.
!
!     ON RETURN
!
!        A       IF DPOCO OR DPOFA WAS USED TO FACTOR  A  THEN
!                DPODI PRODUCES THE UPPER HALF OF INVERSE(A) .
!                IF DQRDC WAS USED TO DECOMPOSE  X  THEN
!                DPODI PRODUCES THE UPPER HALF OF INVERSE(TRANS(X)*X)
!                WHERE TRANS(X) IS THE TRANSPOSE.
!                ELEMENTS OF  A  BELOW THE DIAGONAL ARE UNCHANGED.
!                IF THE UNITS DIGIT OF JOB IS ZERO,  A  IS UNCHANGED.
!
!        DET     REAL*8(2)
!                DETERMINANT OF  A  OR OF  TRANS(X)*X  IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DET(1) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
!        AND IF DPOCO OR DPOFA HAS SET INFO .EQ. 0 .
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL
!     FORTRAN MOD
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      REAL*8 S
      INTEGER I,J,JM1,K,KP1
!
!     COMPUTE DETERMINANT
!
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         S = 10.0D0
         DO 50 I = 1, N
            DET(1) = A(I,I)**2*DET(1)
!        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DET(1) .GE. 1.0D0) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     COMPUTE INVERSE(R)
!
      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM  INVERSE(R) * TRANS(INVERSE(R))
!
         DO 130 J = 1, N
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = A(K,J)
               CALL DAXPY(K,T,A(1,J),1,A(1,K),1)
  110       CONTINUE
  120       CONTINUE
            T = A(J,J)
            CALL DSCAL(J,T,A(1,J),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPPCO(AP,N,RCOND,Z,INFO)
      INTEGER N,INFO
      REAL*8 AP(1),Z(1)
      REAL*8 RCOND
!
!     DPPCO FACTORS A REAL*8 SYMMETRIC POSITIVE DEFINITE
!     MATRIX STORED IN PACKED FORM
!     AND ESTIMATES THE CONDITION OF THE MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DPPFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DPPCO BY DPPSL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DPPCO BY DPPSL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DPPCO BY DPPDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW DPPCO BY DPPDI.
!
!     ON ENTRY
!
!        AP      REAL*8 (N*(N+1)/2)
!                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
!                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
!                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
!                SEE COMMENTS BELOW FOR DETAILS.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED
!                FORM, SO THAT  A = TRANS(R)*R .
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
!
!        RCOND   REAL*8
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.
!
!        Z       REAL*8(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!                IF  INFO .NE. 0 , Z  IS UNCHANGED.
!
!        INFO    INTEGER
!                = 0  FOR NORMAL RETURN.
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
!
!     PACKED STORAGE
!
!          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
!          TRIANGLE OF A SYMMETRIC MATRIX.
!
!                K = 0
!                DO 20 J = 1, N
!                   DO 10 I = 1, J
!                      K = K + 1
!                      AP(K) = A(I,J)
!             10    CONTINUE
!             20 CONTINUE
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DPPFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,DREAL,DSIGN
!
!     INTERNAL VARIABLES
!
      REAL*8 EK,T,WK,WKM
      REAL*8 ANORM,S,SM,YNORM
      INTEGER I,IJ,J,JM1,J1,K,KB,KJ,KK,KP1
!
!
!     FIND NORM OF A
!
      J1 = 1
      DO 30 J = 1, N
         Z(J) = DASUM(J,AP(J1),1)
         IJ = J1
         J1 = J1 + J
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + DABS(AP(IJ))
            IJ = IJ + 1
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
!
!     FACTOR
!
      CALL DPPFA(AP,N,INFO)
      IF (INFO .NE. 0) GO TO 180
!
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!        SOLVE TRANS(R)*W = E
!
         EK = 1.0D0
         DO 50 J = 1, N
            Z(J) = 0.0D0
   50    CONTINUE
         KK = 0
         DO 110 K = 1, N
            KK = KK + K
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
            IF (DABS(EK-Z(K)) .LE. AP(KK)) GO TO 60
               S = AP(KK)/DABS(EK-Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = DABS(WK)
            SM = DABS(WKM)
            WK = WK/AP(KK)
            WKM = WKM/AP(KK)
            KP1 = K + 1
            KJ = KK + K
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + DABS(Z(J)+WKM*AP(KJ))
                  Z(J) = Z(J) + WK*AP(KJ)
                  S = S + DABS(Z(J))
                  KJ = KJ + J
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  KJ = KK + K
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*AP(KJ)
                     KJ = KJ + J
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
!
!        SOLVE R*Y = W
!
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. AP(KK)) GO TO 120
               S = AP(KK)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/AP(KK)
            KK = KK - K
            T = -Z(K)
            CALL DAXPY(K-1,T,AP(KK+1),1,Z(1),1)
  130    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
!
         YNORM = 1.0D0
!
!        SOLVE TRANS(R)*V = Y
!
         DO 150 K = 1, N
            Z(K) = Z(K) - DDOT(K-1,AP(KK+1),1,Z(1),1)
            KK = KK + K
            IF (DABS(Z(K)) .LE. AP(KK)) GO TO 140
               S = AP(KK)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/AP(KK)
  150    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
!
!        SOLVE R*Z = V
!
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. AP(KK)) GO TO 160
               S = AP(KK)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/AP(KK)
            KK = KK - K
            T = -Z(K)
            CALL DAXPY(K-1,T,AP(KK+1),1,Z(1),1)
  170    CONTINUE
!        MAKE ZNORM = 1.0
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
!
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
  180 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPPFA(AP,N,INFO)
      INTEGER N,INFO
      REAL*8 AP(1)
!
!     DPPFA FACTORS A REAL*8 SYMMETRIC POSITIVE DEFINITE
!     MATRIX STORED IN PACKED FORM.
!
!     DPPFA IS USUALLY CALLED BY DPPCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR DPPCO) = (1 + 18/N)*(TIME FOR DPPFA) .
!
!     ON ENTRY
!
!        AP      REAL*8 (N*(N+1)/2)
!                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
!                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
!                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
!                SEE COMMENTS BELOW FOR DETAILS.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED
!                FORM, SO THAT  A = TRANS(R)*R .
!
!        INFO    INTEGER
!                = 0  FOR NORMAL RETURN.
!                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT
!                     POSITIVE DEFINITE.
!
!
!     PACKED STORAGE
!
!          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
!          TRIANGLE OF A SYMMETRIC MATRIX.
!
!                K = 0
!                DO 20 J = 1, N
!                   DO 10 I = 1, J
!                      K = K + 1
!                      AP(K) = A(I,J)
!             10    CONTINUE
!             20 CONTINUE
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DDOT
!     FORTRAN DSQRT
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      REAL*8 S
      INTEGER J,JJ,JM1,K,KJ,KK
!     BEGIN BLOCK WITH ...EXITS TO 40
!
!
         JJ = 0
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            JM1 = J - 1
            KJ = JJ
            KK = 0
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               KJ = KJ + 1
               T = AP(KJ) - DDOT(K-1,AP(KK+1),1,AP(JJ+1),1)
               KK = KK + K
               T = T/AP(KK)
               AP(KJ) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            JJ = JJ + J
            S = AP(JJ) - S
!     ......EXIT
            IF (S .LE. 0.0D0) GO TO 40
            AP(JJ) = DSQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPPSL(AP,N,B)
      INTEGER N
      REAL*8 AP(1),B(1)
!
!     DPPSL SOLVES THE REAL*8 SYMMETRIC POSITIVE DEFINITE
!     SYSTEM A * X = B
!     USING THE FACTORS COMPUTED BY DPPCO OR DPPFA.
!
!     ON ENTRY
!
!        AP      REAL*8 (N*(N+1)/2)
!                THE OUTPUT FROM DPPCO OR DPPFA.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        B       REAL*8(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
!        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
!        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
!        CORRECTLY AND  INFO .EQ. 0 .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL DPPCO(AP,N,RCOND,Z,INFO)
!           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
!           DO 10 J = 1, P
!              CALL DPPSL(AP,N,C(1,J))
!        10 CONTINUE
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DDOT
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      INTEGER K,KB,KK
!
      KK = 0
      DO 10 K = 1, N
         T = DDOT(K-1,AP(KK+1),1,B(1),1)
         KK = KK + K
         B(K) = (B(K) - T)/AP(KK)
   10 CONTINUE
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/AP(KK)
         KK = KK - K
         T = -B(K)
         CALL DAXPY(K-1,T,AP(KK+1),1,B(1),1)
   20 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPPDI(AP,N,DET,JOB)
      INTEGER N,JOB
      REAL*8 AP(1)
      REAL*8 DET(2)
!
!     DPPDI COMPUTES THE DETERMINANT AND INVERSE
!     OF A REAL*8 SYMMETRIC POSITIVE DEFINITE MATRIX
!     USING THE FACTORS COMPUTED BY DPPCO OR DPPFA .
!
!     ON ENTRY
!
!        AP      REAL*8 (N*(N+1)/2)
!                THE OUTPUT FROM DPPCO OR DPPFA.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        JOB     INTEGER
!                = 11   BOTH DETERMINANT AND INVERSE.
!                = 01   INVERSE ONLY.
!                = 10   DETERMINANT ONLY.
!
!     ON RETURN
!
!        AP      THE UPPER TRIANGULAR HALF OF THE INVERSE .
!                THE STRICT LOWER TRIANGLE IS UNALTERED.
!
!        DET     REAL*8(2)
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DET(1) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
!        AND IF DPOCO OR DPOFA HAS SET INFO .EQ. 0 .
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL
!     FORTRAN MOD
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      REAL*8 S
      INTEGER I,II,J,JJ,JM1,J1,K,KJ,KK,KP1,K1
!
!     COMPUTE DETERMINANT
!
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         S = 10.0D0
         II = 0
         DO 50 I = 1, N
            II = II + I
            DET(1) = AP(II)**2*DET(1)
!        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DET(1) .GE. 1.0D0) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     COMPUTE INVERSE(R)
!
      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         KK = 0
         DO 100 K = 1, N
            K1 = KK + 1
            KK = KK + K
            AP(KK) = 1.0D0/AP(KK)
            T = -AP(KK)
            CALL DSCAL(K-1,T,AP(K1),1)
            KP1 = K + 1
            J1 = KK + 1
            KJ = KK + K
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = AP(KJ)
               AP(KJ) = 0.0D0
               CALL DAXPY(K,T,AP(K1),1,AP(J1),1)
               J1 = J1 + J
               KJ = KJ + J
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM  INVERSE(R) * TRANS(INVERSE(R))
!
         JJ = 0
         DO 130 J = 1, N
            J1 = JJ + 1
            JJ = JJ + J
            JM1 = J - 1
            K1 = 1
            KJ = J1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = AP(KJ)
               CALL DAXPY(K,T,AP(J1),1,AP(K1),1)
               K1 = K1 + K
               KJ = KJ + 1
  110       CONTINUE
  120       CONTINUE
            T = AP(JJ)
            CALL DSCAL(J,T,AP(J1),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPBCO(ABD,LDA,N,M,RCOND,Z,INFO)
      INTEGER LDA,N,M,INFO
      REAL*8 ABD(LDA,1),Z(1)
      REAL*8 RCOND
!
!     DPBCO FACTORS A REAL*8 SYMMETRIC POSITIVE DEFINITE
!     MATRIX STORED IN BAND FORM AND ESTIMATES THE CONDITION OF THE
!     MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DPBFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DPBCO BY DPBSL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DPBCO BY DPBSL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DPBCO BY DPBDI.
!
!     ON ENTRY
!
!        ABD     REAL*8(LDA, N)
!                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER
!                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE
!                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE
!                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!                LDA MUST BE .GE. M + 1 .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        M       INTEGER
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!                0 .LE. M .LT. N .
!
!     ON RETURN
!
!        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND
!                FORM, SO THAT  A = TRANS(R)*R .
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
!
!        RCOND   REAL*8
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.
!
!        Z       REAL*8(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!                IF  INFO .NE. 0 , Z  IS UNCHANGED.
!
!        INFO    INTEGER
!                = 0  FOR NORMAL RETURN.
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
!
!     BAND STORAGE
!
!           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,
!           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.
!
!                   M = (BAND WIDTH ABOVE DIAGONAL)
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-M)
!                      DO 10 I = I1, J
!                         K = I-J+M+1
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           THIS USES  M + 1  ROWS OF  A , EXCEPT FOR THE  M BY M
!           UPPER LEFT TRIANGLE, WHICH IS IGNORED.
!
!     EXAMPLE..  IF THE ORIGINAL MATRIX IS
!
!           11 12 13  0  0  0
!           12 22 23 24  0  0
!           13 23 33 34 35  0
!            0 24 34 44 45 46
!            0  0 35 45 55 56
!            0  0  0 46 56 66
!
!     THEN  N = 6 , M = 2  AND  ABD  SHOULD CONTAIN
!
!            *  * 13 24 35 46
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DPBFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,MAX0,MIN0,DREAL,DSIGN
!
!     INTERNAL VARIABLES
!
      REAL*8 EK,T,WK,WKM
      REAL*8 ANORM,S,SM,YNORM
      INTEGER I,J,J2,K,KB,KP1,L,LA,LB,LM,MU
!
!
!     FIND NORM OF A
!
      DO 30 J = 1, N
         L = MIN0(J,M+1)
         MU = MAX0(M+2-J,1)
         Z(J) = DASUM(L,ABD(MU,J),1)
         K = J - L
         IF (M .LT. MU) GO TO 20
         DO 10 I = MU, M
            K = K + 1
            Z(K) = Z(K) + DABS(ABD(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
!
!     FACTOR
!
      CALL DPBFA(ABD,LDA,N,M,INFO)
      IF (INFO .NE. 0) GO TO 180
!
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!        SOLVE TRANS(R)*W = E
!
         EK = 1.0D0
         DO 50 J = 1, N
            Z(J) = 0.0D0
   50    CONTINUE
         DO 110 K = 1, N
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
            IF (DABS(EK-Z(K)) .LE. ABD(M+1,K)) GO TO 60
               S = ABD(M+1,K)/DABS(EK-Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = DABS(WK)
            SM = DABS(WKM)
            WK = WK/ABD(M+1,K)
            WKM = WKM/ABD(M+1,K)
            KP1 = K + 1
            J2 = MIN0(K+M,N)
            I = M + 1
            IF (KP1 .GT. J2) GO TO 100
               DO 70 J = KP1, J2
                  I = I - 1
                  SM = SM + DABS(Z(J)+WKM*ABD(I,J))
                  Z(J) = Z(J) + WK*ABD(I,J)
                  S = S + DABS(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  I = M + 1
                  DO 80 J = KP1, J2
                     I = I - 1
                     Z(J) = Z(J) + T*ABD(I,J)
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
!
!        SOLVE  R*Y = W
!
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 120
               S = ABD(M+1,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
            LM = MIN0(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            T = -Z(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  130    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
!
         YNORM = 1.0D0
!
!        SOLVE TRANS(R)*V = Y
!
         DO 150 K = 1, N
            LM = MIN0(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            Z(K) = Z(K) - DDOT(LM,ABD(LA,K),1,Z(LB),1)
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 140
               S = ABD(M+1,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
  150    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
!
!        SOLVE  R*Z = W
!
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 160
               S = ABD(M+1,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
            LM = MIN0(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            T = -Z(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  170    CONTINUE
!        MAKE ZNORM = 1.0
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
!
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
  180 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPBFA(ABD,LDA,N,M,INFO)
      INTEGER LDA,N,M,INFO
      REAL*8 ABD(LDA,1)
!
!     DPBFA FACTORS A REAL*8 SYMMETRIC POSITIVE DEFINITE
!     MATRIX STORED IN BAND FORM.
!
!     DPBFA IS USUALLY CALLED BY DPBCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!
!     ON ENTRY
!
!        ABD     REAL*8(LDA, N)
!                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER
!                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE
!                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE
!                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!                LDA MUST BE .GE. M + 1 .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        M       INTEGER
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!                0 .LE. M .LT. N .
!
!     ON RETURN
!
!        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND
!                FORM, SO THAT  A = TRANS(R)*R .
!
!        INFO    INTEGER
!                = 0  FOR NORMAL RETURN.
!                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT
!                     POSITIVE DEFINITE.
!
!     BAND STORAGE
!
!           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,
!           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.
!
!                   M = (BAND WIDTH ABOVE DIAGONAL)
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-M)
!                      DO 10 I = I1, J
!                         K = I-J+M+1
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DDOT
!     FORTRAN MAX0,DSQRT
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      REAL*8 S
      INTEGER IK,J,JK,K,MU
!     BEGIN BLOCK WITH ...EXITS TO 40
!
!
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            IK = M + 1
            JK = MAX0(J-M,1)
            MU = MAX0(M+2-J,1)
            IF (M .LT. MU) GO TO 20
            DO 10 K = MU, M
               T = ABD(K,J) - DDOT(K-MU,ABD(IK,JK),1,ABD(MU,J),1)
               T = T/ABD(M+1,JK)
               ABD(K,J) = T
               S = S + T*T
               IK = IK - 1
               JK = JK + 1
   10       CONTINUE
   20       CONTINUE
            S = ABD(M+1,J) - S
!     ......EXIT
            IF (S .LE. 0.0D0) GO TO 40
            ABD(M+1,J) = DSQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPBSL(ABD,LDA,N,M,B)
      INTEGER LDA,N,M
      REAL*8 ABD(LDA,1),B(1)
!
!     DPBSL SOLVES THE REAL*8 SYMMETRIC POSITIVE DEFINITE
!     BAND SYSTEM  A*X = B
!     USING THE FACTORS COMPUTED BY DPBCO OR DPBFA.
!
!     ON ENTRY
!
!        ABD     REAL*8(LDA, N)
!                THE OUTPUT FROM DPBCO OR DPBFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        M       INTEGER
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!        B       REAL*8(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
!        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
!        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
!        CORRECTLY AND  INFO .EQ. 0 .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL DPBCO(ABD,LDA,N,RCOND,Z,INFO)
!           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
!           DO 10 J = 1, P
!              CALL DPBSL(ABD,LDA,N,C(1,J))
!        10 CONTINUE
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DDOT
!     FORTRAN MIN0
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      INTEGER K,KB,LA,LB,LM
!
!     SOLVE TRANS(R)*Y = B
!
      DO 10 K = 1, N
         LM = MIN0(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         T = DDOT(LM,ABD(LA,K),1,B(LB),1)
         B(K) = (B(K) - T)/ABD(M+1,K)
   10 CONTINUE
!
!     SOLVE R*X = Y
!
      DO 20 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         B(K) = B(K)/ABD(M+1,K)
         T = -B(K)
         CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   20 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DPBDI(ABD,LDA,N,M,DET)
      INTEGER LDA,N,M
      REAL*8 ABD(LDA,1)
      REAL*8 DET(2)
!
!     DPBDI COMPUTES THE DETERMINANT
!     OF A REAL*8 SYMMETRIC POSITIVE DEFINITE BAND MATRIX
!     USING THE FACTORS COMPUTED BY DPBCO OR DPBFA.
!     IF THE INVERSE IS NEEDED, USE DPBSL  N  TIMES.
!
!     ON ENTRY
!
!        ABD     REAL*8(LDA, N)
!                THE OUTPUT FROM DPBCO OR DPBFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        M       INTEGER
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!     ON RETURN
!
!        DET     REAL*8(2)
!                DETERMINANT OF ORIGINAL MATRIX IN THE FORM
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DET(1) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!
!     INTERNAL VARIABLES
!
      REAL*8 S
      INTEGER I
!
!     COMPUTE DETERMINANT
!
      DET(1) = 1.0D0
      DET(2) = 0.0D0
      S = 10.0D0
      DO 50 I = 1, N
         DET(1) = ABD(M+1,I)**2*DET(1)
!     ...EXIT
         IF (DET(1) .EQ. 0.0D0) GO TO 60
   10    IF (DET(1) .GE. 1.0D0) GO TO 20
            DET(1) = S*DET(1)
            DET(2) = DET(2) - 1.0D0
         GO TO 10
   20    CONTINUE
   30    IF (DET(1) .LT. S) GO TO 40
            DET(1) = DET(1)/S
            DET(2) = DET(2) + 1.0D0
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DSICO(A,LDA,N,KPVT,RCOND,Z)
      INTEGER LDA,N,KPVT(1)
      REAL*8 A(LDA,1),Z(1)
      REAL*8 RCOND
!
!     DSICO FACTORS A REAL*8 SYMMETRIC MATRIX BY ELIMINATION
!     WITH SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE
!     MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DSIFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DSICO BY DSISL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DSICO BY DSISL.
!     TO COMPUTE  INVERSE(A) , FOLLOW DSICO BY DSIDI.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DSICO BY DSIDI.
!     TO COMPUTE  INERTIA(A), FOLLOW DSICO BY DSIDI.
!
!     ON ENTRY
!
!        A       REAL*8(LDA, N)
!                THE SYMMETRIC MATRIX TO BE FACTORED.
!                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     OUTPUT
!
!        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
!                WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
!
!        KPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        RCOND   REAL*8
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!
!        Z       REAL*8(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DSIFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,IABS,DSIGN
!
!     INTERNAL VARIABLES
!
      REAL*8 AK,AKM1,BK,BKM1,DENOM,EK,T
      REAL*8 ANORM,S,YNORM
      INTEGER I,INFO,J,JM1,K,KP,KPS,KS
!
!
!     FIND NORM OF A USING ONLY UPPER HALF
!
      DO 30 J = 1, N
         Z(J) = DASUM(J,A(1,J),1)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + DABS(A(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
!
!     FACTOR
!
      CALL DSIFA(A,LDA,N,KPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!     SOLVE U*D*W = E
!
      EK = 1.0D0
      DO 50 J = 1, N
         Z(J) = 0.0D0
   50 CONTINUE
      K = N
   60 IF (K .EQ. 0) GO TO 120
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = IABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 90
               S = DABS(A(K,K))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   90       CONTINUE
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 110
  100    CONTINUE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
      GO TO 60
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!
!     SOLVE TRANS(U)*Y = W
!
      K = 1
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2) Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!
      YNORM = 1.0D0
!
!     SOLVE U*D*V = Y
!
      K = N
  170 IF (K .EQ. 0) GO TO 230
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = IABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
            IF (KS .EQ. 2) CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 200
               S = DABS(A(K,K))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 220
  210    CONTINUE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
      GO TO 170
  230 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
!     SOLVE TRANS(U)*Z = V
!
      K = 1
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2) Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         K = K + KS
      GO TO 240
  270 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
  END subroutine

  SUBROUTINE DSIFA(A,LDA,N,KPVT,INFO)
      INTEGER LDA,N,KPVT(1),INFO
      REAL*8 A(LDA,1)
!
!     DSIFA FACTORS A REAL*8 SYMMETRIC MATRIX BY ELIMINATION
!     WITH SYMMETRIC PIVOTING.
!
!     TO SOLVE  A*X = B , FOLLOW DSIFA BY DSISL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DSIFA BY DSISL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DSIFA BY DSIDI.
!     TO COMPUTE  INERTIA(A) , FOLLOW DSIFA BY DSIDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW DSIFA BY DSIDI.
!
!     ON ENTRY
!
!        A       REAL*8(LDA,N)
!                THE SYMMETRIC MATRIX TO BE FACTORED.
!                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
!                WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
!
!        KPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
!                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
!                     BUT IT DOES INDICATE THAT DSISL OR DSIDI MAY
!                     DIVIDE BY ZERO IF CALLED.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSWAP,IDAMAX
!     FORTRAN DABS,DMAX1,DSQRT
!
!     INTERNAL VARIABLES
!
      REAL*8 AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      REAL*8 ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP
      LOGICAL SWAP
!
!
!     INITIALIZE
!
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0
!
      INFO = 0
!
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
!
      K = N
   10 CONTINUE
!
!        LEAVE THE LOOP IF K=0 OR K=1.
!
!     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (A(1,1) .EQ. 0.0D0) INFO = 1
!     ......EXIT
            GO TO 200
   20    CONTINUE
!
!        THIS SECTION OF CODE DETERMINES THE KIND OF
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
!        REQUIRED.
!
         KM1 = K - 1
         ABSAKK = DABS(A(K,K))
!
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!        COLUMN K.
!
         IMAX = IDAMAX(K-1,A(1,K),1)
         COLMAX = DABS(A(IMAX,K))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
!
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!           ROW IMAX.
!
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            DO 40 J = IMAXP1, K
               ROWMAX = DMAX1(ROWMAX,DABS(A(IMAX,J)))
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = IDAMAX(IMAX-1,A(1,IMAX),1)
               ROWMAX = DMAX1(ROWMAX,DABS(A(JMAX,IMAX)))
   50       CONTINUE
            IF (DABS(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100
!
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
!
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
!
!           1 X 1 PIVOT BLOCK.
!
            IF (.NOT.SWAP) GO TO 120
!
!              PERFORM AN INTERCHANGE.
!
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  T = A(J,K)
                  A(J,K) = A(IMAX,J)
                  A(IMAX,J) = T
  110          CONTINUE
  120       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
            DO 130 JJ = 1, KM1
               J = K - JJ
               MULK = -A(J,K)/A(K,K)
               T = MULK
               CALL DAXPY(J,T,A(1,K),1,A(1,J),1)
               A(J,K) = MULK
  130       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
            IF (.NOT.SWAP) GO TO 160
!
!              PERFORM AN INTERCHANGE.
!
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  T = A(J,K-1)
                  A(J,K-1) = A(IMAX,J)
                  A(IMAX,J) = T
  150          CONTINUE
               T = A(K-1,K)
               A(K-1,K) = A(IMAX,K)
               A(IMAX,K) = T
  160       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = A(K,K)/A(K-1,K)
               AKM1 = A(K-1,K-1)/A(K-1,K)
               DENOM = 1.0D0 - AK*AKM1
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  BK = A(J,K)/A(K-1,K)
                  BKM1 = A(J,K-1)/A(K-1,K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL DAXPY(J,T,A(1,K),1,A(1,J),1)
                  T = MULKM1
                  CALL DAXPY(J,T,A(1,K-1),1,A(1,J),1)
                  A(J,K) = MULK
                  A(J,K-1) = MULKM1
  170          CONTINUE
  180       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DSISL(A,LDA,N,KPVT,B)
      INTEGER LDA,N,KPVT(1)
      REAL*8 A(LDA,1),B(1)
!
!     DSISL SOLVES THE REAL*8 SYMMETRIC SYSTEM
!     A * X = B
!     USING THE FACTORS COMPUTED BY DSIFA.
!
!     ON ENTRY
!
!        A       REAL*8(LDA,N)
!                THE OUTPUT FROM DSIFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        KPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DSIFA.
!
!        B       REAL*8(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO MAY OCCUR IF  DSICO  HAS SET RCOND .EQ. 0.0
!        OR  DSIFA  HAS SET INFO .NE. 0  .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL DSIFA(A,LDA,N,KPVT,INFO)
!           IF (INFO .NE. 0) GO TO ...
!           DO 10 J = 1, P
!              CALL DSISL(A,LDA,N,KPVT,C(1,J))
!        10 CONTINUE
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DDOT
!     FORTRAN IABS
!
!     INTERNAL VARIABLES.
!
      REAL*8 AK,AKM1,BK,BKM1,DENOM,TEMP
      INTEGER K,KP
!
!     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
!     D INVERSE TO B.
!
      K = N
   10 IF (K .EQ. 0) GO TO 80
         IF (KPVT(K) .LT. 0) GO TO 40
!
!           1 X 1 PIVOT BLOCK.
!
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
!
!                 INTERCHANGE.
!
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
!
!              APPLY THE TRANSFORMATION.
!
               CALL DAXPY(K-1,B(K),A(1,K),1,B(1),1)
   30       CONTINUE
!
!           APPLY D INVERSE.
!
            B(K) = B(K)/A(K,K)
            K = K - 1
         GO TO 70
   40    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
!
!                 INTERCHANGE.
!
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
!
!              APPLY THE TRANSFORMATION.
!
               CALL DAXPY(K-2,B(K),A(1,K),1,B(1),1)
               CALL DAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)
   60       CONTINUE
!
!           APPLY D INVERSE.
!
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = B(K)/A(K-1,K)
            BKM1 = B(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
   70    CONTINUE
      GO TO 10
   80 CONTINUE
!
!     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
!
      K = 1
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
!
!           1 X 1 PIVOT BLOCK.
!
            IF (K .EQ. 1) GO TO 110
!
!              APPLY THE TRANSFORMATION.
!
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
!
!                 INTERCHANGE.
!
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            K = K + 1
         GO TO 150
  120    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
            IF (K .EQ. 1) GO TO 140
!
!              APPLY THE TRANSFORMATION.
!
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)
               B(K+1) = B(K+1) + DDOT(K-1,A(1,K+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
!
!                 INTERCHANGE.
!
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DSIDI(A,LDA,N,KPVT,DET,INERT,WORK,JOB)
      INTEGER LDA,N,JOB
      REAL*8 A(LDA,1),WORK(1)
      REAL*8 DET(2)
      INTEGER KPVT(1),INERT(3)
!
!     DSIDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE
!     OF A REAL*8 SYMMETRIC MATRIX USING THE FACTORS FROM
!     DSIFA.
!
!     ON ENTRY
!
!        A       REAL*8(LDA,N)
!                THE OUTPUT FROM DSIFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY A.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX A.
!
!        KPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DSIFA.
!
!        WORK    REAL*8(N)
!                WORK VECTOR.  CONTENTS DESTROYED.
!
!        JOB     INTEGER
!                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE
!                   IF  C .NE. 0, THE INVERSE IS COMPUTED,
!                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,
!                   IF  A .NE. 0, THE INERTIA IS COMPUTED.
!
!                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.
!
!     ON RETURN
!
!        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
!
!        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
!               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE
!               IS NEVER REFERENCED.
!
!        DET    REAL*8(2)
!               DETERMINANT OF ORIGINAL MATRIX.
!               DETERMINANT = DET(1) * 10.0**DET(2)
!               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0
!               OR DET(1) = 0.0.
!
!        INERT  INTEGER(3)
!               THE INERTIA OF THE ORIGINAL MATRIX.
!               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.
!               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.
!               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED
!        AND  DSICO  HAS SET RCOND .EQ. 0.0
!        OR  DSIFA  HAS SET  INFO .NE. 0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DCOPY,DDOT,DSWAP
!     FORTRAN DABS,IABS,MOD
!
!     INTERNAL VARIABLES.
!
      REAL*8 AKKP1,TEMP
      REAL*8 TEN,D,T,AK,AKP1
      INTEGER J,JB,K,KM1,KS,KSTEP
      LOGICAL NOINV,NODET,NOERT
!
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
      NOERT = MOD(JOB,1000)/100 .EQ. 0
!
      IF (NODET .AND. NOERT) GO TO 140
         IF (NOERT) GO TO 10
            INERT(1) = 0
            INERT(2) = 0
            INERT(3) = 0
   10    CONTINUE
         IF (NODET) GO TO 20
            DET(1) = 1.0D0
            DET(2) = 0.0D0
            TEN = 10.0D0
   20    CONTINUE
         T = 0.0D0
         DO 130 K = 1, N
            D = A(K,K)
!
!           CHECK IF 1 BY 1
!
            IF (KPVT(K) .GT. 0) GO TO 50
!
!              2 BY 2 BLOCK
!              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)
!                      (S  C)
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
!
               IF (T .NE. 0.0D0) GO TO 30
                  T = DABS(A(K,K+1))
                  D = (D/T)*A(K+1,K+1) - T
               GO TO 40
   30          CONTINUE
                  D = T
                  T = 0.0D0
   40          CONTINUE
   50       CONTINUE
!
            IF (NOERT) GO TO 60
               IF (D .GT. 0.0D0) INERT(1) = INERT(1) + 1
               IF (D .LT. 0.0D0) INERT(2) = INERT(2) + 1
               IF (D .EQ. 0.0D0) INERT(3) = INERT(3) + 1
   60       CONTINUE
!
            IF (NODET) GO TO 120
               DET(1) = D*DET(1)
               IF (DET(1) .EQ. 0.0D0) GO TO 110
   70             IF (DABS(DET(1)) .GE. 1.0D0) GO TO 80
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - 1.0D0
                  GO TO 70
   80             CONTINUE
   90             IF (DABS(DET(1)) .LT. TEN) GO TO 100
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + 1.0D0
                  GO TO 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
!
!     COMPUTE INVERSE(A)
!
      IF (NOINV) GO TO 270
         K = 1
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            IF (KPVT(K) .LT. 0) GO TO 180
!
!              1 BY 1
!
               A(K,K) = 1.0D0/A(K,K)
               IF (KM1 .LT. 1) GO TO 170
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)
                  DO 160 J = 1, KM1
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
!
!              2 BY 2
!
               T = DABS(A(K,K+1))
               AK = A(K,K)/T
               AKP1 = A(K+1,K+1)/T
               AKKP1 = A(K,K+1)/T
               D = T*(AK*AKP1 - 1.0D0)
               A(K,K) = AKP1/D
               A(K+1,K+1) = AK/D
               A(K,K+1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL DCOPY(KM1,A(1,K+1),1,WORK,1)
                  DO 190 J = 1, KM1
                     A(J,K+1) = DDOT(J,A(1,J),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  190             CONTINUE
                  A(K+1,K+1) = A(K+1,K+1) + DDOT(KM1,WORK,1,A(1,K+1),1)
                  A(K,K+1) = A(K,K+1) + DDOT(KM1,A(1,K),1,A(1,K+1),1)
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)
                  DO 200 J = 1, KM1
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  200             CONTINUE
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
!
!           SWAP
!
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               CALL DSWAP(KS,A(1,KS),1,A(1,K),1)
               DO 230 JB = KS, K
                  J = K + KS - JB
                  TEMP = A(J,K)
                  A(J,K) = A(KS,J)
                  A(KS,J) = TEMP
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  TEMP = A(KS,K+1)
                  A(KS,K+1) = A(K,K+1)
                  A(K,K+1) = TEMP
  240          CONTINUE
  250       CONTINUE
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DSPCO(AP,N,KPVT,RCOND,Z)
      INTEGER N,KPVT(1)
      REAL*8 AP(1),Z(1)
      REAL*8 RCOND
!
!     DSPCO FACTORS A REAL*8 SYMMETRIC MATRIX STORED IN
!     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING AND ESTIMATES
!     THE CONDITION OF THE MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DSPFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DSPCO BY DSPSL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPCO BY DSPSL.
!     TO COMPUTE  INVERSE(A) , FOLLOW DSPCO BY DSPDI.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPCO BY DSPDI.
!     TO COMPUTE  INERTIA(A), FOLLOW DSPCO BY DSPDI.
!
!     ON ENTRY
!
!        AP      REAL*8 (N*(N+1)/2)
!                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
!                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
!                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
!                SEE COMMENTS BELOW FOR DETAILS.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     OUTPUT
!
!        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
!                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
!
!        KPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        RCOND   REAL*8
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!
!        Z       REAL*8(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     PACKED STORAGE
!
!          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
!          TRIANGLE OF A SYMMETRIC MATRIX.
!
!                K = 0
!                DO 20 J = 1, N
!                   DO 10 I = 1, J
!                      K = K + 1
!                      AP(K) = A(I,J)
!             10    CONTINUE
!             20 CONTINUE
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DSPFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,IABS,DSIGN
!
!     INTERNAL VARIABLES
!
      REAL*8 AK,AKM1,BK,BKM1,DENOM,EK,T
      REAL*8 ANORM,S,YNORM
      INTEGER I,IJ,IK,IKM1,IKP1,INFO,J,JM1,J1
      INTEGER K,KK,KM1K,KM1KM1,KP,KPS,KS
!
!
!     FIND NORM OF A USING ONLY UPPER HALF
!
      J1 = 1
      DO 30 J = 1, N
         Z(J) = DASUM(J,AP(J1),1)
         IJ = J1
         J1 = J1 + J
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + DABS(AP(IJ))
            IJ = IJ + 1
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
!
!     FACTOR
!
      CALL DSPFA(AP,N,KPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!     SOLVE U*D*W = E
!
      EK = 1.0D0
      DO 50 J = 1, N
         Z(J) = 0.0D0
   50 CONTINUE
      K = N
      IK = (N*(N - 1))/2
   60 IF (K .EQ. 0) GO TO 120
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = IABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL DAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL DAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (DABS(Z(K)) .LE. DABS(AP(KK))) GO TO 90
               S = DABS(AP(KK))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   90       CONTINUE
            IF (AP(KK) .NE. 0.0D0) Z(K) = Z(K)/AP(KK)
            IF (AP(KK) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 110
  100    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/AP(KM1K)
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/AP(KM1K)
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 60
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!
!     SOLVE TRANS(U)*Y = W
!
      K = 1
      IK = 0
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + DDOT(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2) Z(K+1) = Z(K+1) + DDOT(K-1,AP(IKP1+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!
      YNORM = 1.0D0
!
!     SOLVE U*D*V = Y
!
      K = N
      IK = N*(N - 1)/2
  170 IF (K .EQ. 0) GO TO 230
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = IABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL DAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
            IF (KS .EQ. 2) CALL DAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (DABS(Z(K)) .LE. DABS(AP(KK))) GO TO 200
               S = DABS(AP(KK))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (AP(KK) .NE. 0.0D0) Z(K) = Z(K)/AP(KK)
            IF (AP(KK) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 220
  210    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/AP(KM1K)
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/AP(KM1K)
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 170
  230 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
!     SOLVE TRANS(U)*Z = V
!
      K = 1
      IK = 0
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + DDOT(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2) Z(K+1) = Z(K+1) + DDOT(K-1,AP(IKP1+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 240
  270 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
  END subroutine

  SUBROUTINE DSPFA(AP,N,KPVT,INFO)
      INTEGER N,KPVT(1),INFO
      REAL*8 AP(1)
!
!     DSPFA FACTORS A REAL*8 SYMMETRIC MATRIX STORED IN
!     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.
!
!     TO SOLVE  A*X = B , FOLLOW DSPFA BY DSPSL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPFA BY DSPSL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPFA BY DSPDI.
!     TO COMPUTE  INERTIA(A) , FOLLOW DSPFA BY DSPDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW DSPFA BY DSPDI.
!
!     ON ENTRY
!
!        AP      REAL*8 (N*(N+1)/2)
!                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
!                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
!                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
!                SEE COMMENTS BELOW FOR DETAILS.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     OUTPUT
!
!        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
!                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
!
!        KPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
!                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
!                     BUT IT DOES INDICATE THAT DSPSL OR DSPDI MAY
!                     DIVIDE BY ZERO IF CALLED.
!
!     PACKED STORAGE
!
!          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
!          TRIANGLE OF A SYMMETRIC MATRIX.
!
!                K = 0
!                DO 20 J = 1, N
!                   DO 10 I = 1, J
!                      K = K + 1
!                      AP(K)  = A(I,J)
!             10    CONTINUE
!             20 CONTINUE
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSWAP,IDAMAX
!     FORTRAN DABS,DMAX1,DSQRT
!
!     INTERNAL VARIABLES
!
      REAL*8 AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      REAL*8 ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
      LOGICAL SWAP
!
!
!     INITIALIZE
!
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0
!
      INFO = 0
!
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
!
      K = N
      IK = (N*(N - 1))/2
   10 CONTINUE
!
!        LEAVE THE LOOP IF K=0 OR K=1.
!
!     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (AP(1) .EQ. 0.0D0) INFO = 1
!     ......EXIT
            GO TO 200
   20    CONTINUE
!
!        THIS SECTION OF CODE DETERMINES THE KIND OF
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
!        REQUIRED.
!
         KM1 = K - 1
         KK = IK + K
         ABSAKK = DABS(AP(KK))
!
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!        COLUMN K.
!
         IMAX = IDAMAX(K-1,AP(IK+1),1)
         IMK = IK + IMAX
         COLMAX = DABS(AP(IMK))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
!
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!           ROW IMAX.
!
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            IM = IMAX*(IMAX - 1)/2
            IMJ = IM + 2*IMAX
            DO 40 J = IMAXP1, K
               ROWMAX = DMAX1(ROWMAX,DABS(AP(IMJ)))
               IMJ = IMJ + J
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = IDAMAX(IMAX-1,AP(IM+1),1)
               JMIM = JMAX + IM
               ROWMAX = DMAX1(ROWMAX,DABS(AP(JMIM)))
   50       CONTINUE
            IMIM = IMAX + IM
            IF (DABS(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100
!
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
!
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
!
!           1 X 1 PIVOT BLOCK.
!
            IF (.NOT.SWAP) GO TO 120
!
!              PERFORM AN INTERCHANGE.
!
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)
               IMJ = IK + IMAX
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  JK = IK + J
                  T = AP(JK)
                  AP(JK) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  110          CONTINUE
  120       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
            IJ = IK - (K - 1)
            DO 130 JJ = 1, KM1
               J = K - JJ
               JK = IK + J
               MULK = -AP(JK)/AP(KK)
               T = MULK
               CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
               IJJ = IJ + J
               AP(JK) = MULK
               IJ = IJ - (J - 1)
  130       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
            KM1K = IK + K - 1
            IKM1 = IK - (K - 1)
            IF (.NOT.SWAP) GO TO 160
!
!              PERFORM AN INTERCHANGE.
!
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)
               IMJ = IKM1 + IMAX
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  JKM1 = IKM1 + J
                  T = AP(JKM1)
                  AP(JKM1) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  150          CONTINUE
               T = AP(KM1K)
               AP(KM1K) = AP(IMK)
               AP(IMK) = T
  160       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = AP(KK)/AP(KM1K)
               KM1KM1 = IKM1 + K - 1
               AKM1 = AP(KM1KM1)/AP(KM1K)
               DENOM = 1.0D0 - AK*AKM1
               IJ = IK - (K - 1) - (K - 2)
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  JK = IK + J
                  BK = AP(JK)/AP(KM1K)
                  JKM1 = IKM1 + J
                  BKM1 = AP(JKM1)/AP(KM1K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
                  T = MULKM1
                  CALL DAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)
                  AP(JK) = MULK
                  AP(JKM1) = MULKM1
                  IJJ = IJ + J
                  IJ = IJ - (J - 1)
  170          CONTINUE
  180       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         IK = IK - (K - 1)
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DSPSL(AP,N,KPVT,B)
      INTEGER N,KPVT(1)
      REAL*8 AP(1),B(1)
!
!     DSISL SOLVES THE REAL*8 SYMMETRIC SYSTEM
!     A * X = B
!     USING THE FACTORS COMPUTED BY DSPFA.
!
!     ON ENTRY
!
!        AP      REAL*8(N*(N+1)/2)
!                THE OUTPUT FROM DSPFA.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        KPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DSPFA.
!
!        B       REAL*8(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO MAY OCCUR IF  DSPCO  HAS SET RCOND .EQ. 0.0
!        OR  DSPFA  HAS SET INFO .NE. 0  .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL DSPFA(AP,N,KPVT,INFO)
!           IF (INFO .NE. 0) GO TO ...
!           DO 10 J = 1, P
!              CALL DSPSL(AP,N,KPVT,C(1,J))
!        10 CONTINUE
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DDOT
!     FORTRAN IABS
!
!     INTERNAL VARIABLES.
!
      REAL*8 AK,AKM1,BK,BKM1,DENOM,TEMP
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
!
!     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
!     D INVERSE TO B.
!
      K = N
      IK = (N*(N - 1))/2
   10 IF (K .EQ. 0) GO TO 80
         KK = IK + K
         IF (KPVT(K) .LT. 0) GO TO 40
!
!           1 X 1 PIVOT BLOCK.
!
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
!
!                 INTERCHANGE.
!
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
!
!              APPLY THE TRANSFORMATION.
!
               CALL DAXPY(K-1,B(K),AP(IK+1),1,B(1),1)
   30       CONTINUE
!
!           APPLY D INVERSE.
!
            B(K) = B(K)/AP(KK)
            K = K - 1
            IK = IK - K
         GO TO 70
   40    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
            IKM1 = IK - (K - 1)
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
!
!                 INTERCHANGE.
!
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
!
!              APPLY THE TRANSFORMATION.
!
               CALL DAXPY(K-2,B(K),AP(IK+1),1,B(1),1)
               CALL DAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
   60       CONTINUE
!
!           APPLY D INVERSE.
!
            KM1K = IK + K - 1
            KK = IK + K
            AK = AP(KK)/AP(KM1K)
            KM1KM1 = IKM1 + K - 1
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = B(K)/AP(KM1K)
            BKM1 = B(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0D0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
            IK = IK - (K + 1) - K
   70    CONTINUE
      GO TO 10
   80 CONTINUE
!
!     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
!
      K = 1
      IK = 0
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
!
!           1 X 1 PIVOT BLOCK.
!
            IF (K .EQ. 1) GO TO 110
!
!              APPLY THE TRANSFORMATION.
!
               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
!
!                 INTERCHANGE.
!
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            IK = IK + K
            K = K + 1
         GO TO 150
  120    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
            IF (K .EQ. 1) GO TO 140
!
!              APPLY THE TRANSFORMATION.
!
               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
               IKP1 = IK + K
               B(K+1) = B(K+1) + DDOT(K-1,AP(IKP1+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
!
!                 INTERCHANGE.
!
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            IK = IK + K + K + 1
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DSPDI(AP,N,KPVT,DET,INERT,WORK,JOB)
      INTEGER N,JOB
      REAL*8 AP(1),WORK(1)
      REAL*8 DET(2)
      INTEGER KPVT(1),INERT(3)
!
!     DSPDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE
!     OF A REAL*8 SYMMETRIC MATRIX USING THE FACTORS FROM
!     DSPFA, WHERE THE MATRIX IS STORED IN PACKED FORM.
!
!     ON ENTRY
!
!        AP      REAL*8 (N*(N+1)/2)
!                THE OUTPUT FROM DSPFA.
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX A.
!
!        KPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DSPFA.
!
!        WORK    REAL*8(N)
!                WORK VECTOR.  CONTENTS IGNORED.
!
!        JOB     INTEGER
!                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE
!                   IF  C .NE. 0, THE INVERSE IS COMPUTED,
!                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,
!                   IF  A .NE. 0, THE INERTIA IS COMPUTED.
!
!                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.
!
!     ON RETURN
!
!        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
!
!        AP     CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
!               THE ORIGINAL MATRIX, STORED IN PACKED FORM.
!               THE COLUMNS OF THE UPPER TRIANGLE ARE STORED
!               SEQUENTIALLY IN A ONE-DIMENSIONAL ARRAY.
!
!        DET    REAL*8(2)
!               DETERMINANT OF ORIGINAL MATRIX.
!               DETERMINANT = DET(1) * 10.0**DET(2)
!               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0
!               OR DET(1) = 0.0.
!
!        INERT  INTEGER(3)
!               THE INERTIA OF THE ORIGINAL MATRIX.
!               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.
!               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.
!               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INVERSE IS REQUESTED
!        AND  DSPCO  HAS SET RCOND .EQ. 0.0
!        OR  DSPFA  HAS SET  INFO .NE. 0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DCOPY,DDOT,DSWAP
!     FORTRAN DABS,IABS,MOD
!
!     INTERNAL VARIABLES.
!
      REAL*8 AKKP1,TEMP
      REAL*8 TEN,D,T,AK,AKP1
      INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1
      INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP
      LOGICAL NOINV,NODET,NOERT
!
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
      NOERT = MOD(JOB,1000)/100 .EQ. 0
!
      IF (NODET .AND. NOERT) GO TO 140
         IF (NOERT) GO TO 10
            INERT(1) = 0
            INERT(2) = 0
            INERT(3) = 0
   10    CONTINUE
         IF (NODET) GO TO 20
            DET(1) = 1.0D0
            DET(2) = 0.0D0
            TEN = 10.0D0
   20    CONTINUE
         T = 0.0D0
         IK = 0
         DO 130 K = 1, N
            KK = IK + K
            D = AP(KK)
!
!           CHECK IF 1 BY 1
!
            IF (KPVT(K) .GT. 0) GO TO 50
!
!              2 BY 2 BLOCK
!              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)
!                      (S  C)
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
!
               IF (T .NE. 0.0D0) GO TO 30
                  IKP1 = IK + K
                  KKP1 = IKP1 + K
                  T = DABS(AP(KKP1))
                  D = (D/T)*AP(KKP1+1) - T
               GO TO 40
   30          CONTINUE
                  D = T
                  T = 0.0D0
   40          CONTINUE
   50       CONTINUE
!
            IF (NOERT) GO TO 60
               IF (D .GT. 0.0D0) INERT(1) = INERT(1) + 1
               IF (D .LT. 0.0D0) INERT(2) = INERT(2) + 1
               IF (D .EQ. 0.0D0) INERT(3) = INERT(3) + 1
   60       CONTINUE
!
            IF (NODET) GO TO 120
               DET(1) = D*DET(1)
               IF (DET(1) .EQ. 0.0D0) GO TO 110
   70             IF (DABS(DET(1)) .GE. 1.0D0) GO TO 80
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - 1.0D0
                  GO TO 70
   80             CONTINUE
   90             IF (DABS(DET(1)) .LT. TEN) GO TO 100
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + 1.0D0
                  GO TO 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
            IK = IK + K
  130    CONTINUE
  140 CONTINUE
!
!     COMPUTE INVERSE(A)
!
      IF (NOINV) GO TO 270
         K = 1
         IK = 0
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            KK = IK + K
            IKP1 = IK + K
            KKP1 = IKP1 + K
            IF (KPVT(K) .LT. 0) GO TO 180
!
!              1 BY 1
!
               AP(KK) = 1.0D0/AP(KK)
               IF (KM1 .LT. 1) GO TO 170
                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 160 J = 1, KM1
                     JK = IK + J
                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  160             CONTINUE
                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
!
!              2 BY 2
!
               T = DABS(AP(KKP1))
               AK = AP(KK)/T
               AKP1 = AP(KKP1+1)/T
               AKKP1 = AP(KKP1)/T
               D = T*(AK*AKP1 - 1.0D0)
               AP(KK) = AKP1/D
               AP(KKP1+1) = AK/D
               AP(KKP1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL DCOPY(KM1,AP(IKP1+1),1,WORK,1)
                  IJ = 0
                  DO 190 J = 1, KM1
                     JKP1 = IKP1 + J
                     AP(JKP1) = DDOT(J,AP(IJ+1),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)
                     IJ = IJ + J
  190             CONTINUE
                  AP(KKP1+1) = AP(KKP1+1)+ DDOT(KM1,WORK,1,AP(IKP1+1),1)
                  AP(KKP1) = AP(KKP1) + DDOT(KM1,AP(IK+1),1,AP(IKP1+1),1)
                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 200 J = 1, KM1
                     JK = IK + J
                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  200             CONTINUE
                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
!
!           SWAP
!
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               IKS = (KS*(KS - 1))/2
               CALL DSWAP(KS,AP(IKS+1),1,AP(IK+1),1)
               KSJ = IK + KS
               DO 230 JB = KS, K
                  J = K + KS - JB
                  JK = IK + J
                  TEMP = AP(JK)
                  AP(JK) = AP(KSJ)
                  AP(KSJ) = TEMP
                  KSJ = KSJ - (J - 1)
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  KSKP1 = IKP1 + KS
                  TEMP = AP(KSKP1)
                  AP(KSKP1) = AP(KKP1)
                  AP(KKP1) = TEMP
  240          CONTINUE
  250       CONTINUE
            IK = IK + K
            IF (KSTEP .EQ. 2) IK = IK + K + 1
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DTRCO(T,LDT,N,RCOND,Z,JOB)
      INTEGER LDT,N,JOB
      REAL*8 T(LDT,1),Z(1)
      REAL*8 RCOND
!
!     DTRCO ESTIMATES THE CONDITION OF A REAL*8 TRIANGULAR
!     MATRIX.
!
!     ON ENTRY
!
!        T       REAL*8(LDT,N)
!                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO
!                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
!                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
!                USED TO STORE OTHER INFORMATION.
!
!        LDT     INTEGER
!                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
!
!        N       INTEGER
!                N IS THE ORDER OF THE SYSTEM.
!
!        JOB     INTEGER
!                = 0         T  IS LOWER TRIANGULAR.
!                = NONZERO   T  IS UPPER TRIANGULAR.
!
!     ON RETURN
!
!        RCOND   REAL*8
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  T .
!                FOR THE SYSTEM  T*X = B , RELATIVE PERTURBATIONS
!                IN  T  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  T  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!
!        Z       REAL*8(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  T  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,DSIGN
!
!     INTERNAL VARIABLES
!
      REAL*8 W,WK,WKM,EK
      REAL*8 TNORM,YNORM,S,SM
      INTEGER I1,J,J1,J2,K,KK,L
      LOGICAL LOWER
!
      LOWER = JOB .EQ. 0
!
!     COMPUTE 1-NORM OF T
!
      TNORM = 0.0D0
      DO 10 J = 1, N
         L = J
         IF (LOWER) L = N + 1 - J
         I1 = 1
         IF (LOWER) I1 = J
         TNORM = DMAX1(TNORM,DASUM(L,T(I1,J),1))
   10 CONTINUE
!
!     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .
!     TRANS(T)  IS THE TRANSPOSE OF T .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF Y .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!     SOLVE TRANS(T)*Y = E
!
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 KK = 1, N
         K = KK
         IF (LOWER) K = N + 1 - KK
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(T(K,K))) GO TO 30
            S = DABS(T(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (T(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/T(K,K)
            WKM = WKM/T(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         IF (KK .EQ. N) GO TO 90
            J1 = K + 1
            IF (LOWER) J1 = 1
            J2 = N
            IF (LOWER) J2 = K - 1
            DO 60 J = J1, J2
               SM = SM + DABS(Z(J)+WKM*T(K,J))
               Z(J) = Z(J) + WK*T(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               W = WKM - WK
               WK = WKM
               DO 70 J = J1, J2
                  Z(J) = Z(J) + W*T(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!
      YNORM = 1.0D0
!
!     SOLVE T*Z = Y
!
      DO 130 KK = 1, N
         K = N + 1 - KK
         IF (LOWER) K = KK
         IF (DABS(Z(K)) .LE. DABS(T(K,K))) GO TO 110
            S = DABS(T(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  110    CONTINUE
         IF (T(K,K) .NE. 0.0D0) Z(K) = Z(K)/T(K,K)
         IF (T(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         I1 = 1
         IF (LOWER) I1 = K + 1
         IF (KK .GE. N) GO TO 120
            W = -Z(K)
            CALL DAXPY(N-KK,W,T(I1,K),1,Z(I1),1)
  120    CONTINUE
  130 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
      IF (TNORM .NE. 0.0D0) RCOND = YNORM/TNORM
      IF (TNORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
  END subroutine

  SUBROUTINE DTRSL(T,LDT,N,B,JOB,INFO)
      INTEGER LDT,N,JOB,INFO
      REAL*8 T(LDT,1),B(1)
!
!
!     DTRSL SOLVES SYSTEMS OF THE FORM
!
!                   T * X = B
!     OR
!                   TRANS(T) * X = B
!
!     WHERE T IS A TRIANGULAR MATRIX OF ORDER N. HERE TRANS(T)
!     DENOTES THE TRANSPOSE OF THE MATRIX T.
!
!     ON ENTRY
!
!         T         REAL*8(LDT,N)
!                   T CONTAINS THE MATRIX OF THE SYSTEM. THE ZERO
!                   ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
!                   THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
!                   USED TO STORE OTHER INFORMATION.
!
!         LDT       INTEGER
!                   LDT IS THE LEADING DIMENSION OF THE ARRAY T.
!
!         N         INTEGER
!                   N IS THE ORDER OF THE SYSTEM.
!
!         B         REAL*8(N).
!                   B CONTAINS THE RIGHT HAND SIDE OF THE SYSTEM.
!
!         JOB       INTEGER
!                   JOB SPECIFIES WHAT KIND OF SYSTEM IS TO BE SOLVED.
!                   IF JOB IS
!
!                        00   SOLVE T*X=B, T LOWER TRIANGULAR,
!                        01   SOLVE T*X=B, T UPPER TRIANGULAR,
!                        10   SOLVE TRANS(T)*X=B, T LOWER TRIANGULAR,
!                        11   SOLVE TRANS(T)*X=B, T UPPER TRIANGULAR.
!
!     ON RETURN
!
!         B         B CONTAINS THE SOLUTION, IF INFO .EQ. 0.
!                   OTHERWISE B IS UNALTERED.
!
!         INFO      INTEGER
!                   INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR.
!                   OTHERWISE INFO CONTAINS THE INDEX OF
!                   THE FIRST ZERO DIAGONAL ELEMENT OF T.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DDOT
!     FORTRAN MOD
!
!     INTERNAL VARIABLES
!
      REAL*8 TEMP
      INTEGER CASE,J,JJ
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 150
!
!        CHECK FOR ZERO DIAGONAL ELEMENTS.
!
         DO 10 INFO = 1, N
!     ......EXIT
            IF (T(INFO,INFO) .EQ. 0.0D0) GO TO 150
   10    CONTINUE
         INFO = 0
!
!        DETERMINE THE TASK AND GO TO IT.
!
         CASE = 1
         IF (MOD(JOB,10) .NE. 0) CASE = 2
         IF (MOD(JOB,100)/10 .NE. 0) CASE = CASE + 2
         GO TO (20,50,80,110), CASE
!
!        SOLVE T*X=B FOR T LOWER TRIANGULAR
!
   20    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 40
            DO 30 J = 2, N
               TEMP = -B(J-1)
               CALL DAXPY(N-J+1,TEMP,T(J,J-1),1,B(J),1)
               B(J) = B(J)/T(J,J)
   30       CONTINUE
   40       CONTINUE
         GO TO 140
!
!        SOLVE T*X=B FOR T UPPER TRIANGULAR.
!
   50    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 70
            DO 60 JJ = 2, N
               J = N - JJ + 1
               TEMP = -B(J+1)
               CALL DAXPY(J,TEMP,T(1,J+1),1,B(1),1)
               B(J) = B(J)/T(J,J)
   60       CONTINUE
   70       CONTINUE
         GO TO 140
!
!        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR.
!
   80    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 100
            DO 90 JJ = 2, N
               J = N - JJ + 1
               B(J) = B(J) - DDOT(JJ-1,T(J+1,J),1,B(J+1),1)
               B(J) = B(J)/T(J,J)
   90       CONTINUE
  100       CONTINUE
         GO TO 140
!
!        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR.
!
  110    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 130
            DO 120 J = 2, N
               B(J) = B(J) - DDOT(J-1,T(1,J),1,B(1),1)
               B(J) = B(J)/T(J,J)
  120       CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DTRDI(T,LDT,N,DET,JOB,INFO)
      INTEGER LDT,N,JOB,INFO
      REAL*8 T(LDT,1),DET(2)
!
!     DTRDI COMPUTES THE DETERMINANT AND INVERSE OF A REAL*8
!     TRIANGULAR MATRIX.
!
!     ON ENTRY
!
!        T       REAL*8(LDT,N)
!                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO
!                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
!                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
!                USED TO STORE OTHER INFORMATION.
!
!        LDT     INTEGER
!                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
!
!        N       INTEGER
!                N IS THE ORDER OF THE SYSTEM.
!
!        JOB     INTEGER
!                = 010       NO DET, INVERSE OF LOWER TRIANGULAR.
!                = 011       NO DET, INVERSE OF UPPER TRIANGULAR.
!                = 100       DET, NO INVERSE.
!                = 110       DET, INVERSE OF LOWER TRIANGULAR.
!                = 111       DET, INVERSE OF UPPER TRIANGULAR.
!
!     ON RETURN
!
!        T       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE UNCHANGED.
!
!        DET     REAL*8(2)
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!        INFO    INTEGER
!                INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR
!                AND THE INVERSE IS REQUESTED.
!                OTHERWISE INFO CONTAINS THE INDEX OF
!                A ZERO DIAGONAL ELEMENT OF T.
!
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL
!     FORTRAN DABS,MOD
!
!     INTERNAL VARIABLES
!
      REAL*8 TEMP
      REAL*8 TEN
      INTEGER I,J,K,KB,KM1,KP1
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 180
!
!        COMPUTE DETERMINANT
!
         IF (JOB/100 .EQ. 0) GO TO 70
            DET(1) = 1.0D0
            DET(2) = 0.0D0
            TEN = 10.0D0
            DO 50 I = 1, N
               DET(1) = T(I,I)*DET(1)
!           ...EXIT
               IF (DET(1) .EQ. 0.0D0) GO TO 60
   10          IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
                  DET(1) = TEN*DET(1)
                  DET(2) = DET(2) - 1.0D0
               GO TO 10
   20          CONTINUE
   30          IF (DABS(DET(1)) .LT. TEN) GO TO 40
                  DET(1) = DET(1)/TEN
                  DET(2) = DET(2) + 1.0D0
               GO TO 30
   40          CONTINUE
   50       CONTINUE
   60       CONTINUE
   70    CONTINUE
!
!        COMPUTE INVERSE OF UPPER TRIANGULAR
!
         IF (MOD(JOB/10,10) .EQ. 0) GO TO 170
            IF (MOD(JOB,10) .EQ. 0) GO TO 120
!              BEGIN BLOCK PERMITTING ...EXITS TO 110
                  DO 100 K = 1, N
                     INFO = K
!              ......EXIT
                     IF (T(K,K) .EQ. 0.0D0) GO TO 110
                     T(K,K) = 1.0D0/T(K,K)
                     TEMP = -T(K,K)
                     CALL DSCAL(K-1,TEMP,T(1,K),1)
                     KP1 = K + 1
                     IF (N .LT. KP1) GO TO 90
                     DO 80 J = KP1, N
                        TEMP = T(K,J)
                        T(K,J) = 0.0D0
                        CALL DAXPY(K,TEMP,T(1,K),1,T(1,J),1)
   80                CONTINUE
   90                CONTINUE
  100             CONTINUE
                  INFO = 0
  110          CONTINUE
            GO TO 160
  120       CONTINUE
!
!              COMPUTE INVERSE OF LOWER TRIANGULAR
!
               DO 150 KB = 1, N
                  K = N + 1 - KB
                  INFO = K
!     ............EXIT
                  IF (T(K,K) .EQ. 0.0D0) GO TO 180
                  T(K,K) = 1.0D0/T(K,K)
                  TEMP = -T(K,K)
                  IF (K .NE. N) CALL DSCAL(N-K,TEMP,T(K+1,K),1)
                  KM1 = K - 1
                  IF (KM1 .LT. 1) GO TO 140
                  DO 130 J = 1, KM1
                     TEMP = T(K,J)
                     T(K,J) = 0.0D0
                     CALL DAXPY(N-K+1,TEMP,T(K,K),1,T(K,J),1)
  130             CONTINUE
  140             CONTINUE
  150          CONTINUE
               INFO = 0
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DGTSL(N,C,D,E,B,INFO)
      INTEGER N,INFO
      REAL*8 C(1),D(1),E(1),B(1)
!
!     DGTSL GIVEN A GENERAL TRIDIAGONAL MATRIX AND A RIGHT HAND
!     SIDE WILL FIND THE SOLUTION.
!
!     ON ENTRY
!
!        N       INTEGER
!                IS THE ORDER OF THE TRIDIAGONAL MATRIX.
!
!        C       REAL*8(N)
!                IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX.
!                C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL.
!                ON OUTPUT C IS DESTROYED.
!
!        D       REAL*8(N)
!                IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
!                ON OUTPUT D IS DESTROYED.
!
!        E       REAL*8(N)
!                IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX.
!                E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL.
!                ON OUTPUT E IS DESTROYED.
!
!        B       REAL*8(N)
!                IS THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B       IS THE SOLUTION VECTOR.
!
!        INFO    INTEGER
!                = 0 NORMAL VALUE.
!                = K IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES
!                    EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN
!                    THIS IS DETECTED.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
!
!     NO EXTERNALS
!     FORTRAN DABS
!
!     INTERNAL VARIABLES
!
      INTEGER K,KB,KP1,NM1,NM2
      REAL*8 T
!     BEGIN BLOCK PERMITTING ...EXITS TO 100
!
         INFO = 0
         C(1) = D(1)
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 40
            D(1) = E(1)
            E(1) = 0.0D0
            E(N) = 0.0D0
!
            DO 30 K = 1, NM1
               KP1 = K + 1
!
!              FIND THE LARGEST OF THE TWO ROWS
!
               IF (DABS(C(KP1)) .LT. DABS(C(K))) GO TO 10
!
!                 INTERCHANGE ROW
!
                  T = C(KP1)
                  C(KP1) = C(K)
                  C(K) = T
                  T = D(KP1)
                  D(KP1) = D(K)
                  D(K) = T
                  T = E(KP1)
                  E(KP1) = E(K)
                  E(K) = T
                  T = B(KP1)
                  B(KP1) = B(K)
                  B(K) = T
   10          CONTINUE
!
!              ZERO ELEMENTS
!
               IF (C(K) .NE. 0.0D0) GO TO 20
                  INFO = K
!     ............EXIT
                  GO TO 100
   20          CONTINUE
               T = -C(KP1)/C(K)
               C(KP1) = D(KP1) + T*D(K)
               D(KP1) = E(KP1) + T*E(K)
               E(KP1) = 0.0D0
               B(KP1) = B(KP1) + T*B(K)
   30       CONTINUE
   40    CONTINUE
         IF (C(N) .NE. 0.0D0) GO TO 50
            INFO = N
         GO TO 90
   50    CONTINUE
!
!           BACK SOLVE
!
            NM2 = N - 2
            B(N) = B(N)/C(N)
            IF (N .EQ. 1) GO TO 80
               B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)
               IF (NM2 .LT. 1) GO TO 70
               DO 60 KB = 1, NM2
                  K = NM2 - KB + 1
                  B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)
   60          CONTINUE
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
!
      RETURN
  END subroutine

  SUBROUTINE DPTSL(N,D,E,B)
      INTEGER N
      REAL*8 D(1400),E(1400),B(1400)
!
!     DPTSL GIVEN A POSITIVE DEFINITE TRIDIAGONAL MATRIX AND A RIGHT
!     HAND SIDE WILL FIND THE SOLUTION.
!
!     ON ENTRY
!
!        N        INTEGER
!                 IS THE ORDER OF THE TRIDIAGONAL MATRIX.
!
!        D        REAL*8(N)
!                 IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
!                 ON OUTPUT D IS DESTROYED.
!
!        E        REAL*8(N)
!                 IS THE OFFDIAGONAL OF THE TRIDIAGONAL MATRIX.
!                 E(1) THROUGH E(N-1) SHOULD CONTAIN THE
!                 OFFDIAGONAL.
!
!        B        REAL*8(N)
!                 IS THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B        CONTAINS THE SOLUTION.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
!
!     NO EXTERNALS
!     FORTRAN MOD
!
!     INTERNAL VARIABLES
!
      INTEGER K,KBM1,KE,KF,KP1,NM1,NM1D2
      REAL*8 T1,T2
!
!     CHECK FOR 1 X 1 CASE
!
      IF (N .NE. 1) GO TO 10
         B(1) = B(1)/D(1)
      GO TO 70
   10 CONTINUE
         NM1 = N - 1
         NM1D2 = NM1/2
         IF (N .EQ. 2) GO TO 30
            KBM1 = N - 1
!
!           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF
!           SUPERDIAGONAL
!
            DO 20 K = 1, NM1D2
               T1 = E(K)/D(K)
               D(K+1) = D(K+1) - T1*E(K)
               B(K+1) = B(K+1) - T1*B(K)
               T2 = E(KBM1)/D(KBM1+1)
               D(KBM1) = D(KBM1) - T2*E(KBM1)
               B(KBM1) = B(KBM1) - T2*B(KBM1+1)
               KBM1 = KBM1 - 1
   20       CONTINUE
   30    CONTINUE
         KP1 = NM1D2 + 1
!
!        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER
!
         IF (MOD(N,2) .NE. 0) GO TO 40
            T1 = E(KP1)/D(KP1)
            D(KP1+1) = D(KP1+1) - T1*E(KP1)
            B(KP1+1) = B(KP1+1) - T1*B(KP1)
            KP1 = KP1 + 1
   40    CONTINUE
!
!        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP
!        AND BOTTOM
!
         B(KP1) = B(KP1)/D(KP1)
         IF (N .EQ. 2) GO TO 60
            K = KP1 - 1
            KE = KP1 + NM1D2 - 1
            DO 50 KF = KP1, KE
               B(K) = (B(K) - E(K)*B(K+1))/D(K)
               B(KF+1) = (B(KF+1) - E(KF)*B(KF))/D(KF+1)
               K = K - 1
   50       CONTINUE
   60    CONTINUE
         IF (MOD(N,2) .EQ. 0) B(1) = (B(1) - E(1)*B(2))/D(1)
   70 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DCHDC(A,LDA,P,WORK,JPVT,JOB,INFO)
      INTEGER LDA,P,JPVT(1),JOB,INFO
      REAL*8 A(LDA,1),WORK(1)
!
!     DCHDC COMPUTES THE CHOLESKY DECOMPOSITION OF A POSITIVE DEFINITE
!     MATRIX.  A PIVOTING OPTION ALLOWS THE USER TO ESTIMATE THE
!     CONDITION OF A POSITIVE DEFINITE MATRIX OR DETERMINE THE RANK
!     OF A POSITIVE SEMIDEFINITE MATRIX.
!
!     ON ENTRY
!
!         A      REAL*8(LDA,P).
!                A CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO
!                BE COMPUTED.  ONLT THE UPPER HALF OF A NEED BE STORED.
!                THE LOWER PART OF THE ARRAY A IS NOT REFERENCED.
!
!         LDA    INTEGER.
!                LDA IS THE LEADING DIMENSION OF THE ARRAY A.
!
!         P      INTEGER.
!                P IS THE ORDER OF THE MATRIX.
!
!         WORK   REAL*8.
!                WORK IS A WORK ARRAY.
!
!         JPVT   INTEGER(P).
!                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
!                OF THE PIVOT ELEMENTS, IF PIVOTING HAS BEEN REQUESTED.
!                EACH DIAGONAL ELEMENT A(K,K)
!                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
!                VALUE OF JPVT(K).
!
!                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
!                                      ELEMENT.
!
!                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE ELEMENT.
!
!                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL ELEMENT.
!
!                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL ELEMENTS
!                ARE MOVED BY SYMMETRIC ROW AND COLUMN INTERCHANGES TO
!                THE BEGINNING OF THE ARRAY A AND FINAL
!                ELEMENTS TO THE END.  BOTH INITIAL AND FINAL ELEMENTS
!                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
!                FREE ELEMENTS ARE MOVED.  AT THE K-TH STAGE OF THE
!                REDUCTION, IF A(K,K) IS OCCUPIED BY A FREE ELEMENT
!                IT IS INTERCHANGED WITH THE LARGEST FREE ELEMENT
!                A(L,L) WITH L .GE. K.  JPVT IS NOT REFERENCED IF
!                JOB .EQ. 0.
!
!        JOB     INTEGER.
!                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
!                IF JOB .EQ. 0, NO PIVOTING IS DONE.
!                IF JOB .NE. 0, PIVOTING IS DONE.
!
!     ON RETURN
!
!         A      A CONTAINS IN ITS UPPER HALF THE CHOLESKY FACTOR
!                OF THE MATRIX A AS IT HAS BEEN PERMUTED BY PIVOTING.
!
!         JPVT   JPVT(J) CONTAINS THE INDEX OF THE DIAGONAL ELEMENT
!                OF A THAT WAS MOVED INTO THE J-TH POSITION,
!                PROVIDED PIVOTING WAS REQUESTED.
!
!         INFO   CONTAINS THE INDEX OF THE LAST POSITIVE DIAGONAL
!                ELEMENT OF THE CHOLESKY FACTOR.
!
!     FOR POSITIVE DEFINITE MATRICES INFO = P IS THE NORMAL RETURN.
!     FOR PIVOTING WITH POSITIVE SEMIDEFINITE MATRICES INFO WILL
!     IN GENERAL BE LESS THAN P.  HOWEVER, INFO MAY BE GREATER THAN
!     THE RANK OF A, SINCE ROUNDING ERROR CAN CAUSE AN OTHERWISE ZERO
!     ELEMENT TO BE POSITIVE. INDEFINITE SYSTEMS WILL ALWAYS CAUSE
!     INFO TO BE LESS THAN P.
!
!     LINPACK. THIS VERSION DATED 03/19/79 .
!     J.J. DONGARRA AND G.W. STEWART, ARGONNE NATIONAL LABORATORY AND
!     UNIVERSITY OF MARYLAND.
!
!
!     BLAS DAXPY,DSWAP
!     FORTRAN DSQRT
!
!     INTERNAL VARIABLES
!
      INTEGER PU,PL,PLP1,I,J,JP,JT,K,KB,KM1,KP1,L,MAXL
      REAL*8 TEMP
      REAL*8 MAXDIA
      LOGICAL SWAPK,NEGK
!
      PL = 1
      PU = 0
      INFO = P
      IF (JOB .EQ. 0) GO TO 160
!
!        PIVOTING HAS BEEN REQUESTED. REARRANGE THE
!        THE ELEMENTS ACCORDING TO JPVT.
!
         DO 70 K = 1, P
            SWAPK = JPVT(K) .GT. 0
            NEGK = JPVT(K) .LT. 0
            JPVT(K) = K
            IF (NEGK) JPVT(K) = -JPVT(K)
            IF (.NOT.SWAPK) GO TO 60
               IF (K .EQ. PL) GO TO 50
                  CALL DSWAP(PL-1,A(1,K),1,A(1,PL),1)
                  TEMP = A(K,K)
                  A(K,K) = A(PL,PL)
                  A(PL,PL) = TEMP
                  PLP1 = PL + 1
                  IF (P .LT. PLP1) GO TO 40
                  DO 30 J = PLP1, P
                     IF (J .GE. K) GO TO 10
                        TEMP = A(PL,J)
                        A(PL,J) = A(J,K)
                        A(J,K) = TEMP
                     GO TO 20
   10                CONTINUE
                     IF (J .EQ. K) GO TO 20
                        TEMP = A(K,J)
                        A(K,J) = A(PL,J)
                        A(PL,J) = TEMP
   20                CONTINUE
   30             CONTINUE
   40             CONTINUE
                  JPVT(K) = JPVT(PL)
                  JPVT(PL) = K
   50          CONTINUE
               PL = PL + 1
   60       CONTINUE
   70    CONTINUE
         PU = P
         IF (P .LT. PL) GO TO 150
         DO 140 KB = PL, P
            K = P - KB + PL
            IF (JPVT(K) .GE. 0) GO TO 130
               JPVT(K) = -JPVT(K)
               IF (PU .EQ. K) GO TO 120
                  CALL DSWAP(K-1,A(1,K),1,A(1,PU),1)
                  TEMP = A(K,K)
                  A(K,K) = A(PU,PU)
                  A(PU,PU) = TEMP
                  KP1 = K + 1
                  IF (P .LT. KP1) GO TO 110
                  DO 100 J = KP1, P
                     IF (J .GE. PU) GO TO 80
                        TEMP = A(K,J)
                        A(K,J) = A(J,PU)
                        A(J,PU) = TEMP
                     GO TO 90
   80                CONTINUE
                     IF (J .EQ. PU) GO TO 90
                        TEMP = A(K,J)
                        A(K,J) = A(PU,J)
                        A(PU,J) = TEMP
   90                CONTINUE
  100             CONTINUE
  110             CONTINUE
                  JT = JPVT(K)
                  JPVT(K) = JPVT(PU)
                  JPVT(PU) = JT
  120          CONTINUE
               PU = PU - 1
  130       CONTINUE
  140    CONTINUE
  150    CONTINUE
  160 CONTINUE
      DO 270 K = 1, P
!
!        REDUCTION LOOP.
!
         MAXDIA = A(K,K)
         KP1 = K + 1
         MAXL = K
!
!        DETERMINE THE PIVOT ELEMENT.
!
         IF (K .LT. PL .OR. K .GE. PU) GO TO 190
            DO 180 L = KP1, PU
               IF (A(L,L) .LE. MAXDIA) GO TO 170
                  MAXDIA = A(L,L)
                  MAXL = L
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
!
!        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE.
!
         IF (MAXDIA .GT. 0.0D0) GO TO 200
            INFO = K - 1
!     ......EXIT
            GO TO 280
  200    CONTINUE
         IF (K .EQ. MAXL) GO TO 210
!
!           START THE PIVOTING AND UPDATE JPVT.
!
            KM1 = K - 1
            CALL DSWAP(KM1,A(1,K),1,A(1,MAXL),1)
            A(MAXL,MAXL) = A(K,K)
            A(K,K) = MAXDIA
            JP = JPVT(MAXL)
            JPVT(MAXL) = JPVT(K)
            JPVT(K) = JP
  210    CONTINUE
!
!        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS.
!
         WORK(K) = DSQRT(A(K,K))
         A(K,K) = WORK(K)
         IF (P .LT. KP1) GO TO 260
         DO 250 J = KP1, P
            IF (K .EQ. MAXL) GO TO 240
               IF (J .GE. MAXL) GO TO 220
                  TEMP = A(K,J)
                  A(K,J) = A(J,MAXL)
                  A(J,MAXL) = TEMP
               GO TO 230
  220          CONTINUE
               IF (J .EQ. MAXL) GO TO 230
                  TEMP = A(K,J)
                  A(K,J) = A(MAXL,J)
                  A(MAXL,J) = TEMP
  230          CONTINUE
  240       CONTINUE
            A(K,J) = A(K,J)/WORK(K)
            WORK(J) = A(K,J)
            TEMP = -A(K,J)
            CALL DAXPY(J-K,TEMP,WORK(KP1),1,A(KP1,J),1)
  250    CONTINUE
  260    CONTINUE
  270 CONTINUE
  280 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DCHUD(R,LDR,P,X,Z,LDZ,NZ,Y,RHO,C,S)
      INTEGER LDR,P,LDZ,NZ
      REAL*8 RHO(1),C(1)
      REAL*8 R(LDR,1),X(1),Z(LDZ,1),Y(1),S(1)
!
!     DCHUD UPDATES AN AUGMENTED CHOLESKY DECOMPOSITION OF THE
!     TRIANGULAR PART OF AN AUGMENTED QR DECOMPOSITION.  SPECIFICALLY,
!     GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P, A ROW VECTOR
!     X, A COLUMN VECTOR Z, AND A SCALAR Y, DCHUD DETERMINES A
!     UNTIARY MATRIX U AND A SCALAR ZETA SUCH THAT
!
!
!                              (R  Z)     (RR   ZZ )
!                         U  * (    )  =  (        ) ,
!                              (X  Y)     ( 0  ZETA)
!
!     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN
!     OBTAINED FROM THE FACTORIZATION OF A LEAST SQUARES
!     PROBLEM, THEN RR AND ZZ ARE THE FACTORS CORRESPONDING TO
!     THE PROBLEM WITH THE OBSERVATION (X,Y) APPENDED.  IN THIS
!     CASE, IF RHO IS THE NORM OF THE RESIDUAL VECTOR, THEN THE
!     NORM OF THE RESIDUAL VECTOR OF THE UPDATED PROBLEM IS
!     DSQRT(RHO**2 + ZETA**2).  DCHUD WILL SIMULTANEOUSLY UPDATE
!     SEVERAL TRIPLETS (Z,Y,RHO).
!     FOR A LESS TERSE DESCRIPTION OF WHAT DCHUD DOES AND HOW
!     IT MAY BE APPLIED, SEE THE LINPACK GUIDE.
!
!     THE MATRIX U IS DETERMINED AS THE PRODUCT U(P)*...*U(1),
!     WHERE U(I) IS A ROTATION IN THE (I,P+1) PLANE OF THE
!     FORM
!
!                       (     C(I)      S(I) )
!                       (                    ) .
!                       (    -S(I)      C(I) )
!
!     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS REAL*8.
!
!     ON ENTRY
!
!         R      REAL*8(LDR,P), WHERE LDR .GE. P.
!                R CONTAINS THE UPPER TRIANGULAR MATRIX
!                THAT IS TO BE UPDATED.  THE PART OF R
!                BELOW THE DIAGONAL IS NOT REFERENCED.
!
!         LDR    INTEGER.
!                LDR IS THE LEADING DIMENSION OF THE ARRAY R.
!
!         P      INTEGER.
!                P IS THE ORDER OF THE MATRIX R.
!
!         X      REAL*8(P).
!                X CONTAINS THE ROW TO BE ADDED TO R.  X IS
!                NOT ALTERED BY DCHUD.
!
!         Z      REAL*8(LDZ,NZ), WHERE LDZ .GE. P.
!                Z IS AN ARRAY CONTAINING NZ P-VECTORS TO
!                BE UPDATED WITH R.
!
!         LDZ    INTEGER.
!                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z.
!
!         NZ     INTEGER.
!                NZ IS THE NUMBER OF VECTORS TO BE UPDATED
!                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO
!                ARE NOT REFERENCED.
!
!         Y      REAL*8(NZ).
!                Y CONTAINS THE SCALARS FOR UPDATING THE VECTORS
!                Z.  Y IS NOT ALTERED BY DCHUD.
!
!         RHO    REAL*8(NZ).
!                RHO CONTAINS THE NORMS OF THE RESIDUAL
!                VECTORS THAT ARE TO BE UPDATED.  IF RHO(J)
!                IS NEGATIVE, IT IS LEFT UNALTERED.
!
!     ON RETURN
!
!         RC
!         RHO    CONTAIN THE UPDATED QUANTITIES.
!         Z
!
!         C      REAL*8(P).
!                C CONTAINS THE COSINES OF THE TRANSFORMING
!                ROTATIONS.
!
!         S      REAL*8(P).
!                S CONTAINS THE SINES OF THE TRANSFORMING
!                ROTATIONS.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     DCHUD USES THE FOLLOWING FUNCTIONS AND SUBROUTINES.
!
!     EXTENDED BLAS DROTG
!     FORTRAN DSQRT
!
      INTEGER I,J,JM1
      REAL*8 AZETA,SCALE
      REAL*8 T,XJ,ZETA
!
!     UPDATE R.
!
      DO 30 J = 1, P
         XJ = X(J)
!
!        APPLY THE PREVIOUS ROTATIONS.
!
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            T = C(I)*R(I,J) + S(I)*XJ
            XJ = C(I)*XJ - S(I)*R(I,J)
            R(I,J) = T
   10    CONTINUE
   20    CONTINUE
!
!        COMPUTE THE NEXT ROTATION.
!
         CALL DROTG(R(J,J),XJ,C(J),S(J))
   30 CONTINUE
!
!     IF REQUIRED, UPDATE Z AND RHO.
!
      IF (NZ .LT. 1) GO TO 70
      DO 60 J = 1, NZ
         ZETA = Y(J)
         DO 40 I = 1, P
            T = C(I)*Z(I,J) + S(I)*ZETA
            ZETA = C(I)*ZETA - S(I)*Z(I,J)
            Z(I,J) = T
   40    CONTINUE
         AZETA = DABS(ZETA)
         IF (AZETA .EQ. 0.0D0 .OR. RHO(J) .LT. 0.0D0) GO TO 50
            SCALE = AZETA + RHO(J)
            RHO(J) = SCALE*DSQRT((AZETA/SCALE)**2+(RHO(J)/SCALE)**2)
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DCHDD(R,LDR,P,X,Z,LDZ,NZ,Y,RHO,C,S,INFO)
      INTEGER LDR,P,LDZ,NZ,INFO
      REAL*8 R(LDR,1),X(1),Z(LDZ,1),Y(1),S(1)
      REAL*8 RHO(1),C(1)
!
!     DCHDD DOWNDATES AN AUGMENTED CHOLESKY DECOMPOSITION OR THE
!     TRIANGULAR FACTOR OF AN AUGMENTED QR DECOMPOSITION.
!     SPECIFICALLY, GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P,  A
!     ROW VECTOR X, A COLUMN VECTOR Z, AND A SCALAR Y, DCHDD
!     DETERMINEDS A ORTHOGONAL MATRIX U AND A SCALAR ZETA SUCH THAT
!
!                        (R   Z )     (RR  ZZ)
!                    U * (      )  =  (      ) ,
!                        (0 ZETA)     ( X   Y)
!
!     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN OBTAINED
!     FROM THE FACTORIZATION OF A LEAST SQUARES PROBLEM, THEN
!     RR AND ZZ ARE THE FACTORS CORRESPONDING TO THE PROBLEM
!     WITH THE OBSERVATION (X,Y) REMOVED.  IN THIS CASE, IF RHO
!     IS THE NORM OF THE RESIDUAL VECTOR, THEN THE NORM OF
!     THE RESIDUAL VECTOR OF THE DOWNDATED PROBLEM IS
!     DSQRT(RHO**2 - ZETA**2). DCHDD WILL SIMULTANEOUSLY DOWNDATE
!     SEVERAL TRIPLETS (Z,Y,RHO) ALONG WITH R.
!     FOR A LESS TERSE DESCRIPTION OF WHAT DCHDD DOES AND HOW
!     IT MAY BE APPLIED, SEE THE LINPACK GUIDE.
!
!     THE MATRIX U IS DETERMINED AS THE PRODUCT U(1)*...*U(P)
!     WHERE U(I) IS A ROTATION IN THE (P+1,I)-PLANE OF THE
!     FORM
!
!                       ( C(I)     -S(I)     )
!                       (                    ) .
!                       ( S(I)       C(I)    )
!
!     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS REAL*8.
!
!     THE USER IS WARNED THAT A GIVEN DOWNDATING PROBLEM MAY
!     BE IMPOSSIBLE TO ACCOMPLISH OR MAY PRODUCE
!     INACCURATE RESULTS.  FOR EXAMPLE, THIS CAN HAPPEN
!     IF X IS NEAR A VECTOR WHOSE REMOVAL WILL REDUCE THE
!     RANK OF R.  BEWARE.
!
!     ON ENTRY
!
!         R      REAL*8(LDR,P), WHERE LDR .GE. P.
!                R CONTAINS THE UPPER TRIANGULAR MATRIX
!                THAT IS TO BE DOWNDATED.  THE PART OF  R
!                BELOW THE DIAGONAL IS NOT REFERENCED.
!
!         LDR    INTEGER.
!                LDR IS THE LEADING DIMENSION FO THE ARRAY R.
!
!         P      INTEGER.
!                P IS THE ORDER OF THE MATRIX R.
!
!         X      REAL*8(P).
!                X CONTAINS THE ROW VECTOR THAT IS TO
!                BE REMOVED FROM R.  X IS NOT ALTERED BY DCHDD.
!
!         Z      REAL*8(LDZ,NZ), WHERE LDZ .GE. P.
!                Z IS AN ARRAY OF NZ P-VECTORS WHICH
!                ARE TO BE DOWNDATED ALONG WITH R.
!
!         LDZ    INTEGER.
!                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z.
!
!         NZ     INTEGER.
!                NZ IS THE NUMBER OF VECTORS TO BE DOWNDATED
!                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO
!                ARE NOT REFERENCED.
!
!         Y      REAL*8(NZ).
!                Y CONTAINS THE SCALARS FOR THE DOWNDATING
!                OF THE VECTORS Z.  Y IS NOT ALTERED BY DCHDD.
!
!         RHO    REAL*8(NZ).
!                RHO CONTAINS THE NORMS OF THE RESIDUAL
!                VECTORS THAT ARE TO BE DOWNDATED.
!
!     ON RETURN
!
!         R
!         Z      CONTAIN THE DOWNDATED QUANTITIES.
!         RHO
!
!         C      REAL*8(P).
!                C CONTAINS THE COSINES OF THE TRANSFORMING
!                ROTATIONS.
!
!         S      REAL*8(P).
!                S CONTAINS THE SINES OF THE TRANSFORMING
!                ROTATIONS.
!
!         INFO   INTEGER.
!                INFO IS SET AS FOLLOWS.
!
!                   INFO = 0  IF THE ENTIRE DOWNDATING
!                             WAS SUCCESSFUL.
!
!                   INFO =-1  IF R COULD NOT BE DOWNDATED.
!                             IN THIS CASE, ALL QUANTITIES
!                             ARE LEFT UNALTERED.
!
!                   INFO = 1  IF SOME RHO COULD NOT BE
!                             DOWNDATED.  THE OFFENDING RHOS ARE
!                             SET TO -1.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     DCHDD USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
!
!     FORTRAN DABS
!     BLAS DDOT, DNRM2
!
      INTEGER I,II,J
      REAL*8 A,ALPHA,AZETA,NORM
      REAL*8 T,ZETA,B,XX
!
!     SOLVE THE SYSTEM TRANS(R)*A = X, PLACING THE RESULT
!     IN THE ARRAY S.
!
      INFO = 0
      S(1) = X(1)/R(1,1)
      IF (P .LT. 2) GO TO 20
      DO 10 J = 2, P
         S(J) = X(J) - DDOT(J-1,R(1,J),1,S,1)
         S(J) = S(J)/R(J,J)
   10 CONTINUE
   20 CONTINUE
      NORM = DNRM2(P,S,1)
      IF (NORM .LT. 1.0D0) GO TO 30
         INFO = -1
      GO TO 120
   30 CONTINUE
         ALPHA = DSQRT(1.0D0-NORM**2)
!
!        DETERMINE THE TRANSFORMATIONS.
!
         DO 40 II = 1, P
            I = P - II + 1
            SCALE = ALPHA + DABS(S(I))
            A = ALPHA/SCALE
            B = S(I)/SCALE
            NORM = DSQRT(A**2+B**2+0.0D0**2)
            C(I) = A/NORM
            S(I) = B/NORM
            ALPHA = SCALE*NORM
   40    CONTINUE
!
!        APPLY THE TRANSFORMATIONS TO R.
!
         DO 60 J = 1, P
            XX = 0.0D0
            DO 50 II = 1, J
               I = J - II + 1
               T = C(I)*XX + S(I)*R(I,J)
               R(I,J) = C(I)*R(I,J) - S(I)*XX
               XX = T
   50       CONTINUE
   60    CONTINUE
!
!        IF REQUIRED, DOWNDATE Z AND RHO.
!
         IF (NZ .LT. 1) GO TO 110
         DO 100 J = 1, NZ
            ZETA = Y(J)
            DO 70 I = 1, P
               Z(I,J) = (Z(I,J) - S(I)*ZETA)/C(I)
               ZETA = C(I)*ZETA - S(I)*Z(I,J)
   70       CONTINUE
            AZETA = DABS(ZETA)
            IF (AZETA .LE. RHO(J)) GO TO 80
               INFO = 1
               RHO(J) = -1.0D0
            GO TO 90
   80       CONTINUE
               RHO(J) = RHO(J)*DSQRT(1.0D0-(AZETA/RHO(J))**2)
   90       CONTINUE
  100    CONTINUE
  110    CONTINUE
  120 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DCHEX(R,LDR,P,K,L,Z,LDZ,NZ,C,S,JOB)
      INTEGER LDR,P,K,L,LDZ,NZ,JOB
      REAL*8 R(LDR,1),Z(LDZ,1),S(1)
      REAL*8 C(1)
!
!     DCHEX UPDATES THE CHOLESKY FACTORIZATION
!
!                   A = TRANS(R)*R
!
!     OF A POSITIVE DEFINITE MATRIX A OF ORDER P UNDER DIAGONAL
!     PERMUTATIONS OF THE FORM
!
!                   TRANS(E)*A*E
!
!     WHERE E IS A PERMUTATION MATRIX.  SPECIFICALLY, GIVEN
!     AN UPPER TRIANGULAR MATRIX R AND A PERMUTATION MATRIX
!     E (WHICH IS SPECIFIED BY K, L, AND JOB), DCHEX DETERMINES
!     A ORTHOGONAL MATRIX U SUCH THAT
!
!                           U*R*E = RR,
!
!     WHERE RR IS UPPER TRIANGULAR.  AT THE USERS OPTION, THE
!     TRANSFORMATION U WILL BE MULTIPLIED INTO THE ARRAY Z.
!     IF A = TRANS(X)*X, SO THAT R IS THE TRIANGULAR PART OF THE
!     QR FACTORIZATION OF X, THEN RR IS THE TRIANGULAR PART OF THE
!     QR FACTORIZATION OF X*E, I.E. X WITH ITS COLUMNS PERMUTED.
!     FOR A LESS TERSE DESCRIPTION OF WHAT DCHEX DOES AND HOW
!     IT MAY BE APPLIED, SEE THE LINPACK GUIDE.
!
!     THE MATRIX Q IS DETERMINED AS THE PRODUCT U(L-K)*...*U(1)
!     OF PLANE ROTATIONS OF THE FORM
!
!                           (    C(I)       S(I) )
!                           (                    ) ,
!                           (    -S(I)      C(I) )
!
!     WHERE C(I) IS REAL*8, THE ROWS THESE ROTATIONS OPERATE
!     ON ARE DESCRIBED BELOW.
!
!     THERE ARE TWO TYPES OF PERMUTATIONS, WHICH ARE DETERMINED
!     BY THE VALUE OF JOB.
!
!     1. RIGHT CIRCULAR SHIFT (JOB = 1).
!
!         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER.
!
!                1,...,K-1,L,K,K+1,...,L-1,L+1,...,P.
!
!         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I)
!         ACTS IN THE (L-I,L-I+1)-PLANE.
!
!     2. LEFT CIRCULAR SHIFT (JOB = 2).
!         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER
!
!                1,...,K-1,K+1,K+2,...,L,K,L+1,...,P.
!
!         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I)
!         ACTS IN THE (K+I-1,K+I)-PLANE.
!
!     ON ENTRY
!
!         R      REAL*8(LDR,P), WHERE LDR.GE.P.
!                R CONTAINS THE UPPER TRIANGULAR FACTOR
!                THAT IS TO BE UPDATED.  ELEMENTS OF R
!                BELOW THE DIAGONAL ARE NOT REFERENCED.
!
!         LDR    INTEGER.
!                LDR IS THE LEADING DIMENSION OF THE ARRAY R.
!
!         P      INTEGER.
!                P IS THE ORDER OF THE MATRIX R.
!
!         K      INTEGER.
!                K IS THE FIRST COLUMN TO BE PERMUTED.
!
!         L      INTEGER.
!                L IS THE LAST COLUMN TO BE PERMUTED.
!                L MUST BE STRICTLY GREATER THAN K.
!
!         Z      REAL*8(LDZ,NZ), WHERE LDZ.GE.P.
!                Z IS AN ARRAY OF NZ P-VECTORS INTO WHICH THE
!                TRANSFORMATION U IS MULTIPLIED.  Z IS
!                NOT REFERENCED IF NZ = 0.
!
!         LDZ    INTEGER.
!                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z.
!
!         NZ     INTEGER.
!                NZ IS THE NUMBER OF COLUMNS OF THE MATRIX Z.
!
!         JOB    INTEGER.
!                JOB DETERMINES THE TYPE OF PERMUTATION.
!                       JOB = 1  RIGHT CIRCULAR SHIFT.
!                       JOB = 2  LEFT CIRCULAR SHIFT.
!
!     ON RETURN
!
!         R      CONTAINS THE UPDATED FACTOR.
!
!         Z      CONTAINS THE UPDATED MATRIX Z.
!
!         C      REAL*8(P).
!                C CONTAINS THE COSINES OF THE TRANSFORMING ROTATIONS.
!
!         S      REAL*8(P).
!                S CONTAINS THE SINES OF THE TRANSFORMING ROTATIONS.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     DCHEX USES THE FOLLOWING FUNCTIONS AND SUBROUTINES.
!
!     BLAS DROTG
!     FORTRAN MIN0
!
      INTEGER I,II,IL,IU,J,JJ,KM1,KP1,LMK,LM1
      REAL*8 RJP1J,T
!
!     INITIALIZE
!
      KM1 = K - 1
      KP1 = K + 1
      LMK = L - K
      LM1 = L - 1
!
!     PERFORM THE APPROPRIATE TASK.
!
      GO TO (10,130), JOB
!
!     RIGHT CIRCULAR SHIFT.
!
   10 CONTINUE
!
!        REORDER THE COLUMNS.
!
         DO 20 I = 1, L
            II = L - I + 1
            S(I) = R(II,L)
   20    CONTINUE
         DO 40 JJ = K, LM1
            J = LM1 - JJ + K
            DO 30 I = 1, J
               R(I,J+1) = R(I,J)
   30       CONTINUE
            R(J+1,J+1) = 0.0D0
   40    CONTINUE
         IF (K .EQ. 1) GO TO 60
            DO 50 I = 1, KM1
               II = L - I + 1
               R(I,K) = S(II)
   50       CONTINUE
   60    CONTINUE
!
!        CALCULATE THE ROTATIONS.
!
         T = S(1)
         DO 70 I = 1, LMK
            CALL DROTG(S(I+1),T,C(I),S(I))
            T = S(I+1)
   70    CONTINUE
         R(K,K) = T
         DO 90 J = KP1, P
            IL = MAX0(1,L-J+1)
            DO 80 II = IL, LMK
               I = L - II
               T = C(II)*R(I,J) + S(II)*R(I+1,J)
               R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)
               R(I,J) = T
   80       CONTINUE
   90    CONTINUE
!
!        IF REQUIRED, APPLY THE TRANSFORMATIONS TO Z.
!
         IF (NZ .LT. 1) GO TO 120
         DO 110 J = 1, NZ
            DO 100 II = 1, LMK
               I = L - II
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
               Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)
               Z(I,J) = T
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      GO TO 260
!
!     LEFT CIRCULAR SHIFT
!
  130 CONTINUE
!
!        REORDER THE COLUMNS
!
         DO 140 I = 1, K
            II = LMK + I
            S(II) = R(I,K)
  140    CONTINUE
         DO 160 J = K, LM1
            DO 150 I = 1, J
               R(I,J) = R(I,J+1)
  150       CONTINUE
            JJ = J - KM1
            S(JJ) = R(J+1,J+1)
  160    CONTINUE
         DO 170 I = 1, K
            II = LMK + I
            R(I,L) = S(II)
  170    CONTINUE
         DO 180 I = KP1, L
            R(I,L) = 0.0D0
  180    CONTINUE
!
!        REDUCTION LOOP.
!
         DO 220 J = K, P
            IF (J .EQ. K) GO TO 200
!
!              APPLY THE ROTATIONS.
!
               IU = MIN0(J-1,L-1)
               DO 190 I = K, IU
                  II = I - K + 1
                  T = C(II)*R(I,J) + S(II)*R(I+1,J)
                  R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)
                  R(I,J) = T
  190          CONTINUE
  200       CONTINUE
            IF (J .GE. L) GO TO 210
               JJ = J - K + 1
               T = S(JJ)
               CALL DROTG(R(J,J),T,C(JJ),S(JJ))
  210       CONTINUE
  220    CONTINUE
!
!        APPLY THE ROTATIONS TO Z.
!
         IF (NZ .LT. 1) GO TO 250
         DO 240 J = 1, NZ
            DO 230 I = K, LM1
               II = I - KM1
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
               Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)
               Z(I,J) = T
  230       CONTINUE
  240    CONTINUE
  250    CONTINUE
  260 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)
      INTEGER LDX,N,P,JOB
      INTEGER JPVT(1)
      REAL*8 X(LDX,1),QRAUX(1),WORK(1)
!
!     DQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR
!     FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING
!     BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE
!     PERFORMED AT THE USERS OPTION.
!
!     ON ENTRY
!
!        X       REAL*8(LDX,P), WHERE LDX .GE. N.
!                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE
!                COMPUTED.
!
!        LDX     INTEGER.
!                LDX IS THE LEADING DIMENSION OF THE ARRAY X.
!
!        N       INTEGER.
!                N IS THE NUMBER OF ROWS OF THE MATRIX X.
!
!        P       INTEGER.
!                P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
!
!        JPVT    INTEGER(P).
!                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
!                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X
!                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
!                VALUE OF JPVT(K).
!
!                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
!                                      COLUMN.
!
!                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN.
!
!                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN.
!
!                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS
!                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL
!                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS
!                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
!                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE
!                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN
!                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST
!                REDUCED NORM.  JPVT IS NOT REFERENCED IF
!                JOB .EQ. 0.
!
!        WORK    REAL*8(P).
!                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF
!                JOB .EQ. 0.
!
!        JOB     INTEGER.
!                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
!                IF JOB .EQ. 0, NO PIVOTING IS DONE.
!                IF JOB .NE. 0, PIVOTING IS DONE.
!
!     ON RETURN
!
!        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER
!                TRIANGULAR MATRIX R OF THE QR FACTORIZATION.
!                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM
!                WHICH THE ORTHOGONAL PART OF THE DECOMPOSITION
!                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS
!                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT
!                OF THE ORIGINAL MATRIX X BUT THAT OF X
!                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT.
!
!        QRAUX   REAL*8(P).
!                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER
!                THE ORTHOGONAL PART OF THE DECOMPOSITION.
!
!        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE
!                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO
!                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     DQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
!
!     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2
!     FORTRAN DABS,DMAX1,MIN0,DSQRT
!
!     INTERNAL VARIABLES
!
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU
      REAL*8 MAXNRM,TT
      REAL*8 NRMXL,T
      LOGICAL NEGJ,SWAPJ
!
!
      PL = 1
      PU = 0
      IF (JOB .EQ. 0) GO TO 60
!
!        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
!        ACCORDING TO JPVT.
!
         DO 20 J = 1, P
            SWAPJ = JPVT(J) .GT. 0
            NEGJ = JPVT(J) .LT. 0
            JPVT(J) = J
            IF (NEGJ) JPVT(J) = -J
            IF (.NOT.SWAPJ) GO TO 10
               IF (J .NE. PL) CALL DSWAP(N,X(1,PL),1,X(1,J),1)
               JPVT(J) = JPVT(PL)
               JPVT(PL) = J
               PL = PL + 1
   10       CONTINUE
   20    CONTINUE
         PU = P
         DO 50 JJ = 1, P
            J = P - JJ + 1
            IF (JPVT(J) .GE. 0) GO TO 40
               JPVT(J) = -JPVT(J)
               IF (J .EQ. PU) GO TO 30
                  CALL DSWAP(N,X(1,PU),1,X(1,J),1)
                  JP = JPVT(PU)
                  JPVT(PU) = JPVT(J)
                  JPVT(J) = JP
   30          CONTINUE
               PU = PU - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
!
!     COMPUTE THE NORMS OF THE FREE COLUMNS.
!
      IF (PU .LT. PL) GO TO 80
      DO 70 J = PL, PU
         QRAUX(J) = DNRM2(N,X(1,J),1)
         WORK(J) = QRAUX(J)
   70 CONTINUE
   80 CONTINUE
!
!     PERFORM THE HOUSEHOLDER REDUCTION OF X.
!
      LUP = MIN0(N,P)
      DO 200 L = 1, LUP
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120
!
!           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
!           INTO THE PIVOT POSITION.
!
            MAXNRM = 0.0D0
            MAXJ = L
            DO 100 J = L, PU
               IF (QRAUX(J) .LE. MAXNRM) GO TO 90
                  MAXNRM = QRAUX(J)
                  MAXJ = J
   90          CONTINUE
  100       CONTINUE
            IF (MAXJ .EQ. L) GO TO 110
               CALL DSWAP(N,X(1,L),1,X(1,MAXJ),1)
               QRAUX(MAXJ) = QRAUX(L)
               WORK(MAXJ) = WORK(L)
               JP = JPVT(MAXJ)
               JPVT(MAXJ) = JPVT(L)
               JPVT(L) = JP
  110       CONTINUE
  120    CONTINUE
         QRAUX(L) = 0.0D0
         IF (L .EQ. N) GO TO 190
!
!           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.
!
            NRMXL = DNRM2(N-L+1,X(L,L),1)
            IF (NRMXL .EQ. 0.0D0) GO TO 180
               IF (X(L,L) .NE. 0.0D0) NRMXL = DSIGN(NRMXL,X(L,L))
               CALL DSCAL(N-L+1,1.0D0/NRMXL,X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
!
!              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
!              UPDATING THE NORMS.
!
               LP1 = L + 1
               IF (P .LT. LP1) GO TO 170
               DO 160 J = LP1, P
                  T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
                  CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150
                  IF (QRAUX(J) .EQ. 0.0D0) GO TO 150
                     TT = 1.0D0 - (DABS(X(L,J))/QRAUX(J))**2
                     TT = DMAX1(TT,0.0D0)
                     T = TT
                     TT = 1.0D0 + 0.05D0*TT*(QRAUX(J)/WORK(J))**2
                     IF (TT .EQ. 1.0D0) GO TO 130
                        QRAUX(J) = QRAUX(J)*DSQRT(T)
                     GO TO 140
  130                CONTINUE
                        QRAUX(J) = DNRM2(N-L,X(L+1,J),1)
                        WORK(J) = QRAUX(J)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE
!
!              SAVE THE TRANSFORMATION.
!
               QRAUX(L) = X(L,L)
               X(L,L) = -NRMXL
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)
      INTEGER LDX,N,K,JOB,INFO
      REAL*8 X(LDX,1),QRAUX(1),Y(1),QY(1),QTY(1),B(1),RSD(1),XB(1)
!
!     DQRSL APPLIES THE OUTPUT OF DQRDC TO COMPUTE COORDINATE
!     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.
!     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX
!
!            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
!
!     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL
!     N X P MATRIX X THAT WAS INPUT TO DQRDC (IF NO PIVOTING WAS
!     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR
!     ORIGINAL ORDER).  DQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q
!     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT
!
!              XK = Q * (R)
!                       (0)
!
!     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS
!     X AND QRAUX.
!
!     ON ENTRY
!
!        X      REAL*8(LDX,P).
!               X CONTAINS THE OUTPUT OF DQRDC.
!
!        LDX    INTEGER.
!               LDX IS THE LEADING DIMENSION OF THE ARRAY X.
!
!        N      INTEGER.
!               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST
!               HAVE THE SAME VALUE AS N IN DQRDC.
!
!        K      INTEGER.
!               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K
!               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE
!               SAME AS IN THE CALLING SEQUENCE TO DQRDC.
!
!        QRAUX  REAL*8(P).
!               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM DQRDC.
!
!        Y      REAL*8(N)
!               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED
!               BY DQRSL.
!
!        JOB    INTEGER.
!               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS
!               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING
!               MEANING.
!
!                    IF A.NE.0, COMPUTE QY.
!                    IF B,C,D, OR E .NE. 0, COMPUTE QTY.
!                    IF C.NE.0, COMPUTE B.
!                    IF D.NE.0, COMPUTE RSD.
!                    IF E.NE.0, COMPUTE XB.
!
!               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB
!               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR
!               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING
!               SEQUENCE.
!
!     ON RETURN
!
!        QY     REAL*8(N).
!               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN
!               REQUESTED.
!
!        QTY    REAL*8(N).
!               QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS
!               BEEN REQUESTED.  HERE TRANS(Q) IS THE
!               TRANSPOSE OF THE MATRIX Q.
!
!        B      REAL*8(K)
!               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM
!
!                    MINIMIZE NORM2(Y - XK*B),
!
!               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT
!               IF PIVOTING WAS REQUESTED IN DQRDC, THE J-TH
!               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)
!               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO DQRDC.)
!
!        RSD    REAL*8(N).
!               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,
!               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS
!               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE
!               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.
!
!        XB     REAL*8(N).
!               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,
!               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO
!               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE
!               OF X.
!
!        INFO   INTEGER.
!               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS
!               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN
!               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO
!               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.
!
!     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED
!     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE
!     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM.
!     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME
!     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A
!     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE
!     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS
!     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE
!     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE
!     COMPUTED.  THUS THE CALLING SEQUENCE
!
!          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
!
!     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD
!     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING
!     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR
!     A SINGLE CALLINNG SEQUENCE.
!
!          1. (Y,QTY,B) (RSD) (XB) (QY)
!
!          2. (Y,QTY,RSD) (B) (XB) (QY)
!
!          3. (Y,QTY,XB) (B) (RSD) (QY)
!
!          4. (Y,QY) (QTY,B) (RSD) (XB)
!
!          5. (Y,QY) (QTY,RSD) (B) (XB)
!
!          6. (Y,QY) (QTY,XB) (B) (RSD)
!
!     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO
!     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     DQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
!
!     BLAS DAXPY,DCOPY,DDOT
!     FORTRAN DABS,MIN0,MOD
!
!     INTERNAL VARIABLES
!
      INTEGER I,J,JJ,JU,KP1
      REAL*8 T,TEMP
      LOGICAL CB,CQY,CQTY,CR,CXB
!
!
!     SET INFO FLAG.
!
      INFO = 0
!
!     DETERMINE WHAT IS TO BE COMPUTED.
!
      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN0(K,N-1)
!
!     SPECIAL ACTION WHEN N=1.
!
      IF (JU .NE. 0) GO TO 40
         IF (CQY) QY(1) = Y(1)
         IF (CQTY) QTY(1) = Y(1)
         IF (CXB) XB(1) = Y(1)
         IF (.NOT.CB) GO TO 30
            IF (X(1,1) .NE. 0.0D0) GO TO 10
               INFO = 1
            GO TO 20
   10       CONTINUE
               B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
         IF (CR) RSD(1) = 0.0D0
      GO TO 250
   40 CONTINUE
!
!        SET UP TO COMPUTE QY OR QTY.
!
         IF (CQY) CALL DCOPY(N,Y,1,QY,1)
         IF (CQTY) CALL DCOPY(N,Y,1,QTY,1)
         IF (.NOT.CQY) GO TO 70
!
!           COMPUTE QY.
!
            DO 60 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 50
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QY(J),1)
                  X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
         IF (.NOT.CQTY) GO TO 100
!
!           COMPUTE TRANS(Q)*Y.
!
            DO 90 J = 1, JU
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 80
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
                  X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        SET UP TO COMPUTE B, RSD, OR XB.
!
         IF (CB) CALL DCOPY(K,QTY,1,B,1)
         KP1 = K + 1
         IF (CXB) CALL DCOPY(K,QTY,1,XB,1)
         IF (CR .AND. K .LT. N) CALL DCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120
            DO 110 I = KP1, N
               XB(I) = 0.0D0
  110       CONTINUE
  120    CONTINUE
         IF (.NOT.CR) GO TO 140
            DO 130 I = 1, K
               RSD(I) = 0.0D0
  130       CONTINUE
  140    CONTINUE
         IF (.NOT.CB) GO TO 190
!
!           COMPUTE B.
!
            DO 170 JJ = 1, K
               J = K - JJ + 1
               IF (X(J,J) .NE. 0.0D0) GO TO 150
                  INFO = J
!           ......EXIT
                  GO TO 180
  150          CONTINUE
               B(J) = B(J)/X(J,J)
               IF (J .EQ. 1) GO TO 160
                  T = -B(J)
                  CALL DAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240
!
!           COMPUTE RSD OR XB AS REQUIRED.
!
            DO 230 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 220
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  IF (.NOT.CR) GO TO 200
                     T = -DDOT(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
                  IF (.NOT.CXB) GO TO 210
                     T = -DDOT(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
                  X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      REAL*8 X(LDX,1),S(1),E(1),U(LDU,1),V(LDV,1),WORK(1)
!
!
!     DSVDC IS A SUBROUTINE TO REDUCE A REAL*8 NXP MATRIX X
!     BY ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
!     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
!     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
!     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
!
!     ON ENTRY
!
!         X         REAL*8(LDX,P), WHERE LDX.GE.N.
!                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
!                   DECOMPOSITION IS TO BE COMPUTED.  X IS
!                   DESTROYED BY DSVDC.
!
!         LDX       INTEGER.
!                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.
!
!         N         INTEGER.
!                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X.
!
!         P         INTEGER.
!                   P IS THE NUMBER OF ROWS OF THE MATRIX X.
!
!         LDU       INTEGER.
!                   LDU IS THE LEADING DIMENSION OF THE ARRAY U.
!                   (SEE BELOW).
!
!         LDV       INTEGER.
!                   LDV IS THE LEADING DIMENSION OF THE ARRAY V.
!                   (SEE BELOW).
!
!         WORK      REAL*8(N).
!                   WORK IS A SCRATCH ARRAY.
!
!         JOB       INTEGER.
!                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR
!                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB
!                   WITH THE FOLLOWING MEANING
!
!                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR
!                                  VECTORS.
!                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS
!                                  IN U.
!                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR
!                                  VECTORS IN U.
!                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
!                                  VECTORS.
!                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
!                                  IN V.
!
!     ON RETURN
!
!         S         REAL*8(MM), WHERE MM=MIN(N+1,P).
!                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
!                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
!                   ORDER OF MAGNITUDE.
!
!         E         REAL*8(P).
!                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
!                   DISCUSSION OF INFO FOR EXCEPTIONS.
!
!         U         REAL*8(LDU,K), WHERE LDU.GE.N.  IF
!                                   JOBA.EQ.1 THEN K.EQ.N, IF JOBA.GE.2
!                                   THEN K.EQ.MIN(N,P).
!                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
!                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
!                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X
!                   IN THE SUBROUTINE CALL.
!
!         V         REAL*8(LDV,P), WHERE LDV.GE.P.
!                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
!                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N,
!                   THEN V MAY BE IDENTIFIED WITH X IN THE
!                   SUBROUTINE CALL.
!
!         INFO      INTEGER.
!                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
!                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
!                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
!                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
!                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
!                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX
!                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
!                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U)
!                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR
!                   VALUES OF X AND B ARE THE SAME.
!
!     LINPACK. THIS VERSION DATED 03/19/79 .
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     DSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
!
!     EXTERNAL DROT
!     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2,DROTG
!     FORTRAN DABS,DMAX1,MAX0,MIN0,MOD,DSQRT
!
!     INTERNAL VARIABLES
!
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,&
         MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      REAL*8 T,R
      REAL*8 B,C,CS,EL,EMM1,F,G,SCALE,SHIFT,SL,SM,SN,&
         SMM1,T1,TEST,ZTEST
      LOGICAL WANTU,WANTV
!
!
!     SET THE MAXIMUM NUMBER OF ITERATIONS.
!
      MAXIT = 30
!
!     DETERMINE WHAT IS TO BE COMPUTED.
!
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN0(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
!
!     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
!     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
!
      INFO = 0
      NCT = MIN0(N-1,P)
      NRT = MAX0(0,MIN0(P-2,N))
      LU = MAX0(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
!
!           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
!           PLACE THE L-TH DIAGONAL IN S(L).
!
            S(L) = DNRM2(N-L+1,X(L,L),1)
            IF (S(L) .EQ. 0.0D0) GO TO 10
               IF (X(L,L) .NE. 0.0D0) S(L) = DSIGN(S(L),X(L,L))
               CALL DSCAL(N-L+1,1.0D0/S(L),X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .EQ. 0.0D0) GO TO 30
!
!              APPLY THE TRANSFORMATION.
!
               T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
!
!           PLACE THE L-TH ROW OF X INTO  E FOR THE
!           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
!
            E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
!
!           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
!           MULTIPLICATION.
!
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
!
!           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
!           L-TH SUPER-DIAGONAL IN E(L).
!
            E(L) = DNRM2(P-L,E(LP1),1)
            IF (E(L) .EQ. 0.0D0) GO TO 80
               IF (E(LP1) .NE. 0.0D0) E(L) = DSIGN(E(L),E(LP1))
               CALL DSCAL(P-L,1.0D0/E(L),E(LP1),1)
               E(LP1) = 1.0D0 + E(LP1)
   80       CONTINUE
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0D0) GO TO 120
!
!              APPLY THE TRANSFORMATION.
!
               DO 90 I = LP1, N
                  WORK(I) = 0.0D0
   90          CONTINUE
               DO 100 J = LP1, P
                  CALL DAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL DAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
!
!              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
!              BACK MULTIPLICATION.
!
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
!
!     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
!
      M = MIN0(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0D0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0D0
!
!     IF REQUIRED, GENERATE U.
!
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = 0.0D0
  180       CONTINUE
            U(J,J) = 1.0D0
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0D0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
               DO 210 J = LP1, NCU
                  T = -DDOT(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL DAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL DSCAL(N-L+1,-1.0D0,U(L,L),1)
               U(L,L) = 1.0D0 + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = 0.0D0
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = 0.0D0
  260          CONTINUE
               U(L,L) = 1.0D0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
!
!     IF IT IS REQUIRED, GENERATE V.
!
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0D0) GO TO 320
               DO 310 J = LP1, P
                  T = -DDOT(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL DAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = 0.0D0
  330       CONTINUE
            V(L,L) = 1.0D0
  340    CONTINUE
  350 CONTINUE
!
!     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
!
      MM = M
      ITER = 0
  360 CONTINUE
!
!        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
!
!     ...EXIT
         IF (M .EQ. 0) GO TO 620
!
!        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
!        FLAG AND RETURN.
!
         IF (ITER .LT. MAXIT) GO TO 370
            INFO = M
!     ......EXIT
            GO TO 620
  370    CONTINUE
!
!        THIS SECTION OF THE PROGRAM INSPECTS FOR
!        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
!        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
!
!           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
!           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
!           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
!                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
!           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
!
         DO 390 LL = 1, M
            L = M - LL
!        ...EXIT
            IF (L .EQ. 0) GO TO 400
            TEST = DABS(S(L)) + DABS(S(L+1))
            ZTEST = TEST + DABS(E(L))
            IF (ZTEST .NE. TEST) GO TO 380
               E(L) = 0.0D0
!        ......EXIT
               GO TO 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
         IF (L .NE. M - 1) GO TO 410
            KASE = 4
         GO TO 480
  410    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 430 LLS = LP1, MP1
               LS = M - LLS + LP1
!           ...EXIT
               IF (LS .EQ. L) GO TO 440
               TEST = 0.0D0
               IF (LS .NE. M) TEST = TEST + DABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + DABS(E(LS-1))
               ZTEST = TEST + DABS(S(LS))
               IF (ZTEST .NE. TEST) GO TO 420
                  S(LS) = 0.0D0
!           ......EXIT
                  GO TO 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
            IF (LS .NE. L) GO TO 450
               KASE = 3
            GO TO 470
  450       CONTINUE
            IF (LS .NE. M) GO TO 460
               KASE = 1
            GO TO 470
  460       CONTINUE
               KASE = 2
               L = LS
  470       CONTINUE
  480    CONTINUE
         L = L + 1
!
!        PERFORM THE TASK INDICATED BY KASE.
!
         GO TO (490,520,540,570), KASE
!
!        DEFLATE NEGLIGIBLE S(M).
!
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0D0
            DO 510 KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               IF (K .EQ. L) GO TO 500
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
  500          CONTINUE
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
         GO TO 610
!
!        SPLIT AT NEGLIGIBLE S(L).
!
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0D0
            DO 530 K = L, M
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL DROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
         GO TO 610
!
!        PERFORM ONE QR STEP.
!
  540    CONTINUE
!
!           CALCULATE THE SHIFT.
!
            SCALE = DMAX1(DABS(S(M)),DABS(S(M-1)),DABS(E(M-1)),DABS(S(L)),DABS(E(L)))
            SM = S(M)/SCALE
            SMM1 = S(M-1)/SCALE
            EMM1 = E(M-1)/SCALE
            SL = S(L)/SCALE
            EL = E(L)/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0
            C = (SM*EMM1)**2
            SHIFT = 0.0D0
            IF (B .EQ. 0.0D0 .AND. C .EQ. 0.0D0) GO TO 550
               SHIFT = DSQRT(B**2+C)
               IF (B .LT. 0.0D0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  550       CONTINUE
            F = (SL + SM)*(SL - SM) + SHIFT
            G = SL*EL
!
!           CHASE ZEROS.
!
            MM1 = M - 1
            DO 560 K = L, MM1
               CALL DROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL DROTG(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N) CALL DROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
!
!        CONVERGENCE.
!
  570    CONTINUE
!
!           MAKE THE SINGULAR VALUE  POSITIVE.
!
            IF (S(L) .GE. 0.0D0) GO TO 580
               S(L) = -S(L)
               IF (WANTV) CALL DSCAL(P,-1.0D0,V(1,L),1)
  580       CONTINUE
!
!           ORDER THE SINGULAR VALUE.
!
  590       IF (L .EQ. MM) GO TO 600
!           ...EXIT
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P) CALL DSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N) CALL DSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            ITER = 0
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
  END subroutine

  REAL*8 FUNCTION DASUM(N,DX,INCX)
!
!     TAKES THE SUM OF THE ABSOLUTE VALUES.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL*8 DX(1),DTEMP
      INTEGER I,INCX,M,MP1,N,NINCX
!
      DASUM = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DTEMP = DTEMP + DABS(DX(I))
   10 CONTINUE
      DASUM = DTEMP
      RETURN
!
!        CODE FOR INCREMENT EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I + 1)) + DABS(DX(I + 2))&
           + DABS(DX(I + 3)) + DABS(DX(I + 4)) + DABS(DX(I + 5))
   50 CONTINUE
   60 DASUM = DTEMP
      RETURN
  END function

  SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL*8 DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
!
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE  DCOPY(N,DX,INCX,DY,INCY)
!
!     COPIES A VECTOR, X, TO A VECTOR, Y.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL*8 DX(1),DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
  END subroutine

  REAL*8 FUNCTION DDOT(N,DX,INCX,DY,INCY)
!
!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL*8 DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
!
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +&
           DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
  END function

  REAL*8 FUNCTION DMACH(JOB)
      INTEGER JOB
!
!     SMACH COMPUTES MACHINE PARAMETERS OF FLOATING POINT
!     ARITHMETIC FOR USE IN TESTING ONLY.  NOT REQUIRED BY
!     LINPACK PROPER.
!
!     IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES,
!     THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS.
!     ASSUME THE COMPUTER HAS
!
!        B = BASE OF ARITHMETIC
!        T = NUMBER OF BASE  B  DIGITS
!        L = SMALLEST POSSIBLE EXPONENT
!        U = LARGEST POSSIBLE EXPONENT
!
!     THEN
!
!        EPS = B**(1-T)
!        TINY = 100.0*B**(-L+T)
!        HUGE = 0.01*B**(U-T)
!
!     DMACH SAME AS SMACH EXCEPT T, L, U APPLY TO
!     REAL*8.
!
!     CMACH SAME AS SMACH EXCEPT IF COMPLEX DIVISION
!     IS DONE BY
!
!        1/(X+I*Y) = (X-I*Y)/(X**2+Y**2)
!
!     THEN
!
!        TINY = SQRT(TINY)
!        HUGE = SQRT(HUGE)
!
!
!     JOB IS 1, 2 OR 3 FOR EPSILON, TINY AND HUGE, RESPECTIVELY.
!
      REAL*8 EPS,TINY,HUGE,S
!
      EPS = 1.0D0
   10 EPS = EPS/2.0D0
      S = 1.0D0 + EPS
      IF (S .GT. 1.0D0) GO TO 10
      EPS = 2.0D0*EPS
!
      S = 1.0D0
   20 TINY = S
      S = S/16.0D0
      IF (S*1.0 .NE. 0.0D0) GO TO 20
      TINY = (TINY/EPS)*100.0
      HUGE = 1.0D0/TINY
!
      IF (JOB .EQ. 1) DMACH = EPS
      IF (JOB .EQ. 2) DMACH = TINY
      IF (JOB .EQ. 3) DMACH = HUGE
      RETURN
  END function

  REAL*8 FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER          NEXT
      REAL*8   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
!
!     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
!     INCREMENT INCX .
!     IF    N .LE. 0 RETURN WITH RESULT = 0.
!     IF N .GE. 1 THEN INCX MUST BE .GE. 1
!
!           C.L.LAWSON, 1978 JAN 08
!
!     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
!     HOPEFULLY APPLICABLE TO ALL MACHINES.
!         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
!         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
!     WHERE
!         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
!         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
!         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
!
!     BRIEF OUTLINE OF ALGORITHM..
!
!     PHASE 1    SCANS ZERO COMPONENTS.
!     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
!     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
!     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
!     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
!
!     VALUES FOR CUTLO AND CUTHI..
!     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
!     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
!                   UNIVAC AND DEC AT 2**(-103)
!                   THUS CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
!                   THUS CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
!                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
!     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
!     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
!
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
!
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
!                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
!
!                        PHASE 1.  SUM IS ZERO
!
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
!
!                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
!
!                                PREPARE FOR PHASE 4.
!
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
!
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
!
!
!                  PREPARE FOR PHASE 3.
!
   75 SUM = (SUM * XMAX) * XMAX
!
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
   85 HITEST = CUTHI/FLOAT( N )
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
!
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
  END function

  SUBROUTINE  DROT (N,DX,INCX,DY,INCY,C,S)
!
!     APPLIES A PLANE ROTATION.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL*8 DX(1),DY(1),DTEMP,C,S
      INTEGER I,INCX,INCY,IX,IY,N
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = C*DX(IX) + S*DY(IY)
        DY(IY) = C*DY(IY) - S*DX(IX)
        DX(IX) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
   20 DO 30 I = 1,N
        DTEMP = C*DX(I) + S*DY(I)
        DY(I) = C*DY(I) - S*DX(I)
        DX(I) = DTEMP
   30 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE DROTG(DA,DB,C,S)
!
!     CONSTRUCT GIVENS PLANE ROTATION.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL*8 DA,DB,C,S,ROE,SCALE,R,Z
!
      ROE = DB
      IF( DABS(DA) .GT. DABS(DB) ) ROE = DA
      SCALE = DABS(DA) + DABS(DB)
      IF( SCALE .NE. 0.0D0 ) GO TO 10
         C = 1.0D0
         S = 0.0D0
         R = 0.0D0
         GO TO 20
   10 R = SCALE*DSQRT((DA/SCALE)**2 + (DB/SCALE)**2)
      R = DSIGN(1.0D0,ROE)*R
      C = DA/R
      S = DB/R
   20 Z = 1.0D0
      IF( DABS(DA) .GT. DABS(DB) ) Z = S
      IF( DABS(DB) .GE. DABS(DA) .AND. C .NE. 0.0D0 ) Z = 1.0D0/C
      DA = R
      DB = Z
      RETURN
  END subroutine

  SUBROUTINE  DSCAL(N,DA,DX,INCX)
!
!     SCALES A VECTOR BY A CONSTANT.
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL*8 DA,DX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
!
!        CODE FOR INCREMENT EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
  END subroutine

  SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)
!
!     INTERCHANGES TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL*8 DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!       CLEAN-UP LOOP
!
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
  END subroutine

  INTEGER FUNCTION IDAMAX(N,DX,INCX)
!
!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL*8 DX(1),DMAX
      INTEGER I,INCX,IX,N
!
      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
!
!        CODE FOR INCREMENT EQUAL TO 1
!
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
  END function

  subroutine gen_oh(code, num, x, y, z, w, a, b, v)
       implicit logical(a-z)
       REAL*8 x(*),y(*),z(*),w(*)
       REAL*8 a,b,v
       integer code
       integer num
       REAL*8 c
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated from C to fortran77 by hand.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
!vw
!vw    Given a point on a sphere (specified by a and b), generate all
!vw    the equivalent points under Oh symmetry, making grid points with
!vw    weight v.
!vw    The variable num is increased by the number of different points
!vw    generated.
!vw
!vw    Depending on code, there are 6...48 different but equivalent
!vw    points.
!vw
!vw    code=1:   (0,0,1) etc                                (  6 points)
!vw    code=2:   (0,a,a) etc, a=1/sqrt(2)                   ( 12 points)
!vw    code=3:   (a,a,a) etc, a=1/sqrt(3)                   (  8 points)
!vw    code=4:   (a,a,b) etc, b=sqrt(1-2 a^2)               ( 24 points)
!vw    code=5:   (a,b,0) etc, b=sqrt(1-a^2), a input        ( 24 points)
!vw    code=6:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input  ( 48 points)
!vw
       goto (1,2,3,4,5,6) code
       write (6,*) 'Gen_Oh: Invalid Code'
       stop 
    1  continue
       a=1.0d0
       x(1) =  a
       y(1) =  0.0d0
       z(1) =  0.0d0
       w(1) =  v
       x(2) = -a
       y(2) =  0.0d0
       z(2) =  0.0d0
       w(2) =  v
       x(3) =  0.0d0
       y(3) =  a
       z(3) =  0.0d0
       w(3) =  v
       x(4) =  0.0d0
       y(4) = -a
       z(4) =  0.0d0
       w(4) =  v
       x(5) =  0.0d0
       y(5) =  0.0d0
       z(5) =  a
       w(5) =  v
       x(6) =  0.0d0
       y(6) =  0.0d0
       z(6) = -a
       w(6) =  v
       num=num+6
       return
!vw
    2  continue
       a=sqrt(0.5d0)
       x( 1) =  0d0
       y( 1) =  a
       z( 1) =  a
       w( 1) =  v
       x( 2) =  0d0
       y( 2) = -a
       z( 2) =  a
       w( 2) =  v
       x( 3) =  0d0
       y( 3) =  a
       z( 3) = -a
       w( 3) =  v
       x( 4) =  0d0
       y( 4) = -a
       z( 4) = -a
       w( 4) =  v
       x( 5) =  a
       y( 5) =  0d0
       z( 5) =  a
       w( 5) =  v
       x( 6) = -a
       y( 6) =  0d0
       z( 6) =  a
       w( 6) =  v
       x( 7) =  a
       y( 7) =  0d0
       z( 7) = -a
       w( 7) =  v
       x( 8) = -a
       y( 8) =  0d0
       z( 8) = -a
       w( 8) =  v
       x( 9) =  a
       y( 9) =  a
       z( 9) =  0d0
       w( 9) =  v
       x(10) = -a
       y(10) =  a
       z(10) =  0d0
       w(10) =  v
       x(11) =  a
       y(11) = -a
       z(11) =  0d0
       w(11) =  v
       x(12) = -a
       y(12) = -a
       z(12) =  0d0
       w(12) =  v
       num=num+12
       return
!vw
    3  continue
       a = sqrt(1d0/3d0)
       x(1) =  a
       y(1) =  a
       z(1) =  a
       w(1) =  v
       x(2) = -a
       y(2) =  a
       z(2) =  a
       w(2) =  v
       x(3) =  a
       y(3) = -a
       z(3) =  a
       w(3) =  v
       x(4) = -a
       y(4) = -a
       z(4) =  a
       w(4) =  v
       x(5) =  a
       y(5) =  a
       z(5) = -a
       w(5) =  v
       x(6) = -a
       y(6) =  a
       z(6) = -a
       w(6) =  v
       x(7) =  a
       y(7) = -a
       z(7) = -a
       w(7) =  v
       x(8) = -a
       y(8) = -a
       z(8) = -a
       w(8) =  v
       num=num+8
       return
!vw
    4  continue
       b = sqrt(1d0 - 2d0*a*a)
       x( 1) =  a
       y( 1) =  a
       z( 1) =  b
       w( 1) =  v
       x( 2) = -a
       y( 2) =  a
       z( 2) =  b
       w( 2) =  v
       x( 3) =  a
       y( 3) = -a
       z( 3) =  b
       w( 3) =  v
       x( 4) = -a
       y( 4) = -a
       z( 4) =  b
       w( 4) =  v
       x( 5) =  a
       y( 5) =  a
       z( 5) = -b
       w( 5) =  v
       x( 6) = -a
       y( 6) =  a
       z( 6) = -b
       w( 6) =  v
       x( 7) =  a
       y( 7) = -a
       z( 7) = -b
       w( 7) =  v
       x( 8) = -a
       y( 8) = -a
       z( 8) = -b
       w( 8) =  v
       x( 9) =  a
       y( 9) =  b
       z( 9) =  a
       w( 9) =  v
       x(10) = -a
       y(10) =  b
       z(10) =  a
       w(10) =  v
       x(11) =  a
       y(11) = -b
       z(11) =  a
       w(11) =  v
       x(12) = -a
       y(12) = -b
       z(12) =  a
       w(12) =  v
       x(13) =  a
       y(13) =  b
       z(13) = -a
       w(13) =  v
       x(14) = -a
       y(14) =  b
       z(14) = -a
       w(14) =  v
       x(15) =  a
       y(15) = -b
       z(15) = -a
       w(15) =  v
       x(16) = -a
       y(16) = -b
       z(16) = -a
       w(16) =  v
       x(17) =  b
       y(17) =  a
       z(17) =  a
       w(17) =  v
       x(18) = -b
       y(18) =  a
       z(18) =  a
       w(18) =  v
       x(19) =  b
       y(19) = -a
       z(19) =  a
       w(19) =  v
       x(20) = -b
       y(20) = -a
       z(20) =  a
       w(20) =  v
       x(21) =  b
       y(21) =  a
       z(21) = -a
       w(21) =  v
       x(22) = -b
       y(22) =  a
       z(22) = -a
       w(22) =  v
       x(23) =  b
       y(23) = -a
       z(23) = -a
       w(23) =  v
       x(24) = -b
       y(24) = -a
       z(24) = -a
       w(24) =  v
       num=num+24
       return
!vw
    5  continue
       b=sqrt(1d0-a*a)
       x( 1) =  a
       y( 1) =  b
       z( 1) =  0d0
       w( 1) =  v
       x( 2) = -a
       y( 2) =  b
       z( 2) =  0d0
       w( 2) =  v
       x( 3) =  a
       y( 3) = -b
       z( 3) =  0d0
       w( 3) =  v
       x( 4) = -a
       y( 4) = -b
       z( 4) =  0d0
       w( 4) =  v
       x( 5) =  b
       y( 5) =  a
       z( 5) =  0d0
       w( 5) =  v
       x( 6) = -b
       y( 6) =  a
       z( 6) =  0d0
       w( 6) =  v
       x( 7) =  b
       y( 7) = -a
       z( 7) =  0d0
       w( 7) =  v
       x( 8) = -b
       y( 8) = -a
       z( 8) =  0d0
       w( 8) =  v
       x( 9) =  a
       y( 9) =  0d0
       z( 9) =  b
       w( 9) =  v
       x(10) = -a
       y(10) =  0d0
       z(10) =  b
       w(10) =  v
       x(11) =  a
       y(11) =  0d0
       z(11) = -b
       w(11) =  v
       x(12) = -a
       y(12) =  0d0
       z(12) = -b
       w(12) =  v
       x(13) =  b
       y(13) =  0d0
       z(13) =  a
       w(13) =  v
       x(14) = -b
       y(14) =  0d0
       z(14) =  a
       w(14) =  v
       x(15) =  b
       y(15) =  0d0
       z(15) = -a
       w(15) =  v
       x(16) = -b
       y(16) =  0d0
       z(16) = -a
       w(16) =  v
       x(17) =  0d0
       y(17) =  a
       z(17) =  b
       w(17) =  v
       x(18) =  0d0
       y(18) = -a
       z(18) =  b
       w(18) =  v
       x(19) =  0d0
       y(19) =  a
       z(19) = -b
       w(19) =  v
       x(20) =  0d0
       y(20) = -a
       z(20) = -b
       w(20) =  v
       x(21) =  0d0
       y(21) =  b
       z(21) =  a
       w(21) =  v
       x(22) =  0d0
       y(22) = -b
       z(22) =  a
       w(22) =  v
       x(23) =  0d0
       y(23) =  b
       z(23) = -a
       w(23) =  v
       x(24) =  0d0
       y(24) = -b
       z(24) = -a
       w(24) =  v
       num=num+24
       return
!vw
    6  continue
       c=sqrt(1d0 - a*a - b*b)
       x( 1) =  a
       y( 1) =  b
       z( 1) =  c
       w( 1) =  v
       x( 2) = -a
       y( 2) =  b
       z( 2) =  c
       w( 2) =  v
       x( 3) =  a
       y( 3) = -b
       z( 3) =  c
       w( 3) =  v
       x( 4) = -a
       y( 4) = -b
       z( 4) =  c
       w( 4) =  v
       x( 5) =  a
       y( 5) =  b
       z( 5) = -c
       w( 5) =  v
       x( 6) = -a
       y( 6) =  b
       z( 6) = -c
       w( 6) =  v
       x( 7) =  a
       y( 7) = -b
       z( 7) = -c
       w( 7) =  v
       x( 8) = -a
       y( 8) = -b
       z( 8) = -c
       w( 8) =  v
       x( 9) =  a
       y( 9) =  c
       z( 9) =  b
       w( 9) =  v
       x(10) = -a
       y(10) =  c
       z(10) =  b
       w(10) =  v
       x(11) =  a
       y(11) = -c
       z(11) =  b
       w(11) =  v
       x(12) = -a
       y(12) = -c
       z(12) =  b
       w(12) =  v
       x(13) =  a
       y(13) =  c
       z(13) = -b
       w(13) =  v
       x(14) = -a
       y(14) =  c
       z(14) = -b
       w(14) =  v
       x(15) =  a
       y(15) = -c
       z(15) = -b
       w(15) =  v
       x(16) = -a
       y(16) = -c
       z(16) = -b
       w(16) =  v
       x(17) =  b
       y(17) =  a
       z(17) =  c
       w(17) =  v
       x(18) = -b
       y(18) =  a
       z(18) =  c
       w(18) =  v
       x(19) =  b
       y(19) = -a
       z(19) =  c
       w(19) =  v
       x(20) = -b
       y(20) = -a
       z(20) =  c
       w(20) =  v
       x(21) =  b
       y(21) =  a
       z(21) = -c
       w(21) =  v
       x(22) = -b
       y(22) =  a
       z(22) = -c
       w(22) =  v
       x(23) =  b
       y(23) = -a
       z(23) = -c
       w(23) =  v
       x(24) = -b
       y(24) = -a
       z(24) = -c
       w(24) =  v
       x(25) =  b
       y(25) =  c
       z(25) =  a
       w(25) =  v
       x(26) = -b
       y(26) =  c
       z(26) =  a
       w(26) =  v
       x(27) =  b
       y(27) = -c
       z(27) =  a
       w(27) =  v
       x(28) = -b
       y(28) = -c
       z(28) =  a
       w(28) =  v
       x(29) =  b
       y(29) =  c
       z(29) = -a
       w(29) =  v
       x(30) = -b
       y(30) =  c
       z(30) = -a
       w(30) =  v
       x(31) =  b
       y(31) = -c
       z(31) = -a
       w(31) =  v
       x(32) = -b
       y(32) = -c
       z(32) = -a
       w(32) =  v
       x(33) =  c
       y(33) =  a
       z(33) =  b
       w(33) =  v
       x(34) = -c
       y(34) =  a
       z(34) =  b
       w(34) =  v
       x(35) =  c
       y(35) = -a
       z(35) =  b
       w(35) =  v
       x(36) = -c
       y(36) = -a
       z(36) =  b
       w(36) =  v
       x(37) =  c
       y(37) =  a
       z(37) = -b
       w(37) =  v
       x(38) = -c
       y(38) =  a
       z(38) = -b
       w(38) =  v
       x(39) =  c
       y(39) = -a
       z(39) = -b
       w(39) =  v
       x(40) = -c
       y(40) = -a
       z(40) = -b
       w(40) =  v
       x(41) =  c
       y(41) =  b
       z(41) =  a
       w(41) =  v
       x(42) = -c
       y(42) =  b
       z(42) =  a
       w(42) =  v
       x(43) =  c
       y(43) = -b
       z(43) =  a
       w(43) =  v
       x(44) = -c
       y(44) = -b
       z(44) =  a
       w(44) =  v
       x(45) =  c
       y(45) =  b
       z(45) = -a
       w(45) =  v
       x(46) = -c
       y(46) =  b
       z(46) = -a
       w(46) =  v
       x(47) =  c
       y(47) = -b
       z(47) = -a
       w(47) =  v
       x(48) = -c
       y(48) = -b
       z(48) = -a
       w(48) =  v
       num=num+48
       return
  end subroutine

  SUBROUTINE LD0006(X,Y,Z,W,N)
       REAL*8 X(   6)
       REAL*8 Y(   6)
       REAL*8 Z(   6)
       REAL*8 W(   6)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV    6-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.1666666666666667D+0
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0014(X,Y,Z,W,N)
       REAL*8 X(  14)
       REAL*8 Y(  14)
       REAL*8 Z(  14)
       REAL*8 W(  14)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV   14-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.6666666666666667D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7500000000000000D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0026(X,Y,Z,W,N)
       REAL*8 X(  26)
       REAL*8 Y(  26)
       REAL*8 Z(  26)
       REAL*8 W(  26)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV   26-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.4761904761904762D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3809523809523810D-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0038(X,Y,Z,W,N)
       REAL*8 X(  38)
       REAL*8 Y(  38)
       REAL*8 Z(  38)
       REAL*8 W(  38)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV   38-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.9523809523809524D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4597008433809831D+0
       V=0.2857142857142857D-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0050(X,Y,Z,W,N)
       REAL*8 X(  50)
       REAL*8 Y(  50)
       REAL*8 Z(  50)
       REAL*8 W(  50)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV   50-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.1269841269841270D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2257495590828924D-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2109375000000000D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3015113445777636D+0
       V=0.2017333553791887D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0074(X,Y,Z,W,N)
       REAL*8 X(  74)
       REAL*8 Y(  74)
       REAL*8 Z(  74)
       REAL*8 W(  74)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV   74-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.5130671797338464D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1660406956574204D-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2958603896103896D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4803844614152614D+0
       V=0.2657620708215946D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3207726489807764D+0
       V=0.1652217099371571D-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0086(X,Y,Z,W,N)
       REAL*8 X(  86)
       REAL*8 Y(  86)
       REAL*8 Z(  86)
       REAL*8 W(  86)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV   86-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.1154401154401154D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1194390908585628D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3696028464541502D+0
       V=0.1111055571060340D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6943540066026664D+0
       V=0.1187650129453714D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3742430390903412D+0
       V=0.1181230374690448D-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0110(X,Y,Z,W,N)
       REAL*8 X( 110)
       REAL*8 Y( 110)
       REAL*8 Z( 110)
       REAL*8 W( 110)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  110-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.3828270494937162D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9793737512487512D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1851156353447362D+0
       V=0.8211737283191111D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6904210483822922D+0
       V=0.9942814891178103D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3956894730559419D+0
       V=0.9595471336070963D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4783690288121502D+0
       V=0.9694996361663028D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0146(X,Y,Z,W,N)
       REAL*8 X( 146)
       REAL*8 Y( 146)
       REAL*8 Z( 146)
       REAL*8 W( 146)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  146-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.5996313688621381D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7372999718620756D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7210515360144488D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6764410400114264D+0
       V=0.7116355493117555D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4174961227965453D+0
       V=0.6753829486314477D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1574676672039082D+0
       V=0.7574394159054034D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1403553811713183D+0
       B=0.4493328323269557D+0
       V=0.6991087353303262D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0170(X,Y,Z,W,N)
       REAL*8 X( 170)
       REAL*8 Y( 170)
       REAL*8 Z( 170)
       REAL*8 W( 170)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  170-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.5544842902037365D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6071332770670752D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6383674773515093D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2551252621114134D+0
       V=0.5183387587747790D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6743601460362766D+0
       V=0.6317929009813725D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318910696719410D+0
       V=0.6201670006589077D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2613931360335988D+0
       V=0.5477143385137348D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4990453161796037D+0
       B=0.1446630744325115D+0
       V=0.5968383987681156D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0194(X,Y,Z,W,N)
       REAL*8 X( 194)
       REAL*8 Y( 194)
       REAL*8 Z( 194)
       REAL*8 W( 194)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  194-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.1782340447244611D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5716905949977102D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5573383178848738D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6712973442695226D+0
       V=0.5608704082587997D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2892465627575439D+0
       V=0.5158237711805383D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4446933178717437D+0
       V=0.5518771467273614D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1299335447650067D+0
       V=0.4106777028169394D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3457702197611283D+0
       V=0.5051846064614808D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1590417105383530D+0
       B=0.8360360154824589D+0
       V=0.5530248916233094D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0230(X,Y,Z,W,N)
       REAL*8 X( 230)
       REAL*8 Y( 230)
       REAL*8 Z( 230)
       REAL*8 W( 230)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  230-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=-0.5522639919727325D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4450274607445226D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4492044687397611D+0
       V=0.4496841067921404D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2520419490210201D+0
       V=0.5049153450478750D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6981906658447242D+0
       V=0.3976408018051883D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6587405243460960D+0
       V=0.4401400650381014D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4038544050097660D-1
       V=0.1724544350544401D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5823842309715585D+0
       V=0.4231083095357343D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3545877390518688D+0
       V=0.5198069864064399D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2272181808998187D+0
       B=0.4864661535886647D+0
       V=0.4695720972568883D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0266(X,Y,Z,W,N)
       REAL*8 X( 266)
       REAL*8 Y( 266)
       REAL*8 Z( 266)
       REAL*8 W( 266)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  266-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=-0.1313769127326952D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2522728704859336D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4186853881700583D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7039373391585475D+0
       V=0.5315167977810885D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1012526248572414D+0
       V=0.4047142377086219D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4647448726420539D+0
       V=0.4112482394406990D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3277420654971629D+0
       V=0.3595584899758782D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6620338663699974D+0
       V=0.4256131351428158D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8506508083520399D+0
       V=0.4229582700647240D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3233484542692899D+0
       B=0.1153112011009701D+0
       V=0.4080914225780505D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2314790158712601D+0
       B=0.5244939240922365D+0
       V=0.4071467593830964D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0302(X,Y,Z,W,N)
       REAL*8 X( 302)
       REAL*8 Y( 302)
       REAL*8 Z( 302)
       REAL*8 W( 302)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  302-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.8545911725128148D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3599119285025571D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3515640345570105D+0
       V=0.3449788424305883D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6566329410219612D+0
       V=0.3604822601419882D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4729054132581005D+0
       V=0.3576729661743367D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9618308522614784D-1
       V=0.2352101413689164D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2219645236294178D+0
       V=0.3108953122413675D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7011766416089545D+0
       V=0.3650045807677255D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2644152887060663D+0
       V=0.2982344963171804D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5718955891878961D+0
       V=0.3600820932216460D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2510034751770465D+0
       B=0.8000727494073952D+0
       V=0.3571540554273387D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1233548532583327D+0
       B=0.4127724083168531D+0
       V=0.3392312205006170D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0350(X,Y,Z,W,N)
       REAL*8 X( 350)
       REAL*8 Y( 350)
       REAL*8 Z( 350)
       REAL*8 W( 350)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  350-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.3006796749453936D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3050627745650771D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7068965463912316D+0
       V=0.1621104600288991D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4794682625712025D+0
       V=0.3005701484901752D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1927533154878019D+0
       V=0.2990992529653774D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6930357961327123D+0
       V=0.2982170644107595D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3608302115520091D+0
       V=0.2721564237310992D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6498486161496169D+0
       V=0.3033513795811141D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1932945013230339D+0
       V=0.3007949555218533D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3800494919899303D+0
       V=0.2881964603055307D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2899558825499574D+0
       B=0.7934537856582316D+0
       V=0.2958357626535696D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9684121455103957D-1
       B=0.8280801506686862D+0
       V=0.3036020026407088D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1833434647041659D+0
       B=0.9074658265305127D+0
       V=0.2832187403926303D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0434(X,Y,Z,W,N)
       REAL*8 X( 434)
       REAL*8 Y( 434)
       REAL*8 Z( 434)
       REAL*8 W( 434)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  434-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.5265897968224436D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2548219972002607D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2512317418927307D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6909346307509111D+0
       V=0.2530403801186355D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1774836054609158D+0
       V=0.2014279020918528D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4914342637784746D+0
       V=0.2501725168402936D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6456664707424256D+0
       V=0.2513267174597564D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2861289010307638D+0
       V=0.2302694782227416D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7568084367178018D-1
       V=0.1462495621594614D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3927259763368002D+0
       V=0.2445373437312980D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8818132877794288D+0
       V=0.2417442375638981D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9776428111182649D+0
       V=0.1910951282179532D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2054823696403044D+0
       B=0.8689460322872412D+0
       V=0.2416930044324775D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5905157048925271D+0
       B=0.7999278543857286D+0
       V=0.2512236854563495D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5550152361076807D+0
       B=0.7717462626915901D+0
       V=0.2496644054553086D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9371809858553722D+0
       B=0.3344363145343455D+0
       V=0.2236607760437849D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0590(X,Y,Z,W,N)
       REAL*8 X( 590)
       REAL*8 Y( 590)
       REAL*8 Z( 590)
       REAL*8 W( 590)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  590-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.3095121295306187D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1852379698597489D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7040954938227469D+0
       V=0.1871790639277744D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6807744066455243D+0
       V=0.1858812585438317D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6372546939258752D+0
       V=0.1852028828296213D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5044419707800358D+0
       V=0.1846715956151242D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4215761784010967D+0
       V=0.1818471778162769D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3317920736472123D+0
       V=0.1749564657281154D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2384736701421887D+0
       V=0.1617210647254411D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1459036449157763D+0
       V=0.1384737234851692D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6095034115507196D-1
       V=0.9764331165051050D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6116843442009876D+0
       V=0.1857161196774078D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3964755348199858D+0
       V=0.1705153996395864D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1724782009907724D+0
       V=0.1300321685886048D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5610263808622060D+0
       B=0.3518280927733519D+0
       V=0.1842866472905286D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4742392842551980D+0
       B=0.2634716655937950D+0
       V=0.1802658934377451D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5984126497885380D+0
       B=0.1816640840360209D+0
       V=0.1849830560443660D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3791035407695563D+0
       B=0.1720795225656878D+0
       V=0.1713904507106709D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2778673190586244D+0
       B=0.8213021581932511D-1
       V=0.1555213603396808D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5033564271075117D+0
       B=0.8999205842074875D-1
       V=0.1802239128008525D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0770(X,Y,Z,W,N)
       REAL*8 X( 770)
       REAL*8 Y( 770)
       REAL*8 Z( 770)
       REAL*8 W( 770)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  770-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.2192942088181184D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1436433617319080D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1421940344335877D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5087204410502360D-1
       V=0.6798123511050502D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1228198790178831D+0
       V=0.9913184235294912D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2026890814408786D+0
       V=0.1180207833238949D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2847745156464294D+0
       V=0.1296599602080921D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3656719078978026D+0
       V=0.1365871427428316D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4428264886713469D+0
       V=0.1402988604775325D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5140619627249735D+0
       V=0.1418645563595609D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6306401219166803D+0
       V=0.1421376741851662D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6716883332022612D+0
       V=0.1423996475490962D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6979792685336881D+0
       V=0.1431554042178567D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1446865674195309D+0
       V=0.9254401499865368D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3390263475411216D+0
       V=0.1250239995053509D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5335804651263506D+0
       V=0.1394365843329230D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6944024393349413D-1
       B=0.2355187894242326D+0
       V=0.1127089094671749D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2269004109529460D+0
       B=0.4102182474045730D+0
       V=0.1345753760910670D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8025574607775339D-1
       B=0.6214302417481605D+0
       V=0.1424957283316783D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1467999527896572D+0
       B=0.3245284345717394D+0
       V=0.1261523341237750D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1571507769824727D+0
       B=0.5224482189696630D+0
       V=0.1392547106052696D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2365702993157246D+0
       B=0.6017546634089558D+0
       V=0.1418761677877656D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7714815866765732D-1
       B=0.4346575516141163D+0
       V=0.1338366684479554D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3062936666210730D+0
       B=0.4908826589037616D+0
       V=0.1393700862676131D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3822477379524787D+0
       B=0.5648768149099500D+0
       V=0.1415914757466932D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD0974(X,Y,Z,W,N)
       REAL*8 X( 974)
       REAL*8 Y( 974)
       REAL*8 Z( 974)
       REAL*8 W( 974)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV  974-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.1438294190527431D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1125772288287004D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4292963545341347D-1
       V=0.4948029341949241D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1051426854086404D+0
       V=0.7357990109125470D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1750024867623087D+0
       V=0.8889132771304384D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2477653379650257D+0
       V=0.9888347838921435D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3206567123955957D+0
       V=0.1053299681709471D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3916520749849983D+0
       V=0.1092778807014578D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4590825874187624D+0
       V=0.1114389394063227D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5214563888415861D+0
       V=0.1123724788051555D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6253170244654199D+0
       V=0.1125239325243814D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6637926744523170D+0
       V=0.1126153271815905D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6910410398498301D+0
       V=0.1130286931123841D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7052907007457760D+0
       V=0.1134986534363955D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1236686762657990D+0
       V=0.6823367927109931D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2940777114468387D+0
       V=0.9454158160447096D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4697753849207649D+0
       V=0.1074429975385679D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6334563241139567D+0
       V=0.1129300086569132D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5974048614181342D-1
       B=0.2029128752777523D+0
       V=0.8436884500901954D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1375760408473636D+0
       B=0.4602621942484054D+0
       V=0.1075255720448885D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3391016526336286D+0
       B=0.5030673999662036D+0
       V=0.1108577236864462D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1271675191439820D+0
       B=0.2817606422442134D+0
       V=0.9566475323783357D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2693120740413512D+0
       B=0.4331561291720157D+0
       V=0.1080663250717391D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1419786452601918D+0
       B=0.6256167358580814D+0
       V=0.1126797131196295D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6709284600738255D-1
       B=0.3798395216859157D+0
       V=0.1022568715358061D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7057738183256172D-1
       B=0.5517505421423520D+0
       V=0.1108960267713108D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2783888477882155D+0
       B=0.6029619156159187D+0
       V=0.1122790653435766D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1979578938917407D+0
       B=0.3589606329589096D+0
       V=0.1032401847117460D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2087307061103274D+0
       B=0.5348666438135476D+0
       V=0.1107249382283854D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4055122137872836D+0
       B=0.5674997546074373D+0
       V=0.1121780048519972D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD1202(X,Y,Z,W,N)
       REAL*8 X(1202)
       REAL*8 Y(1202)
       REAL*8 Z(1202)
       REAL*8 W(1202)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 1202-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.1105189233267572D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9205232738090741D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9133159786443561D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3712636449657089D-1
       V=0.3690421898017899D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9140060412262223D-1
       V=0.5603990928680660D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1531077852469906D+0
       V=0.6865297629282609D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2180928891660612D+0
       V=0.7720338551145630D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2839874532200175D+0
       V=0.8301545958894795D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3491177600963764D+0
       V=0.8686692550179628D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4121431461444309D+0
       V=0.8927076285846890D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4718993627149127D+0
       V=0.9060820238568219D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5273145452842337D+0
       V=0.9119777254940867D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6209475332444019D+0
       V=0.9128720138604181D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6569722711857291D+0
       V=0.9130714935691735D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6841788309070143D+0
       V=0.9152873784554116D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012604330123631D+0
       V=0.9187436274321654D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1072382215478166D+0
       V=0.5176977312965694D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2582068959496968D+0
       V=0.7331143682101417D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4172752955306717D+0
       V=0.8463232836379928D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5700366911792503D+0
       V=0.9031122694253992D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9827986018263947D+0
       B=0.1771774022615325D+0
       V=0.6485778453163257D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9624249230326228D+0
       B=0.2475716463426288D+0
       V=0.7435030910982369D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9402007994128811D+0
       B=0.3354616289066489D+0
       V=0.7998527891839054D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9320822040143202D+0
       B=0.3173615246611977D+0
       V=0.8101731497468018D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9043674199393299D+0
       B=0.4090268427085357D+0
       V=0.8483389574594331D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8912407560074747D+0
       B=0.3854291150669224D+0
       V=0.8556299257311812D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8676435628462708D+0
       B=0.4932221184851285D+0
       V=0.8803208679738260D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8581979986041619D+0
       B=0.4785320675922435D+0
       V=0.8811048182425720D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8396753624049856D+0
       B=0.4507422593157064D+0
       V=0.8850282341265444D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8165288564022188D+0
       B=0.5632123020762100D+0
       V=0.9021342299040653D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8015469370783529D+0
       B=0.5434303569693900D+0
       V=0.9010091677105086D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7773563069070351D+0
       B=0.5123518486419871D+0
       V=0.9022692938426915D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7661621213900394D+0
       B=0.6394279634749102D+0
       V=0.9158016174693465D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7553584143533510D+0
       B=0.6269805509024392D+0
       V=0.9131578003189435D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7344305757559503D+0
       B=0.6031161693096310D+0
       V=0.9107813579482705D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7043837184021765D+0
       B=0.5693702498468441D+0
       V=0.9105760258970126D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD1454(X,Y,Z,W,N)
       REAL*8 X(1454)
       REAL*8 Y(1454)
       REAL*8 Z(1454)
       REAL*8 W(1454)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 1454-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.7777160743261247D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7557646413004701D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3229290663413854D-1
       V=0.2841633806090617D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8036733271462222D-1
       V=0.4374419127053555D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1354289960531653D+0
       V=0.5417174740872172D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1938963861114426D+0
       V=0.6148000891358593D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2537343715011275D+0
       V=0.6664394485800705D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3135251434752570D+0
       V=0.7025039356923220D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3721558339375338D+0
       V=0.7268511789249627D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4286809575195696D+0
       V=0.7422637534208629D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4822510128282994D+0
       V=0.7509545035841214D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5320679333566263D+0
       V=0.7548535057718401D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6172998195394274D+0
       V=0.7554088969774001D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6510679849127481D+0
       V=0.7553147174442808D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6777315251687360D+0
       V=0.7564767653292297D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6963109410648741D+0
       V=0.7587991808518730D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7058935009831749D+0
       V=0.7608261832033027D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9955546194091857D+0
       V=0.4021680447874916D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9734115901794209D+0
       V=0.5804871793945964D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9275693732388626D+0
       V=0.6792151955945159D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8568022422795103D+0
       V=0.7336741211286294D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7623495553719372D+0
       V=0.7581866300989608D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5707522908892223D+0
       B=0.4387028039889501D+0
       V=0.7538257859800743D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5196463388403083D+0
       B=0.3858908414762617D+0
       V=0.7483517247053123D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4646337531215351D+0
       B=0.3301937372343854D+0
       V=0.7371763661112059D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4063901697557691D+0
       B=0.2725423573563777D+0
       V=0.7183448895756934D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3456329466643087D+0
       B=0.2139510237495250D+0
       V=0.6895815529822191D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2831395121050332D+0
       B=0.1555922309786647D+0
       V=0.6480105801792886D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2197682022925330D+0
       B=0.9892878979686097D-1
       V=0.5897558896594636D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1564696098650355D+0
       B=0.4598642910675510D-1
       V=0.5095708849247346D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6027356673721295D+0
       B=0.3376625140173426D+0
       V=0.7536906428909755D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5496032320255096D+0
       B=0.2822301309727988D+0
       V=0.7472505965575118D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4921707755234567D+0
       B=0.2248632342592540D+0
       V=0.7343017132279698D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4309422998598483D+0
       B=0.1666224723456479D+0
       V=0.7130871582177445D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3664108182313672D+0
       B=0.1086964901822169D+0
       V=0.6817022032112776D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2990189057758436D+0
       B=0.5251989784120085D-1
       V=0.6380941145604121D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6268724013144998D+0
       B=0.2297523657550023D+0
       V=0.7550381377920310D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5707324144834607D+0
       B=0.1723080607093800D+0
       V=0.7478646640144802D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5096360901960365D+0
       B=0.1140238465390513D+0
       V=0.7335918720601220D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4438729938312456D+0
       B=0.5611522095882537D-1
       V=0.7110120527658118D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6419978471082389D+0
       B=0.1164174423140873D+0
       V=0.7571363978689501D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5817218061802611D+0
       B=0.5797589531445219D-1
       V=0.7489908329079234D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD1730(X,Y,Z,W,N)
       REAL*8 X(1730)
       REAL*8 Y(1730)
       REAL*8 Z(1730)
       REAL*8 W(1730)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 1730-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.6309049437420976D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6398287705571748D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6357185073530720D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2860923126194662D-1
       V=0.2221207162188168D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7142556767711522D-1
       V=0.3475784022286848D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1209199540995559D+0
       V=0.4350742443589804D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1738673106594379D+0
       V=0.4978569136522127D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2284645438467734D+0
       V=0.5435036221998053D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2834807671701512D+0
       V=0.5765913388219542D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3379680145467339D+0
       V=0.6001200359226003D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3911355454819537D+0
       V=0.6162178172717512D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4422860353001403D+0
       V=0.6265218152438485D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4907781568726057D+0
       V=0.6323987160974212D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5360006153211468D+0
       V=0.6350767851540569D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6142105973596603D+0
       V=0.6354362775297107D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6459300387977504D+0
       V=0.6352302462706235D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718056125089225D+0
       V=0.6358117881417972D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6910888533186254D+0
       V=0.6373101590310117D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7030467416823252D+0
       V=0.6390428961368665D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8354951166354646D-1
       V=0.3186913449946576D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2050143009099486D+0
       V=0.4678028558591711D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3370208290706637D+0
       V=0.5538829697598626D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4689051484233963D+0
       V=0.6044475907190476D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5939400424557334D+0
       V=0.6313575103509012D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1394983311832261D+0
       B=0.4097581162050343D-1
       V=0.4078626431855630D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1967999180485014D+0
       B=0.8851987391293348D-1
       V=0.4759933057812725D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2546183732548967D+0
       B=0.1397680182969819D+0
       V=0.5268151186413440D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3121281074713875D+0
       B=0.1929452542226526D+0
       V=0.5643048560507316D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3685981078502492D+0
       B=0.2467898337061562D+0
       V=0.5914501076613073D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4233760321547856D+0
       B=0.3003104124785409D+0
       V=0.6104561257874195D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4758671236059246D+0
       B=0.3526684328175033D+0
       V=0.6230252860707806D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5255178579796463D+0
       B=0.4031134861145713D+0
       V=0.6305618761760796D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5718025633734589D+0
       B=0.4509426448342351D+0
       V=0.6343092767597889D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2686927772723415D+0
       B=0.4711322502423248D-1
       V=0.5176268945737826D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3306006819904809D+0
       B=0.9784487303942695D-1
       V=0.5564840313313692D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3904906850594983D+0
       B=0.1505395810025273D+0
       V=0.5856426671038980D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4479957951904390D+0
       B=0.2039728156296050D+0
       V=0.6066386925777091D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5027076848919780D+0
       B=0.2571529941121107D+0
       V=0.6208824962234458D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5542087392260217D+0
       B=0.3092191375815670D+0
       V=0.6296314297822907D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6020850887375187D+0
       B=0.3593807506130276D+0
       V=0.6340423756791859D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4019851409179594D+0
       B=0.5063389934378671D-1
       V=0.5829627677107342D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4635614567449800D+0
       B=0.1032422269160612D+0
       V=0.6048693376081110D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5215860931591575D+0
       B=0.1566322094006254D+0
       V=0.6202362317732461D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5758202499099271D+0
       B=0.2098082827491099D+0
       V=0.6299005328403779D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6259893683876795D+0
       B=0.2618824114553391D+0
       V=0.6347722390609353D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5313795124811891D+0
       B=0.5263245019338556D-1
       V=0.6203778981238834D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5893317955931995D+0
       B=0.1061059730982005D+0
       V=0.6308414671239979D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6426246321215801D+0
       B=0.1594171564034221D+0
       V=0.6362706466959498D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6511904367376113D+0
       B=0.5354789536565540D-1
       V=0.6375414170333233D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD2030(X,Y,Z,W,N)
       REAL*8 X(2030)
       REAL*8 Y(2030)
       REAL*8 Z(2030)
       REAL*8 W(2030)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 2030-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.4656031899197431D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5421549195295507D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2540835336814348D-1
       V=0.1778522133346553D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6399322800504915D-1
       V=0.2811325405682796D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1088269469804125D+0
       V=0.3548896312631459D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1570670798818287D+0
       V=0.4090310897173364D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2071163932282514D+0
       V=0.4493286134169965D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2578914044450844D+0
       V=0.4793728447962723D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3085687558169623D+0
       V=0.5015415319164265D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3584719706267024D+0
       V=0.5175127372677937D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4070135594428709D+0
       V=0.5285522262081019D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4536618626222638D+0
       V=0.5356832703713962D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4979195686463577D+0
       V=0.5397914736175170D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5393075111126999D+0
       V=0.5416899441599930D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6115617676843916D+0
       V=0.5419308476889938D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6414308435160159D+0
       V=0.5416936902030596D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6664099412721607D+0
       V=0.5419544338703164D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6859161771214913D+0
       V=0.5428983656630975D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6993625593503890D+0
       V=0.5442286500098193D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7062393387719380D+0
       V=0.5452250345057301D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7479028168349763D-1
       V=0.2568002497728530D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1848951153969366D+0
       V=0.3827211700292145D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3059529066581305D+0
       V=0.4579491561917824D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4285556101021362D+0
       V=0.5042003969083574D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5468758653496526D+0
       V=0.5312708889976025D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6565821978343439D+0
       V=0.5438401790747117D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1253901572367117D+0
       B=0.3681917226439641D-1
       V=0.3316041873197344D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1775721510383941D+0
       B=0.7982487607213301D-1
       V=0.3899113567153771D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2305693358216114D+0
       B=0.1264640966592335D+0
       V=0.4343343327201309D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2836502845992063D+0
       B=0.1751585683418957D+0
       V=0.4679415262318919D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3361794746232590D+0
       B=0.2247995907632670D+0
       V=0.4930847981631031D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3875979172264824D+0
       B=0.2745299257422246D+0
       V=0.5115031867540091D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4374019316999074D+0
       B=0.3236373482441118D+0
       V=0.5245217148457367D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4851275843340022D+0
       B=0.3714967859436741D+0
       V=0.5332041499895321D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5303391803806868D+0
       B=0.4175353646321745D+0
       V=0.5384583126021542D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726197380596287D+0
       B=0.4612084406355461D+0
       V=0.5411067210798852D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2431520732564863D+0
       B=0.4258040133043952D-1
       V=0.4259797391468714D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3002096800895869D+0
       B=0.8869424306722721D-1
       V=0.4604931368460021D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3558554457457432D+0
       B=0.1368811706510655D+0
       V=0.4871814878255202D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4097782537048887D+0
       B=0.1860739985015033D+0
       V=0.5072242910074885D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4616337666067458D+0
       B=0.2354235077395853D+0
       V=0.5217069845235350D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5110707008417874D+0
       B=0.2842074921347011D+0
       V=0.5315785966280310D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5577415286163795D+0
       B=0.3317784414984102D+0
       V=0.5376833708758905D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6013060431366950D+0
       B=0.3775299002040700D+0
       V=0.5408032092069521D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3661596767261781D+0
       B=0.4599367887164592D-1
       V=0.4842744917904866D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4237633153506581D+0
       B=0.9404893773654421D-1
       V=0.5048926076188130D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4786328454658452D+0
       B=0.1431377109091971D+0
       V=0.5202607980478373D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5305702076789774D+0
       B=0.1924186388843570D+0
       V=0.5309932388325743D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5793436224231788D+0
       B=0.2411590944775190D+0
       V=0.5377419770895208D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6247069017094747D+0
       B=0.2886871491583605D+0
       V=0.5411696331677717D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4874315552535204D+0
       B=0.4804978774953206D-1
       V=0.5197996293282420D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5427337322059053D+0
       B=0.9716857199366665D-1
       V=0.5311120836622945D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5943493747246700D+0
       B=0.1465205839795055D+0
       V=0.5384309319956951D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6421314033564943D+0
       B=0.1953579449803574D+0
       V=0.5421859504051886D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6020628374713980D+0
       B=0.4916375015738108D-1
       V=0.5390948355046314D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6529222529856881D+0
       B=0.9861621540127005D-1
       V=0.5433312705027845D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD2354(X,Y,Z,W,N)
       REAL*8 X(2354)
       REAL*8 Y(2354)
       REAL*8 Z(2354)
       REAL*8 W(2354)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 2354-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.3922616270665292D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4703831750854424D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4678202801282136D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2290024646530589D-1
       V=0.1437832228979900D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5779086652271284D-1
       V=0.2303572493577644D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9863103576375984D-1
       V=0.2933110752447454D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1428155792982185D+0
       V=0.3402905998359838D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1888978116601463D+0
       V=0.3759138466870372D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2359091682970210D+0
       V=0.4030638447899798D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2831228833706171D+0
       V=0.4236591432242211D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3299495857966693D+0
       V=0.4390522656946746D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3758840802660796D+0
       V=0.4502523466626247D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4204751831009480D+0
       V=0.4580577727783541D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4633068518751051D+0
       V=0.4631391616615899D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5039849474507313D+0
       V=0.4660928953698676D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5421265793440747D+0
       V=0.4674751807936953D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6092660230557310D+0
       V=0.4676414903932920D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6374654204984869D+0
       V=0.4674086492347870D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6615136472609892D+0
       V=0.4674928539483207D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6809487285958127D+0
       V=0.4680748979686447D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6952980021665196D+0
       V=0.4690449806389040D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7041245497695400D+0
       V=0.4699877075860818D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6744033088306065D-1
       V=0.2099942281069176D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1678684485334166D+0
       V=0.3172269150712804D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2793559049539613D+0
       V=0.3832051358546523D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3935264218057639D+0
       V=0.4252193818146985D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5052629268232558D+0
       V=0.4513807963755000D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6107905315437531D+0
       V=0.4657797469114178D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1135081039843524D+0
       B=0.3331954884662588D-1
       V=0.2733362800522836D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1612866626099378D+0
       B=0.7247167465436538D-1
       V=0.3235485368463559D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2100786550168205D+0
       B=0.1151539110849745D+0
       V=0.3624908726013453D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2592282009459942D+0
       B=0.1599491097143677D+0
       V=0.3925540070712828D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3081740561320203D+0
       B=0.2058699956028027D+0
       V=0.4156129781116235D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3564289781578164D+0
       B=0.2521624953502911D+0
       V=0.4330644984623263D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4035587288240703D+0
       B=0.2982090785797674D+0
       V=0.4459677725921312D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4491671196373903D+0
       B=0.3434762087235733D+0
       V=0.4551593004456795D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4928854782917489D+0
       B=0.3874831357203437D+0
       V=0.4613341462749918D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5343646791958988D+0
       B=0.4297814821746926D+0
       V=0.4651019618269806D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5732683216530990D+0
       B=0.4699402260943537D+0
       V=0.4670249536100625D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2214131583218986D+0
       B=0.3873602040643895D-1
       V=0.3549555576441708D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2741796504750071D+0
       B=0.8089496256902013D-1
       V=0.3856108245249010D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3259797439149485D+0
       B=0.1251732177620872D+0
       V=0.4098622845756882D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3765441148826891D+0
       B=0.1706260286403185D+0
       V=0.4286328604268950D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4255773574530558D+0
       B=0.2165115147300408D+0
       V=0.4427802198993945D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4727795117058430D+0
       B=0.2622089812225259D+0
       V=0.4530473511488561D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5178546895819012D+0
       B=0.3071721431296201D+0
       V=0.4600805475703138D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5605141192097460D+0
       B=0.3508998998801138D+0
       V=0.4644599059958017D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6004763319352512D+0
       B=0.3929160876166931D+0
       V=0.4667274455712508D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3352842634946949D+0
       B=0.4202563457288019D-1
       V=0.4069360518020356D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3891971629814670D+0
       B=0.8614309758870850D-1
       V=0.4260442819919195D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4409875565542281D+0
       B=0.1314500879380001D+0
       V=0.4408678508029063D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4904893058592484D+0
       B=0.1772189657383859D+0
       V=0.4518748115548597D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5375056138769549D+0
       B=0.2228277110050294D+0
       V=0.4595564875375116D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5818255708669969D+0
       B=0.2677179935014386D+0
       V=0.4643988774315846D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6232334858144959D+0
       B=0.3113675035544165D+0
       V=0.4668827491646946D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4489485354492058D+0
       B=0.4409162378368174D-1
       V=0.4400541823741973D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5015136875933150D+0
       B=0.8939009917748489D-1
       V=0.4514512890193797D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5511300550512623D+0
       B=0.1351806029383365D+0
       V=0.4596198627347549D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5976720409858000D+0
       B=0.1808370355053196D+0
       V=0.4648659016801781D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6409956378989354D+0
       B=0.2257852192301602D+0
       V=0.4675502017157673D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5581222330827514D+0
       B=0.4532173421637160D-1
       V=0.4598494476455523D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6074705984161695D+0
       B=0.9117488031840314D-1
       V=0.4654916955152048D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6532272537379033D+0
       B=0.1369294213140155D+0
       V=0.4684709779505137D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6594761494500487D+0
       B=0.4589901487275583D-1
       V=0.4691445539106986D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD2702(X,Y,Z,W,N)
       REAL*8 X(2702)
       REAL*8 Y(2702)
       REAL*8 Z(2702)
       REAL*8 W(2702)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 2702-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.2998675149888161D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4077860529495355D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2065562538818703D-1
       V=0.1185349192520667D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5250918173022379D-1
       V=0.1913408643425751D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8993480082038376D-1
       V=0.2452886577209897D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1306023924436019D+0
       V=0.2862408183288702D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1732060388531418D+0
       V=0.3178032258257357D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2168727084820249D+0
       V=0.3422945667633690D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2609528309173586D+0
       V=0.3612790520235922D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3049252927938952D+0
       V=0.3758638229818521D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3483484138084404D+0
       V=0.3868711798859953D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3908321549106406D+0
       V=0.3949429933189938D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4320210071894814D+0
       V=0.4006068107541156D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4715824795890053D+0
       V=0.4043192149672723D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5091984794078453D+0
       V=0.4064947495808078D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5445580145650803D+0
       V=0.4075245619813152D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6072575796841768D+0
       V=0.4076423540893566D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6339484505755803D+0
       V=0.4074280862251555D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6570718257486958D+0
       V=0.4074163756012244D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6762557330090709D+0
       V=0.4077647795071246D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6911161696923790D+0
       V=0.4084517552782530D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012841911659961D+0
       V=0.4092468459224052D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7064559272410020D+0
       V=0.4097872687240906D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6123554989894765D-1
       V=0.1738986811745028D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1533070348312393D+0
       V=0.2659616045280191D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2563902605244206D+0
       V=0.3240596008171533D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3629346991663361D+0
       V=0.3621195964432943D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4683949968987538D+0
       V=0.3868838330760539D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5694479240657952D+0
       V=0.4018911532693111D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6634465430993955D+0
       V=0.4089929432983252D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1033958573552305D+0
       B=0.3034544009063584D-1
       V=0.2279907527706409D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1473521412414395D+0
       B=0.6618803044247135D-1
       V=0.2715205490578897D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1924552158705967D+0
       B=0.1054431128987715D+0
       V=0.3057917896703976D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2381094362890328D+0
       B=0.1468263551238858D+0
       V=0.3326913052452555D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2838121707936760D+0
       B=0.1894486108187886D+0
       V=0.3537334711890037D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3291323133373415D+0
       B=0.2326374238761579D+0
       V=0.3700567500783129D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3736896978741460D+0
       B=0.2758485808485768D+0
       V=0.3825245372589122D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4171406040760013D+0
       B=0.3186179331996921D+0
       V=0.3918125171518296D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4591677985256915D+0
       B=0.3605329796303794D+0
       V=0.3984720419937579D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4994733831718418D+0
       B=0.4012147253586509D+0
       V=0.4029746003338211D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5377731830445096D+0
       B=0.4403050025570692D+0
       V=0.4057428632156627D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5737917830001331D+0
       B=0.4774565904277483D+0
       V=0.4071719274114857D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2027323586271389D+0
       B=0.3544122504976147D-1
       V=0.2990236950664119D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2516942375187273D+0
       B=0.7418304388646328D-1
       V=0.3262951734212878D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3000227995257181D+0
       B=0.1150502745727186D+0
       V=0.3482634608242413D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3474806691046342D+0
       B=0.1571963371209364D+0
       V=0.3656596681700892D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3938103180359209D+0
       B=0.1999631877247100D+0
       V=0.3791740467794218D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4387519590455703D+0
       B=0.2428073457846535D+0
       V=0.3894034450156905D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4820503960077787D+0
       B=0.2852575132906155D+0
       V=0.3968600245508371D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5234573778475101D+0
       B=0.3268884208674639D+0
       V=0.4019931351420050D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5627318647235282D+0
       B=0.3673033321675939D+0
       V=0.4052108801278599D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5996390607156954D+0
       B=0.4061211551830290D+0
       V=0.4068978613940934D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3084780753791947D+0
       B=0.3860125523100059D-1
       V=0.3454275351319704D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3589988275920223D+0
       B=0.7928938987104867D-1
       V=0.3629963537007920D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4078628415881973D+0
       B=0.1212614643030087D+0
       V=0.3770187233889873D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4549287258889735D+0
       B=0.1638770827382693D+0
       V=0.3878608613694378D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5000278512957279D+0
       B=0.2065965798260176D+0
       V=0.3959065270221274D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5429785044928199D+0
       B=0.2489436378852235D+0
       V=0.4015286975463570D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5835939850491711D+0
       B=0.2904811368946891D+0
       V=0.4050866785614717D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6216870353444856D+0
       B=0.3307941957666609D+0
       V=0.4069320185051913D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4151104662709091D+0
       B=0.4064829146052554D-1
       V=0.3760120964062763D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4649804275009218D+0
       B=0.8258424547294755D-1
       V=0.3870969564418064D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5124695757009662D+0
       B=0.1251841962027289D+0
       V=0.3955287790534055D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5574711100606224D+0
       B=0.1679107505976331D+0
       V=0.4015361911302668D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5998597333287227D+0
       B=0.2102805057358715D+0
       V=0.4053836986719548D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6395007148516600D+0
       B=0.2518418087774107D+0
       V=0.4073578673299117D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5188456224746252D+0
       B=0.4194321676077518D-1
       V=0.3954628379231406D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5664190707942778D+0
       B=0.8457661551921499D-1
       V=0.4017645508847530D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6110464353283153D+0
       B=0.1273652932519396D+0
       V=0.4059030348651293D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6526430302051563D+0
       B=0.1698173239076354D+0
       V=0.4080565809484880D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6167551880377548D+0
       B=0.4266398851548864D-1
       V=0.4063018753664651D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6607195418355383D+0
       B=0.8551925814238349D-1
       V=0.4087191292799671D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD3074(X,Y,Z,W,N)
       REAL*8 X(3074)
       REAL*8 Y(3074)
       REAL*8 Z(3074)
       REAL*8 W(3074)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 3074-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.2599095953754734D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3603134089687541D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3586067974412447D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1886108518723392D-1
       V=0.9831528474385880D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4800217244625303D-1
       V=0.1605023107954450D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8244922058397242D-1
       V=0.2072200131464099D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1200408362484023D+0
       V=0.2431297618814187D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1595773530809965D+0
       V=0.2711819064496707D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2002635973434064D+0
       V=0.2932762038321116D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2415127590139982D+0
       V=0.3107032514197368D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2828584158458477D+0
       V=0.3243808058921213D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3239091015338138D+0
       V=0.3349899091374030D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3643225097962194D+0
       V=0.3430580688505218D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4037897083691802D+0
       V=0.3490124109290343D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4420247515194127D+0
       V=0.3532148948561955D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4787572538464938D+0
       V=0.3559862669062833D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5137265251275234D+0
       V=0.3576224317551411D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5466764056654611D+0
       V=0.3584050533086076D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6054859420813535D+0
       V=0.3584903581373224D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6308106701764562D+0
       V=0.3582991879040586D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6530369230179584D+0
       V=0.3582371187963125D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718609524611158D+0
       V=0.3584353631122350D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6869676499894013D+0
       V=0.3589120166517785D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6980467077240748D+0
       V=0.3595445704531601D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7048241721250522D+0
       V=0.3600943557111074D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5591105222058232D-1
       V=0.1456447096742039D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1407384078513916D+0
       V=0.2252370188283782D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2364035438976309D+0
       V=0.2766135443474897D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3360602737818170D+0
       V=0.3110729491500851D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4356292630054665D+0
       V=0.3342506712303391D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5321569415256174D+0
       V=0.3491981834026860D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6232956305040554D+0
       V=0.3576003604348932D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9469870086838469D-1
       B=0.2778748387309470D-1
       V=0.1921921305788564D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1353170300568141D+0
       B=0.6076569878628364D-1
       V=0.2301458216495632D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1771679481726077D+0
       B=0.9703072762711040D-1
       V=0.2604248549522893D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2197066664231751D+0
       B=0.1354112458524762D+0
       V=0.2845275425870697D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2624783557374927D+0
       B=0.1750996479744100D+0
       V=0.3036870897974840D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3050969521214442D+0
       B=0.2154896907449802D+0
       V=0.3188414832298066D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3472252637196021D+0
       B=0.2560954625740152D+0
       V=0.3307046414722089D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3885610219026360D+0
       B=0.2965070050624096D+0
       V=0.3398330969031360D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4288273776062765D+0
       B=0.3363641488734497D+0
       V=0.3466757899705373D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4677662471302948D+0
       B=0.3753400029836788D+0
       V=0.3516095923230054D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5051333589553359D+0
       B=0.4131297522144286D+0
       V=0.3549645184048486D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5406942145810492D+0
       B=0.4494423776081795D+0
       V=0.3570415969441392D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5742204122576457D+0
       B=0.4839938958841502D+0
       V=0.3581251798496118D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1865407027225188D+0
       B=0.3259144851070796D-1
       V=0.2543491329913348D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2321186453689432D+0
       B=0.6835679505297343D-1
       V=0.2786711051330776D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2773159142523882D+0
       B=0.1062284864451989D+0
       V=0.2985552361083679D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3219200192237254D+0
       B=0.1454404409323047D+0
       V=0.3145867929154039D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3657032593944029D+0
       B=0.1854018282582510D+0
       V=0.3273290662067609D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4084376778363622D+0
       B=0.2256297412014750D+0
       V=0.3372705511943501D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4499004945751427D+0
       B=0.2657104425000896D+0
       V=0.3448274437851510D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4898758141326335D+0
       B=0.3052755487631557D+0
       V=0.3503592783048583D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5281547442266309D+0
       B=0.3439863920645423D+0
       V=0.3541854792663162D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5645346989813992D+0
       B=0.3815229456121914D+0
       V=0.3565995517909428D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5988181252159848D+0
       B=0.4175752420966734D+0
       V=0.3578802078302898D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2850425424471603D+0
       B=0.3562149509862536D-1
       V=0.2958644592860982D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3324619433027876D+0
       B=0.7330318886871096D-1
       V=0.3119548129116835D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3785848333076282D+0
       B=0.1123226296008472D+0
       V=0.3250745225005984D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4232891028562115D+0
       B=0.1521084193337708D+0
       V=0.3355153415935208D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4664287050829722D+0
       B=0.1921844459223610D+0
       V=0.3435847568549328D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5078458493735726D+0
       B=0.2321360989678303D+0
       V=0.3495786831622488D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5473779816204180D+0
       B=0.2715886486360520D+0
       V=0.3537767805534621D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5848617133811376D+0
       B=0.3101924707571355D+0
       V=0.3564459815421428D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6201348281584888D+0
       B=0.3476121052890973D+0
       V=0.3578464061225468D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3852191185387871D+0
       B=0.3763224880035108D-1
       V=0.3239748762836212D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4325025061073423D+0
       B=0.7659581935637135D-1
       V=0.3345491784174287D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4778486229734490D+0
       B=0.1163381306083900D+0
       V=0.3429126177301782D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5211663693009000D+0
       B=0.1563890598752899D+0
       V=0.3492420343097421D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5623469504853703D+0
       B=0.1963320810149200D+0
       V=0.3537399050235257D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6012718188659246D+0
       B=0.2357847407258738D+0
       V=0.3566209152659172D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6378179206390117D+0
       B=0.2743846121244060D+0
       V=0.3581084321919782D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4836936460214534D+0
       B=0.3895902610739024D-1
       V=0.3426522117591512D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5293792562683797D+0
       B=0.7871246819312640D-1
       V=0.3491848770121379D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726281253100033D+0
       B=0.1187963808202981D+0
       V=0.3539318235231476D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6133658776169068D+0
       B=0.1587914708061787D+0
       V=0.3570231438458694D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6515085491865307D+0
       B=0.1983058575227646D+0
       V=0.3586207335051714D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5778692716064976D+0
       B=0.3977209689791542D-1
       V=0.3541196205164025D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6207904288086192D+0
       B=0.7990157592981152D-1
       V=0.3574296911573953D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6608688171046802D+0
       B=0.1199671308754309D+0
       V=0.3591993279818963D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6656263089489130D+0
       B=0.4015955957805969D-1
       V=0.3595855034661997D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD3470(X,Y,Z,W,N)
       REAL*8 X(3470)
       REAL*8 Y(3470)
       REAL*8 Z(3470)
       REAL*8 W(3470)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 3470-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.2040382730826330D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3178149703889544D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1721420832906233D-1
       V=0.8288115128076110D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4408875374981770D-1
       V=0.1360883192522954D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7594680813878681D-1
       V=0.1766854454542662D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1108335359204799D+0
       V=0.2083153161230153D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1476517054388567D+0
       V=0.2333279544657158D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1856731870860615D+0
       V=0.2532809539930247D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2243634099428821D+0
       V=0.2692472184211158D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2633006881662727D+0
       V=0.2819949946811885D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3021340904916283D+0
       V=0.2920953593973030D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3405594048030089D+0
       V=0.2999889782948352D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3783044434007372D+0
       V=0.3060292120496902D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4151194767407910D+0
       V=0.3105109167522192D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4507705766443257D+0
       V=0.3136902387550312D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4850346056573187D+0
       V=0.3157984652454632D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5176950817792470D+0
       V=0.3170516518425422D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5485384240820989D+0
       V=0.3176568425633755D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6039117238943308D+0
       V=0.3177198411207062D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6279956655573113D+0
       V=0.3175519492394733D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6493636169568952D+0
       V=0.3174654952634756D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6677644117704504D+0
       V=0.3175676415467654D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6829368572115624D+0
       V=0.3178923417835410D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6946195818184121D+0
       V=0.3183788287531909D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7025711542057026D+0
       V=0.3188755151918807D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7066004767140119D+0
       V=0.3191916889313849D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5132537689946062D-1
       V=0.1231779611744508D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1297994661331225D+0
       V=0.1924661373839880D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2188852049401307D+0
       V=0.2380881867403424D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3123174824903457D+0
       V=0.2693100663037885D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4064037620738195D+0
       V=0.2908673382834366D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4984958396944782D+0
       V=0.3053914619381535D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5864975046021365D+0
       V=0.3143916684147777D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6686711634580175D+0
       V=0.3187042244055363D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8715738780835950D-1
       B=0.2557175233367578D-1
       V=0.1635219535869790D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1248383123134007D+0
       B=0.5604823383376681D-1
       V=0.1968109917696070D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1638062693383378D+0
       B=0.8968568601900765D-1
       V=0.2236754342249974D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2035586203373176D+0
       B=0.1254086651976279D+0
       V=0.2453186687017181D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2436798975293774D+0
       B=0.1624780150162012D+0
       V=0.2627551791580541D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2838207507773806D+0
       B=0.2003422342683208D+0
       V=0.2767654860152220D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3236787502217692D+0
       B=0.2385628026255263D+0
       V=0.2879467027765895D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3629849554840691D+0
       B=0.2767731148783578D+0
       V=0.2967639918918702D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4014948081992087D+0
       B=0.3146542308245309D+0
       V=0.3035900684660351D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4389818379260225D+0
       B=0.3519196415895088D+0
       V=0.3087338237298308D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4752331143674377D+0
       B=0.3883050984023654D+0
       V=0.3124608838860167D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5100457318374018D+0
       B=0.4235613423908649D+0
       V=0.3150084294226743D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5432238388954868D+0
       B=0.4574484717196220D+0
       V=0.3165958398598402D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5745758685072442D+0
       B=0.4897311639255524D+0
       V=0.3174320440957372D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1723981437592809D+0
       B=0.3010630597881105D-1
       V=0.2182188909812599D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2149553257844597D+0
       B=0.6326031554204694D-1
       V=0.2399727933921445D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2573256081247422D+0
       B=0.9848566980258631D-1
       V=0.2579796133514652D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2993163751238106D+0
       B=0.1350835952384266D+0
       V=0.2727114052623535D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3407238005148000D+0
       B=0.1725184055442181D+0
       V=0.2846327656281355D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3813454978483264D+0
       B=0.2103559279730725D+0
       V=0.2941491102051334D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4209848104423343D+0
       B=0.2482278774554860D+0
       V=0.3016049492136107D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4594519699996300D+0
       B=0.2858099509982883D+0
       V=0.3072949726175648D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4965640166185930D+0
       B=0.3228075659915428D+0
       V=0.3114768142886460D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5321441655571562D+0
       B=0.3589459907204151D+0
       V=0.3143823673666223D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5660208438582166D+0
       B=0.3939630088864310D+0
       V=0.3162269764661535D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5980264315964364D+0
       B=0.4276029922949089D+0
       V=0.3172164663759821D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2644215852350733D+0
       B=0.3300939429072552D-1
       V=0.2554575398967435D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3090113743443063D+0
       B=0.6803887650078501D-1
       V=0.2701704069135677D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3525871079197808D+0
       B=0.1044326136206709D+0
       V=0.2823693413468940D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3950418005354029D+0
       B=0.1416751597517679D+0
       V=0.2922898463214289D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4362475663430163D+0
       B=0.1793408610504821D+0
       V=0.3001829062162428D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4760661812145854D+0
       B=0.2170630750175722D+0
       V=0.3062890864542953D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5143551042512103D+0
       B=0.2545145157815807D+0
       V=0.3108328279264746D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5509709026935597D+0
       B=0.2913940101706601D+0
       V=0.3140243146201245D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5857711030329428D+0
       B=0.3274169910910705D+0
       V=0.3160638030977130D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6186149917404392D+0
       B=0.3623081329317265D+0
       V=0.3171462882206275D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3586894569557064D+0
       B=0.3497354386450040D-1
       V=0.2812388416031796D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4035266610019441D+0
       B=0.7129736739757095D-1
       V=0.2912137500288045D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4467775312332510D+0
       B=0.1084758620193165D+0
       V=0.2993241256502206D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4883638346608543D+0
       B=0.1460915689241772D+0
       V=0.3057101738983822D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5281908348434601D+0
       B=0.1837790832369980D+0
       V=0.3105319326251432D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5661542687149311D+0
       B=0.2212075390874021D+0
       V=0.3139565514428167D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6021450102031452D+0
       B=0.2580682841160985D+0
       V=0.3161543006806366D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6360520783610050D+0
       B=0.2940656362094121D+0
       V=0.3172985960613294D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4521611065087196D+0
       B=0.3631055365867002D-1
       V=0.2989400336901431D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4959365651560963D+0
       B=0.7348318468484350D-1
       V=0.3054555883947677D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5376815804038283D+0
       B=0.1111087643812648D+0
       V=0.3104764960807702D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5773314480243768D+0
       B=0.1488226085145408D+0
       V=0.3141015825977616D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6148113245575056D+0
       B=0.1862892274135151D+0
       V=0.3164520621159896D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6500407462842380D+0
       B=0.2231909701714456D+0
       V=0.3176652305912204D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5425151448707213D+0
       B=0.3718201306118944D-1
       V=0.3105097161023939D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5841860556907931D+0
       B=0.7483616335067346D-1
       V=0.3143014117890550D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6234632186851500D+0
       B=0.1125990834266120D+0
       V=0.3168172866287200D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6602934551848843D+0
       B=0.1501303813157619D+0
       V=0.3181401865570968D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6278573968375105D+0
       B=0.3767559930245720D-1
       V=0.3170663659156037D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6665611711264577D+0
       B=0.7548443301360158D-1
       V=0.3185447944625510D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD3890(X,Y,Z,W,N)
       REAL*8 X(3890)
       REAL*8 Y(3890)
       REAL*8 Z(3890)
       REAL*8 W(3890)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 3890-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.1807395252196920D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2848008782238827D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2836065837530581D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1587876419858352D-1
       V=0.7013149266673816D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4069193593751206D-1
       V=0.1162798021956766D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7025888115257997D-1
       V=0.1518728583972105D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1027495450028704D+0
       V=0.1798796108216934D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1371457730893426D+0
       V=0.2022593385972785D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1727758532671953D+0
       V=0.2203093105575464D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2091492038929037D+0
       V=0.2349294234299855D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2458813281751915D+0
       V=0.2467682058747003D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2826545859450066D+0
       V=0.2563092683572224D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3191957291799622D+0
       V=0.2639253896763318D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3552621469299578D+0
       V=0.2699137479265108D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3906329503406230D+0
       V=0.2745196420166739D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4251028614093031D+0
       V=0.2779529197397593D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4584777520111870D+0
       V=0.2803996086684265D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4905711358710193D+0
       V=0.2820302356715842D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5212011669847385D+0
       V=0.2830056747491068D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5501878488737995D+0
       V=0.2834808950776839D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6025037877479342D+0
       V=0.2835282339078929D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6254572689549016D+0
       V=0.2833819267065800D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6460107179528248D+0
       V=0.2832858336906784D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6639541138154251D+0
       V=0.2833268235451244D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6790688515667495D+0
       V=0.2835432677029253D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6911338580371512D+0
       V=0.2839091722743049D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6999385956126490D+0
       V=0.2843308178875841D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7053037748656896D+0
       V=0.2846703550533846D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4732224387180115D-1
       V=0.1051193406971900D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1202100529326803D+0
       V=0.1657871838796974D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2034304820664855D+0
       V=0.2064648113714232D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2912285643573002D+0
       V=0.2347942745819741D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3802361792726768D+0
       V=0.2547775326597726D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4680598511056146D+0
       V=0.2686876684847025D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5528151052155599D+0
       V=0.2778665755515867D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6329386307803041D+0
       V=0.2830996616782929D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8056516651369069D-1
       B=0.2363454684003124D-1
       V=0.1403063340168372D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1156476077139389D+0
       B=0.5191291632545936D-1
       V=0.1696504125939477D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1520473382760421D+0
       B=0.8322715736994519D-1
       V=0.1935787242745390D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1892986699745931D+0
       B=0.1165855667993712D+0
       V=0.2130614510521968D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2270194446777792D+0
       B=0.1513077167409504D+0
       V=0.2289381265931048D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2648908185093273D+0
       B=0.1868882025807859D+0
       V=0.2418630292816186D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3026389259574136D+0
       B=0.2229277629776224D+0
       V=0.2523400495631193D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3400220296151384D+0
       B=0.2590951840746235D+0
       V=0.2607623973449605D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3768217953335510D+0
       B=0.2951047291750847D+0
       V=0.2674441032689209D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4128372900921884D+0
       B=0.3307019714169930D+0
       V=0.2726432360343356D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4478807131815630D+0
       B=0.3656544101087634D+0
       V=0.2765787685924545D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4817742034089257D+0
       B=0.3997448951939695D+0
       V=0.2794428690642224D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5143472814653344D+0
       B=0.4327667110812024D+0
       V=0.2814099002062895D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5454346213905650D+0
       B=0.4645196123532293D+0
       V=0.2826429531578994D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5748739313170252D+0
       B=0.4948063555703345D+0
       V=0.2832983542550884D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1599598738286342D+0
       B=0.2792357590048985D-1
       V=0.1886695565284976D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1998097412500951D+0
       B=0.5877141038139065D-1
       V=0.2081867882748234D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2396228952566202D+0
       B=0.9164573914691377D-1
       V=0.2245148680600796D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2792228341097746D+0
       B=0.1259049641962687D+0
       V=0.2380370491511872D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3184251107546741D+0
       B=0.1610594823400863D+0
       V=0.2491398041852455D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3570481164426244D+0
       B=0.1967151653460898D+0
       V=0.2581632405881230D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3949164710492144D+0
       B=0.2325404606175168D+0
       V=0.2653965506227417D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318617293970503D+0
       B=0.2682461141151439D+0
       V=0.2710857216747087D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4677221009931678D+0
       B=0.3035720116011973D+0
       V=0.2754434093903659D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5023417939270955D+0
       B=0.3382781859197439D+0
       V=0.2786579932519380D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5355701836636128D+0
       B=0.3721383065625942D+0
       V=0.2809011080679474D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5672608451328771D+0
       B=0.4049346360466055D+0
       V=0.2823336184560987D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5972704202540162D+0
       B=0.4364538098633802D+0
       V=0.2831101175806309D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2461687022333596D+0
       B=0.3070423166833368D-1
       V=0.2221679970354546D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2881774566286831D+0
       B=0.6338034669281885D-1
       V=0.2356185734270703D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3293963604116978D+0
       B=0.9742862487067941D-1
       V=0.2469228344805590D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3697303822241377D+0
       B=0.1323799532282290D+0
       V=0.2562726348642046D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4090663023135127D+0
       B=0.1678497018129336D+0
       V=0.2638756726753028D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4472819355411712D+0
       B=0.2035095105326114D+0
       V=0.2699311157390862D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4842513377231437D+0
       B=0.2390692566672091D+0
       V=0.2746233268403837D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5198477629962928D+0
       B=0.2742649818076149D+0
       V=0.2781225674454771D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5539453011883145D+0
       B=0.3088503806580094D+0
       V=0.2805881254045684D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5864196762401251D+0
       B=0.3425904245906614D+0
       V=0.2821719877004913D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6171484466668390D+0
       B=0.3752562294789468D+0
       V=0.2830222502333124D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3350337830565727D+0
       B=0.3261589934634747D-1
       V=0.2457995956744870D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3775773224758284D+0
       B=0.6658438928081572D-1
       V=0.2551474407503706D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4188155229848973D+0
       B=0.1014565797157954D+0
       V=0.2629065335195311D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4586805892009344D+0
       B=0.1368573320843822D+0
       V=0.2691900449925075D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4970895714224235D+0
       B=0.1724614851951608D+0
       V=0.2741275485754276D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5339505133960747D+0
       B=0.2079779381416412D+0
       V=0.2778530970122595D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5691665792531440D+0
       B=0.2431385788322288D+0
       V=0.2805010567646741D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6026387682680377D+0
       B=0.2776901883049853D+0
       V=0.2822055834031040D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6342676150163307D+0
       B=0.3113881356386632D+0
       V=0.2831016901243473D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4237951119537067D+0
       B=0.3394877848664351D-1
       V=0.2624474901131803D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4656918683234929D+0
       B=0.6880219556291447D-1
       V=0.2688034163039377D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5058857069185980D+0
       B=0.1041946859721635D+0
       V=0.2738932751287636D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5443204666713996D+0
       B=0.1398039738736393D+0
       V=0.2777944791242523D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5809298813759742D+0
       B=0.1753373381196155D+0
       V=0.2806011661660987D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6156416039447128D+0
       B=0.2105215793514010D+0
       V=0.2824181456597460D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6483801351066604D+0
       B=0.2450953312157051D+0
       V=0.2833585216577828D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5103616577251688D+0
       B=0.3485560643800719D-1
       V=0.2738165236962878D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5506738792580681D+0
       B=0.7026308631512033D-1
       V=0.2778365208203180D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5889573040995292D+0
       B=0.1059035061296403D+0
       V=0.2807852940418966D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6251641589516930D+0
       B=0.1414823925236026D+0
       V=0.2827245949674705D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6592414921570178D+0
       B=0.1767207908214530D+0
       V=0.2837342344829828D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5930314017533384D+0
       B=0.3542189339561672D-1
       V=0.2809233907610981D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6309812253390175D+0
       B=0.7109574040369549D-1
       V=0.2829930809742694D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6666296011353230D+0
       B=0.1067259792282730D+0
       V=0.2841097874111479D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6703715271049922D+0
       B=0.3569455268820809D-1
       V=0.2843455206008783D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD4334(X,Y,Z,W,N)
       REAL*8 X(4334)
       REAL*8 Y(4334)
       REAL*8 Z(4334)
       REAL*8 W(4334)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 4334-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.1449063022537883D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2546377329828424D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1462896151831013D-1
       V=0.6018432961087496D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3769840812493139D-1
       V=0.1002286583263673D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6524701904096891D-1
       V=0.1315222931028093D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9560543416134648D-1
       V=0.1564213746876724D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1278335898929198D+0
       V=0.1765118841507736D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1613096104466031D+0
       V=0.1928737099311080D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1955806225745371D+0
       V=0.2062658534263270D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2302935218498028D+0
       V=0.2172395445953787D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2651584344113027D+0
       V=0.2262076188876047D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2999276825183209D+0
       V=0.2334885699462397D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3343828669718798D+0
       V=0.2393355273179203D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3683265013750518D+0
       V=0.2439559200468863D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4015763206518108D+0
       V=0.2475251866060002D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4339612026399770D+0
       V=0.2501965558158773D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4653180651114582D+0
       V=0.2521081407925925D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4954893331080803D+0
       V=0.2533881002388081D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5243207068924930D+0
       V=0.2541582900848261D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5516590479041704D+0
       V=0.2545365737525860D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6012371927804176D+0
       V=0.2545726993066799D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6231574466449819D+0
       V=0.2544456197465555D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6429416514181271D+0
       V=0.2543481596881064D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6604124272943595D+0
       V=0.2543506451429194D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6753851470408250D+0
       V=0.2544905675493763D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6876717970626160D+0
       V=0.2547611407344429D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6970895061319234D+0
       V=0.2551060375448869D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7034746912553310D+0
       V=0.2554291933816039D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7067017217542295D+0
       V=0.2556255710686343D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4382223501131123D-1
       V=0.9041339695118195D-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1117474077400006D+0
       V=0.1438426330079022D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1897153252911440D+0
       V=0.1802523089820518D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2724023009910331D+0
       V=0.2060052290565496D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3567163308709902D+0
       V=0.2245002248967466D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4404784483028087D+0
       V=0.2377059847731150D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5219833154161411D+0
       V=0.2468118955882525D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5998179868977553D+0
       V=0.2525410872966528D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6727803154548222D+0
       V=0.2553101409933397D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7476563943166086D-1
       B=0.2193168509461185D-1
       V=0.1212879733668632D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1075341482001416D+0
       B=0.4826419281533887D-1
       V=0.1472872881270931D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1416344885203259D+0
       B=0.7751191883575742D-1
       V=0.1686846601010828D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1766325315388586D+0
       B=0.1087558139247680D+0
       V=0.1862698414660208D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2121744174481514D+0
       B=0.1413661374253096D+0
       V=0.2007430956991861D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2479669443408145D+0
       B=0.1748768214258880D+0
       V=0.2126568125394796D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2837600452294113D+0
       B=0.2089216406612073D+0
       V=0.2224394603372113D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3193344933193984D+0
       B=0.2431987685545972D+0
       V=0.2304264522673135D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3544935442438745D+0
       B=0.2774497054377770D+0
       V=0.2368854288424087D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3890571932288154D+0
       B=0.3114460356156915D+0
       V=0.2420352089461772D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4228581214259090D+0
       B=0.3449806851913012D+0
       V=0.2460597113081295D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4557387211304052D+0
       B=0.3778618641248256D+0
       V=0.2491181912257687D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4875487950541643D+0
       B=0.4099086391698978D+0
       V=0.2513528194205857D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5181436529962997D+0
       B=0.4409474925853973D+0
       V=0.2528943096693220D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5473824095600661D+0
       B=0.4708094517711291D+0
       V=0.2538660368488136D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5751263398976174D+0
       B=0.4993275140354637D+0
       V=0.2543868648299022D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1489515746840028D+0
       B=0.2599381993267017D-1
       V=0.1642595537825183D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1863656444351767D+0
       B=0.5479286532462190D-1
       V=0.1818246659849308D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2238602880356348D+0
       B=0.8556763251425254D-1
       V=0.1966565649492420D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2612723375728160D+0
       B=0.1177257802267011D+0
       V=0.2090677905657991D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2984332990206190D+0
       B=0.1508168456192700D+0
       V=0.2193820409510504D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3351786584663333D+0
       B=0.1844801892177727D+0
       V=0.2278870827661928D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3713505522209120D+0
       B=0.2184145236087598D+0
       V=0.2348283192282090D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4067981098954663D+0
       B=0.2523590641486229D+0
       V=0.2404139755581477D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4413769993687534D+0
       B=0.2860812976901373D+0
       V=0.2448227407760734D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4749487182516394D+0
       B=0.3193686757808996D+0
       V=0.2482110455592573D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5073798105075426D+0
       B=0.3520226949547602D+0
       V=0.2507192397774103D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5385410448878654D+0
       B=0.3838544395667890D+0
       V=0.2524765968534880D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5683065353670530D+0
       B=0.4146810037640963D+0
       V=0.2536052388539425D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5965527620663510D+0
       B=0.4443224094681121D+0
       V=0.2542230588033068D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2299227700856157D+0
       B=0.2865757664057584D-1
       V=0.1944817013047896D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2695752998553267D+0
       B=0.5923421684485993D-1
       V=0.2067862362746635D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3086178716611389D+0
       B=0.9117817776057715D-1
       V=0.2172440734649114D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3469649871659077D+0
       B=0.1240593814082605D+0
       V=0.2260125991723423D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3845153566319655D+0
       B=0.1575272058259175D+0
       V=0.2332655008689523D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4211600033403215D+0
       B=0.1912845163525413D+0
       V=0.2391699681532458D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4567867834329882D+0
       B=0.2250710177858171D+0
       V=0.2438801528273928D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4912829319232061D+0
       B=0.2586521303440910D+0
       V=0.2475370504260665D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5245364793303812D+0
       B=0.2918112242865407D+0
       V=0.2502707235640574D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5564369788915756D+0
       B=0.3243439239067890D+0
       V=0.2522031701054241D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5868757697775287D+0
       B=0.3560536787835351D+0
       V=0.2534511269978784D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6157458853519617D+0
       B=0.3867480821242581D+0
       V=0.2541284914955151D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3138461110672113D+0
       B=0.3051374637507278D-1
       V=0.2161509250688394D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3542495872050569D+0
       B=0.6237111233730755D-1
       V=0.2248778513437852D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3935751553120181D+0
       B=0.9516223952401907D-1
       V=0.2322388803404617D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4317634668111147D+0
       B=0.1285467341508517D+0
       V=0.2383265471001355D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4687413842250821D+0
       B=0.1622318931656033D+0
       V=0.2432476675019525D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5044274237060283D+0
       B=0.1959581153836453D+0
       V=0.2471122223750674D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5387354077925727D+0
       B=0.2294888081183837D+0
       V=0.2500291752486870D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5715768898356105D+0
       B=0.2626031152713945D+0
       V=0.2521055942764682D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6028627200136111D+0
       B=0.2950904075286713D+0
       V=0.2534472785575503D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6325039812653463D+0
       B=0.3267458451113286D+0
       V=0.2541599713080121D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3981986708423407D+0
       B=0.3183291458749821D-1
       V=0.2317380975862936D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4382791182133300D+0
       B=0.6459548193880908D-1
       V=0.2378550733719775D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4769233057218166D+0
       B=0.9795757037087952D-1
       V=0.2428884456739118D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5140823911194238D+0
       B=0.1316307235126655D+0
       V=0.2469002655757292D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5496977833862983D+0
       B=0.1653556486358704D+0
       V=0.2499657574265851D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5837047306512727D+0
       B=0.1988931724126510D+0
       V=0.2521676168486082D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6160349566926879D+0
       B=0.2320174581438950D+0
       V=0.2535935662645334D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6466185353209440D+0
       B=0.2645106562168662D+0
       V=0.2543356743363214D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4810835158795404D+0
       B=0.3275917807743992D-1
       V=0.2427353285201535D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5199925041324341D+0
       B=0.6612546183967181D-1
       V=0.2468258039744386D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5571717692207494D+0
       B=0.9981498331474143D-1
       V=0.2500060956440310D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5925789250836378D+0
       B=0.1335687001410374D+0
       V=0.2523238365420979D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6261658523859670D+0
       B=0.1671444402896463D+0
       V=0.2538399260252846D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6578811126669331D+0
       B=0.2003106382156076D+0
       V=0.2546255927268069D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5609624612998100D+0
       B=0.3337500940231335D-1
       V=0.2500583360048449D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5979959659984670D+0
       B=0.6708750335901803D-1
       V=0.2524777638260203D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6330523711054002D+0
       B=0.1008792126424850D+0
       V=0.2540951193860656D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6660960998103972D+0
       B=0.1345050343171794D+0
       V=0.2549524085027472D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6365384364585819D+0
       B=0.3372799460737052D-1
       V=0.2542569507009158D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6710994302899275D+0
       B=0.6755249309678028D-1
       V=0.2552114127580376D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD4802(X,Y,Z,W,N)
       REAL*8 X(4802)
       REAL*8 Y(4802)
       REAL*8 Z(4802)
       REAL*8 W(4802)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 4802-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.9687521879420705D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2307897895367918D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2297310852498558D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2335728608887064D-1
       V=0.7386265944001919D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4352987836550653D-1
       V=0.8257977698542210D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6439200521088801D-1
       V=0.9706044762057630D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9003943631993181D-1
       V=0.1302393847117003D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1196706615548473D+0
       V=0.1541957004600968D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1511715412838134D+0
       V=0.1704459770092199D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1835982828503801D+0
       V=0.1827374890942906D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2165081259155405D+0
       V=0.1926360817436107D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2496208720417563D+0
       V=0.2008010239494833D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2827200673567900D+0
       V=0.2075635983209175D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3156190823994346D+0
       V=0.2131306638690909D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3481476793749115D+0
       V=0.2176562329937335D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3801466086947226D+0
       V=0.2212682262991018D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4114652119634011D+0
       V=0.2240799515668565D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4419598786519751D+0
       V=0.2261959816187525D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4714925949329543D+0
       V=0.2277156368808855D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4999293972879466D+0
       V=0.2287351772128336D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5271387221431248D+0
       V=0.2293490814084085D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5529896780837761D+0
       V=0.2296505312376273D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6000856099481712D+0
       V=0.2296793832318756D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6210562192785175D+0
       V=0.2295785443842974D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6401165879934240D+0
       V=0.2295017931529102D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6571144029244334D+0
       V=0.2295059638184868D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718910821718863D+0
       V=0.2296232343237362D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6842845591099010D+0
       V=0.2298530178740771D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6941353476269816D+0
       V=0.2301579790280501D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012965242212991D+0
       V=0.2304690404996513D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7056471428242644D+0
       V=0.2307027995907102D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4595557643585895D-1
       V=0.9312274696671092D-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1049316742435023D+0
       V=0.1199919385876926D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1773548879549274D+0
       V=0.1598039138877690D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2559071411236127D+0
       V=0.1822253763574900D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3358156837985898D+0
       V=0.1988579593655040D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4155835743763893D+0
       V=0.2112620102533307D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4937894296167472D+0
       V=0.2201594887699007D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5691569694793316D+0
       V=0.2261622590895036D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6405840854894251D+0
       V=0.2296458453435705D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7345133894143348D-1
       B=0.2177844081486067D-1
       V=0.1006006990267000D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1009859834044931D+0
       B=0.4590362185775188D-1
       V=0.1227676689635876D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1324289619748758D+0
       B=0.7255063095690877D-1
       V=0.1467864280270117D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1654272109607127D+0
       B=0.1017825451960684D+0
       V=0.1644178912101232D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1990767186776461D+0
       B=0.1325652320980364D+0
       V=0.1777664890718961D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2330125945523278D+0
       B=0.1642765374496765D+0
       V=0.1884825664516690D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2670080611108287D+0
       B=0.1965360374337889D+0
       V=0.1973269246453848D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3008753376294316D+0
       B=0.2290726770542238D+0
       V=0.2046767775855328D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3344475596167860D+0
       B=0.2616645495370823D+0
       V=0.2107600125918040D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3675709724070786D+0
       B=0.2941150728843141D+0
       V=0.2157416362266829D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4001000887587812D+0
       B=0.3262440400919066D+0
       V=0.2197557816920721D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318956350436028D+0
       B=0.3578835350611916D+0
       V=0.2229192611835437D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4628239056795531D+0
       B=0.3888751854043678D+0
       V=0.2253385110212775D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4927563229773636D+0
       B=0.4190678003222840D+0
       V=0.2271137107548774D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5215687136707969D+0
       B=0.4483151836883852D+0
       V=0.2283414092917525D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5491402346984905D+0
       B=0.4764740676087880D+0
       V=0.2291161673130077D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5753520160126075D+0
       B=0.5034021310998277D+0
       V=0.2295313908576598D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1388326356417754D+0
       B=0.2435436510372806D-1
       V=0.1438204721359031D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1743686900537244D+0
       B=0.5118897057342652D-1
       V=0.1607738025495257D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2099737037950268D+0
       B=0.8014695048539634D-1
       V=0.1741483853528379D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2454492590908548D+0
       B=0.1105117874155699D+0
       V=0.1851918467519151D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2807219257864278D+0
       B=0.1417950531570966D+0
       V=0.1944628638070613D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3156842271975842D+0
       B=0.1736604945719597D+0
       V=0.2022495446275152D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3502090945177752D+0
       B=0.2058466324693981D+0
       V=0.2087462382438514D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3841684849519686D+0
       B=0.2381284261195919D+0
       V=0.2141074754818308D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4174372367906016D+0
       B=0.2703031270422569D+0
       V=0.2184640913748162D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4498926465011892D+0
       B=0.3021845683091309D+0
       V=0.2219309165220329D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4814146229807701D+0
       B=0.3335993355165720D+0
       V=0.2246123118340624D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5118863625734701D+0
       B=0.3643833735518232D+0
       V=0.2266062766915125D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5411947455119144D+0
       B=0.3943789541958179D+0
       V=0.2280072952230796D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5692301500357246D+0
       B=0.4234320144403542D+0
       V=0.2289082025202583D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5958857204139576D+0
       B=0.4513897947419260D+0
       V=0.2294012695120025D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2156270284785766D+0
       B=0.2681225755444491D-1
       V=0.1722434488736947D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2532385054909710D+0
       B=0.5557495747805614D-1
       V=0.1830237421455091D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2902564617771537D+0
       B=0.8569368062950249D-1
       V=0.1923855349997633D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3266979823143256D+0
       B=0.1167367450324135D+0
       V=0.2004067861936271D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3625039627493614D+0
       B=0.1483861994003304D+0
       V=0.2071817297354263D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3975838937548699D+0
       B=0.1803821503011405D+0
       V=0.2128250834102103D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318396099009774D+0
       B=0.2124962965666424D+0
       V=0.2174513719440102D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4651706555732742D+0
       B=0.2445221837805913D+0
       V=0.2211661839150214D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4974752649620969D+0
       B=0.2762701224322987D+0
       V=0.2240665257813102D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5286517579627517D+0
       B=0.3075627775211328D+0
       V=0.2262439516632620D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5586001195731895D+0
       B=0.3382311089826877D+0
       V=0.2277874557231869D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5872229902021319D+0
       B=0.3681108834741399D+0
       V=0.2287854314454994D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6144258616235123D+0
       B=0.3970397446872839D+0
       V=0.2293268499615575D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2951676508064861D+0
       B=0.2867499538750441D-1
       V=0.1912628201529828D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3335085485472725D+0
       B=0.5867879341903510D-1
       V=0.1992499672238701D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3709561760636381D+0
       B=0.8961099205022284D-1
       V=0.2061275533454027D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4074722861667498D+0
       B=0.1211627927626297D+0
       V=0.2119318215968572D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4429923648839117D+0
       B=0.1530748903554898D+0
       V=0.2167416581882652D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4774428052721736D+0
       B=0.1851176436721877D+0
       V=0.2206430730516600D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5107446539535904D+0
       B=0.2170829107658179D+0
       V=0.2237186938699523D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5428151370542935D+0
       B=0.2487786689026271D+0
       V=0.2260480075032884D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5735699292556964D+0
       B=0.2800239952795016D+0
       V=0.2277098884558542D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6029253794562866D+0
       B=0.3106445702878119D+0
       V=0.2287845715109671D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6307998987073145D+0
       B=0.3404689500841194D+0
       V=0.2293547268236294D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3752652273692719D+0
       B=0.2997145098184479D-1
       V=0.2056073839852528D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4135383879344028D+0
       B=0.6086725898678011D-1
       V=0.2114235865831876D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4506113885153907D+0
       B=0.9238849548435643D-1
       V=0.2163175629770551D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4864401554606072D+0
       B=0.1242786603851851D+0
       V=0.2203392158111650D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5209708076611709D+0
       B=0.1563086731483386D+0
       V=0.2235473176847839D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5541422135830122D+0
       B=0.1882696509388506D+0
       V=0.2260024141501235D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5858880915113817D+0
       B=0.2199672979126059D+0
       V=0.2277675929329182D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6161399390603444D+0
       B=0.2512165482924867D+0
       V=0.2289102112284834D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6448296482255090D+0
       B=0.2818368701871888D+0
       V=0.2295027954625118D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4544796274917948D+0
       B=0.3088970405060312D-1
       V=0.2161281589879992D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4919389072146628D+0
       B=0.6240947677636835D-1
       V=0.2201980477395102D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5279313026985183D+0
       B=0.9430706144280313D-1
       V=0.2234952066593166D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5624169925571135D+0
       B=0.1263547818770374D+0
       V=0.2260540098520838D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5953484627093287D+0
       B=0.1583430788822594D+0
       V=0.2279157981899988D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6266730715339185D+0
       B=0.1900748462555988D+0
       V=0.2291296918565571D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6563363204278871D+0
       B=0.2213599519592567D+0
       V=0.2297533752536649D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5314574716585696D+0
       B=0.3152508811515374D-1
       V=0.2234927356465995D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5674614932298185D+0
       B=0.6343865291465561D-1
       V=0.2261288012985219D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6017706004970264D+0
       B=0.9551503504223951D-1
       V=0.2280818160923688D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6343471270264178D+0
       B=0.1275440099801196D+0
       V=0.2293773295180159D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6651494599127802D+0
       B=0.1593252037671960D+0
       V=0.2300528767338634D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6050184986005704D+0
       B=0.3192538338496105D-1
       V=0.2281893855065666D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6390163550880400D+0
       B=0.6402824353962306D-1
       V=0.2295720444840727D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6711199107088448D+0
       B=0.9609805077002909D-1
       V=0.2303227649026753D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6741354429572275D+0
       B=0.3211853196273233D-1
       V=0.2304831913227114D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD5294(X,Y,Z,W,N)
       REAL*8 X(5294)
       REAL*8 Y(5294)
       REAL*8 Z(5294)
       REAL*8 W(5294)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 5294-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.9080510764308163D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2084824361987793D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2303261686261450D-1
       V=0.5011105657239616D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3757208620162394D-1
       V=0.5942520409683854D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5821912033821852D-1
       V=0.9564394826109721D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8403127529194872D-1
       V=0.1185530657126338D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1122927798060578D+0
       V=0.1364510114230331D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1420125319192987D+0
       V=0.1505828825605415D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1726396437341978D+0
       V=0.1619298749867023D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2038170058115696D+0
       V=0.1712450504267789D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2352849892876508D+0
       V=0.1789891098164999D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2668363354312461D+0
       V=0.1854474955629795D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2982941279900452D+0
       V=0.1908148636673661D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3295002922087076D+0
       V=0.1952377405281833D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3603094918363593D+0
       V=0.1988349254282232D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3905857895173920D+0
       V=0.2017079807160050D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4202005758160837D+0
       V=0.2039473082709094D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4490310061597227D+0
       V=0.2056360279288953D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4769586160311491D+0
       V=0.2068525823066865D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5038679887049750D+0
       V=0.2076724877534488D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5296454286519961D+0
       V=0.2081694278237885D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5541776207164850D+0
       V=0.2084157631219326D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5990467321921213D+0
       V=0.2084381531128593D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6191467096294587D+0
       V=0.2083476277129307D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6375251212901849D+0
       V=0.2082686194459732D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6540514381131168D+0
       V=0.2082475686112415D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6685899064391510D+0
       V=0.2083139860289915D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6810013009681648D+0
       V=0.2084745561831237D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6911469578730340D+0
       V=0.2087091313375890D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6988956915141736D+0
       V=0.2089718413297697D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7041335794868720D+0
       V=0.2092003303479793D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7067754398018567D+0
       V=0.2093336148263241D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3840368707853623D-1
       V=0.7591708117365267D-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9835485954117399D-1
       V=0.1083383968169186D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1665774947612998D+0
       V=0.1403019395292510D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2405702335362910D+0
       V=0.1615970179286436D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3165270770189046D+0
       V=0.1771144187504911D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3927386145645443D+0
       V=0.1887760022988168D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4678825918374656D+0
       V=0.1973474670768214D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5408022024266935D+0
       V=0.2033787661234659D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6104967445752438D+0
       V=0.2072343626517331D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6760910702685738D+0
       V=0.2091177834226918D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6655644120217392D-1
       B=0.1936508874588424D-1
       V=0.9316684484675566D-4
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9446246161270182D-1
       B=0.4252442002115869D-1
       V=0.1116193688682976D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1242651925452509D+0
       B=0.6806529315354374D-1
       V=0.1298623551559414D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1553438064846751D+0
       B=0.9560957491205369D-1
       V=0.1450236832456426D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1871137110542670D+0
       B=0.1245931657452888D+0
       V=0.1572719958149914D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2192612628836257D+0
       B=0.1545385828778978D+0
       V=0.1673234785867195D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2515682807206955D+0
       B=0.1851004249723368D+0
       V=0.1756860118725188D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2838535866287290D+0
       B=0.2160182608272384D+0
       V=0.1826776290439367D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3159578817528521D+0
       B=0.2470799012277111D+0
       V=0.1885116347992865D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3477370882791392D+0
       B=0.2781014208986402D+0
       V=0.1933457860170574D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3790576960890540D+0
       B=0.3089172523515731D+0
       V=0.1973060671902064D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4097938317810200D+0
       B=0.3393750055472244D+0
       V=0.2004987099616311D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4398256572859637D+0
       B=0.3693322470987730D+0
       V=0.2030170909281499D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4690384114718480D+0
       B=0.3986541005609877D+0
       V=0.2049461460119080D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4973216048301053D+0
       B=0.4272112491408562D+0
       V=0.2063653565200186D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5245681526132446D+0
       B=0.4548781735309936D+0
       V=0.2073507927381027D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5506733911803888D+0
       B=0.4815315355023251D+0
       V=0.2079764593256122D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5755339829522475D+0
       B=0.5070486445801855D+0
       V=0.2083150534968778D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1305472386056362D+0
       B=0.2284970375722366D-1
       V=0.1262715121590664D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1637327908216477D+0
       B=0.4812254338288384D-1
       V=0.1414386128545972D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1972734634149637D+0
       B=0.7531734457511935D-1
       V=0.1538740401313898D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2308694653110130D+0
       B=0.1039043639882017D+0
       V=0.1642434942331432D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2643899218338160D+0
       B=0.1334526587117626D+0
       V=0.1729790609237496D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2977171599622171D+0
       B=0.1636414868936382D+0
       V=0.1803505190260828D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3307293903032310D+0
       B=0.1942195406166568D+0
       V=0.1865475350079657D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3633069198219073D+0
       B=0.2249752879943753D+0
       V=0.1917182669679069D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3953346955922727D+0
       B=0.2557218821820032D+0
       V=0.1959851709034382D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4267018394184914D+0
       B=0.2862897925213193D+0
       V=0.1994529548117882D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4573009622571704D+0
       B=0.3165224536636518D+0
       V=0.2022138911146548D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4870279559856109D+0
       B=0.3462730221636496D+0
       V=0.2043518024208592D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5157819581450322D+0
       B=0.3754016870282835D+0
       V=0.2059450313018110D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5434651666465393D+0
       B=0.4037733784993613D+0
       V=0.2070685715318472D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5699823887764627D+0
       B=0.4312557784139123D+0
       V=0.2077955310694373D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5952403350947741D+0
       B=0.4577175367122110D+0
       V=0.2081980387824712D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2025152599210369D+0
       B=0.2520253617719557D-1
       V=0.1521318610377956D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2381066653274425D+0
       B=0.5223254506119000D-1
       V=0.1622772720185755D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2732823383651612D+0
       B=0.8060669688588620D-1
       V=0.1710498139420709D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3080137692611118D+0
       B=0.1099335754081255D+0
       V=0.1785911149448736D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3422405614587601D+0
       B=0.1399120955959857D+0
       V=0.1850125313687736D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3758808773890420D+0
       B=0.1702977801651705D+0
       V=0.1904229703933298D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4088458383438932D+0
       B=0.2008799256601680D+0
       V=0.1949259956121987D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4410450550841152D+0
       B=0.2314703052180836D+0
       V=0.1986161545363960D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4723879420561312D+0
       B=0.2618972111375892D+0
       V=0.2015790585641370D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5027843561874343D+0
       B=0.2920013195600270D+0
       V=0.2038934198707418D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5321453674452458D+0
       B=0.3216322555190551D+0
       V=0.2056334060538251D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5603839113834030D+0
       B=0.3506456615934198D+0
       V=0.2068705959462289D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5874150706875146D+0
       B=0.3789007181306267D+0
       V=0.2076753906106002D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6131559381660038D+0
       B=0.4062580170572782D+0
       V=0.2081179391734803D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2778497016394506D+0
       B=0.2696271276876226D-1
       V=0.1700345216228943D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3143733562261912D+0
       B=0.5523469316960465D-1
       V=0.1774906779990410D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3501485810261827D+0
       B=0.8445193201626464D-1
       V=0.1839659377002642D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3851430322303653D+0
       B=0.1143263119336083D+0
       V=0.1894987462975169D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4193013979470415D+0
       B=0.1446177898344475D+0
       V=0.1941548809452595D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4525585960458567D+0
       B=0.1751165438438091D+0
       V=0.1980078427252384D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4848447779622947D+0
       B=0.2056338306745660D+0
       V=0.2011296284744488D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5160871208276894D+0
       B=0.2359965487229226D+0
       V=0.2035888456966776D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5462112185696926D+0
       B=0.2660430223139146D+0
       V=0.2054516325352142D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5751425068101757D+0
       B=0.2956193664498032D+0
       V=0.2067831033092635D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6028073872853596D+0
       B=0.3245763905312779D+0
       V=0.2076485320284876D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6291338275278409D+0
       B=0.3527670026206972D+0
       V=0.2081141439525255D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3541797528439391D+0
       B=0.2823853479435550D-1
       V=0.1834383015469222D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3908234972074657D+0
       B=0.5741296374713106D-1
       V=0.1889540591777677D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4264408450107590D+0
       B=0.8724646633650199D-1
       V=0.1936677023597375D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4609949666553286D+0
       B=0.1175034422915616D+0
       V=0.1976176495066504D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4944389496536006D+0
       B=0.1479755652628428D+0
       V=0.2008536004560983D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5267194884346086D+0
       B=0.1784740659484352D+0
       V=0.2034280351712291D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5577787810220990D+0
       B=0.2088245700431244D+0
       V=0.2053944466027758D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5875563763536670D+0
       B=0.2388628136570763D+0
       V=0.2068077642882360D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6159910016391269D+0
       B=0.2684308928769185D+0
       V=0.2077250949661599D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6430219602956268D+0
       B=0.2973740761960252D+0
       V=0.2082062440705320D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4300647036213646D+0
       B=0.2916399920493977D-1
       V=0.1934374486546626D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4661486308935531D+0
       B=0.5898803024755659D-1
       V=0.1974107010484300D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5009658555287261D+0
       B=0.8924162698525409D-1
       V=0.2007129290388658D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5344824270447704D+0
       B=0.1197185199637321D+0
       V=0.2033736947471293D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5666575997416371D+0
       B=0.1502300756161382D+0
       V=0.2054287125902493D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5974457471404752D+0
       B=0.1806004191913564D+0
       V=0.2069184936818894D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6267984444116886D+0
       B=0.2106621764786252D+0
       V=0.2078883689808782D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6546664713575417D+0
       B=0.2402526932671914D+0
       V=0.2083886366116359D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5042711004437253D+0
       B=0.2982529203607657D-1
       V=0.2006593275470817D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5392127456774380D+0
       B=0.6008728062339922D-1
       V=0.2033728426135397D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726819437668618D+0
       B=0.9058227674571398D-1
       V=0.2055008781377608D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6046469254207278D+0
       B=0.1211219235803400D+0
       V=0.2070651783518502D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6350716157434952D+0
       B=0.1515286404791580D+0
       V=0.2080953335094320D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6639177679185454D+0
       B=0.1816314681255552D+0
       V=0.2086284998988521D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5757276040972253D+0
       B=0.3026991752575440D-1
       V=0.2055549387644668D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6090265823139755D+0
       B=0.6078402297870770D-1
       V=0.2071871850267654D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6406735344387661D+0
       B=0.9135459984176636D-1
       V=0.2082856600431965D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6706397927793709D+0
       B=0.1218024155966590D+0
       V=0.2088705858819358D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6435019674426665D+0
       B=0.3052608357660639D-1
       V=0.2083995867536322D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6747218676375681D+0
       B=0.6112185773983089D-1
       V=0.2090509712889637D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine

  SUBROUTINE LD5810(X,Y,Z,W,N)
       REAL*8 X(5810)
       REAL*8 Y(5810)
       REAL*8 Z(5810)
       REAL*8 W(5810)
       INTEGER N
       REAL*8 A,B,V
!VW
!VW    LEBEDEV 5810-POINT ANGULAR GRID
!VW
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated using a C to fortran77 conversion
!hvd   tool written by Dr. Christoph van Wuellen.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
       N=1
       V=0.9735347946175486D-5
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1907581241803167D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1901059546737578D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1182361662400277D-1
       V=0.3926424538919212D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3062145009138958D-1
       V=0.6667905467294382D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5329794036834243D-1
       V=0.8868891315019135D-4
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7848165532862220D-1
       V=0.1066306000958872D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1054038157636201D+0
       V=0.1214506743336128D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1335577797766211D+0
       V=0.1338054681640871D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1625769955502252D+0
       V=0.1441677023628504D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1921787193412792D+0
       V=0.1528880200826557D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2221340534690548D+0
       V=0.1602330623773609D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2522504912791132D+0
       V=0.1664102653445244D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2823610860679697D+0
       V=0.1715845854011323D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3123173966267560D+0
       V=0.1758901000133069D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3419847036953789D+0
       V=0.1794382485256736D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3712386456999758D+0
       V=0.1823238106757407D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3999627649876828D+0
       V=0.1846293252959976D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4280466458648093D+0
       V=0.1864284079323098D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4553844360185711D+0
       V=0.1877882694626914D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4818736094437834D+0
       V=0.1887716321852025D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5074138709260629D+0
       V=0.1894381638175673D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5319061304570707D+0
       V=0.1898454899533629D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5552514978677286D+0
       V=0.1900497929577815D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5981009025246183D+0
       V=0.1900671501924092D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6173990192228116D+0
       V=0.1899837555533510D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6351365239411131D+0
       V=0.1899014113156229D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6512010228227200D+0
       V=0.1898581257705106D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6654758363948120D+0
       V=0.1898804756095753D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6778410414853370D+0
       V=0.1899793610426402D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6881760887484110D+0
       V=0.1901464554844117D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6963645267094598D+0
       V=0.1903533246259542D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7023010617153579D+0
       V=0.1905556158463228D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7059004636628753D+0
       V=0.1907037155663528D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3552470312472575D-1
       V=0.5992997844249967D-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9151176620841283D-1
       V=0.9749059382456978D-4
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1566197930068980D+0
       V=0.1241680804599158D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2265467599271907D+0
       V=0.1437626154299360D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2988242318581361D+0
       V=0.1584200054793902D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3717482419703886D+0
       V=0.1694436550982744D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4440094491758889D+0
       V=0.1776617014018108D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5145337096756642D+0
       V=0.1836132434440077D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5824053672860230D+0
       V=0.1876494727075983D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6468283961043370D+0
       V=0.1899906535336482D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6095964259104373D-1
       B=0.1787828275342931D-1
       V=0.8143252820767350D-4
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8811962270959388D-1
       B=0.3953888740792096D-1
       V=0.9998859890887728D-4
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1165936722428831D+0
       B=0.6378121797722990D-1
       V=0.1156199403068359D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1460232857031785D+0
       B=0.8985890813745037D-1
       V=0.1287632092635513D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1761197110181755D+0
       B=0.1172606510576162D+0
       V=0.1398378643365139D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2066471190463718D+0
       B=0.1456102876970995D+0
       V=0.1491876468417391D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2374076026328152D+0
       B=0.1746153823011775D+0
       V=0.1570855679175456D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2682305474337051D+0
       B=0.2040383070295584D+0
       V=0.1637483948103775D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2989653312142369D+0
       B=0.2336788634003698D+0
       V=0.1693500566632843D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3294762752772209D+0
       B=0.2633632752654219D+0
       V=0.1740322769393633D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3596390887276086D+0
       B=0.2929369098051601D+0
       V=0.1779126637278296D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3893383046398812D+0
       B=0.3222592785275512D+0
       V=0.1810908108835412D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4184653789358347D+0
       B=0.3512004791195743D+0
       V=0.1836529132600190D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4469172319076166D+0
       B=0.3796385677684537D+0
       V=0.1856752841777379D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4745950813276976D+0
       B=0.4074575378263879D+0
       V=0.1872270566606832D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5014034601410262D+0
       B=0.4345456906027828D+0
       V=0.1883722645591307D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5272493404551239D+0
       B=0.4607942515205134D+0
       V=0.1891714324525297D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5520413051846366D+0
       B=0.4860961284181720D+0
       V=0.1896827480450146D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5756887237503077D+0
       B=0.5103447395342790D+0
       V=0.1899628417059528D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1225039430588352D+0
       B=0.2136455922655793D-1
       V=0.1123301829001669D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1539113217321372D+0
       B=0.4520926166137188D-1
       V=0.1253698826711277D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1856213098637712D+0
       B=0.7086468177864818D-1
       V=0.1366266117678531D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2174998728035131D+0
       B=0.9785239488772918D-1
       V=0.1462736856106918D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2494128336938330D+0
       B=0.1258106396267210D+0
       V=0.1545076466685412D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2812321562143480D+0
       B=0.1544529125047001D+0
       V=0.1615096280814007D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3128372276456111D+0
       B=0.1835433512202753D+0
       V=0.1674366639741759D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3441145160177973D+0
       B=0.2128813258619585D+0
       V=0.1724225002437900D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3749567714853510D+0
       B=0.2422913734880829D+0
       V=0.1765810822987288D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4052621732015610D+0
       B=0.2716163748391453D+0
       V=0.1800104126010751D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4349335453522385D+0
       B=0.3007127671240280D+0
       V=0.1827960437331284D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4638776641524965D+0
       B=0.3294470677216479D+0
       V=0.1850140300716308D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4920046410462687D+0
       B=0.3576932543699155D+0
       V=0.1867333507394938D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5192273554861704D+0
       B=0.3853307059757764D+0
       V=0.1880178688638289D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5454609081136522D+0
       B=0.4122425044452694D+0
       V=0.1889278925654758D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5706220661424140D+0
       B=0.4383139587781027D+0
       V=0.1895213832507346D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5946286755181518D+0
       B=0.4634312536300553D+0
       V=0.1898548277397420D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1905370790924295D+0
       B=0.2371311537781979D-1
       V=0.1349105935937341D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2242518717748009D+0
       B=0.4917878059254806D-1
       V=0.1444060068369326D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2577190808025936D+0
       B=0.7595498960495142D-1
       V=0.1526797390930008D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2908724534927187D+0
       B=0.1036991083191100D+0
       V=0.1598208771406474D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3236354020056219D+0
       B=0.1321348584450234D+0
       V=0.1659354368615331D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3559267359304543D+0
       B=0.1610316571314789D+0
       V=0.1711279910946440D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3876637123676956D+0
       B=0.1901912080395707D+0
       V=0.1754952725601440D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4187636705218842D+0
       B=0.2194384950137950D+0
       V=0.1791247850802529D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4491449019883107D+0
       B=0.2486155334763858D+0
       V=0.1820954300877716D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4787270932425445D+0
       B=0.2775768931812335D+0
       V=0.1844788524548449D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5074315153055574D+0
       B=0.3061863786591120D+0
       V=0.1863409481706220D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5351810507738336D+0
       B=0.3343144718152556D+0
       V=0.1877433008795068D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5619001025975381D+0
       B=0.3618362729028427D+0
       V=0.1887444543705232D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5875144035268046D+0
       B=0.3886297583620408D+0
       V=0.1894009829375006D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6119507308734495D+0
       B=0.4145742277792031D+0
       V=0.1897683345035198D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2619733870119463D+0
       B=0.2540047186389353D-1
       V=0.1517327037467653D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2968149743237949D+0
       B=0.5208107018543989D-1
       V=0.1587740557483543D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3310451504860488D+0
       B=0.7971828470885599D-1
       V=0.1649093382274097D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3646215567376676D+0
       B=0.1080465999177927D+0
       V=0.1701915216193265D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3974916785279360D+0
       B=0.1368413849366629D+0
       V=0.1746847753144065D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4295967403772029D+0
       B=0.1659073184763559D+0
       V=0.1784555512007570D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4608742854473447D+0
       B=0.1950703730454614D+0
       V=0.1815687562112174D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4912598858949903D+0
       B=0.2241721144376724D+0
       V=0.1840864370663302D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5206882758945558D+0
       B=0.2530655255406489D+0
       V=0.1860676785390006D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5490940914019819D+0
       B=0.2816118409731066D+0
       V=0.1875690583743703D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5764123302025542D+0
       B=0.3096780504593238D+0
       V=0.1886453236347225D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6025786004213506D+0
       B=0.3371348366394987D+0
       V=0.1893501123329645D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6275291964794956D+0
       B=0.3638547827694396D+0
       V=0.1897366184519868D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3348189479861771D+0
       B=0.2664841935537443D-1
       V=0.1643908815152736D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3699515545855295D+0
       B=0.5424000066843495D-1
       V=0.1696300350907768D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4042003071474669D+0
       B=0.8251992715430854D-1
       V=0.1741553103844483D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4375320100182624D+0
       B=0.1112695182483710D+0
       V=0.1780015282386092D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4699054490335947D+0
       B=0.1402964116467816D+0
       V=0.1812116787077125D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5012739879431952D+0
       B=0.1694275117584291D+0
       V=0.1838323158085421D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5315874883754966D+0
       B=0.1985038235312689D+0
       V=0.1859113119837737D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5607937109622117D+0
       B=0.2273765660020893D+0
       V=0.1874969220221698D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5888393223495521D+0
       B=0.2559041492849764D+0
       V=0.1886375612681076D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6156705979160163D+0
       B=0.2839497251976899D+0
       V=0.1893819575809276D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6412338809078123D+0
       B=0.3113791060500690D+0
       V=0.1897794748256767D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4076051259257167D+0
       B=0.2757792290858463D-1
       V=0.1738963926584846D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4423788125791520D+0
       B=0.5584136834984293D-1
       V=0.1777442359873466D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4760480917328258D+0
       B=0.8457772087727143D-1
       V=0.1810010815068719D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5085838725946297D+0
       B=0.1135975846359248D+0
       V=0.1836920318248129D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5399513637391218D+0
       B=0.1427286904765053D+0
       V=0.1858489473214328D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5701118433636380D+0
       B=0.1718112740057635D+0
       V=0.1875079342496592D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5990240530606021D+0
       B=0.2006944855985351D+0
       V=0.1887080239102310D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6266452685139695D+0
       B=0.2292335090598907D+0
       V=0.1894905752176822D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6529320971415942D+0
       B=0.2572871512353714D+0
       V=0.1898991061200695D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4791583834610126D+0
       B=0.2826094197735932D-1
       V=0.1809065016458791D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5130373952796940D+0
       B=0.5699871359683649D-1
       V=0.1836297121596799D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5456252429628476D+0
       B=0.8602712528554394D-1
       V=0.1858426916241869D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5768956329682385D+0
       B=0.1151748137221281D+0
       V=0.1875654101134641D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6068186944699046D+0
       B=0.1442811654136362D+0
       V=0.1888240751833503D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6353622248024907D+0
       B=0.1731930321657680D+0
       V=0.1896497383866979D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6624927035731797D+0
       B=0.2017619958756061D+0
       V=0.1900775530219121D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5484933508028488D+0
       B=0.2874219755907391D-1
       V=0.1858525041478814D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5810207682142106D+0
       B=0.5778312123713695D-1
       V=0.1876248690077947D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6120955197181352D+0
       B=0.8695262371439526D-1
       V=0.1889404439064607D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6416944284294319D+0
       B=0.1160893767057166D+0
       V=0.1898168539265290D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6697926391731260D+0
       B=0.1450378826743251D+0
       V=0.1902779940661772D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6147594390585488D+0
       B=0.2904957622341456D-1
       V=0.1890125641731815D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6455390026356783D+0
       B=0.5823809152617197D-1
       V=0.1899434637795751D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6747258588365477D+0
       B=0.8740384899884715D-1
       V=0.1904520856831751D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6772135750395347D+0
       B=0.2919946135808105D-1
       V=0.1905534498734563D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
  END subroutine
end module tools_math
