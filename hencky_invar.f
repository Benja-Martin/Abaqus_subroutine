      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA) 
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
	  
      DOUBLE PRECISION K2, K3
      REAL HENCK(3,3), HENCKSQ(3,3)

!!! This subroutine is used to compute the invariants of the Hencky tensor


        CALL GETVRM('LE',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)

      HENCK(1,1) = ARRAY(1) ! Symmetric tensor
      HENCK(2,2) = ARRAY(2)
      HENCK(3,3) = ARRAY(3)
      HENCK(1,2) = ARRAY(4)
      HENCK(1,3) = ARRAY(5)
      HENCK(2,3) = ARRAY(6)
      HENCK(2,1) = ARRAY(4)
      HENCK(3,1) = ARRAY(5)
      HENCK(3,2) = ARRAY(6)
      CALL MULTI33(HENCK,HENCK,HENCKSQ)
	  
      K2 = 0.5*(
     1      (HENCK(1,1) + HENCK(2,2) + HENCK(3,3))**2 + 
     2      (HENCKSQ(1,1) + HENCKSQ(2,2) + HENCKSQ(3,3)))
	  
      K3 = 7.348*HENCK(1,1)*HENCK(2,2)*HENCK(3,3)
	  
	  
C         CALL GETVRM('LEP',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
C      1 MATLAYO,LACCFLA)

C       LEPRINCIP(1) = ARRAY(1) ! Minimum principal value : ln(lambda_1)
C       LEPRINCIP(2) = ARRAY(2)
C       LEPRINCIP(3) = ARRAY(3) ! Maximum principal value ln : (lambda_3)

C       K2 = SQRT((LOG(LEPRINCIP(1)))**2 + (LOG(LEPRINCIP(2)))**2 + 
C      1     (LOG(LEPRINCIP(2)))**2)
	  
C       K3 = (3*SQRT(6)/(K2**3))*LOG(LEPRINCIP(1))*LOG(LEPRINCIP(2))*
C      1      LOG(LEPRINCIP(3))
	  
	  
	  
c 	  HENCKY strain
      UVAR(1) = HENCK(1,1)
      UVAR(2) = HENCK(2,2)
      UVAR(3) = HENCK(3,3)
      UVAR(4) = HENCK(1,2)
      UVAR(5) = HENCK(1,3)
      UVAR(6) = HENCK(2,3)
	  
      UVAR(7) = K2
      UVAR(8) = K3
C       UVAR(9) = LEPRINCIP(3)
	  
C       UVAR(10) = PS(1)
C       UVAR(11) = PS(2)
C       UVAR(12) = PS(3)
	  
C       UVAR(13) = PE(1)
C       UVAR(14) = PE(2)
C       UVAR(15) = PE(3)
	  
C       UVAR(16) = ANPE(1,1)
C       UVAR(17) = ANPE(2,2)
C       UVAR(18) = ANPE(3,3)
C       UVAR(19) = ANPE(1,2)
C       UVAR(20) = ANPE(1,3)
C       UVAR(21) = ANPE(2,3)
C       UVAR(22) = ANPE(2,1)
C       UVAR(23) = ANPE(3,1)
C       UVAR(24) = ANPE(3,2)



      RETURN
      END


C **********************************************************************************
C **********************************************************************************
C **********************************************************************************
C 						SUBROUTINES
C **********************************************************************************
C **********************************************************************************
C **********************************************************************************


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!    Subroutine to compute the inverse of a 3x3 matrix.
!!!!	It uses the determinant and inverse of the transpose of the comatrix.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       SUBROUTINE M33INV(A, AINV, OK_FLAG, DET)


C       REAL A(3,3), AINV(3,3), COFACTOR(3,3)
C       LOGICAL OK_FLAG
C       PARAMETER EPS = 1.0D-10
C       DOUBLE PRECISION DET


C       DET =   A(1,1)*A(2,2)*A(3,3)
C      1       - A(1,1)*A(2,3)*A(3,2)
C      2       - A(1,2)*A(2,1)*A(3,3)
C      3       + A(1,2)*A(2,3)*A(3,1)
C      4       + A(1,3)*A(2,1)*A(3,2)
C      5       - A(1,3)*A(2,2)*A(3,1)

C       IF (ABS(DET) .LE. EPS) THEN
C          AINV = 0.0D0
C          OK_FLAG = .FALSE.
C          RETURN
C       END IF

C       COFACTOR(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))
C       COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
C       COFACTOR(1,3) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))
C       COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
C       COFACTOR(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))
C       COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
C       COFACTOR(3,1) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))
C       COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
C       COFACTOR(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))

C       AINV = TRANSPOSE(COFACTOR) / DET
C       OK_FLAG = .TRUE.

C       RETURN
C       END
	  
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C !!!!    Subroutine to compute the multiplication of 2 3x3 matrices.

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE MULTI33(A,B,RES)

      REAL A(3,3), B(3,3), RES(3,3)
	  
      DO I=1, 3
         DO J=1, 3
             RES(I,J) = 0
             DO P=1, 3
                 RES(I,J) = RES(I,J) + (A(I,P)*B(P,J))
             END DO
         END DO
      END DO
	  
      RETURN
      END
	  
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C !!!!    Subroutine to compute the addition of 2 3x3 matrices.

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
C       SUBROUTINE ADD33(A,B,RES)
	  
C       REAL A(3,3), B(3,3), RES(3,3)
	  
C       DO I=1, 3
C          DO J=1, 3
C              RES(I,J) = A(I,J) + B(I,J)
C          END DO
C       END DO
	  
C       RETURN
C       END
	  
	  





