!DEC$ ATTRIBUTES ALIAS:"uvarm"::UVARM
!!! Subroutine Abaqus calculant le critere de cavitation de Lopez-Pamies
!!! dans le cas d'un matériau neo-hookéen. Pour une définition plus 
!!! complète, l'utilisateur est renvoyé vers "Lopez-Pamies et al. 2011
!!! cavitation in elastomeric solids: II" notamment les équations 31 et 46
!!! ainsi que "Lefevre et al., 2015, Cavitation in rubber: an elastic
!!! instability or a fracture phenomenon?" et l'equation 9
!!! La valeur de UVAR(1) est 1 si les conditions sont reunies et 0 sinon
!!!
!!! Abaqus subroutine that compute the Lopez-Pamies cavitation criterion
!!! applied to a neo-hookean material 
!!! cf. Lopez-Pamies et al., 2011, Cavitation in elastomeric solids: II—Onset-of-cavitation surfaces for Neo-Hookean materials
!!! cf. Lefevre et al., 2015, Cavitation in rubber: an elastic instability or a fracture phenomenon?
!!! The value of UVAR(1) is 1 if the conditions are met and 0 otherwise.
C
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,
     1 DTIME,CMNAME,ORNAME, NUVARM,NOEL,
     2 NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,
     3 COORD,JMAC,JMATYP,MATLAYO,LACCFLA) 
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3)
      DIMENSION T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15)
      DIMENSION JMAC(*),JMATYP(*),COORD(*)
      DOUBLE PRECISION T1, T2, T3
      DOUBLE PRECISION MU
C
!!! Valeur du module de cisaillement à fixer.
!!! Shear modulus value 
      MU = 2.D0
C
      CALL GETVRM('SP',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
      T1 = ARRAY(1) ! S_min
      T2 = ARRAY(2) ! S_mid
      T3 = ARRAY(3) ! S_max
C
      CRIT = 8.D0*(T1*T2*T3) - 12.D0*MU*(T1*T2+T2*T3+T3*T1) + 18.D0*MU**2*(T1+T2+T3) - 35.D0*MU**3
	IF ( (CRIT>0) .AND. (T1>3/2*MU) .AND. (T2>3/2*MU) .AND. (T3>3/2*MU) ) THEN
      		UVAR(1) = 1.D0
	ELSE 
		UVAR(1) = 0.D0
	END IF 
C
      RETURN
      END