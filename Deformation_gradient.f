!DEC$ ATTRIBUTES ALIAS:"uvarm"::UVARM
!!! Subroutine Abaqus permettant de calculer les composantes 
!!! du tenseur gradient de la deformation. Attention en utilisant
!!! ces donnees, elles sont exprimees avec un repere qui n'est
!!! pas forcement celui dans lequel les chargements ont ete donnes
!!! 
!!! Abaqus subroutine for calculating the components of the 
!!! deformation gradient tensor. Be careful when using these
!!! data, they are expressed with a reference frame that is 
!!! not necessarily the one in which the loads were given.
C
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA) 
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*), COORD(*)
C
C
        CALL GETVRM('DG',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)
C
      UVAR(1) = ARRAY(1) ! Non-symmetric tensor
      UVAR(2) = ARRAY(2)
      UVAR(3) = ARRAY(3)
      UVAR(4) = ARRAY(4)
      UVAR(5) = ARRAY(5)
      UVAR(6) = ARRAY(6)
      UVAR(7) = ARRAY(7)
      UVAR(8) = ARRAY(8)
      UVAR(9) = ARRAY(9)
C
      RETURN
      END