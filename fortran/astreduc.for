*
*   -------------------------------------------------------------------
*
*                              ASTREDUC.FOR
*
*   This file contains Astrodynamic subroutines and functions to
*   implement reduction calculations. These routines are described in Ch3.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                   2007
*                             by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*            31 may 07  david vallado
*                           3rd edition baseline
*    changes :
*            28 Jan 04  David Vallado
*                         Fix headers
*             3 Apr 03  David Vallado
*                         New baseline
*            14 Feb 03  David Vallado
*                         Misc updates
*             2 Oct 02  David Vallado
*                         Reform to match Matlab scripts
*            12 Mar 02  David Vallado
*                         Add time record generation and searches
*                         Note - ER and km both work here except sidereal time
*                         Update documentation
*                         Fix AST for 1 Jan 97 change
*                         Add combined routine GCRF_ITRF
*            14 May 01  David Vallado
*                         2nd edition baseline
*            23 Nov 87  David Vallado
*                         Original baseline
*
*     *****************************************************************
*
*     Uses object files:
*       astutil,
*       astmath,
*       asttime
*     Uses common files:
*       astmath.cmn
*       astreduc.cmn
*
*
*  The routines for this module are listed below. They are broken into 
*  distinct areas. Notice that the arguments for the conversion routines
*  include the YMD HMS in UTC, while the analysis routines input the TTDB
*  term. This is intended to facilitate both functions. The CONVTIME routine
*  is used to aid in the conversion of time for any analysis functions.
*

*  ------------------- Conversion of coordinates -----------------------
*      SUBROUTINE InitReduc   ( FileN1 )
*
*
*      SUBROUTINE GCRF_ITRF     ( rj2000,vj2000, Direction, rITRF,vITRF,
*     &                        TTT, JDUT1, LOD, xp, yp, terms )
*
*      SUBROUTINE FK4         ( rJ2000, vJ2000, Direction, rFK4, vFK4 )
*
*  ------------ Individual routines for testing and analysis -----------
*
*      SUBROUTINE Precession  ( TTT, Prec )
*
*      SUBROUTINE GCRF_MOD    ( rJ2000, vJ2000, Direction, rMOD, vMOD,
*     &                         TTT )
*
*      SUBROUTINE Nutation    ( TTT,
*     &                         DeltaPsi, TrueEps, MeanEps, Omega, Nut )
*
*      SUBROUTINE GCRF_TOD    ( rj2000, vj2000, Direction, rTOD, vTOD,
*     &                         TTT )
*
*      SUBROUTINE SIDEREAL    ( JDUT1,DeltaPsi,MeanEps,Omega,LOD,
*     &                         ST,STDot,OmegaEarth,terms )
*
*      SUBROUTINE GCRF_PEF    ( rj2000, vj2000, Direction, rPEF, vPEF,
*     &                         TTT, JDUT1, LOD, terms )
*
*      SUBROUTINE POLARM      ( xp, yp,  PM )
*
*
*      SUBROUTINE TrueMean    ( Order, TTT, Terms, Nutteme )
*
*
*      SUBROUTINE GCRF_TEME   ( rMOD, vMOD, Direction, rTM, vTM,
*     &                         order, TTT, Terms )
*
*

*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE INITREDUC
*
*  This procedure initializes the nutation matricies needed for reduction
*    calculations. The routine needs the filename of the files as input.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    FileN1      - Name for nutation file         nutation.dat
*    FileN2      - Name for planetary nutation    nut85.dat
*
*  Outputs       :
*    IAr80       - Integers for FK5 1980
*    RAr80       - Reals for FK5 1980
*    IAr00       - Integers for IAU 2000
*    RAr00       - Reals for IAU 2000
*    IAr96       - Integers for IAU 1996
*    RAr96       - Reals for IAU 1996
*    pIAr96      - Integers for planetary IAU 1996
*    pRAr96      - Reals for planetary IAU 1996
*
*  Locals        :
*    convrt      - conversion factor to degrees
*    i           - Index
*
*  Coupling      :
*    None        -
*
*  References    :
*
* -------------------------------------------------------------------------------

      SUBROUTINE InitReduc   ( FileN1 )
c      SUBROUTINE InitReduc   ( FileN1, FileN2, IAr00, RAr00, IAr80,
c     &                         RAr80, IAr96, RAr96, pIAr96, pRAr96 )
       IMPLICIT NONE
       CHARACTER*64 FileN1
c       INTEGER IAr00(5,106), IAr96(5,263), pIAr96(10,112)
c       REAL*8  RAr00(4,106), RAr96(4,263), pRAr96(4,112)

        INCLUDE 'astreduc.cmn'

* ----------------------------  Locals  -------------------------------
       INTEGER i
       REAL*8 Convrt

        ! --------------------  Implementation   ----------------------
c       Convrt= 0.0001D0/3600.0D0 ! 0.0001" to deg

       ! Open this for future use in the program
       OPEN( UNIT=44, FILE='timerec.rec',FORM='UNFORMATTED',
     &       ACCESS='DIRECT',RECL=80,STATUS='UNKNOWN')

       OPEN( 42, FILE = FileN1, STATUS = 'OLD' )
       READ( 42,*)
       READ( 42,*)
       ! ---------------- Read in IAU 2000 Theroy coefficients --------
       DO i = 1, 106
         Read( 42,*)
c         Read( 42,*) IAr00(1,i),IAr00(2,i),IAr00(3,i),IAr00(4,i),
c     &               IAr00(5,i),RAr00(1,i),RAr00(2,i),RAr00(3,i),
c     &               RAr00(4,i),RAr00(5,i),RAr00(6,i)
c         RAr00(1,i)= RAr00(1,i) * Convrt
c         RAr00(2,i)= RAr00(2,i) * Convrt
c         RAr00(3,i)= RAr00(3,i) * Convrt
c         RAr00(4,i)= RAr00(4,i) * Convrt
c         RAr00(5,i)= RAr00(5,i) * Convrt
c         RAr00(6,i)= RAr00(6,i) * Convrt
       ENDDO

       ! ---------------- Read in 1980 IAU Theory of Nutation ---------
       READ( 42,*)
       READ( 42,*)
       Convrt= 0.0001D0/3600.0D0 ! 0.0001" to deg
       DO i = 1, 106
         Read( 42,*) IAr80(1,i),IAr80(2,i),IAr80(3,i),IAr80(4,i),
     &               IAr80(5,i),RAr80(1,i),RAr80(2,i),RAr80(3,i),
     &               RAr80(4,i)
         RAr80(1,i)= RAr80(1,i) * Convrt
         RAr80(2,i)= RAr80(2,i) * Convrt
         RAr80(3,i)= RAr80(3,i) * Convrt
         RAr80(4,i)= RAr80(4,i) * Convrt
       ENDDO

       ! ---------------- Read in 1996 IAU Theory of Nutation ---------
c       READ( 42,*)
c       READ( 42,*)
c       Convrt= 0.000001D0/3600.0D0 ! 0.000001" to deg
c       DO i = 1, 263
c         Read( 42,*) IAr96(1,i),IAr96(2,i),IAr96(3,i),IAr96(4,i),
c     &               IAr96(5,i),RAr96(1,i),RAr96(2,i),RAr96(3,i),
c     &               RAr96(4,i),RAr96(5,i),RAr96(6,i)
c         RAr96(1,i)= RAr96(1,i) * Convrt
c         RAr96(2,i)= RAr96(2,i) * Convrt
c         RAr96(3,i)= RAr96(3,i) * Convrt
c         RAr96(4,i)= RAr96(4,i) * Convrt
c         RAr96(5,i)= RAr96(5,i) * Convrt
c         RAr96(6,i)= RAr96(6,i) * Convrt
c       ENDDO
c
c       CLOSE( 42 )

       ! ------------- Read in 1996 planetary correction terms --------
c       OPEN( 42, FILE = FileN2, STATUS = 'OLD' )
c       READ( 42,*)
c       READ( 42,*)
c       Convrt= 0.000001D0/3600.0D0 ! 0.000001" to deg
c       DO i= 1, 112
c         Read( 42,*) pIAr96(1,i),pIAr96(2,i),pIAr96(3,i),pIAr96(4,i),
c     &               pIAr96(5,i),pIAr96(6,i),pIAr96(7,i),pIAr96(8,i),
c     &               pIAr96(9,i),pIAr96(10,i),pRAr96(1,i),pRAr96(2,i),
c     &               pRAr96(3,i),pRAr96(4,i)
c         pRAr96(1,i)= pRAr96(1,i) * Convrt
c         pRAr96(2,i)= pRAr96(2,i) * Convrt
c         pRAr96(3,i)= pRAr96(3,i) * Convrt
c         pRAr96(4,i)= pRAr96(4,i) * Convrt
c        ENDDO
c
       CLOSE( 42 )

      RETURN
      END ! Subroutine InitReduc
*

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE FK4
*
*  this subroutine converts vectors between the B1950 and J2000 epochs.  Be
*    aware that this process is not exact. There are different secular rates
*    for each system, and there are differences in the central location. The
*    matrices are multiplied directly for speed.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position vec ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity vec
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    rFK4        - Position vector FK4
*    vFK4        - Velocity vector FK4
*
*  Locals        :
*    r11,r12,r13 - Components of rotation matrix
*    r21,r22,r23 - Components of rotation matrix
*    r31,r32,r33 - Components of rotation matrix
*
*  Coupling      :
*
*
*  References    :
*    Vallado       2007, 240
*
* ----------------------------------------------------------------------------

      SUBROUTINE FK4         ( rJ2000, vJ2000, Direction, rFK4, vFK4 )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rJ2000(3), vJ2000(3), rFK4(3), vFK4(3)

* ----------------------------  Locals  -------------------------------
        REAL*8 r11, r12, r13, r21, r22, r23, r31, r32, r33

        r11=  0.9999256794956877D0
        r12= -0.0111814832204662D0 
        r13= -0.0048590038153592D0 

        r21=  0.0111814832391717D0
        r22=  0.9999374848933135D0 
        r23= -0.0000271625947142D0 

        r31=  0.0048590037723143D0 
        r32= -0.0000271702937440D0 
        r33=  0.9999881946043742D0 

        IF ( Direction .eq. 'TOO ' ) THEN
            rFK4(1) = r11*rJ2000(1) + r21*rJ2000(2) + r31*rJ2000(3)
            rFK4(2) = r12*rJ2000(1) + r22*rJ2000(2) + r32*rJ2000(3)
            rFK4(3) = r13*rJ2000(1) + r23*rJ2000(2) + r33*rJ2000(3)
            vFK4(1) = r11*vJ2000(1) + r21*vJ2000(2) + r31*vJ2000(3)
            vFK4(2) = r12*vJ2000(1) + r22*vJ2000(2) + r32*vJ2000(3)
            vFK4(3) = r13*vJ2000(1) + r23*vJ2000(2) + r33*vJ2000(3)
          ELSE
            rJ2000(1) = r11*rFK4(1) + r12*rFK4(2) + r13*rFK4(3)
            rJ2000(2) = r21*rFK4(1) + r22*rFK4(2) + r23*rFK4(3)
            rJ2000(3) = r31*rFK4(1) + r32*rFK4(2) + r33*rFK4(3)
            vJ2000(1) = r11*vFK4(1) + r12*vFK4(2) + r13*vFK4(3)
            vJ2000(2) = r21*vFK4(1) + r22*vFK4(2) + r23*vFK4(3)
            vJ2000(3) = r31*vFK4(1) + r32*vFK4(2) + r33*vFK4(3)
          ENDIF

      RETURN
      END   ! SUBROUTINE FK4
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE PRECESSION
*
*  this function calulates the transformation matrix for precession.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    TTT         - Julian Centuries of TT         centuries
*
*  Outputs       :
*    Prec        - Precession transformation (eci-mod)
*
*  Locals        :
*    TTT2        - TTT squared
*    TTT3        - TTT cubed
*    Zeta        - PRECESSION ANGLE               rad
*    z           - PRECESSION ANGLE               rad
*    Theta       - PRECESSION ANGLE               rad
*
*  Coupling      :
*    none
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE PRECESSION  ( TTT, Prec )
        IMPLICIT NONE
        REAL*8 TTT, Prec(3,3)

* ----------------------------  Locals  -------------------------------
        REAL*8 zeta, z, theta
        REAL*8 TTT2, TTT3
        REAL*8 coszeta, sinzeta, costheta, sintheta, cosz, sinz


        INCLUDE 'astmath.cmn'

        ! --------------------- PRECESSION angles ---------------------
        TTT2= TTT * TTT
        TTT3= TTT2 * TTT
        Zeta = 2306.2181D0*TTT + 0.30188D0*TTT2 + 0.017998D0*TTT3
        Theta= 2004.3109D0*TTT - 0.42665D0*TTT2 - 0.041833D0*TTT3
        z    = 2306.2181D0*TTT + 1.09468D0*TTT2 + 0.018203D0*TTT3

        Zeta = Zeta  * Deg2Rad / 3600.0D0
        Theta= Theta * Deg2Rad / 3600.0D0
        Z    = Z     * Deg2Rad / 3600.0D0

        coszeta   = DCOS(zeta)
        sinzeta   = DSIN(zeta)
        costheta = DCOS(theta)
        sintheta = DSIN(theta)
        cosz    = DCOS(z)
        sinz    = DSIN(z)

        ! ----------------- form matrix  J2000 to MOD -----------------
        prec(1,1) =  coszeta * costheta * cosz - sinzeta * sinz
        prec(1,2) = -sinzeta * costheta * cosz - coszeta * sinz
        prec(1,3) = -sintheta * cosz

        prec(2,1) =  coszeta * costheta * sinz + sinzeta * cosz
        prec(2,2) = -sinzeta * costheta * sinz + coszeta * cosz
        prec(2,3) = -sintheta * sinz

        prec(3,1) =  coszeta * sintheta
        prec(3,2) = -sinzeta * sintheta
        prec(3,3) =  costheta

      RETURN
      END   ! SUBROUTINE PRECESSION
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_MOD
*
*  this function transfroms a vector between the mean equator mean equinox of
*    epoch (eci) and the mean equator mean equinox of date (mod).
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    TTT         - Julian Centuries of TT         centuries
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    rMOD        - Position vector of date        ER, km, etc
*                    Mean Equator, Mean Equinox
*    vMOD        - Velocity vector of date
*                    Mean Equator, Mean Equinox
*
*  Locals        :
*    Prec        - Precession transformation (eci-mod)
*    tmt         - transformation matrix transpose
*    Zeta        - PRECESSION ANGLE               rad
*    z           - PRECESSION ANGLE               rad
*    Theta       - PRECESSION ANGLE               rad
*
*  Coupling      :
*    PRECESSION  - Find precession matrix
*    ROT2        - Rotation about the second axis
*    ROT3        - Rotation about the third axis
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE GCRF_MOD  ( rJ2000, vJ2000, Direction, rMOD, vMOD, TTT)
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rJ2000(3), vJ2000(3), rMOD(3), vMOD(3), TTT

* ----------------------------  Locals  -------------------------------
        REAL*8 Prec(3,3),tmt(3,3)

        ! ---------------- get transformation matrices ----------------
        CALL Precession ( TTT,Prec)

        ! ------------------- Perform matrix mmpy ---------------------
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( Prec, rJ2000, 3,3,1,3,3,3, rmod )
            CALL MATMULT     ( Prec, vJ2000, 3,3,1,3,3,3, vmod )
          ELSE
            CALL MATTRANS( prec, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt, rmod  , 3,3,1,3,3,3, rj2000 )
            CALL MATMULT     ( tmt, vmod  , 3,3,1,3,3,3, vj2000 )
          ENDIF

c       ! --------------------- Perform rotations ---------------------
c        IF ( Direction .eq. 'TOO ' ) THEN
c            CALL ROT3( rJ2000, -Zeta, Temp1 )
c            CALL ROT2( Temp1 , Theta, Temp  )
c            CALL ROT3( Temp  ,  -z  , rMOD  )
c            CALL ROT3( vJ2000, -Zeta, Temp1 )
c            CALL ROT2( Temp1 , Theta, Temp  )
c            CALL ROT3( Temp  ,  -z  , vMOD  )
c          ELSE
c            CALL ROT3( rMOD  ,   z  , Temp   )
c            CALL ROT2( Temp  ,-Theta, Temp1  )
c            CALL ROT3( Temp1 ,  Zeta, rJ2000 )
c            CALL ROT3( vMOD  ,   z  , Temp   )
c            CALL ROT2( Temp  ,-Theta, Temp1  )
c            CALL ROT3( Temp1 ,  Zeta, vJ2000 )
c          ENDIF

      RETURN
      END   ! SUBROUTINE GCRF_MOD
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE NUTATION
*
*  this function calulates the transformation matrix for nutation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    TTT         - Julian Centuries of TT
*
*  Outputs       :
*    Nut         - Nutation transformation (mod-tod)
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    Omega       -                                rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*
*  Locals        :
*    TTT2        - TTT squared
*    TTT3        - TTT cubed
*    MeanEps     - Mean obliquity of the ecliptic rad
*    l           -                                rad
*    ll          -                                rad
*    F           -                                rad
*    D           -                                rad
*    DeltaEps    - Change in obliquity            rad
*
*  Coupling      :
*    none
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE Nutation    ( TTT,
     &                         DeltaPsi, TrueEps, MeanEps, Omega, Nut )
        IMPLICIT NONE
        REAL*8  Omega,
     &          DeltaPsi, TrueEps, TTT, Nut(3,3),meaneps

* ----------------------------  Locals  -------------------------------
        INTEGER i
        REAL*8 cospsi,sinpsi,coseps,sineps,costrueeps,sintrueeps
        REAL*8  DeltaEps, Tempval, TTT4, TTT2, TTT3, l, l1, f,d

        INCLUDE 'astreduc.cmn'
!         INCLUDE 'astmath.cmn'
!         INCLUDE 'astconst.cmn'
        REAL*8 Small,    Rad2Deg,  Deg2Rad,  HalfPi,
     &         Pi,       TwoPi,    Infinite, Undefined

        REAL*8     rekm,     mu,     omegaearth, flat,     EESqrd, auer

        DATA rekm       /6378.137D0/
c        DATA mu         /398600.4418D0/
        DATA mu         /1.32712440018D11/       ! Sun
        DATA omegaearth /7.2921158553D-5/
        DATA flat       /0.003352810665D0/       ! f = 1.0/298.257223563
        DATA EESqrd     /0.006694379990D0/       ! 2f - f**2
        DATA AUER       /23454.79095228D0/       ! 149597870.0/6378.137

        DATA Small      /0.00000001D0/
        DATA Infinite   /999999.9D0/
        DATA Undefined  /999999.1D0/
c        DATA Halfpi     /1.57079632679489662D0/
c        DATA Pi         /3.14159265358979324D0/
c        DATA TwoPi      /6.28318530717958648D0/
c        DATA Rad2Deg    /57.2957795130823208D0/
c        DATA Deg2Rad    /0.01745329251994329D0/

c       use machine precision instead
        Pi      =  4.0D0 * DATAN(1.0D0)
        HalfPi  =  0.5D0*pi
        TwoPi   =  2.0D0*pi
        Rad2Deg = 180.0D0/pi
        Deg2Rad = pi/180.0D0

        ! ---- Determine coefficients for IAU 1980 NUTATION Theory ----
        TTT2= TTT*TTT
        TTT3= TTT2*TTT
        TTT4= TTT2*TTT2

        MeanEps = -46.8150D0*TTT - 0.00059D0*TTT2 + 0.001813D0*TTT3 +
     &            84381.448D0
        MeanEps = DMOD( MeanEps/3600.0D0,360.0D0 )
        MeanEps = MeanEps * Deg2Rad

c     ---- Old values ----
c        rr   = 360.0
c        l    =  134.9629814 + (1325*rr  + 198.8673981)*TTT + 0.0086972 *TTT2 + 0.00001778*TTT3
c        l1   =  357.5277233 + (  99*rr  + 359.05034  )*TTT - 0.00016028*TTT2 - 0.00000333*TTT3
c        F    =   93.2719103 + (1342*rr  +  82.0175381)*TTT - 0.0036825 *TTT2 + 0.00000306*TTT3
c        D    =  297.8503631 + (1236*rr  + 307.111480 )*TTT - 0.00191417*TTT2 + 0.00000528*TTT3
c        Omega=  125.0445222 - (   5*rr  + 134.1362608)*TTT + 0.0020708 *TTT2 + 0.00000222*TTT3

        l    =  134.96340251D0 + ( 1717915923.2178D0*TTT +
     &          31.8792D0*TTT2 + 0.051635D0*TTT3 - 0.00024470D0*TTT4 )
     &          / 3600.0D0
        l1   =  357.52910918D0 + (  129596581.0481D0*TTT -
     &           0.5532D0*TTT2 - 0.000136D0*TTT3 - 0.00001149*TTT4 )
     &          / 3600.0D0
        F    =   93.27209062D0 + ( 1739527262.8478D0*TTT -
     &          12.7512D0*TTT2 + 0.001037D0*TTT3 + 0.00000417*TTT4 )
     &          / 3600.0D0
        D    =  297.85019547D0 + ( 1602961601.2090D0*TTT -
     &           6.3706D0*TTT2 + 0.006593D0*TTT3 - 0.00003169*TTT4 )
     &          / 3600.0D0
        Omega=  125.04455501D0 + (   -6962890.2665D0*TTT +
     &           7.4722D0*TTT2 + 0.007702D0*TTT3 - 0.00005939*TTT4 )
     &          / 3600.0D0

        l    = DMOD( l,360.0D0 )     * Deg2Rad
        l1   = DMOD( l1,360.0D0 )    * Deg2Rad
        F    = DMOD( F,360.0D0 )     * Deg2Rad
        D    = DMOD( D,360.0D0 )     * Deg2Rad
        Omega= DMOD( Omega,360.0D0 ) * Deg2Rad

        DeltaPsi= 0.0D0
        DeltaEps= 0.0D0

        DO i= 106, 1, -1
            Tempval= IAr80(1,i)*l + IAr80(2,i)*l1 + IAr80(3,i)*F +
     &               IAr80(4,i)*D + IAr80(5,i)*Omega
            DeltaPsi= DeltaPsi + (RAr80(1,i)+RAr80(2,i)*TTT) *
     &                 DSIN( TempVal )
            DeltaEps= DeltaEps + (RAr80(3,i)+RAr80(4,i)*TTT) *
     &                 DCOS( TempVal )
          ENDDO

        ! --------------- Find NUTATION Parameters --------------------
        DeltaPsi = DMOD( DeltaPsi,360.0D0 ) * Deg2Rad
        DeltaEps = DMOD( DeltaEps,360.0D0 ) * Deg2Rad
        TrueEps  = MeanEps + DeltaEps

        cospsi  = DCOS(deltapsi)
        sinpsi  = DSIN(deltapsi)
        coseps  = DCOS(meaneps)
        sineps  = DSIN(meaneps)
        costrueeps = DCOS(trueeps)
        sintrueeps = DSIN(trueeps)

        nut(1,1) =  cospsi
        nut(1,2) = -coseps * sinpsi
        nut(1,3) = -sineps * sinpsi

        nut(2,1) =  costrueeps * sinpsi
        nut(2,2) =  costrueeps * coseps * cospsi + sintrueeps * sineps
        nut(2,3) =  costrueeps * sineps * cospsi - sintrueeps * coseps

        nut(3,1) =  sintrueeps * sinpsi
        nut(3,2) =  sintrueeps * coseps * cospsi - sineps * costrueeps
        nut(3,3) =  sintrueeps * sineps * cospsi + costrueeps * coseps

      RETURN
      END   ! SUBROUTINE NUTATION
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_TOD
*
*  this function transforms vectors between the mean equator mean equinox of
*    date (j2000) and the true equator true equinox of date (tod).
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    TTT         - Julian Centuries of TT
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    rTOD        - Position vector of date
*                    True Equator, True Equinox
*    vTOD        - Velocity vector of date
*                    True Equator, True Equinox
*
*  Locals        :
*    Prec        - Precession transformation (eci-mod)
*    Nut         - Nutation transformation (mod-tod)
*    tmt         - transformation matrix transpose
*    Omega       -                                rad
*    DeltaEps    - Change in obliquity            rad
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*
*  Coupling      :
*    PRECESSION  - Find precession matrix
*    NUTATION    - Find nutation matrix
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------
      SUBROUTINE GCRF_TOD    ( rj2000, vj2000, Direction, rTOD, vTOD,
     &                        TTT )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rj2000(3),vj2000(3),rtod(3),vtod(3), TTT

        INCLUDE 'astreduc.cmn'

* ----------------------------  Locals  -------------------------------
        REAL*8  meaneps,omega,DeltaPsi, TrueEps,Prec(3,3),Nut(3,3),
     &          tm(3,3),tmt(3,3)

        ! ---------------- get transformation matrices ----------------
        CALL Precession ( TTT,Prec)

        CALL nutation( ttt,  deltapsi,trueeps,meaneps,omega,nut)

        ! ------------------- Perform matrix mmpy ---------------------
        CALL MATMULT     ( Nut , Prec  , 3,3,3,3,3,3, tm )
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( tm  , rJ2000, 3,3,1,3,3,3, rtod )
            CALL MATMULT     ( tm  , vJ2000, 3,3,1,3,3,3, vtod )
          ELSE
            CALL MATTRANS( tm, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt , rtod  , 3,3,1,3,3,3, rj2000 )
            CALL MATMULT     ( tmt , vtod  , 3,3,1,3,3,3, vj2000 )
          ENDIF

c        ! --------------------- Perform rotations ---------------------
c        IF ( Direction .eq. 'TOO ' ) THEN
c            CALL ROT1( rMOD ,   Eps       , Temp  )
c            CALL ROT3( Temp , -DeltaPsi, Temp1 )
c            CALL ROT1( Temp1, -TrueEps , rTOD  )
c            CALL ROT1( vMOD ,   Eps       , Temp  )
c            CALL ROT3( Temp , -DeltaPsi, Temp1 )
c            CALL ROT1( Temp1, -TrueEps , vTOD  )
c          ELSE
c            CALL ROT1( rTOD ,  TrueEps , Temp  )
c            CALL ROT3( Temp ,  DeltaPsi, Temp1 )
c            CALL ROT1( Temp1,  -Eps       , rMOD  )
c            CALL ROT1( vTOD ,  TrueEps , Temp  )
c            CALL ROT3( Temp ,  DeltaPsi, Temp1 )
c            CALL ROT1( Temp1,  -Eps       , vMOD  )
c          ENDIF

      RETURN
      END   ! SUBROUTINE GCRF_TOD
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE SIDEREAL
*
*  this function calulates the transformation matrix that accounts for the
*    effects of nutation. Notice that deltaspi should not be moded to a
*    positive number because it is multiplied rather than used in a
*    trigonometric argument.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JDUT1       - Julian Date of UT1             days from 4713 BC
*    DeltaPsi    - NUTATION ANGLE                 rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*    Omega       -                                rad
*    LOD         - Excess length of day           sec
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    St          - Sidereal Time transformation (tod-pef)
*    StDot       - Sidereal Time rate transformation (tod-pef)
*    Omegaearth  - rotation of the earth          rad
*
*  Locals        :
*    GST         - Mean Greenwich SIDEREAL Time   0 to 2Pi rad
*    AST         - Apparent GST                   0 to 2Pi rad
*    Hr          - hour                           hr
*    minute         - minutes                        minute
*    SEC         - seconds                        SEC
*
*  Coupling      :
*    none
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE SIDEREAL    ( JDUT1,DeltaPsi,MeanEps,Omega,LOD,
     &                         ST,STDot,OmegaEarth,terms )
        IMPLICIT NONE
        REAL*8 JDUT1, DeltaPsi, GSTIME, Omega, LOD,
     &         ThetaSa,St(3,3),stdot(3,3),Omegaearth(3)
        EXTERNAL GSTIME
        INTEGER terms

* ----------------------------  Locals  -------------------------------
        REAL*8 GST, AST, Conv1, meaneps

        INCLUDE 'astmath.cmn'

        OmegaEarth(1) = 0.0D0
        OmegaEarth(2) = 0.0D0
        OmegaEarth(3) = 0.0D0

        Conv1 = pi / (180.0D0*3600.0D0)

        ! ------------------------ Find Mean GST ----------------------
        GST= GSTIME( JDUT1 )

        ! ------------------------ Find Mean AST ----------------------
        IF ((terms.gt.0) .and. (JDUT1.gt.2450449.5D0)) THEN
            AST= GST + DeltaPsi* DCOS(MeanEps)
     &           + 0.00264D0*Conv1*DSIN(Omega)
     &           + 0.000063D0*Conv1*DSIN(2.0D0*Omega)
          ELSE
            AST= GST + DeltaPsi* DCOS(MeanEps)
          ENDIF

        st(1,1) =  DCOS(ast)
        st(1,2) =  DSIN(ast)
        st(1,3) =  0.0

        st(2,1) = -DSIN(ast)
        st(2,2) =  DCOS(ast)
        st(2,3) =  0.0

        st(3,1) =  0.0
        st(3,2) =  0.0
        st(3,3) =  1.0

        ! ------------ compute sidereal time rate matrix --------------
        ThetaSA   =  7.29211514670698D-05 * (1.0D0 - LOD/86400.0D0)
        omegaearth(3) = ThetaSa

        stdot(1,1) = -omegaearth(3) * DSIN(ast)
        stdot(1,2) =  omegaearth(3) * DCOS(ast)
        stdot(1,3) =  0.0

        stdot(2,1) = -omegaearth(3) * DCOS(ast)
        stdot(2,2) = -omegaearth(3) * DSIN(ast)
        stdot(2,3) =  0.0

        stdot(3,1) =  0.0
        stdot(3,2) =  0.0
        stdot(3,3) =  0.0

      RETURN
      END   ! SUBROUTINE SIDEREAL Time
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_PEF
*
*  this function transforms a vector between the mean equator, mean equinox of epoch
*    frame (j2000), and the pseudo earth fixed frame (pef).
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    Direction   - Which set of vars to output    FROM  TOO
*    TTT         - Julian Centuries of TT
*    JDUT1       - Julian Date of UT1             days from 4713 BC
*    LOD         - Excess length of day           sec
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    rPEF        - Position vector of date        ER
*                    True Equator, True Equinox
*    vPEF        - Velocity vector of date        ER/TU
*                    True Equator, True Equinox
*
*  Locals        :
*    Prec        - Precession transformation (eci-mod)
*    Nut         - Nutation transformation (mod-tod)
*    St          - Sidereal Time transformation (tod-pef)
*    StDot       - Sidereal Time rate transformation (tod-pef)
*    tmt         - transformation matrix transpose
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*    Omegaearth  - rotation of the earth          rad
*    Omega       -                                rad
*    GST         - Mean Greenwich SIDEREAL Time   0 to 2Pi rad
*    AST         - Apparent GST                   0 to 2Pi rad
*
*  Coupling      :
*    PRECESSION  - Find precession matrix
*    NUTATION    - Find nutation matrix
*    SIDEREAL    - Find sidereal time matrix
*
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE GCRF_PEF    ( rj2000, vj2000, Direction, rPEF, vPEF,
     &                        TTT, JDUT1, LOD, terms )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rj2000(3), vj2000(3), rPEF(3), vPEF(3), JDUT1, LOD
        INTEGER Terms
        EXTERNAL GSTIME

        INCLUDE 'astreduc.cmn'

* ----------------------------  Locals  -------------------------------
        REAL*8 ttt,meaneps,omegaearth(3),wcrossr(3),
     &         TrueEps, GSTIME, Omega, tm(3,3),tmt(3,3),
     &         Prec(3,3),Nut(3,3),st(3,3), DeltaPsi,
     &         stdot(3,3),tm1(3,3)

        ! ---------------- get transformation matrices ----------------
        CALL Precession ( TTT,Prec)

        CALL nutation( ttt,  deltapsi,trueeps,meaneps,omega,nut)

        CALL sidereal( jdut1,deltapsi,meaneps,omega,lod,  st,stdot,
     &                 OmegaEarth,terms )

        ! ------------------- Perform matrix mmpy ---------------------
        CALL MATMULT     ( Nut , Prec  , 3,3,3,3,3,3, tm1 )
        CALL MATMULT     ( St  , tm1    , 3,3,3,3,3,3, tm )
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( tm  , rJ2000, 3,3,1,3,3,3, rPEF )
            CALL MATMULT     ( tm  , vJ2000, 3,3,1,3,3,3, vPEF )
            CALL CROSS( OmegaEarth,rPEF,  wcrossr )
            CALL SUBVEC( vPEF, wcrossr,  vPEF )
          ELSE
            CALL MATTRANS( tm, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt , rPEF  , 3,3,1,3,3,3, rj2000 )
            CALL CROSS( OmegaEarth,rPEF,  wcrossr )
            CALL ADDVEC( vPEF, wcrossr,  vPEF )
            CALL MATMULT     ( tmt , vPEF  , 3,3,1,3,3,3, vj2000 )
          ENDIF

c       ! --------------------- Perform rotations ---------------------
c       IF ( Direction .eq. 'TOO ' ) THEN
c           CALL ROT3( rTOD,   AST  , rPEF )
c           CALL ROT3( vTOD,   AST  , Temp  )
c           vPEF(1)= Temp(1) + rPEF(2)*OmegaEarth(3)
c           vPEF(2)= Temp(2) - rPEF(1)*OmegaEarth(3)
c           vPEF(3)= Temp(3)
c         ELSE
c           CALL ROT3( rPEF, -AST  , rTOD )
c           DO i = 1,4
c               Temp(i)   = vPEF(i)
c             ENDDO
c           Temp(1)= vPEF(1) - rPEF(2)*OmegaEarth(3)
c           Temp(2)= vPEF(2) + rPEF(1)*OmegaEarth(3)
c           CALL ROT3( Temp,  -AST  , vTOD )
c         ENDIF

      RETURN
      END   ! SUBROUTINE GCRF_PEF
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE POLARM
*
*  this function calulates the transformation matrix for polar motion.
*    the units for polar motion are input in rad because it's more
*    efficient to do this in the main routine.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    xp          - Polar motion coefficient       rad
*    yp          - Polar motion coefficient       rad
*
*  Outputs       :
*    PM          - Polar motion transformation (pef-ecef)
*
*  Locals        :
*    None.
*
*  Coupling      :
*    ROT2        - Rotation about the second axis
*    ROT3        - Rotation about the third axis
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE POLARM      ( xp, yp,  PM )
        IMPLICIT NONE
        REAL*8 xp, yp, PM(3,3)

* ----------------------------  Locals  -------------------------------
        REAL*8 cosxp,cosyp,sinxp,sinyp

        cosxp = DCOS(xp)
        sinxp = DSIN(xp)
        cosyp = DCOS(yp)
        sinyp = DSIN(yp)

        pm(1,1) =  cosxp
        pm(1,2) =  sinxp * sinyp
        pm(1,3) =  sinxp * cosyp

        pm(2,1) =  0.0
        pm(2,2) =  cosyp
        pm(2,3) = -sinyp

        pm(3,1) = -sinxp
        pm(3,2) =  cosxp * sinyp
        pm(3,3) =  cosxp * cosyp

        ! Approximate matrix using small angle approximations
c       pm(1,1) =  1.0
c       pm(1,2) =  0.0
c       pm(1,3) =  xp

c       pm(2,1) =  0.0
c       pm(2,2) =  1.0
c       pm(2,3) = -yp

c       pm(3,1) = -xp
c       pm(3,2) =  yp
c       pm(3,3) =  1.0

      RETURN
      END   ! SUBROUTINE PolarM
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_ITRF
*
*  this function calulates the reduction formula do a vector going from the
*    earth centered inertial frame (eci, eci) to the earth fixed (ITRF)
*    frame.  the results take into account the effects of precession, nutation,
*    sidereal time, and polar motion.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    Direction   - Which set of vars to output    FROM  TOO
*    TTT         - Julian Centuries of TT         centuries
*    JDUT1       - Julian Date of UT1             days from 4713 BC
*    LOD         - Excess length of day           sec
*    xp          - Polar motion coefficient       rad
*    yp          - Polar motion coefficient       rad
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    rITRF       - Position vector of date
*                    True Equator, True Equinox
*    vITRF       - Velocity vector of date
*                    True Equator, True Equinox
*
*  Locals        :
*    Prec        - Precession transformation (eci-mod)
*    Nut         - Nutation transformation (mod-tod)
*    St          - Sidereal Time transformation (tod-pef)
*    StDot       - Sidereal Time rate transformation (tod-pef)
*    PM          - Polar motion transformation (pef-ecef)
*    tmt         - transformation matrix transpose
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*    Omegaearth  - rotation of the earth          rad
*
*  Coupling      :
*    PRECESSION  - Find precession matrix
*    NUTATION    - Find nutation matrix
*    SIDEREAL    - Find sidereal time matrix
*    POLARM      - Find polar motion matrix
*    ROT2        - Rotation about the second axis
*    ROT3        - Rotation about the third axis
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE GCRF_ITRF   ( rj2000,vj2000, Direction, rITRF,vITRF,
     &                        TTT, JDUT1, LOD, xp, yp, terms )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rPEF(3), vPEF(3), rITRF(3), vITRF(3), xp, yp,
     &         TTT, Jdut1, LOD
        INTEGER terms

        INCLUDE 'astreduc.cmn'

* ----------------------------  Locals  -------------------------------
        REAL*8 Deltapsi,meaneps,trueeps,omega,
     &         nut(3,3), st(3,3), stdot(3,3), pm(3,3), tm(3,3),
     &         tmt(3,3),  rj2000(3),vj2000(3),prec(3,3),
     &         omegaearth(3),wcrossr(3),tm1(3,3)

        CALL precession(ttt,  prec)

        CALL nutation( ttt,  deltapsi,trueeps,meaneps,omega,nut)

        CALL sidereal(jdut1,deltapsi,meaneps,omega,lod,  st,stdot,
     &                 OmegaEarth,terms )

        CALL polarm(xp,yp, pm)

        ! ------------------- Perform matrix mmpy ---------------------
        CALL MATMULT     ( Nut , Prec  , 3,3,3,3,3,3, tm1  )
        CALL MATMULT     ( St  , tm1    , 3,3,3,3,3,3, tm  )
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( tm  , rJ2000, 3,3,1,3,3,3, rPEF  )
            CALL MATMULT     ( pm  , rPEF  , 3,3,1,3,3,3, ritrf )

            CALL MATMULT     ( tm  , vJ2000, 3,3,1,3,3,3, vPEF  )
            CALL CROSS       ( omegaearth, rPEF, wcrossr )
            CALL SUBVEC      ( vPEF, wcrossr, vPEF )
            CALL MATMULT     ( pm  , vPEF  , 3,3,1,3,3,3, vitrf )
          ELSE
            CALL MATTRANS( tm, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt , ritrf  , 3,3,1,3,3,3, rj2000 )
            CALL MATMULT     ( tmt , vitrf  , 3,3,1,3,3,3, vj2000 )
          ENDIF

chk these accel
c        CALL cross(omegaearth,rPEF,  wcrossr)
c        aecef = pm*(st*nut*prec*aeci - stdot*rPEF - cross(omegaearth,temp) ...
c                - 2.0*cross(omegaearth,vPEF))

c        Write(20,*) 'xy,yp  ',xp*3600*57.29577,' ',yp*3600*57.29577
c        IF ( Direction .eq. 'TOO ' ) THEN
c            rITRF(1)= rPEF(1) + xp*rPEF(3)
c            rITRF(2)= rPEF(2) - yp*rPEF(3)
c            rITRF(3)= rPEF(3) - xp*rPEF(1) + yp*rPEF(2)
c            vITRF(1)= vPEF(1) + xp*vPEF(3)
c            vITRF(2)= vPEF(2) - yp*vPEF(3)
c            vITRF(3)= vPEF(3) - xp*vPEF(1) + yp*vPEF(2)
c          ELSE
c            rPEF(1)= rITRF(1) - xp*rITRF(3)
c            rPEF(2)= rITRF(2) + yp*rITRF(3)
c            rPEF(3)= rITRF(3) + xp*rITRF(1) - yp*rITRF(2)
c            vPEF(1)= vITRF(1) - xp*vITRF(3)
c            vPEF(2)= vITRF(2) + yp*vITRF(3)
c            vPEF(3)= vITRF(3) + xp*vITRF(1) - yp*vITRF(2)
c          ENDIF

      RETURN
      END   ! SUBROUTINE GCRF_ITRF
*
* ----------------------------------------------------------------------------
*
*                           PROCEDURE TRUEMEAN
*
*  this subroutine calculates the transformation matrix to teme.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    TTT         - Julian Centuries of TT
*    Order       - Number of coefficients used in nutation
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    rTm         - Position vector of date
*                    Mean Equator, True Equinox
*    vTm         - Velocity vector of date
*                    Mean Equator, True Equinox
*    DeltaPsiSp  - NUTATION ANGLE                 rad
*    TrueEpsSp   - True obliquity of the ecliptic rad
*
*  Locals        :
*    TTT2        - TTT squared
*    TTT3        - TTT cubed
*    Eps         - Mean obliquity of the ecliptic rad
*    l           -                                rad
*    ll          -                                rad
*    F           -                                rad
*    D           -                                rad
*    Omega       -                                rad
*    DeltaEps    - Change in obliquity            rad
*
*  Coupling      :
*    ROT2        - Rotation about the second axis
*    ROT3        - Rotation about the third axis
*
*  References    :
*    Vallado       2007, 236
*
* ----------------------------------------------------------------------------

      SUBROUTINE TrueMean    ( Order, TTT, Terms, tm )
        IMPLICIT NONE
        REAL*8  TTT, Nutteme(3,3)
        INTEGER Order, Terms

* ----------------------------  Locals  -------------------------------
        INTEGER ii, i
        REAL*8  DeltaEps,Tempval, TTT2, TTT3, rr, l, l1, f, d,
     &          Omega, DeltaPsi,TrueEps,Conv1,
     &          Prec(3,3),meaneps,st(3,3),eqe,tm(3,3),tm1(3,3)
        REAL*8 cospsi,sinpsi,coseps,sineps,costrueeps,sintrueeps

        INCLUDE 'astreduc.cmn'
!         INCLUDE 'astmath.cmn'
!         INCLUDE 'astconst.cmn'
        REAL*8 Small,    Rad2Deg,  Deg2Rad,  HalfPi,
     &         Pi,       TwoPi,    Infinite, Undefined

        REAL*8     rekm,     mu,     omegaearth, flat,     EESqrd, auer

        DATA rekm       /6378.137D0/
c        DATA mu         /398600.4418D0/
        DATA mu         /1.32712440018D11/       ! Sun
        DATA omegaearth /7.2921158553D-5/
        DATA flat       /0.003352810665D0/       ! f = 1.0/298.257223563
        DATA EESqrd     /0.006694379990D0/       ! 2f - f**2
        DATA AUER       /23454.79095228D0/       ! 149597870.0/6378.137

        DATA Small      /0.00000001D0/
        DATA Infinite   /999999.9D0/
        DATA Undefined  /999999.1D0/
c        DATA Halfpi     /1.57079632679489662D0/
c        DATA Pi         /3.14159265358979324D0/
c        DATA TwoPi      /6.28318530717958648D0/
c        DATA Rad2Deg    /57.2957795130823208D0/
c        DATA Deg2Rad    /0.01745329251994329D0/

c       use machine precision instead
        Pi      =  4.0D0 * DATAN(1.0D0)
        HalfPi  =  0.5D0*pi
        TwoPi   =  2.0D0*pi
        Rad2Deg = 180.0D0/pi
        Deg2Rad = pi/180.0D0

        ! ----- Determine coefficients for IAU 1980 NUTATION Theory ----
        Conv1 = pi / (3600.0D0*180.0D0)

        ! ---------------- get transformation matrices ----------------
        CALL Precession ( TTT,Prec)

       TTT2= TTT*TTT
       TTT3= TTT2*TTT
       MeanEps  = 23.439291 - 0.0130042*TTT - 0.000000164*TTT2 +
     &             0.000000504*TTT3
       MeanEps = DMOD( MeanEps,360.0D0 )
       MeanEps = MeanEps * Deg2Rad

       rr   = 360.0  !deg
       l    =  134.9629814 + (1325*rr  + 198.8673981)*TTT +
     &           0.0086972 *TTT2 + 0.00001778*TTT3
       l1   =  357.5277233 + (  99*rr  + 359.05034  )*TTT -
     &           0.00016028*TTT2 - 0.00000333*TTT3
       F    =   93.2719103 + (1342*rr  +  82.0175381)*TTT -
     &           0.0036825 *TTT2 + 0.00000306*TTT3
       D    =  297.8503631 + (1236*rr  + 307.111480 )*TTT -
     &           0.00191417*TTT2 + 0.00000528*TTT3
       Omega=  125.0445222 - (   5*rr  + 134.1362608)*TTT +
     &           0.0020708 *TTT2 + 0.00000222*TTT3
       l    = DMOD( l,360.0D0 )     * Deg2Rad
       l1   = DMOD( l1,360.0D0 )    * Deg2Rad
       F    = DMOD( F,360.0D0 )     * Deg2Rad
       D    = DMOD( D,360.0D0 )     * Deg2Rad
       Omega= DMOD( Omega,360.0D0 ) * Deg2Rad

       DeltaPsi= 0.0
       DeltaEps= 0.0
       DO ii= 1, Order ! make sure the datafile is in the correct order
          i = ii

         Tempval= IAr80(1,i)*l + IAr80(2,i)*l1 + IAr80(3,i)*F +
     &            IAr80(4,i)*D + IAr80(5,i)*Omega
         DeltaPsi= DeltaPsi + (RAr80(1,i)+RAr80(2,i)*TTT) *
     &             DSIN( TempVal )
         DeltaEps= DeltaEps + (RAr80(3,i)+RAr80(4,i)*TTT) *
     &             DCOS( TempVal )
       ENDDO

       ! --------------- Find Approx Nutation Parameters --------------
       DeltaPsi = DMOD( DeltaPsi,360.0D0 ) * Deg2Rad
       DeltaEps = DMOD( DeltaEps,360.0D0 ) * Deg2Rad
       TrueEps  = MeanEps + DeltaEps

        cospsi  = DCOS(deltapsi)
        sinpsi  = DSIN(deltapsi)
        coseps  = DCOS(meaneps)
        sineps  = DSIN(meaneps)
        costrueeps = DCOS(trueeps)
        sintrueeps = DSIN(trueeps)

        nutteme(1,1) =  cospsi
        nutteme(1,2) = -coseps * sinpsi
        nutteme(1,3) = -sineps * sinpsi

        nutteme(2,1) =  costrueeps * sinpsi
        nutteme(2,2) =  costrueeps * coseps * cospsi +
     &                   sintrueeps * sineps
        nutteme(2,3) =  costrueeps * sineps * cospsi -
     &                   sintrueeps * coseps

        nutteme(3,1) =  sintrueeps * sinpsi
        nutteme(3,2) =  sintrueeps * coseps * cospsi -
     &                   sineps * costrueeps
        nutteme(3,3) =  sintrueeps * sineps * cospsi +
     &                   costrueeps * coseps

        IF (terms.gt.0) THEN
            eqe= DeltaPsi* DCOS(MeanEps)
     &           + 0.00264D0*Conv1*DSIN(Omega)
     &           + 0.000063D0*Conv1*DSIN(2.0D0*Omega)
          ELSE
            eqe= DeltaPsi* DCOS(MeanEps)
          ENDIF

        st(1,1) =  DCOS(eqe)
        st(1,2) =  DSIN(eqe)
        st(1,3) =  0.0

        st(2,1) = -DSIN(eqe)
        st(2,2) =  DCOS(eqe)
        st(2,3) =  0.0

        st(3,1) =  0.0
        st(3,2) =  0.0
        st(3,3) =  1.0

        CALL MATMULT     ( nutteme, prec , 3,3,3,3,3,3, tm1 )
        CALL MATMULT     ( st , tm1      , 3,3,3,3,3,3, tm )

      RETURN
      END   ! Subroutine TrueMean
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_TEME
*
*  this function transforms a vector from the mean equator, mean equinox of date
*    frame (j2000), to the true equator mean equinox frame (teme).
*
*  Author        : David Vallado                  719-573-2600   26 Sep 2002
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    Direction   - Which set of vars to output    FROM  TOO
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    Omega       -                                rad
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    rTm         - Position vector of date
*                    True Equator, Mean Equinox
*    vTm         - Velocity vector of date
*                    True Equator, Mean Equinox
*
*  Locals        :
*    GST         - Mean Greenwich SIDEREAL Time   0 to 2Pi rad
*    AST         - Apparent GST                   0 to 2Pi rad
*    Hr          - hour                           hr
*    minute         - minutes                        minute
*    SEC         - seconds                        SEC
*    Temp        - Temporary vector
*    TempVal     - Temporary variable
*
*  Coupling      :
*
*
*  References    :
*    Vallado       2007, 236
*
* ----------------------------------------------------------------------------

      SUBROUTINE GCRF_TEME    ( rj2000, vj2000, Direction, rTM, vTM,
     &                         Order, TTT, Terms )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rj2000(3), vj2000(3), rTM(3), vTM(3), TTT
        INTEGER Order, Terms

        INCLUDE 'astreduc.cmn'

* ----------------------------  Locals  -------------------------------
        REAL*8  nutteme(3,3), tmt(3,3)

        CALL TrueMean ( Order, TTT, Terms, Nutteme )

        ! ------------------- Perform matrix mmpy ---------------------
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( nutteme , rJ2000, 3,3,1,3,3,3, rtm )
            CALL MATMULT     ( nutteme , vJ2000, 3,3,1,3,3,3, vtm )
          ELSE
            CALL MATTRANS( nutteme, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt , rtm  , 3,3,1,3,3,3, rj2000 )
            CALL MATMULT     ( tmt , vtm  , 3,3,1,3,3,3, vj2000 )
          ENDIF


c       IF (Direction .eq. 'TOO ') THEN
c           rtm(1)= rmod(1) - DeltaPsi*DSIN(eps)*rmod(3)
c           rtm(2)= rmod(2) - DeltaEps*rmod(3)
c           rtm(3)= rmod(3) +DeltaPsi*DSIN(eps)*rmod(1) +
c     &             Deltaeps*rmod(2)
c           vtm(1)= vmod(1) - DeltaPsi*DSIN(eps)*vmod(3)
c           vtm(2)= vmod(2) - DeltaEps*vmod(3)
c           vtm(3)= vmod(3) +DeltaPsi*DSIN(eps)*vmod(1) +
c     &             Deltaeps*vmod(2)
c         ELSE
c           rmod(1)= rtm(1) + DeltaPsi*DSIN(eps)*rtm(3)
c           rmod(2)= rtm(2) + Deltaeps*rtm(3)
c           rmod(3)= rtm(3) - DeltaPsi*DSIN(eps)*rtm(1) -
c     &              Deltaeps*rtm(2)
c           vmod(1)= vtm(1) + DeltaPsi*DSIN(eps)*vtm(3)
c           vmod(2)= vtm(2) + Deltaeps*vtm(3)
c           vmod(3)= vtm(3) - DeltaPsi*DSIN(eps)*vtm(1) -
c     &              Deltaeps*vtm(2)
c         ENDIF

      RETURN
      END   ! Subroutine GCRF_TEME

