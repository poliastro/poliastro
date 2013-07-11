*      ----------------------------------------------------------------
*
*                              ASTINTER.FOR
*
*   this file contains fundamental astrodynamic procedures and functions
*   relating to interplanetary calculations.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                   2007
*                             by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              30 may 07  david vallado
*                           3rd edition baseline
*    changes :
*              21 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
*
*     *****************************************************************
*
*     Uses object files:
*         Astmath,
*         asttime,
*         Ast2body
*     Uses common files:
*         Astmath.cmn
*
*
*       SUBROUTINE PlanetRV ( PlanetNum,WhichEpoch,WhichCoord,JD, R,V )
*
*       SUBROUTINE Interplanetary ( R1,R2,Rbo,Rimpact,Mu1,Mut,Mu2,
*     &                             Deltav1,Deltav2,Vbo,Vretro
*

*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE PLANETRV
*
*  this subroutine calculate the planetary ephemerides using the Epoch J2000.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    PlanetNum   - Number of planet               1..9D0
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    R           - XYZ position vector            AU
*    V           - XYZ velocity vector            AU / day
*
*  Locals        :
*    ArgLat      -
*    TrueLon     -
*    LonPer      -
*    TUUT1       - Julian cenuries from Jan 1, 2000
*    TU2         - Tu squared
*    TU3         - TU Cubed
*    N           -
*    obliquity   - angle between ecliptic .and.
*                    Earth equator                rad
*    a           - Semi .or. axis
*    ecc         - eccentricity
*    p           - semi-parameter
*    incl        - inclination
*    omega       - ecliptic long of ascending node
*    argp        - ecliptic arg of perihelion
*    nu          - ecliptic true anomaly
*    m           - ecliptic mean anomaly
*    LLong       - True longitude
*    LongP       - longitude of perihelion
*    e0          -
*
*  Coupling      :
*    LnCom1      -
*    NewtonM     -
*    RandV       -
*
*  References    :
*    Vallado       2007, 995
*
* ----------------------------------------------------------------------------

      SUBROUTINE PlanetRV ( PlanetNum,WhichEpoch,WhichCoord,JD, R,V )
        IMPLICIT NONE
        REAL*8 JD, R(4), V(4)
        INTEGER PlanetNum
        CHARACTER*3 WhichEpoch, WhichCoord

        INCLUDE 'astmath.cmn'

* ----------------------------  Locals  -------------------------------
        REAL*8 llong,longp,m,e0, TUDaySun,ArgLat,
     &         TrueLon,LonPer,Tut1,tut12,Tut13,eps, a,ecc,p,
     &         incl,omega,argp,nu

       ! ---------------- Approximate TTDB with TUT1 -----------------
        Tut1 = ( JD - 2451545.0D0 ) / 36525.0D0
        Tut12= Tut1*Tut1
        Tut13= Tut12*Tut1

        IF ( WhichEpoch .eq. 'J20' ) THEN
           ! --------- Mean equinox of J2000 in radians (XYZ) --------
            IF (PlanetNum.eq.1) THEN
                   a    =  0.387098310D0
                   ecc  =  0.20563175D0 +    0.000020406D0*Tut1
     &                  - 0.0000000284D0*Tut12 - 0.00000000017D0*Tut13
                   incl =  0.12226007D0 -    0.000103875D0*TUT1
     &                  + 0.0000000141D0*TUT12 + 0.00000000072D0*TUT13
                   longp=  1.35186430D0 +    0.002772705D0*TUT1
     &                  - 0.0000002344D0*TUT12 + 0.00000000068D0*TUT13
                   Omega=  0.84353321D0 -    0.002189043D0*TUT1
     &                  - 0.0000015416D0*TUT12 - 0.00000000342D0*TUT13
                   LLong=  4.40260885D0 + 2608.790314157D0*TUT1
     &                  - 0.0000000934D0*TUT12 + 0.00000000003D0*TUT13
              ELSEIF (PlanetNum.eq.2) THEN
                   a    =  0.723329820D0
                   ecc  =  0.00677188D0 -    0.000047766D0*Tut1
     &                  + 0.0000000975D0*Tut12 + 0.00000000044D0*TUT13
                   incl =  0.05924803D0 -    0.000014954D0*TUT1
     &                  - 0.0000005662D0*TUT12 + 0.00000000017D0*TUT13
                   longp=  2.29621986D0 +    0.000084903D0*TUT1
     &                  - 0.0000241260D0*TUT12 - 0.00000009306D0*TUT13
                   Omega=  1.33831707D0 -    0.004852155D0*TUT1
     &                  - 0.0000024881D0*TUT12 - 0.00000000346D0*TUT13
                   LLong=  3.17614670D0 + 1021.328554621D0*TUT1
     &                  + 0.0000000288D0*TUT12 - 0.00000000003D0*TUT13
              ELSEIF (PlanetNum.eq.3) THEN
                   a    =  1.000001018D0
                   ecc  =  0.01670862D0 -    0.000042037D0*Tut1
     &                  - 0.0000001236D0*TUT12 + 0.00000000004D0*TUT13
                   incl =  0.00000000D0 +    0.000227846D0*TUT1
     &                  - 0.0000001625D0*TUT12 - 0.00000000059D0*TUT13
                   LongP=  1.79659565D0 +    0.005629659D0*TUT1
     &                  + 0.0000026225D0*TUT12 + 0.00000000834D0*TUT13
                   Omega=  0.0D0 
                   LLong=  1.75347032D0 +  628.307584919D0*TUT1
     &                  - 0.0000000991D0*TUT12 + 0.00000000000D0*TUT13
              ELSEIF (PlanetNum.eq.4) THEN
                   a    =  1.523679342D0
                   ecc  =  0.09340062D0 +    0.000090483D0*Tut1
     &                  - 0.0000000806D0*Tut12 - 0.00000000035D0*TUT13
                   incl =  0.03228381D0 -    0.000142208D0*TUT1
     &                  - 0.0000003936D0*TUT12 - 0.00000000047D0*TUT13
                   longp=  5.86535757D0 +    0.007747339D0*TUT1
     &                  - 0.0000030231D0*TUT12 + 0.00000000524D0*TUT13
                   Omega=  0.86495189D0 -    0.005148453D0*TUT1
     &                  - 0.0000111689D0*TUT12 - 0.00000003740D0*TUT13
                   LLong=  6.20348092D0 +  334.061243148D0*TUT1
     &                  + 0.0000000456D0*TUT12 - 0.00000000005D0*TUT13
              ELSEIF (PlanetNum.eq.5) THEN
                   a    =  5.202603191D0 +   0.0000001913D0*TUT1
                   ecc  =  0.04849485D0 +    0.000163244D0*Tut1
     &                  - 0.0000004719D0*Tut12 - 0.00000000197D0*TUT13
                   incl =  0.02274635D0 -    0.000034683D0*TUT1
     &                  + 0.0000005791D0*TUT12 + 0.00000000161D0*TUT13
                   longp=  0.25012853D0 +    0.003762101D0*TUT1
     &                  + 0.0000126104D0*TUT12 - 0.00000008011D0*TUT13
                   Omega=  1.75343528D0 +    0.003083697D0*TUT1
     &                  + 0.0000157755D0*TUT12 - 0.00000012273D0*TUT13
                   LLong=  0.59954650D0 +   52.969096509D0*TUT1
     &                  - 0.0000014837D0*TUT12 + 0.00000000007D0*TUT13
              ELSEIF (PlanetNum.eq.6) THEN
                   a    =  9.554909596D0 -   0.0000021389D0*TUT1
                   ecc  =  0.05550862D0 -    0.000346818D0*Tut1
     &                  - 0.0000006456D0*Tut12 + 0.00000000338D0*TUT13
                   incl =  0.04343912D0 +    0.000044532D0*TUT1
     &                  - 0.0000008557D0*TUT12 + 0.00000000031D0*TUT13
                   longp=  1.62414732D0 +    0.009888156D0*TUT1
     &                  + 0.0000092169D0*TUT12 + 0.00000008521D0*TUT13
                   Omega=  1.98383764D0 -    0.004479648D0*TUT1
     &                  - 0.0000032018D0*TUT12 + 0.00000000623D0*TUT13
                   LLong=  0.87401675D0 +   21.329909545D0*TUT1
     &                  + 0.0000036659D0*TUT12 - 0.00000000033D0*TUT13
              ELSEIF (PlanetNum.eq.7) THEN
                   a    = 19.218446062D0 -   0.0000000372D0*TUT1
     &                  + 0.00000000098D0*TUT12
                   ecc  =  0.04629590D0 -    0.000027337D0*Tut1
     &                  + 0.0000000790D0*Tut12 + 0.00000000025D0*TUT13
                   incl =  0.01349482D0 -    0.000029442D0*TUT1
     &                  + 0.0000000609D0*TUT12 + 0.00000000028D0*TUT13
                   longp=  3.01950965D0 +    0.001558939D0*TUT1
     &                  - 0.0000016528D0*TUT12 + 0.00000000721D0*TUT13
                   Omega=  1.29164744D0 +    0.001294094D0*TUT1
     &                  + 0.0000070756D0*TUT12 + 0.00000000182D0*TUT13
                   LLong=  5.48129387D0 +    7.478159856D0*TUT1
     &                  - 0.0000000848D0*TUT12 + 0.00000000010D0*TUT13
              ELSEIF (PlanetNum.eq.8) THEN
                   a    = 30.110386869D0 -   0.0000001663D0*TUT1
     &                  + 0.00000000069D0*TUT12
                   ecc  =  0.00898809D0 +    0.000006408D0*Tut1
     &                  - 0.0000000008D0*TUT12
                   incl =  0.03089149D0 +    0.000003939D0*TUT1
     &                  + 0.0000000040D0*TUT12 - 0.00000000000D0*TUT13
                   longp=  0.83991686D0 +    0.000508915D0*TUT1
     &                  + 0.0000012306D0*TUT12 - 0.00000000000D0*TUT13
                   Omega=  2.30006570D0 -    0.000107601D0*TUT1
     &                  - 0.0000000382D0*TUT12 - 0.00000000136D0*TUT13
                   LLong=  5.31188628D0 +    3.813303564D0*TUT1
     &                  + 0.0000000103D0*TUT12 - 0.00000000003D0*TUT13
              ELSEIF (PlanetNum.eq.9) THEN
                   a    = 39.53758D0
                   ecc  =  0.250877D0 
                   incl =  0.2990156D0 
                   LongP=  3.920268D0 
                   Omega=  1.926957D0
                   LLong=  3.8203049D0
              ENDIF
          ENDIF

        IF ( WhichEpoch .eq. 'ODA' ) THEN
           ! ----------- Mean equinox of date in radians (XYZ) -----------
            IF (PlanetNum.eq.1) THEN
                   a    =  0.387098310D0
                   ecc  =  0.20563175D0 +    0.000020406D0*Tut1
     &                  - 0.0000000284D0*Tut12 - 0.00000000017D0*Tut13
                   incl =  0.12226007D0 +    0.000031791D0*TUT1
     &                  - 0.0000003157D0*TUT12 + 0.00000000093D0*TUT13
                   longp=  1.35186430D0 +    0.027165657D0*TUT1
     &                  + 0.0000051643D0*TUT12 + 0.00000000098D0*TUT13
                   Omega=  0.84353321D0 +    0.020702904D0*TUT1
     &                  + 0.0000030695D0*TUT12 + 0.00000000368D0*TUT13
                   LLong=  4.40260885D0 + 2608.814707111D0*TUT1
     &                  + 0.0000053053D0*TUT12 + 0.00000000031D0*TUT13
              ELSEIF (PlanetNum.eq.2) THEN
                   a    =  0.723329820D0
                   ecc  =  0.00677188D0 -    0.000047766D0*Tut1
     &                  + 0.0000000975D0*TUT12 + 0.00000000044D0*TUT13
                   incl =  0.05924803D0 +    0.000017518D0*TUT1
     &                  - 0.0000000154D0*TUT12 - 0.00000000012D0*TUT13
                   longp=  2.29621986D0 +    0.024473335D0*TUT1
     &                  - 0.0000187338D0*TUT12 - 0.00000009276D0*TUT13
                   Omega=  1.33831707D0 +    0.015727494D0*TUT1
     &                  + 0.0000070974D0*TUT12 - 0.00000000140D0*TUT13
                   LLong=  3.17614670D0 + 1021.352943053D0*TUT1
     &                  + 0.0000054210D0*TUT12 + 0.00000000026D0*TUT13
              ELSEIF (PlanetNum.eq.3) THEN
                   a    =  1.000001018D0
                   ecc  =  0.01670862D0 -    0.000042037D0*Tut1
     &                  - 0.0000001236D0*TUT12 + 0.00000000004D0*TUT13
                   incl =  0.000D0 
                   Omega=  0.0D0 
                   Longp=  1.79659565D0 +    0.030011406D0*TUT1
     &                  + 0.0000080219D0*TUT12 + 0.00000000871D0*TUT13
                   LLong=  1.75347032D0 +  628.331966666D0*TUT1
     &                  + 0.0000053002D0*TUT12 + 0.00000000037D0*TUT13
              ELSEIF (PlanetNum.eq.4) THEN
                   a    =  1.523679342D0
                   ecc  =  0.09340062D0 +    0.000090483D0*Tut1
     &                  - 0.0000000806D0*Tut12 - 0.00000000035D0*TUT13
                   incl =  0.03228381D0 -    0.000010489D0*TUT1
     &                  + 0.0000002227D0*TUT12 - 0.00000000010D0*TUT13
                   longp=  5.86535757D0 +    0.032132089D0*TUT1
     &                  + 0.0000023588D0*TUT12 + 0.00000000555D0*TUT13
                   Omega=  0.86495189D0 +    0.013475553D0*TUT1
     &                  + 0.0000002801D0*TUT12 + 0.00000004058D0*TUT13
                   LLong=  6.20348092D0 +  334.085627899D0*TUT1
     &                  + 0.0000054275D0*TUT12 + 0.00000000026D0*TUT13
              ELSEIF (PlanetNum.eq.5) THEN
                   a    =  5.202603191D0 +   0.0000001913D0*TUT1
                   ecc  =  0.04849485D0 +    0.000163244D0*Tut1
     &                  - 0.0000004719D0*Tut12 - 0.00000000197D0*TUT13
                   incl =  0.02274635D0 -    0.000095934D0*TUT1
     &                  + 0.0000000812D0*TUT12 - 0.00000000007D0*TUT13
                   longp=  0.25012853D0 +    0.028146345D0*TUT1
     &                  + 0.0000179991D0*TUT12 - 0.00000007974D0*TUT13
                   Omega=  1.75343528D0 +    0.017819026D0*TUT1
     &                  + 0.0000070017D0*TUT12 + 0.00000000993D0*TUT13
                   LLong=  0.59954650D0 +   52.993480754D0*TUT1
     &                  + 0.0000039050D0*TUT12 + 0.00000000044D0*TUT13
              ELSEIF (PlanetNum.eq.6) THEN
                   a    =  9.554909596D0 -   0.0000021389D0*Tut1
                   ecc  =  0.05550862D0 -    0.000346818D0*Tut1
     &                  - 0.0000006456D0*Tut12 + 0.00000000338D0*TUT13
                   incl =  0.04343912D0 -    0.000065211D0*TUT1
     &                  - 0.0000002646D0*TUT12 + 0.00000000155D0*TUT13
                   longp=  1.62414732D0 +    0.034274242D0*TUT1
     &                  + 0.0000146184D0*TUT12 + 0.00000008550D0*TUT13
                   Omega=  1.98383764D0 +    0.015308246D0*TUT1
     &                  - 0.0000021061D0*TUT12 - 0.00000004154D0*TUT13
                   LLong=  0.87401675D0 +   21.354295630D0*TUT1
     &                  + 0.0000090673D0*TUT12 - 0.00000000005D0*TUT13
              ELSEIF (PlanetNum.eq.7) THEN
                   a    = 19.218446062D0 -   0.0000000372D0*Tut1
     &                   + 0.00000000098D0*TUT12
                   ecc  =  0.04629590D0 -    0.000027337D0*Tut1
     &                  + 0.0000000790D0*Tut12 + 0.00000000025D0*TUT13
                   incl =  0.01349482D0 +    0.000013516D0*TUT1
     &                  + 0.0000006543D0*TUT12 - 0.00000000161D0*TUT13
                   longp=  3.01950965D0 +    0.025942197D0*TUT1
     &                  + 0.0000037437D0*TUT12 + 0.00000000756D0*TUT13
                   Omega=  1.29164744D0 +    0.009095361D0*TUT1
     &                  + 0.0000233843D0*TUT12 + 0.00000032317D0*TUT13
                   LLong=  5.48129387D0 +    7.502543115D0*TUT1
     &                  + 0.0000053117D0*TUT12 + 0.00000000045D0*TUT13
              ELSEIF (PlanetNum.eq.8) THEN
                   a    = 30.110386869D0 -   0.0000001663D0*Tut1
     &                  + 0.00000000069D0*TUT12
                   ecc  =  0.00898809D0 +    0.000006408D0*Tut1
     &                  - 0.0000000008D0*TUT12
                   incl =  0.03089149D0 -    0.000162459D0*TUT1
     &                  - 0.0000001236D0*TUT12 + 0.00000000049D0*TUT13
                   longp=  0.83991686D0 +    0.024893067D0*TUT1
     &                  + 0.0000066179D0*TUT12 - 0.00000000005D0*TUT13
                   Omega=  2.30006570D0 +    0.019237118D0*TUT1
     &                  + 0.0000045389D0*TUT12 - 0.00000001110D0*TUT13
                   LLong=  5.31188628D0 +    3.837687716D0*TUT1
     &                  + 0.0000053976D0*TUT12 + 0.00000000031D0*TUT13
              ELSEIF (PlanetNum.eq.9) THEN
                   a    = 39.53758D0
                   ecc  =  0.250877D0 
                   incl =  0.2990156D0 
                   LongP=  3.920268D0
                   Omega=  1.926957D0 
                   LLong=  3.8203049D0 
                 ENDIF 
          ENDIF 

        LLong= DMOD( LLong,TwoPI ) 
        LongP= DMOD( LongP,TwoPI ) 
        Omega= DMOD( Omega,TwoPI ) 
        Argp = LongP - Omega 
        M    = LLong - LongP 

        CALL NewTonM( ecc,M,  E0,Nu )

       ! ------------ Find Heliocentric ecliptic r .and. v -------------  
        p      = a*(1.0D0-ecc*ecc)
        ArgLat = Argp+Nu
        TrueLon= Omega+Argp+Nu
        LonPer = Omega+Argp+Pi
        CALL COE2RV(P,ecc,Incl,Omega,Argp,Nu,ArgLat,TrueLon,LonPer,R,V)

       ! --- Correct the velocity because we used TTdb - days!! ------
        TUDaySun= 1.0D0/58.1324409D0   ! 1.0D0 / days per sun TU
        CALL LNCOM1( TUDaySun,V,  V )

        IF ( WhichCoord .eq. 'GEO' ) THEN
           ! ----------- Find obliquity of the ecliptic angle --------
            Eps = 23.439291D0 - 0.0130042D0*TUT1 - 0.000000164D0*TUT12
     &              + 0.000000504D0*TUT13
            Eps = DMOD( Eps,360.0D0 ) 
            Eps = Eps * Deg2Rad

           ! ------------- Rotate to Geocentric coordinates ----------  
            CALL ROT1( R ,-eps, R )
            CALL ROT1( V ,-eps, V )
          ENDIF

       RETURN
       END
*
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE INTERPLANETARY
*
*  this subroutine calculates the delta v's for an interplanetary mission.
*    The transfer assumes circular orbits for each of the planets.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - Radius of planet 1 from sun    km
*    R2          - Radius of planet 2 from sun    km
*    Rbo         - Radius at burnout at planet 1  km
*    Rimpact     - Radius at impact on planet 2   km
*    Mu1         - Grav parameter of planet 1     km3/s2
*    Mut         - Grav parameter of planet Sun   km3/s2
*    Mu2         - Grav parameter of planet 2     km3/s2
*
*  OutPuts       :
*    DeltaV1     - Hyperb Exc vel at planet 1 SOI km/s
*    DeltaV2     - Hyperb Exc vel at planet 2 SOI km/s
*    Vbo         - Burnout vel at planet 1        km/s
*    Vretro      - Retro vel at surface planet 2  km/s
*
*  Locals        :
*    SME1        - Spec Mech Energy of 1st orbit  Km2/s
*    SMEt        - Spec Mech Energy of trans orbitKm2/s
*    SME2        - Spec Mech Energy of 2nd orbit  Km2/s
*    Vcs1        - Vel of 1st orbit at dv 1 point Km/s
*    Vcs2        - Vel of 2nd orbit at dv 2 point Km/s
*    Vt1         - Vel of Trans orbit at dv 1 pnt Km/s
*    Vt2         - Vel of Trans orbit at dv 2 pnt Km/s
*    A           - Semimajor Axis of Trans orbit  Km
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 47-48
*
* ----------------------------------------------------------------------------  

       SUBROUTINE Interplanetary ( R1,R2,Rbo,Rimpact,Mu1,Mut,Mu2,
     &                             Deltav1,Deltav2,Vbo,Vretro )
        IMPLICIT NONE
        REAL*8 R1,R2,Rbo,Rimpact,Mu1,Mut,Mu2,Deltav1,Deltav2,vbo,
     &         vretro

        INCLUDE 'astmath.cmn'

* ----------------------------  Locals  -------------------------------
        REAL*8  SME1,SME2,SMEt, Vcs1, Vcs2, Vt1, Vt2, A

       !  Find a, SME, apogee .and. perigee velocities of trans orbit --
        A   = (R1+R2) * 0.5D0
        SMEt= -Mut/ (2.0D0*A)
        Vt1 = DSQRT( 2.0D0*( (Mut/R1) + SMEt ) )
        Vt2 = DSQRT( 2.0D0*( (Mut/R2) + SMEt ) )

       ! --- Find circular velocities of launch .and. target planet ----
        Vcs1= DSQRT( Mut/R1 )
        Vcs2= DSQRT( Mut/R2 )

       ! --- Find delta velocities DO Hohmann transfer portion  -----
        DeltaV1= DABS( Vt1 - Vcs1 ) 
        DeltaV2= DABS( Vcs2 - Vt2 ) 

       !  Find SME .and. burnout/impact vel of launch / target planets -  
        SME1  = Deltav1*DeltaV1 * 0.5D0 
        SME2  = Deltav2*DeltaV2 * 0.5D0 
        Vbo   = DSQRT( 2.0D0*( (Mu1/Rbo) + SME1 ) ) 
        Vretro= DSQRT( 2.0D0*( (Mu2/Rimpact) + SME2 ) ) 

c        TP= Pi*DSQRT( a*a*a/Mut )   ! Transfer Period in secs

       RETURN
       END
*

