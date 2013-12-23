*   -------------------------------------------------------------------
*
*                              ASTIOD.FOR
*
*   this file contains fundamental astrodynamic procedures and functions
*   relating to the initial orbit determination techniques. see ch 7 for
*   a complete discussion of these routines.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                   2007
*                             by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              15 mar 07  david vallado
*                           3rd edition baseline
*    changes :
*              21 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
*
*     Uses object files:
*         Astmath
*         Asttime
*         Ast2body
*         Astreduc
*     Uses common files:
*         Astmath.cmn
*
*
*      SUBROUTINE SITE        ( Latgd,Alt,Lon, RSecef,VSecef )
*
*      ------------------ Angles-only techniques ----------------------
*
*      SUBROUTINE ANGLESLAPLACE( Delta1,Delta2,Delta3,Alpha1,Alpha2,
*     &                          Alpha3,JD1,JD2,JD3,RS1,RS2,RS3, r2,v2)
*
*      SUBROUTINE ANGLESGAUSS ( Delta1,Delta2,Delta3,Alpha1,Alpha2,
*     &                         Alpha3,JD1,JD2,JD3,RS1,RS2,RS3, r2,v2 )
*
*      ------------------- Conversion techniques ----------------------
*
*      SUBROUTINE RV_RADEC    ( Rijk,Vijk, Direction, rr,RtAsc,Decl,
*     &                         DRr,DRtAsc,DDecl )
*
*      SUBROUTINE RV_TRADEC   ( Rijk,Vijk,RSecef, Direction, Rho,TRtAsc,
*     &                         TDecl,DRho,DTRtAsc,DTDecl )
*
*      SUBROUTINE RV_RAZEL    ( Reci,Veci,Latgd,Lon,alt,TTT,jdut1,lod,
*     &                         xp,yp,terms, Direction,
*     &                         Rho,Az,El,DRho,DAz,DEl )
*
*      SUBROUTINE RV_ELATLON  ( Rijk,Vijk, Direction, rr,EclLat,EclLon,
*     &                         DRr,DEclLat,DEclLon )
*
*      SUBROUTINE RVSEZ_RAZEL ( Rhosez,DRhosez,Direction, Rho,Az,El,
*     &                         DRho,DAz,DEl )
*
*      SUBROUTINE RADEC_ELATLON ( RtAsc,Decl,Direction, EclLat, EclLon)
*
*      SUBROUTINE RADEC_AZEL  ( RtAsc,Decl,LST,Latgd, Direction, Az,El)
*
*      ------------------- Three vector techniques --------------------
*
*      SUBROUTINE GIBBS       ( R1,R2,R3, V2, Theta,Theta1,Copa, Error)
*
*      SUBROUTINE HERRGIBBS   ( R1,R2,R3,JD1,JD2,JD3, V2, Theta,Theta1,
*     &                         Copa, Error )
*
*      ------------------------ Lambert techniques --------------------
*
*      SUBROUTINE LAMBERTUNIV ( ro,r, dm,OverRev, Dtsec, vo,v, Error )
*
*      SUBROUTINE LAMBERTBATTIN ( ro,r, dm,OverRev, Dtsec, vo,v, Error )
*
*      SUBROUTINE TARGET      ( RInt,VInt,RTgt,VTgt, Dm,Kind, Dtsec,
*     &                         V1t,V2t,DV1,DV2, Error  )
*
*
* ---------------------------------------------------------------------------
*
*                           SUBROUTINE SITE
*
*  this subroutine finds the position and velocity vectors for a SITE.  The
*    answer is returned in the Geocentric Equatorial (IJK) coordinate system.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
*    Alt         - Altitude                       km
*    LST         - Local SIDEREAL Time            -2Pi to 2Pi rad
*
*  OutPuts       :
*    RSecef      - ecef SITE position vector      km
*    VSecef      - ecef SITE velocity vector      km/s
*
*  Locals        :
*    EarthRate   - IJK Earth's rotation rate      rad/s
*    SinLat      - Variable containing  DSIN(Lat) rad
*    Temp        - Temporary Real value
*    Rdel        - Rdel component of SITE vector  km
*    Rk          - Rk component of SITE vector    km
*    CEarth      -
*
*  Coupling      :
*    CROSS       - CROSS product of two vectors
*
*  References    :
*    Vallado       2001, 404-407, Alg 47, Ex 7-1
*
* -----------------------------------------------------------------------------  

      SUBROUTINE SITE               ( Latgd,Alt,Lon, RSecef,VSecef )
        IMPLICIT NONE
        REAL*8 LatGd,Alt,Lon,RSecef(3),VSecef(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 SinLat, CEarth, Rdel, Rk

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        SinLat      = DSIN( Latgd ) 

        ! ------  Find Rdel and Rk components of SITE vector  ---------
        CEarth= rekm / DSQRT( 1.0D0 - ( EESqrd*SinLat*SinLat ) )
        Rdel  = ( CEarth + Alt )*DCOS( Latgd ) 
        Rk    = ( (1.0D0-EESqrd)*CEarth + Alt )*SinLat

        ! ---------------  Find SITE position vector  -----------------
        RSecef(1) = Rdel * DCOS( Lon )
        RSecef(2) = Rdel * DSIN( Lon )
        RSecef(3) = Rk

        ! ---------------  Find SITE velocity vector  ------------------
        VSecef(1) = 0.0D0
        VSecef(2) = 0.0D0
        VSecef(3) = 0.0D0
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ANGLESLAPLACE
*
*  this subroutine solves the problem of orbit determination using three
*    optical sightings and the method of Laplace.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Alpha1       - Right Ascension #1            rad
*    Alpha2       - Right Ascension #2            rad
*    Alpha3       - Right Ascension #3            rad
*    Delta1       - Declination #1                rad
*    Delta2       - Declination #2                rad
*    Delta3       - Declination #3                rad
*    JD1          - Julian Date of 1st sighting   Days from 4713 BC
*    JD2          - Julian Date of 2nd sighting   Days from 4713 BC
*    JD3          - Julian Date of 3rd sighting   Days from 4713 BC
*    RS1          - IJK SITE position vector #1   km
*    RS2          - IJK SITE position vector #2   km
*    RS3          - IJK SITE position vector #3   km
*
*  OutPuts        :
*    R            - IJK position vector           km
*    V            - IJK velocity vector           km / s
*
*  Locals         :
*    L1           - Line of SIGHT vector for 1st
*    L2           - Line of SIGHT vector for 2nd
*    L3           - Line of SIGHT vector for 3rd
*    LDot         - 1st derivative of L2
*    LDDot        - 2nd derivative of L2
*    RS2Dot       - 1st Derivative of RS2 - vel
*    RS2DDot      - 2nd Derivative of RS2
*    t12t13       - (t1-t2) * (t1-t3)
*    t21t23       - (t2-t1) * (t2-t3)
*    t31t32       - (t3-t1) * (t3-t2)
*    i            - index
*    D            -
*    D1           -
*    D2           -
*    D3           -
*    D4           -
*    OldR         - Previous iteration on r
*    Rho          - Range from SITE to satellite at t2
*    RhoDot       -
*    DMat         -
*    D1Mat        -
*    D2Mat        -
*    D3Mat        -
*    D4Mat        -
*    EarthRate    - Angular rotation of the earth
*    L2DotRS      - Vector L2 Dotted with RSecef
*    Temp         - Temporary vector
*    Temp1        - Temporary vector
*    Small        - Tolerance
*    Roots        -
*
*  Coupling       :
*    MAG          - Magnitude of a vector
*    DETERMINANT  - Evaluate the determinant of a matrix
*    CROSS        - CROSS product of two vectors
*    NORM         - Normlize a matrix
*    FACTOR       - Find the roots of a polynomial
*
*  References     :
*    Vallado       2001, 413-417
*
* ------------------------------------------------------------------------------  
*
      SUBROUTINE ANGLESLAPLACE ( Delta1,Delta2,Delta3,Alpha1,Alpha2,
     &                      Alpha3,JD1,JD2,JD3,RS1,RS2,RS3, r2,v2 )
        IMPLICIT NONE
        REAL*8 Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,JD2,JD3,
     &         RS1(3),RS2(3),RS3(3),r2(3),v2(3)
        EXTERNAL DETERMINANT, Dot, Mag
* -----------------------------  Locals  ------------------------------
        INTEGER i, j, k
        REAL*8 Small, Poly(16),Roots(15,2), MAG
        REAL*8 DMat(3,3), DMat1(3,3), DMat2(3,3), DMat3(3,3),DMat4(3,3)
        REAL*8 L1(3), L2(3), L3(3), LDot(3), LDDot(3), RS2Dot(3),
     &         RS2DDot(3), EarthRate(3), Temp(3), Temp1(3),magr2,
     &         D, D1, D2, D3, D4, Rho, RhoDot, t1t13, t1t3, t31t3,
     &         tau1, tau3, BigR2, L2DotRS, Determinant, Dot, magtemp,
     &         magtemp1, magrs2

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
c        TUDay        =     0.00933809017716D0
*        TUDay        =    58.132440906D0
        Small        =     0.0000001D0
        EarthRate(1)= 0.0D0
        EarthRate(2)= 0.0D0
        EarthRate(3)= OmegaEarth

        JD1= JD1*86400.0D0    ! days to sec
        JD2= JD2*86400.0D0
        JD3= JD3*86400.0D0

        ! ---------- set middle to 0, find deltas to others -----------
        tau1= JD1-JD2 
        tau3= JD3-JD2 

        ! --------------- Find Line of SIGHT vectors ------------------
        L1(1)= DCOS(Delta1)*DCOS(Alpha1)
        L1(2)= DCOS(Delta1)*DSIN(Alpha1)
        L1(3)= DSIN(Delta1)
        L2(1)= DCOS(Delta2)*DCOS(Alpha2)
        L2(2)= DCOS(Delta2)*DSIN(Alpha2)
        L2(3)= DSIN(Delta2)
        L3(1)= DCOS(Delta3)*DCOS(Alpha3)
        L3(2)= DCOS(Delta3)*DSIN(Alpha3)
        L3(3)= DSIN(Delta3)

        ! -------------------------------------------------------------
*       Using Lagrange Interpolation formula to derive an expression
*       for L(t), substitute t=t2 and differentiate to obtain the
*       derivatives of L.
        ! -------------------------------------------------------------
        t1t13= 1.0D0 / (tau1*(tau1-tau3)) 
        t1t3 = 1.0D0 / (tau1*tau3) 
        t31t3= 1.0D0 / ((tau3-tau1)*tau3) 
        DO i= 1 , 3
            LDot(i)=      ( -tau3 * t1t13 )*L1(i) +
     &               ( (-tau1-tau3) * t1t3  )*L2(i) +
     &                      ( -tau1 * t31t3 )*L3(i)
            LDDot(i)= ( 2.0D0 * t1t13 )*L1(i) +
     &                  ( 2.0D0 * t1t3  )*L2(i) +
     &                  ( 2.0D0 * t31t3 )*L3(i)
          ENDDO
        CALL NORM( LDot,  LDot )
        CALL NORM( LDDot, LDDot )

        ! ------------------- Find 2nd derivative of RSecef ---------------
        CALL CROSS( RS1,RS2, Temp )
        magtemp = MAG(Temp)
        CALL CROSS( RS2,RS3, Temp1 )
        magtemp1 = MAG(Temp1)
*
*      needs a different test xxxx!!  
        IF ( ( DABS(magtemp) .gt. Small ) .and.
     &     ( DABS( magtemp1) .gt. Small )  ) THEN
           ! ------------ All sightings from one SITE -----------------
*          fix this testhere  
            DO i= 1 , 3
                RS2Dot(i)=      ( -tau3 * t1t13 )*RS1(i) +
     &                     ( (-tau1-tau3) * t1t3  )*RS2(i) +
     &                            ( -tau1 * t31t3 )*RS3(i)
                RS2DDot(i)= ( 2.0D0 * t1t13 )*RS1(i) +
     &                        ( 2.0D0 * t1t3  )*RS2(i) +
     &                        ( 2.0D0 * t31t3 )*RS3(i)
              ENDDO

            CALL CROSS( EarthRate,RS2,     RS2Dot )
            CALL CROSS( EarthRate,RS2Dot,  RS2DDot )
          ELSE
            ! ---------- Each sighting from a different SITE ----------
            DO i= 1 , 3
                RS2Dot(i)=      ( -tau3 * t1t13 )*RS1(i) +
     &                     ( (-tau1-tau3) * t1t3  )*RS2(i) +
     &                            ( -tau1 * t31t3 )*RS3(i)
                RS2DDot(i)= ( 2.0D0 * t1t13 )*RS1(i) +
     &                        ( 2.0D0 * t1t3  )*RS2(i) +
     &                        ( 2.0D0 * t31t3 )*RS3(i)
              ENDDO
          ENDIF 

        DO i= 1 , 3
            DMat(i,1) =2.0D0 * L2(i)
            DMat(i,2) =2.0D0 * LDot(i)
            DMat(i,3) =2.0D0 * LDDot(i)

            ! ----------------  Position determinants -----------------
            DMat1(i,1) =L2(i)
            DMat1(i,2) =LDot(i)
            DMat1(i,3) =RS2DDot(i)
            DMat2(i,1) =L2(i)
            DMat2(i,2) =LDot(i)
            DMat2(i,3) =RS2(i)

            ! ------------  Velocity determinants ---------------------
            DMat3(i,1) =L2(i)
            DMat3(i,2) =RS2DDot(i)
            DMat3(i,3) =LDDot(i)
            DMat4(i,1) =L2(i)
            DMat4(i,2) =RS2(i)
            DMat4(i,3) =LDDot(i)
          ENDDO

        D = DETERMINANT(DMat,3) 
        D1= DETERMINANT(DMat1,3) 
        D2= DETERMINANT(DMat2,3) 
        D3= DETERMINANT(DMat3,3) 
        D4= DETERMINANT(DMat4,3) 
* 
      ! ---------------  Iterate to find Rho magnitude ----------------
*     magr= 1.5D0   ! First Guess
*     Write( 'Input initial guess for magr ' )
*     Read( magr )
*     i= 1 
*     REPEAT
*         OldR= magr
*         Rho= -2.0D0*D1/D - 2.0D0*D2/(magr**3*D)
*         magr= DSQRT( Rho*Rho + 2.0D0*Rho*L2DotRS + magRS2**2 )
*         INC(i) 
*         magr= (OldR - magr ) / 2.0D0             ! Simple bissection
*         WriteLn( FileOut,'Rho guesses ',i:2,'Rho ',Rho:14:7,' magr ',magr:14:7,oldr:14:7 )
*! seems to converge, but wrong Numbers
*         INC(i) 
*     UNTIL ( DABS( OldR-magR ) .lt. Small ) .or. ( i .ge. 30 )
   
*
        IF ( DABS(D) .gt. 0.000001D0 ) THEN
            ! --------------- Solve eighth order poly -----------------
            L2DotRS= DOT( L2,RS2 ) 
            magrs2 = MAG(rs2)
            Poly( 1)=  1.0D0  ! r2^8th variable!!!!!!!!!!!!!!
            Poly( 2)=  0.0D0
            Poly( 3)=  (L2DotRS*4.0D0*D1/D - 4.0D0*D1*d1/(D*D)
     &                 - magRS2**2 )
            Poly( 4)=  0.0D0
            Poly( 5)=  0.0D0
            Poly( 6)=  Mu*(L2DotRS*4.0D0*D2/D - 8.0D0*D1*D2/(D*D) )
            Poly( 7)=  0.0D0
            Poly( 8)=  0.0D0
            Poly( 9)=  -4.0D0*Mu*D2*D2/(D*D)
            Poly(10)=  0.0D0
            Poly(11)=  0.0D0
            Poly(12)=  0.0D0
            Poly(13)=  0.0D0
            Poly(14)=  0.0D0
            Poly(15)=  0.0D0
            Poly(16)=  0.0D0
            CALL FACTOR( Poly,8,  Roots )

           ! ------------------ Find correct (xx) root ----------------
            BigR2= 0.0D0 
            DO j= 1 , 8
*                IF ( DABS( Roots(j,2) ) .lt. Small ) THEN
*                    WriteLn( 'Root ',j,Roots(j,1),' + ',Roots(j,2),'j' )
*        temproot= roots(j,1)*roots(j,1)
*        temproot= Temproot*TempRoot*TempRoot*TempRoot +
*                  Poly(3)*TempRoot*TempRoot*TempRoot + Poly(6)*roots(j,1)*Temproot + Poly(9)
*                    WriteLn( FileOut,'Root ',j,Roots(j,1),' + ',Roots(j,2),'j  value = ',temproot )
                    IF ( Roots(j,1) .gt. BigR2 ) THEN
                        BigR2= Roots(j,1)
                      ENDIF
*                  ENDIF
              ENDDO
        Write(*,*) 'BigR2 ',BigR2
        Write(*,*) 'Keep this root ? '
        READ(*,*) BigR2

            Rho= -2.0D0*D1/D - 2.0D0*Mu*D2 / (BigR2*BigR2*BigR2*D) 

            ! --------- Find the middle position vector ---------------
            DO k= 1 , 3
                r2(k)= Rho*L2(k) + RS2(k)
              ENDDO
            magr2 = MAG( r2 )
            ! ---------------- Find RhoDot magnitude ------------------
            RhoDot= -D3/D - Mu*D4/(magr2**3*D)
*        WriteLn( FileOut,'Rho ',Rho:14:7 )
*        WriteLn( FileOut,'RhoDot ',RhoDot:14:7 )

            ! -------------- Find middle velocity vector --------------
            DO i= 1 , 3
                V2(i)= RhoDot*L2(i) + Rho*LDot(i) + RS2Dot(i)
              ENDDO
         ELSE
           Write(*,*) 'Determinant value was zero ',D
         ENDIF

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ANGLESGAUSS
*
*  this subroutine solves the problem of orbit determination using three
*    optical sightings.  The solution SUBROUTINE uses the Gaussian technique.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Alpha1       - Right Ascension #1            rad
*    Alpha2       - Right Ascension #2            rad
*    Alpha3       - Right Ascension #3            rad
*    Delta1       - Declination #1                rad
*    Delta2       - Declination #2                rad
*    Delta3       - Declination #3                rad
*    JD1          - Julian Date of 1st sighting   Days from 4713 BC
*    JD2          - Julian Date of 2nd sighting   Days from 4713 BC
*    JD3          - Julian Date of 3rd sighting   Days from 4713 BC
*    RSecef           - IJK SITE position vector      km
*
*  OutPuts        :
*    R            - IJK position vector at t2     km
*    V            - IJK velocity vector at t2     km / s
*
*  Locals         :
*    L1           - Line of SIGHT vector for 1st
*    L2           - Line of SIGHT vector for 2nd
*    L3           - Line of SIGHT vector for 3rd
*    Tau          - Taylor expansion series about
*                   Tau ( t - to )
*    TauSqr       - Tau squared
*    t21t23       - (t2-t1) * (t2-t3)
*    t31t32       - (t3-t1) * (t3-t2)
*    i            - index
*    D            -
*    Rho          - Range from SITE to sat at t2  km
*    RhoDot       -
*    DMat         -
*    RS1          - SITE vectors
*    RS2          -
*    RS3          -
*    EarthRate    - Velocity of Earth rotation
*    P            -
*    Q            -
*    OldR         -
*    OldV         -
*    F1           - F coefficient
*    G1           -
*    F3           -
*    G3           -
*    L2DotRS      -
*
*  Coupling       :
*    Detrminant   - Evaluate the determinant of a matrix
*    FACTOR       - Find roots of a polynomial
*    MATMULT      - Multiply two matrices together
*    GIBBS        - GIBBS method of orbit determination
*    HGIBBS       - Herrick GIBBS method of orbit determination
*    ANGLE        - ANGLE between two vectors
*
*  References     :
*    Vallado       2001, 417-421, Alg 49, Ex 7-2 (425-427)
*
* ------------------------------------------------------------------------------  
*
      SUBROUTINE ANGLESGAUSS ( Delta1,Delta2,Delta3,Alpha1,Alpha2,
     &                      Alpha3,JD1,JD2,JD3,RS1,RS2,RS3, r2,v2 )
        IMPLICIT NONE
        REAL*8 Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,JD2,JD3,
     &         RS1(3),RS2(3),RS3(3),r2(3),v2(3)
        EXTERNAL Determinant, Dot, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i, ll, j
        REAL*8 small, Roots(15,2), Poly(16),
     &         r1(3), r3(3), L1(3), L2(3), L3(3)
        CHARACTER*12 Error
        REAL*8 LMatIi(3,3), CMat(3,1),RhoMat(3,1), LMatI(3,3),
     &         RSMat(3,3),LIR(3,3), Determinant, Dot, magrs2
        REAL*8 rDot, tau1, tau3, u, uDot, p, MAG, magr2, magr1, magr3,
     &         f1, g1, f3, g3, a, ecc, incl, omega, argp,
     &         Nu, m, l, ArgPer, BigR2, a1, a1u, a3, a3u, d, d1,
     &         d2, c1, c3, L2DotRS, rhoold1, rhoold2, rhoold3,
     &         rad,  theta, theta1, copa, TauSqr

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
*        TUDay        =    58.132440906D0
c        TUDay        =     0.00933809017716D0
        Small        =     0.0000001D0
        rad    = 57.29577951308D0 

        ! ---------- set middle to 0, find deltas to others -----------
        tau1= (JD1-JD2)*86400.0D0
        tau3= (JD3-JD2)*86400.0D0

        ! ----------------  Find Line of SIGHT vectors  ---------------
        L1(1)= DCOS(Delta1)*DCOS(Alpha1)
        L1(2)= DCOS(Delta1)*DSIN(Alpha1)
        L1(3)= DSIN(Delta1)

        L2(1)= DCOS(Delta2)*DCOS(Alpha2)
        L2(2)= DCOS(Delta2)*DSIN(Alpha2)
        L2(3)= DSIN(Delta2)

        L3(1)= DCOS(Delta3)*DCOS(Alpha3)
        L3(2)= DCOS(Delta3)*DSIN(Alpha3)
        L3(3)= DSIN(Delta3)

        ! ------------- Find L matrix and determinant -----------------
        ! --------- Called LMatI since it is only used for determ -----

        DO i= 1 , 3
              LMatIi(i,1) =L1(i)
              LMatIi(i,2) =L2(i)
              LMatIi(i,3) =L3(i)
              RSMat(i,1) =RS1(i)
              RSMat(i,2) =RS2(i)
              RSMat(i,3) =RS3(i)
          ENDDO

        D= DETERMINANT(LMatIi,3) 
        ! ------------------ Now assign the inverse -------------------
        LMatI(1,1) = ( L2(2)*L3(3)-L2(3)*L3(2)) / D
        LMatI(2,1) = (-L1(2)*L3(3)+L1(3)*L3(2)) / D
        LMatI(3,1) = ( L1(2)*L2(3)-L1(3)*L2(2)) / D
        LMatI(1,2) = (-L2(1)*L3(3)+L2(3)*L3(1)) / D
        LMatI(2,2) = ( L1(1)*L3(3)-L1(3)*L3(1)) / D
        LMatI(3,2) = (-L1(1)*L2(3)+L1(3)*L2(1)) / D
        LMatI(1,3) = ( L2(1)*L3(2)-L2(2)*L3(1)) / D
        LMatI(2,3) = (-L1(1)*L3(2)+L1(2)*L3(1)) / D
        LMatI(3,3) = ( L1(1)*L2(2)-L1(2)*L2(1)) / D

        CALL MATMULT( LMatI,RSMat,3,3,3, 3,3,3,   LIR )
*
        ! ------------ Find f and g series at 1st and 3rd obs ---------
*      speed by assuming circ sat vel for uDot here ??
*      some similartities in 1/6t3t1 ...  
        ! --- keep separated this time ----
        a1 =  Tau3 / (Tau3 - Tau1) 
        a1u=  (Tau3*((tau3-Tau1)*(Tau3-Tau1) - Tau3*Tau3 )) /
     &        (6.0D0*(Tau3 - Tau1))
        a3 = -Tau1 / (Tau3 - Tau1) 
        a3u= -(Tau1*((tau3-Tau1)*(Tau3-Tau1) - Tau1*Tau1 )) /
     &        (6.0D0*(Tau3 - Tau1))

        ! --- Form initial guess of r2 ----
        d1=  LIR(2,1)*a1 - LIR(2,2) + LIR(2,3)*a3
        d2=  LIR(2,1)*a1u + LIR(2,3)*a3u

        ! ------- Solve eighth order poly NOT same as LAPLACE ---------
        L2DotRS= DOT( L2,RS2 ) 
        magrs2 = MAG(rs2)
        Poly( 1)=  1.0D0  ! r2^8th variable!!!!!!!!!!!!!!
        Poly( 2)=  0.0D0
        Poly( 3)=  -(D1*D1 + 2.0D0*D1*L2DotRS + magRS2**2)
        Poly( 4)=  0.0D0
        Poly( 5)=  0.0D0
        Poly( 6)=  -2.0D0*Mu*(L2DotRS*D2 + D1*D2)
        Poly( 7)=  0.0D0
        Poly( 8)=  0.0D0
        Poly( 9)=  -Mu*Mu*D2*D2
        Poly(10)=  0.0D0
        Poly(11)=  0.0D0
        Poly(12)=  0.0D0
        Poly(13)=  0.0D0
        Poly(14)=  0.0D0
        Poly(15)=  0.0D0
        Poly(16)=  0.0D0
        CALL FACTOR( Poly,8,  Roots )

        ! ------------------ Select the correct root ------------------
        BigR2= 0.0D0 
        DO j= 1 , 8
*            IF ( DABS( Roots(j,2) ) .lt. Small ) THEN
*     temproot= roots(j,1)*roots(j,1)
*     temproot= Temproot*TempRoot*TempRoot*TempRoot +
*              Poly(3)*TempRoot*TempRoot*TempRoot + Poly(6)*roots(j,1)*Temproot + Poly(9)
*                WriteLn( FileOut,'Root ',j,Roots(j,1),' + ',Roots(j,2),'j  value = ',temproot )
                IF ( Roots(j,1) .gt. BigR2 ) THEN
                    BigR2= Roots(j,1)
                ENDIF  ! IF (
*             ENDIF  ! IF (
          ENDDO
        ! ------------ Solve matrix with u2 better known --------------
        u= Mu / ( BigR2*BigR2*BigR2 ) 

        c1= a1+a1u*u 
        c3= a3+a3u*u 
          CMat(1,1)= -c1
          CMat(2,1)= 1.0D0
          CMat(3,1)= -c3
        CALL MATMULT( LIR,CMat,3,3,1, 3,3,1,  RhoMat )

        Rhoold1=  RhoMat(1,1)/c1
        Rhoold2= -RhoMat(2,1)
        Rhoold3=  RhoMat(3,1)/c3

*
      ! -------- Loop through the refining process ------------  for WHILE () DO
      DO ll= 1 , 3
            Write( *,*) ' Iteration # ',ll
            ! ---------- Now form the three position vectors ----------
            DO i= 1 , 3
                R1(i)=  RhoMat(1,1)*L1(i)/c1 + RS1(i)
                R2(i)= -RhoMat(2,1)*L2(i)    + RS2(i)
                R3(i)=  RhoMat(3,1)*L3(i)/c3 + RS3(i)
              ENDDO

            CALL GIBBS(r1,r2,r3,  v2,theta,theta1,copa,error )

            IF ( (Error .ne. 'ok') .and. (copa .lt. 1.0D0/Rad) ) THEN
                ! --- HGibbs to get middle vector ----
                CALL HERRGIBBS(r1,r2,r3,JD1,JD2,JD3,
     &                     v2,theta,theta1,copa,error )
*                WriteLn( FileOut,'hgibbs ' )
              ENDIF 

            CALL rv2coe( r2,v2, p,a,ecc,incl,omega,argp,Nu,m,u,l,ArgPer)
            magr2 = MAG(r2)

        IF ( ll .le. 2 ) THEN
            ! --- Now get an improved estimate of the f and g series --
*       .or. can the analytic functions be found now??  
            u= Mu / ( magr2**3 )
            rDot= DOT(r2,v2)/magr2
            uDot= (-3.0D0*Mu*RDot) / (magr2**4)

            TauSqr= Tau1*Tau1 
            f1=  1.0D0 - 0.5D0*u*TauSqr -(1.0D0/6.0D0)*UDot*TauSqr*Tau1
     &                 + (1.0D0/24.0D0) * u*u*TauSqr*TauSqr
     &                 + (1.0D0/30.0D0)*U*UDot*TauSqr*TauSqr*Tau1
            g1= Tau1 - (1.0D0/6.0D0)*u*Tau1*TauSqr - (1.0D0/12.0D0) *
     &                 UDot*TauSqr*TauSqr
     &                 + (1.0D0/120.0D0)*u*u*TauSqr*TauSqr*Tau1
     &                 + (1.0D0/120.0D0)*u*UDot*TauSqr*TauSqr*TauSqr
            TauSqr= Tau3*Tau3 
            f3=  1.0D0 - 0.5D0*u*TauSqr -(1.0D0/6.0D0)*UDot*TauSqr*Tau3
     &                 + (1.0D0/24.0D0) * u*u*TauSqr*TauSqr
     &                 + (1.0D0/30.0D0)*U*UDot*TauSqr*TauSqr*Tau3
            g3= Tau3 - (1.0D0/6.0D0)*u*Tau3*TauSqr - (1.0D0/12.0D0) *
     &                 UDot*TauSqr*TauSqr
     &                 + (1.0D0/120.0D0)*u*u*TauSqr*TauSqr*Tau3
     &                 + (1.0D0/120.0D0)*u*UDot*TauSqr*TauSqr*TauSqr
          ELSE
            ! -------- Now use exact method to find f and g -----------
            CALL ANGLE( R1,R2, Theta )
            CALL ANGLE( R2,R3, Theta1 )
            magr1 = MAG(r1)
            magr3 = MAG(r3)

            f1= 1.0D0 - ( (magR1*(1.0D0 - DCOS(Theta)) / p ) )
            g1= ( magR1*magR2*DSIN(-theta) ) / DSQRT( p )  ! - ANGLE because backwards!!
            f3= 1.0D0 - ( (magR3*(1.0D0 - DCOS(Theta1)) / p ) )
            g3= ( magR3*magR2*DSIN(theta1) ) / DSQRT( p )

         ENDIF
            c1=  g3 / (f1*g3 - f3*g1)
            c3= -g1 / (f1*g3 - f3*g1) 
            ! ----- Solve for all three ranges via matrix equation ----
            CMat(1,1)= -c1
            CMat(2,1)= 1.0D0
            CMat(3,1)= -c3
            CALL MATMULT( LIR,CMat,3,3,1, 3,3,1,  RhoMat )

            ! ----------------- Check for convergence -----------------

          ENDDO   ! DO WHILE the ranges are still changing

        ! ---------------- Find all three vectors ri ------------------
        DO i= 1 , 3
            R1(i)=  RhoMat(1,1)*L1(i)/c1 + RS1(i)
            R2(i)= -RhoMat(2,1)*L2(i)    + RS2(i)
            R3(i)=  RhoMat(3,1)*L3(i)/c3 + RS3(i)
          ENDDO
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RV_RADEC
*
*  this subroutine converts the right ascension and declination values with
*    position and velocity vectors of a satellite. Uses velocity vector to
*    find the solution of singular cases.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rijk        - IJK position vector            km
*    Vijk        - IJK velocity vector            km/s
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    Rr          - Radius of the satellite        km
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*    DRr         - Radius of the satellite rate   km/s
*    DRtAsc      - Right Ascension rate           rad/s
*    DDecl       - Declination rate               rad/s
*
*  Locals        :
*    Temp        - Temporary position vector
*    Temp1       - Temporary variable
*
*  Coupling      :
*    DOT         - DOT product of two vectors
*
*  References    :
*    Vallado       2001, 246-248, Alg 25
*
* ------------------------------------------------------------------------------  

      SUBROUTINE RV_RADEC    ( Rijk,Vijk, Direction, rr,RtAsc,Decl,
     &                         DRr,DRtAsc,DDecl )
        IMPLICIT NONE
        REAL*8 Rijk(3),Vijk(3),rr,RtAsc,Decl,DRr,DRtAsc,DDecl
        CHARACTER*4 Direction
        EXTERNAL Dot, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER Small
        REAL*8 Temp, Temp1, Dot, MAG

        ! --------------------  Implementation   ----------------------
        Small        = 0.00000001D0
        IF ( Direction .eq. 'FROM' ) THEN
            Rijk(1)= rr*DCOS(Decl)*DCOS(RtAsc)
            Rijk(2)= rr*DCOS(Decl)*DSIN(RtAsc)
            Rijk(3)= rr*DSIN(Decl)
            Vijk(1)= DRr*DCOS(Decl)*DCOS(RtAsc) -
     &               rr*DSIN(Decl)*DCOS(RtAsc)*DDecl
     &               - rr*DCOS(Decl)*DSIN(RtAsc)*DRtAsc
            Vijk(2)= DRr*DCOS(Decl)*DSIN(RtAsc) -
     &               rr*DSIN(Decl)*DSIN(RtAsc)*DDecl
     &               + rr*DCOS(Decl)*DCOS(RtAsc)*DRtAsc
            Vijk(3)= DRr*DSIN(Decl) + rr*DCOS(Decl)*DDecl
          ELSE
            ! ------------- Calculate Angles and Rates ----------------
            rr = MAG(Rijk)
            Temp= DSQRT( Rijk(1)*Rijk(1) + Rijk(2)*Rijk(2) )
            IF ( Temp .lt. Small ) THEN
                RtAsc= DATAN2( Vijk(2), Vijk(1) )
              ELSE
                RtAsc= DATAN2( Rijk(2), Rijk(1) )
              ENDIF
            Decl= DASIN( Rijk(3)/rr )

            Temp1= -Rijk(2)*Rijk(2) - Rijk(1)*Rijk(1)  ! different now
            DRr= DOT(Rijk,Vijk)/rr 
            IF ( DABS(Temp1) .gt. Small ) THEN
                DRtAsc= ( Vijk(1)*Rijk(2) - Vijk(2)*Rijk(1) ) / Temp1
              ELSE
                DRtAsc= 0.0D0
              ENDIF
            IF ( DABS( Temp ) .gt. Small ) THEN
                DDecl= ( Vijk(3) - DRr*DSIN( Decl ) ) / Temp
              ELSE
                DDecl= 0.0D0
              ENDIF
          ENDIF
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RV_TRADEC
*
*  this subroutine converts topocentric right-ascension declination with
*    position and velocity vectors. Uses velocity vector to find the
*    solution of singular cases.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rijk        - IJK position vector            km
*    Vijk        - IJK velocity vector            km/s
*    RSecef          - IJK SITE position vector       km
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    Rho         - Top Radius of the sat          km
*    TRtAsc      - Top Right Ascension            rad
*    TDecl       - Top Declination                rad
*    DRho        - Top Radius of the sat rate     km/s
*    TDRtAsc     - Top Right Ascension rate       rad/s
*    TDDecl      - Top Declination rate           rad/s
*
*  Locals        :
*    RhoV        - IJK Range Vector from SITE     km
*    DRhoV       - IJK Velocity Vector from SITE  km / s
*    Temp        - Temporary REAL*8 value
*    Temp1       - Temporary REAL*8 value
*    i           - Index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    LNCOM2      - Linear combination of 2 vectors
*    ADDVEC      - Add two vectors
*    DOT         - DOT product of two vectors
*
*  References    :
*    Vallado       2001, 248-250, Alg 26
*
* ------------------------------------------------------------------------------  
*
      SUBROUTINE RV_TRADEC   ( Rijk,Vijk,RSecef, Direction, Rho,TRtAsc,
     &                         TDecl,DRho,DTRtAsc,DTDecl )
        IMPLICIT NONE
        REAL*8 Rijk(3),VIjk(3),RSecef(3),Rho,TRtAsc,TDecl,DRho,DTRtAsc,
     &         DTDecl, MAG
        CHARACTER*4 Direction
        EXTERNAL DOT, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 Small, temp, temp1, RhoV(3),DRhoV(3), Dot, Magrhov

        ! --------------------  Implementation   ----------------------
        Small        = 0.00000001D0
        IF ( Direction .eq. 'FROM' ) THEN
            ! --------  Calculate Topocentric Vectors -----------------
            RhoV(1)= Rho*DCOS(TDecl)*DCOS(TRtAsc)
            RhoV(2)= Rho*DCOS(TDecl)*DSIN(TRtAsc)
            RhoV(3)= Rho*DSIN(TDecl)

            DRhoV(1)= DRho*DCOS(TDecl)*DCOS(TRtAsc)
     &                  - Rho*DSIN(TDecl)*DCOS(TRtAsc)*DTDecl
     &                  - Rho*DCOS(TDecl)*DSIN(TRtAsc)*DTRtAsc
            DRhoV(2)= DRho*DCOS(TDecl)*DSIN(TRtAsc)
     &                  - Rho*DSIN(TDecl)*DSIN(TRtAsc)*DTDecl
     &                  + Rho*DCOS(TDecl)*DCOS(TRtAsc)*DTRtAsc
            DRhoV(3)= DRho*DSIN(TDecl) + Rho*DCOS(TDecl)*DTDecl

            ! ------ Find IJK range vector from SITE to satellite -----
            CALL ADDVEC( RhoV,RSecef,  Rijk )
            DO i=1,3
                Vijk(i)= DRhoV(i)
              ENDDO

          ELSE
            ! ------ Find IJK range vector from SITE to satellite -----
            CALL LNCOM2( 1.0D0,-1.0D0, Rijk,RSecef,  RhoV )
            DO i=1,3
                DRhoV(i)= Vijk(i)  ! Same for topocentric
              ENDDO
            magrhov = MAG(rhoV)

            ! ------- Calculate Topocentric ANGLE and Rate Values -----
            Rho= MAG(Rhov)
            Temp= DSQRT( RhoV(1)*RhoV(1) + RhoV(2)*RhoV(2) )
            IF ( Temp .lt. Small ) THEN
                TRtAsc= DATAN2( DRhoV(2), DRhoV(1) )
              ELSE
                TRtAsc= DATAN2( RhoV(2), RhoV(1) )
              ENDIF

            TDecl= DASIN( RhoV(3)/magRhoV )

            Temp1= -RhoV(2)*RhoV(2) - RhoV(1)*RhoV(1)  ! different now
            DRho= DOT(RhoV,DRhoV)/Rho
            IF ( DABS(Temp1) .gt. Small ) THEN
                DTRtAsc= ( DRhoV(1)*RhoV(2) - DRhoV(2)*RhoV(1) ) / Temp1
              ELSE
                DTRtAsc= 0.0D0
              ENDIF
            IF ( DABS( Temp ) .gt. Small ) THEN
                DTDecl= ( DRhoV(3) - DRho*DSIN( TDecl ) ) / Temp
              ELSE
                DTDecl= 0.0D0
              ENDIF

          ENDIF
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RV_RAZEL
*
*  this subroutine converts Range, Azimuth, and Elevation and their rates with
*    the Geocentric Equatorial (IJK) Position and Velocity vectors.  Notice the
*    value of small as it can affect rate term calculations. Uses velocity
*    vector to find the solution of singular cases.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rijk        - IJK position vector            km
*    Vijk        - IJK velocity vector            km/s
*    RSecef          - IJK SITE Position Vector       km
*    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
*    Lon         - Site longitude                 -Pi to Pi rad
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    Rho         - Satellite Range from SITE      km
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*    El          - Elevation                      -Pi/2 to Pi/2 rad
*    DRho        - Range Rate                     km / s
*    DAz         - Azimuth Rate                   rad / s
*    DEl         - Elevation rate                 rad / s
*
*  Locals        :
*    RhoVijk     - IJK Range Vector from SITE     km
*    DRhoVijk    - IJK Velocity Vector from SITE  km / s
*    Rhosez      - SEZ Range vector from SITE     km
*    DRhosez     - SEZ Velocity vector from SITE  km
*    WCrossR     - CALL CROSS product result      km / s
*    EarthRate   - IJK Earth's rotation rate vec  rad / s
*    TempVec     - Temporary vector
*    Temp        - Temporary REAL*8 value
*    Temp1       - Temporary REAL*8 value
*    i           - Index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    ADDVEC      - Add two vectors
*    CROSS       - CROSS product of two vectors
*    ROT3        - Rotation about the 3rd axis
*    ROT2        - Rotation about the 2nd axis
*    DOT         - DOT product of two vectors
*    RVSEZ_RAZEL - Find R and V from SITE in Topocentric Horizon (SEZ) system
*    LNCOM2      - Combine two vectors and constants
*
*  References    :
*    Vallado       2001, 250-255, Alg 27
*
* ------------------------------------------------------------------------------

      SUBROUTINE RV_RAZEL    ( Reci,Veci,Latgd,Lon,alt,TTT,jdut1,lod,
     &                         xp,yp,terms, Direction,
     &                         Rho,Az,El,DRho,DAz,DEl )
        IMPLICIT NONE
        REAL*8 Reci(3),Veci(3),Latgd,Lon,Alt,Rho,Az,El,DRho,DAz,DEl,
     &         TTT,jdut1,lod,xp,yp
        INTEGER terms
        CHARACTER*4 Direction
        EXTERNAL Dot, MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 Temp, Temp1, Rhoecef(3), magrhosez, RSecef(3),VSecef(3),
     &         DRhoecef(3), Rhosez(3), DRhosez(3), WCrossR(3),
     &         TempVec, Dot, MAG, Lat,
     &         recef(3), vecef(3)
        INTEGER i

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

        ! --------------------  Implementation   ----------------------
        CALL SITE ( Latgd,Alt,Lon, RSecef,VSecef )

        IF ( Direction .eq. 'FROM' ) THEN
            ! --------  Find SEZ range and velocity vectors -----------
            CALL RVSEZ_RAZEL( Rhosez,DRhosez, 'FROM',
     &                        Rho,Az,El,DRho,DAz,DEl )

            ! ---------  Perform SEZ to ECEF transformation -----------
            CALL ROT2( Rhosez ,Latgd-HalfPi, TempVec )
            CALL ROT3( TempVec,   -Lon   , Rhoecef )
            CALL ROT2( DRhosez,Latgd-HalfPi, TempVec )
            CALL ROT3( TempVec,   -Lon   , vecef )

            ! -----------  Find range and velocity vectors ------------
            CALL ADDVEC( Rhoecef,RSecef,Recef )

            CALL GCRF_ITRF ( reci,veci, 'FROM', rECEF,vECEF,
     &                      TTT, JDUT1, LOD, xp, yp, terms )
          ELSE
*
            ! ---------------- convert eci to ecef --------------------
            CALL GCRF_ITRF  ( reci,veci, 'TOO ', rECEF,vECEF,
     &                      TTT, JDUT1, LOD, xp, yp, terms )

            ! ----- find ecef range vector from site to satellite -----
            CALL SUBVEC( recef, rsecef,  rhoecef)
            rho = mag(rhoecef)

            ! ----------- Convert to SEZ for calculations -------------
            CALL ROT3( Rhoecef,    Lon   ,  TempVec )
            CALL ROT2( TempVec,HalfPi-Latgd,   Rhosez   )
            CALL ROT3( vecef  ,    Lon   ,  TempVec )
            CALL ROT2( TempVec,HalfPi-Latgd,  DRhosez   )

            ! ----------- Calculate Azimuth and Elevation -------------
            Temp= DSQRT( Rhosez(1)*Rhosez(1) + Rhosez(2)*Rhosez(2) )
            IF ( Temp .lt. Small ) THEN
                Az = DATAN2( DRhosez(2) , -DRhosez(1) )
              ELSE
                Az = DATAN2( Rhosez(2) , -Rhosez(1) )
              ENDIF

            IF ( ( Temp .lt. Small ) ) THEN   ! directly over the north pole
                El= DSIGN(1.0D0, Rhosez(3))*HalfPi ! +- 90
              ELSE
                magrhosez = MAG(rhosez)
                El= DASIN( Rhosez(3) / magRhosez )
              ENDIF

            ! ---- Calculate Range, Azimuth and Elevation rates -------
            DRho= DOT(Rhosez,DRhosez)/Rho
            IF ( DABS( Temp*Temp ) .gt. Small ) THEN
                DAz= ( DRhosez(1)*Rhosez(2) - DRhosez(2)*Rhosez(1) ) /
     &                 ( Temp*Temp )
              ELSE
                DAz= 0.0D0
              ENDIF

            IF ( DABS( Temp ) .gt. 0.00000001D0 ) THEN
                DEl= ( DRhosez(3) - DRho*DSIN( El ) ) / Temp
              ELSE
                DEl= 0.0D0
              ENDIF

          ENDIF  ! If
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RV_ELATLON
*
*  this subroutine converts ecliptic latitude and longitude with position .and.
*    velocity vectors. Uses velocity vector to find the solution of singular
*    cases.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rijk        - IJK position vector            km
*    Vijk        - IJK velocity vector            km/s
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    Rr          - Radius of the sat              km
*    EclLat      - Ecliptic Latitude              -Pi/2 to Pi/2 rad
*    EclLon      - Ecliptic Longitude             -Pi/2 to Pi/2 rad
*    DRr         - Radius of the sat rate         km/s
*    DEclLat     - Ecliptic Latitude rate         -Pi/2 to Pi/2 rad
*    EEclLon     - Ecliptic Longitude rate        -Pi/2 to Pi/2 rad
*
*  Locals        :
*    Obliquity   - Obliquity of the ecliptic      rad
*    Temp        -
*    Temp1       -
*    Re          - Position vec in eclitpic frame
*    Ve          - Velocity vec in ecliptic frame
*
*  Coupling      :
*    ROT1        - Rotation about 1st axis
*    DOT         - DOT product
*
*  References    :
*    Vallado       2001, 257-259, Eq 4-15
*
* ------------------------------------------------------------------------------

      SUBROUTINE RV_ELATLON  ( Rijk,Vijk, Direction, rr,EclLat,EclLon,
     &                         DRr,DEclLat,DEclLon )
        IMPLICIT NONE
        REAL*8 Rijk(3), Vijk(3), rr,EclLat,EclLon,DRr,DEclLat,DEclLon
        CHARACTER*4 Direction
        EXTERNAL DOT, MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 Dot, Small, Re(3), Ve(3), Obliquity, Temp, Temp1, Mag

        ! --------------------  Implementation   ----------------------
        Small    = 0.00000001D0
        Obliquity= 0.40909280D0  !23.439291D0/rad
        IF ( Direction .eq. 'FROM' ) THEN
            Re(1)= rr*DCOS(EclLat)*DCOS(EclLon)
            Re(2)= rr*DCOS(EclLat)*DSIN(EclLon)
            Re(3)= rr*DSIN(EclLat)

            Ve(1)= DRr*DCOS(EclLat)*DCOS(EclLon)
     &               - rr*DSIN(EclLat)*DCOS(EclLon)*DEclLat
     &               - rr*DCOS(EclLat)*DSIN(EclLon)*DEclLon
            Ve(2)= DRr*DCOS(EclLat)*DSIN(EclLon)
     &               - rr*DSIN(EclLat)*DSIN(EclLon)*DEclLat
     &               + rr*DCOS(EclLat)*DCOS(EclLon)*DEclLon
            Ve(3)= DRr*DSIN(EclLat) + rr*DCOS(EclLat)*DEclLat

            CALL ROT1( Re, -Obliquity, Rijk )
            CALL ROT1( Ve, -Obliquity, Vijk )
          ELSE
            CALL ROT1( Rijk, Obliquity, Re )
            CALL ROT1( Vijk, Obliquity, Ve )

            ! ------------- Calculate Angles and Rates ----------------
            rr= MAG(Re)
            Temp= DSQRT( Re(1)*Re(1) + Re(2)*Re(2) )
            IF ( Temp .lt. Small ) THEN
                Temp1= DSQRT( Ve(1)*Ve(1) + Ve(2)*Ve(2) )
                IF ( DABS(Temp1) .gt. Small ) THEN
                    EclLon= DATAN2( Ve(2) , Ve(1) )
                  ELSE
                    EclLon= 0.0D0
                  ENDIF
              ELSE
                EclLon= DATAN2( Re(2) , Re(1) )
              ENDIF
            EclLat= DASIN( Re(3)/rr )

            Temp1= -Re(2)*Re(2) - Re(1)*Re(1)  ! different now
            DRr= DOT(re,Ve)/rr
            IF ( DABS( Temp1 ) .gt. Small ) THEN
                DEclLon= ( Ve(1)*Re(2) - Ve(2)*Re(1) ) / Temp1
              ELSE
                DEclLon= 0.0D0
              ENDIF
            IF ( DABS( Temp ) .gt. Small ) THEN
                DEclLat= ( Ve(3) - DRr*DSIN( EclLat ) ) / Temp
              ELSE
                DEclLat= 0.0D0
              ENDIF
          ENDIF  ! IF (

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RVSEZ_RAZEL
*
*  this subroutine converts range, azimuth, and elevation values with slant
*    range and velocity vectors for a satellite from a radar SITE in the
*    Topocentric Horizon (SEZ) system.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RhoVec      - SEZ Satellite range vector     km
*    DRhoVec     - SEZ Satellite velocity vector  km / s
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    Rho         - Satellite range from SITE      km
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*    El          - Elevation                      -Pi/2 to Pi/2 rad
*    DRho        - Range Rate                     km / s
*    DAz         - Azimuth Rate                   rad / s
*    DEl         - Elevation rate                 rad / s
*
*  Locals        :
*    SinEl       - Variable for DSIN( El )
*    CosEl       - Variable for DCOS( El )
*    SinAz       - Variable for DSIN( Az )
*    CosAz       - Variable for DCOS( Az )
*    Temp        -
*    Temp1       -
*
*  Coupling      :
*    DOT         - DOT product
*
*  References    :
*    Vallado       2001, 250-251, Eq 4-4, Eq 4-5
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE RVSEZ_RAZEL ( Rhosez,DRhosez,Direction, Rho,Az,El,
     &           DRho,DAz,DEl )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 RhoSez(3), DRhoSez(3),Rho,Az,El,DRho,DAz, DEl, MAG
        EXTERNAL Dot, MAG

* -----------------------------  Locals  ------------------------------
        REAL*8 Temp1, Temp, SinEl, CosEl, SinAz,CosAz, Dot, magrhosez

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        IF ( Direction .eq. 'FROM' ) THEN
            ! -------------------- Initialize values ------------------
            SinEl= DSIN(El)
            CosEl= DCOS(El)
            SinAz= DSIN(Az)
            CosAz= DCOS(Az)

            ! ----------------- Form SEZ range vector -----------------
            Rhosez(1) = -Rho*CosEl*CosAz
            Rhosez(2) =  Rho*CosEl*SinAz
            Rhosez(3) =  Rho*SinEl

            ! --------------- Form SEZ velocity vector ----------------
            DRhosez(1) = -DRho*CosEl*CosAz + Rhosez(3)*DEl*CosAz +
     &                     Rhosez(2)*DAz
            DRhosez(2) =  DRho*CosEl*SinAz - Rhosez(3)*DEl*SinAz -
     &                     Rhosez(1)*DAz
            DRhosez(3) =  DRho*SinEl       + Rho*DEl*CosEl
          ELSE
            ! ----------- Calculate Azimuth and Elevation -------------
            Temp= DSQRT( Rhosez(1)*Rhosez(1) + Rhosez(2)*Rhosez(2) )
            IF ( DABS( Rhosez(2) ) .lt. Small ) THEN
                IF ( Temp .lt. Small ) THEN
                    Temp1= DSQRT( DRhosez(1)*DRhosez(1) +
     &                     DRhosez(2)*DRhosez(2) )
                    Az   =  DATAN2( DRhosez(2)/Temp1 ,
     &                     -DRhosez(1)/Temp1 )
                  ELSE
                    IF ( Rhosez(1) .gt. 0.0D0 ) THEN
                        Az= Pi
                      ELSE
                        Az= 0.0D0
                      ENDIF
                  ENDIF
              ELSE
                Az= DATAN2( Rhosez(2)/Temp , -Rhosez(1)/Temp )
              ENDIF

            IF ( ( Temp .lt. Small ) ) THEN   ! directly over the north pole
                El= DSIGN(1.0D0,Rhosez(3))*HalfPi ! +- 90
              ELSE
                El= DASIN( Rhosez(3) / magRhosez )
              ENDIF

            ! -----  Calculate Range, Azimuth and Elevation rates -----
            DRho= DOT(Rhosez,DRhosez)/Rho 
            IF ( DABS( Temp*Temp ) .gt. Small ) THEN
                DAz= ( DRhosez(1)*Rhosez(2) - DRhosez(2)*Rhosez(1) ) /
     &                 ( Temp*Temp )
              ELSE
                DAz= 0.0D0
              ENDIF

            IF ( DABS( Temp ) .gt. Small ) THEN
                DEl= ( DRhosez(3) - DRho*DSIN( El ) ) / Temp
              ELSE
                DEl= 0.0D0 
              ENDIF
          ENDIF
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RADEC_ELATLON
*
*  this subroutine converts right-ascension declination values with ecliptic
*    latitude and longitude values.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    EclLat      - Ecliptic Latitude              -Pi/2 to Pi/2 rad
*    EclLon      - Ecliptic Longitude             -Pi/2 to Pi/2 rad
*
*  Locals        :
*    Obliquity   - Obliquity of the ecliptic      rad
*    Sinv        -
*    Cosv        -
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2001, 259, Eq 4-19, Eq 4-20
*
* ------------------------------------------------------------------------------  

      SUBROUTINE RADEC_ELATLON ( RtAsc,Decl,Direction, EclLat, EclLon )
        IMPLICIT NONE
        REAL*8 RtAsc,Decl,EclLat,EclLon
        CHARACTER*4 Direction
* -----------------------------  Locals  ------------------------------
        REAL*8 Sinv, Cosv, Obliquity

        ! --------------------  Implementation   ----------------------
        Obliquity= 0.40909280D0  !23.439291D0/rad
        IF ( Direction .eq. 'FROM' ) THEN
            Decl = DASIN( DSIN(EclLat)*DCOS(Obliquity)
     &                    + DCOS(EclLat)*DSIN(Obliquity)*DSIN(EclLon) )
            Sinv = ( -DSIN(EclLat)*DSIN(Obliquity)
     &                 + DCOS(EclLat)*DCOS(Obliquity)*DSIN(EclLon) ) /
     &                 DCOS(Decl)
            Cosv = DCOS(EclLat)*DCOS(EclLon) / DCOS(Decl) 
            RtAsc= DATAN2( Sinv,Cosv ) 
          ELSE
            EclLat= DASIN( -DCOS(Decl)*DSIN(RtAsc)*DSIN(Obliquity)
     &                        + DSIN(Decl)*DCOS(Obliquity) )
            Sinv  = ( DCOS(Decl)*DSIN(RtAsc)*DCOS(Obliquity)
     &                  + DSIN(Decl)*DSIN(Obliquity) ) / DCOS(EclLat)
            Cosv  = DCOS(Decl)*DCOS(RtAsc) / DCOS(EclLat) 
            EclLon= DATAN2( Sinv,Cosv ) 
          ENDIF 
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RADEC_AZEL
*
* this subroutine converts right ascension declination values with
*   azimuth, and elevation.  Notice the range is not defined because
*   Right ascension declination only allows a unit vector to be formed.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RtAsc       - Right Ascension                0.0D0 to 2Pi rad
*    Decl        - Declination                    -Pi/2 to Pi/2 rad
*    LST         - Local SIDEREAL Time            -2Pi to 2Pi rad
*    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*    El          - Elevation                      -Pi/2 to Pi/2 rad
*
*  Locals        :
*    LHA         - Local Hour ANGLE               -2Pi to 2Pi rad
*    Sinv        - Sine value
*    Cosv        - Cosine value
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2001, 255-257, Alg 28
*
* ------------------------------------------------------------------------------  

      SUBROUTINE RADEC_AZEL  ( RtAsc,Decl,LST,Latgd, Direction, Az,El )
        IMPLICIT NONE
        REAL*8 RtAsc,Decl,LST,Latgd,Az,El
        CHARACTER*4 Direction
* -----------------------------  Locals  ------------------------------
        REAL*8 Sinv, Cosv, LHA

        ! --------------------  Implementation   ----------------------
        IF ( Direction .eq. 'FROM' ) THEN
            Decl = DASIN( DSIN(El)*DSIN(Latgd) +
     &                 DCOS(el)*DCOS(Latgd)*DCOS(Az) )

            Sinv = -(DSIN(az)*DCOS(el)*DCOS(Latgd)) /
     &              (DCOS(Latgd)*DCOS(Decl))
            Cosv = (DSIN(el) - DSIN(Latgd)*DSIN(decl)) /
     &              (DCOS(Latgd)*DCOS(Decl))
            LHA  = DATAN2( Sinv,Cosv ) 
            RtAsc= LST - LHA 
          ELSE
            LHA = LST - RtAsc

            El  = DASIN( DSIN(Decl)*DSIN(Latgd) +
     &            DCOS(Decl)*DCOS(Latgd)*DCOS(LHA) )

            Sinv= -DSIN(LHA)*DCOS(Decl)*DCOS(Latgd)/
     &                (DCOS(el)*DCOS(Latgd))
            Cosv= ( DSIN(Decl)-DSIN(el)*DSIN(Latgd) )/
     &             (DCOS(el)*DCOS(Latgd))
            Az  = DATAN2( Sinv,Cosv ) 
          ENDIF 

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE GIBBS
*
*  this subroutine performs the GIBBS method of orbit determination.  This
*    method determines the velocity at the middle point of the 3 given position
*    vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - IJK Position vector #1         km
*    R2          - IJK Position vector #2         km
*    R3          - IJK Position vector #3         km
*
*  OutPuts       :
*    V2          - IJK Velocity Vector for R2     km / s
*    Theta       - ANGLE between vectors          rad
*    Error       - Flag indicating success        'ok',...
*
*  Locals        :
*    tover2      -
*    l           -
*    Small       - Tolerance for roundoff errors
*    r1mr2       - Magnitude of r1 - r2
*    r3mr1       - Magnitude of r3 - r1
*    r2mr3       - Magnitude of r2 - r3
*    p           - P Vector     r2 x r3
*    q           - Q Vector     r3 x r1
*    w           - W Vector     r1 x r2
*    d           - D Vector     p + q + w
*    n           - N Vector (r1)p + (r2)q + (r3)w
*    s           - S Vector
*                    (r2-r3)r1+(r3-r1)r2+(r1-r2)r3
*    b           - B Vector     d x r2
*    Theta1      - Temp ANGLE between the vectors rad
*    Pn          - P Unit Vector
*    R1N         - R1 Unit Vector
*    dn          - D Unit Vector
*    Nn          - N Unit Vector
*    i           - index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    CROSS       - CALL CROSS product of two vectors
*    DOT         - DOT product of two vectors
*    ADD3VEC     - Add three vectors
*    LNCOM2      - Multiply two vectors by two constants
*    LNCOM3      - Add three vectors each multiplied by a constant
*    NORM        - Creates a Unit Vector
*    ANGLE       - ANGLE between two vectors
*
*  References    :
*    Vallado       2001, 432-445, Alg 52, Ex 7-5
*
* ------------------------------------------------------------------------------  
*
      SUBROUTINE GIBBS       ( R1,R2,R3, V2, Theta,Theta1,Copa, Error )
        IMPLICIT NONE
        REAL*8 R1(3), R2(3), R3(3), V2(3), Theta, Theta1, Copa
        CHARACTER*12 Error
        EXTERNAL DOT, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 tover2, l, Small, r1mr2, r3mr1, r2mr3, p(3), q(3), w(3),
     &         d(3), n(3), s(3), b(3), Pn(3), R1N(3), Dn(3), Nn(3),Dot,
     &         magr1, magr2, magr3, mag, magd, magn

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Small= 0.000001D0
        Theta= 0.0D0
        Error = 'ok'
        Theta1= 0.0D0
        magr1 = MAG( R1 )
        magr2 = MAG( R2 )
        magr3 = MAG( R3 )
        DO i= 1 , 3
            V2(i)= 0.0D0
          ENDDO

        CALL CROSS( R2,R3,P )
        CALL CROSS( R3,R1,Q )
        CALL CROSS( R1,R2,W )
        CALL NORM( P,Pn )
        CALL NORM( R1,R1N )
        Copa=  DASIN( DOT( Pn,R1n ) ) 

        IF ( DABS( DOT(R1N,Pn) ) .gt. 0.017452406D0 ) THEN
            Error= 'not coplanar'
          ENDIF

        ! --------------- .or. don't contiNue processing --------------
        CALL ADD3VEC( P,Q,W,D )
        CALL LNCOM3( magr1,magr2,magr3,P,Q,W,N )
        CALL NORM( N,Nn )
        CALL NORM( D,DN )
        magd = MAG(d)
        magn = MAG(n)

        ! -------------------------------------------------------------
*       Determine If  the orbit is possible.  Both D and N must be in
*         the same direction, and non-zero.
        ! -------------------------------------------------------------
        IF ( ( DABS(magd).lt.Small ) .or. ( DABS(magn).lt.Small ) .or.
     &      ( DOT(Nn,dn) .lt. Small ) ) THEN
            Error= 'impossible'
          ELSE
              CALL ANGLE( R1,R2, Theta )
              CALL ANGLE( R2,R3, Theta1 )

              ! ----------- Perform GIBBS method to find V2 -----------
              R1mr2= magr1-magr2
              R3mr1= magr3-magr1
              R2mr3= magr2-magr3
              CALL LNCOM3(R1mr2,R3mr1,R2mr3,R3,R2,R1,S)
              CALL CROSS( d,r2,b )
              L    = DSQRT( mu / (magd*magn) )
              Tover2= L / magr2
              CALL LNCOM2(Tover2,L,B,S,V2)
            ENDIF

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HERRGIBBS
*
*  this subroutine implements the Herrick-GIBBS approximation for orbit
*    determination, and finds the middle velocity vector for the 3 given
*    position vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - IJK Position vector #1         km
*    R2          - IJK Position vector #2         km
*    R3          - IJK Position vector #3         km
*    JD1         - Julian Date of 1st sighting    days from 4713 BC
*    JD2         - Julian Date of 2nd sighting    days from 4713 BC
*    JD3         - Julian Date of 3rd sighting    days from 4713 BC
*
*  OutPuts       :
*    V2          - IJK Velocity Vector for R2     km / s
*    Theta       - ANGLE between vectors          rad
*    Error       - Flag indicating success        'ok',...
*
*  Locals        :
*    Dt21        - time delta between r1 and r2   s
*    Dt31        - time delta between r3 and r1   s
*    Dt32        - time delta between r3 and r2   s
*    p           - P vector    r2 x r3
*    Pn          - P Unit Vector
*    R1N         - R1 Unit Vector
*    Theta1      - temporary ANGLE between vec    rad
*    TolAngle    - Tolerance ANGLE  (1 deg)       rad
*    Term1       - 1st Term for HGibbs expansion
*    Term2       - 2nd Term for HGibbs expansion
*    Term3       - 3rd Term for HGibbs expansion
*    i           - Index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    CROSS       - CALL CROSS product of two vectors
*    DOT         - DOT product of two vectors
*    NORM        - Creates a Unit Vector
*    LNCOM3      - Combination of three scalars and three vectors
*    ANGLE       - ANGLE between two vectors
*
*  References    :
*    Vallado       2001, 439-445, Alg 52, Ex 7-4
*
* ------------------------------------------------------------------------------  
*
      SUBROUTINE HERRGIBBS   ( R1,R2,R3,JD1,JD2,JD3, V2, Theta,Theta1,
     &                         Copa, Error )
        IMPLICIT NONE
        REAL*8 R1(3), R2(3), R3(3), JD1, JD2, JD3, V2(3), Theta,
     &         Theta1, Copa
        CHARACTER*12 Error
        EXTERNAL Dot
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 p(3), Pn(3), R1n(3),Dot, magr1, magr2, magr3,
     &         Dt21, Dt31, Dt32, Term1, Term2, Term3, TolAngle, mag

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Error =  'ok'
        Theta = 0.0D0
        Theta1= 0.0D0
        magr1 = MAG( R1 )
        magr2 = MAG( R2 )
        magr3 = MAG( R3 )
        DO i= 1 , 3
            V2(i)= 0.0D0
          ENDDO
        TolAngle= 0.01745329251994D0
        Dt21= (JD2-JD1)*86400.0D0
        Dt31= (JD3-JD1)*86400.0D0   ! differences in times
        Dt32= (JD3-JD2)*86400.0D0

        CALL CROSS( R2,R3,P )
        CALL NORM( P,Pn )
        CALL NORM( R1,R1N )
        Copa=  DASIN( DOT( Pn,R1n ) )
        IF ( DABS( DOT(R1N,Pn) ) .gt. 0.017452406D0 ) THEN
            Error= 'not coplanar'
          ENDIF

        ! --------------------------------------------------------------
*       Check the size of the angles between the three position vectors.
*       Herrick GIBBS only gives "reasonable" answers when the
*       position vectors are reasonably close.  10 deg is only an estimate.
        ! --------------------------------------------------------------
        CALL ANGLE( R1,R2, Theta )
        CALL ANGLE( R2,R3, Theta1 )
        IF ( (Theta .gt. TolAngle) .or. (Theta1 .gt. TolAngle) ) THEN
            Error= 'ANGLE > 1 deg'
          ENDIF

        ! ----------- Perform Herrick-GIBBS method to find V2 ---------
        Term1= -Dt32*( 1.0D0/(Dt21*Dt31) +
     &        mu/(12.0D0*magr1*magr1*magr1) )
        Term2= (Dt32-Dt21)*( 1.0D0/(Dt21*Dt32) +
     &        mu/(12.0D0*magr2*magr2*magr2) )
        Term3=  Dt21*( 1.0D0/(Dt32*Dt31) +
     &         mu/(12.0D0*magr3*magr3*magr3) )
        CALL LNCOM3( Term1,Term2,Term3,R1,R2,R3, V2 )

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LAMBERTUNIV
*
*  this subroutine solves the Lambert problem for orbit determination and returns
*    the velocity vectors at each of two given position vectors.  The solution
*    uses Universal Variables for calculation and a bissection technique
*    updating psi.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - IJK Position vector 1          km
*    R2          - IJK Position vector 2          km
*    DM          - direction of motion            'L','S'
*    Dtsec        - Time between R1 and R2         s
*
*  OutPuts       :
*    V1          - IJK Velocity vector            km / s
*    V2          - IJK Velocity vector            km / s
*    Error       - Error flag                     'ok', ...
*
*  Locals        :
*    VarA        - Variable of the iteration,
*                  NOT the semi .or. axis!
*    Y           - Area between position vectors
*    Upper       - Upper bound for Z
*    Lower       - Lower bound for Z
*    CosDeltaNu  - Cosine of true anomaly change  rad
*    F           - f expression
*    G           - g expression
*    GDot        - g DOT expression
*    XOld        - Old Universal Variable X
*    XOldCubed   - XOld cubed
*    ZOld        - Old value of z
*    ZNew        - New value of z
*    C2New       - C2(z) FUNCTION
*    C3New       - C3(z) FUNCTION
*    TimeNew     - New time                       s
*    Small       - Tolerance for roundoff errors
*    i, j        - index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    DOT         - DOT product of two vectors
*    FINDC2C3    - Find C2 and C3 functions
*
*  References    :
*    Vallado       2001, 459-464, Alg 55, Ex 7-5
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE LAMBERTUNIV ( ro,r, dm,OverRev, Dtsec, mu,vo,v, Error )
        IMPLICIT NONE
        REAL*8 Ro(3), r(3), Dtsec, Vo(3),v(3),mu
        CHARACTER Dm, OverRev
        CHARACTER*12 Error
        EXTERNAL DOT, MAG

* -----------------------------  Locals  ------------------------------
        INTEGER i, Loops, YNegKtr, NumIter
        REAL*8 VarA, Y, Upper, Lower,CosDeltaNu, magro, magr,
     &         F, G, GDot, XOld, XOldCubed, Dot, mag,
     &         PsiOld, PsiNew, C2New, C3New, dtNew

!         INCLUDE 'astmath.cmn'
!         INCLUDE 'astconst.cmn'
        REAL*8 Small,    Rad2Deg,  Deg2Rad,  HalfPi,
     &         Pi,       TwoPi,    Infinite, Undefined

        REAL*8     rekm,     omegaearth, flat,     EESqrd, auer

        DATA rekm       /6378.137D0/
!         DATA mu         /398600.4418D0/
!         DATA mu         /1.32712440018D11/       ! Sun
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

        ! --------------------  Implementation   ----------------------
        NumIter= 40
        Error  = 'ok' 
        PsiNew = 0.0D0 
        magro = MAG(ro)
        magr = MAG(r)
        DO i= 1 , 3
            vo(i)= 0.0D0
            v(i)= 0.0D0
          ENDDO

        CosDeltaNu= DOT(ro,r)/(magro*magr)
        IF ( Dm .eq. 'L' ) THEN
            VarA = -DSQRT( magro*magr*(1.0D0+CosDeltaNu) )
          ELSE
            VarA =  DSQRT( magro*magr*(1.0D0+CosDeltaNu) )
          ENDIF

        ! ---------------  Form Initial guesses   ---------------------
        PsiOld = 0.0D0 
        PsiNew = 0.0D0 
        xOld   = 0.0D0 
        C2New  = 0.5D0 
        C3New  = 1.0D0/6.0D0

        ! --------- Set up initial bounds for the bissection ----------
        IF ( OverRev .eq. 'N' ) THEN
            Upper= TwoPi*TwoPi
            Lower= -4.0D0*TwoPi 
          ELSE
            Upper= -0.001D0+4.0D0*TwoPi*TwoPi ! at 4, not alw work
            Lower=  0.001D0+TwoPi*TwoPi       ! 2.0D0, makes orbit bigger
          ENDIF                               ! how about 2 revs??xx

        ! -------  Determine If  the orbit is possible at all ---------
        IF ( DABS( VarA ) .gt. Small ) THEN
*
            Loops  = 0 
            YNegKtr= 1  ! y neg ktr
            DtNew = -10.0D0
            DO WHILE ((DABS(dtNew-Dtsec) / Dtsec .ge. 0.000001D0) .and.
     &              (Loops .lt. NumIter) .and. (YNegKtr .le. 10))
                IF ( DABS(C2New) .gt. Small ) THEN
                    Y= magro + magr -
     &                      ( VarA*(1.0D0-PsiOld*C3New)/DSQRT(C2New) )
                  ELSE
                    Y= magro + magr
                  ENDIF
                ! ----------- Check for negative values of y ----------
                IF ( ( VarA .gt. 0.0D0 ) .and. ( Y .lt. 0.0D0 ) ) THEN
                    YNegKtr= 1
                    DO WHILE (( Y.lt.0.0D0 ) .and. ( YNegKtr .lt. 10 ))
                        PsiNew= 0.8D0*(1.0D0/C3New)*( 1.0D0
     &                           - (magro+magr)*DSQRT(C2New)/VarA  )
                        ! -------- Find C2 and C3 functions -----------
                        CALL FINDC2C3( PsiNew, C2New,C3New )
                        PsiOld= PsiNew 
                        Lower= PsiOld 
                        IF ( DABS(C2New) .gt. Small ) THEN
                            Y= magro + magr - ( VarA*(1.0D0-
     &                               PsiOld*C3New)/DSQRT(C2New) )
                          ELSE
                            Y= magro + magr
                          ENDIF
                        YNegKtr = YNegKtr + 1
                      ENDDO ! while
                  ENDIF  ! If  y neg

                IF ( YNegKtr .lt. 10 ) THEN
                    IF ( DABS(C2New) .gt. Small ) THEN
                        XOld= DSQRT( Y/C2New )
                      ELSE
                        XOld= 0.0D0
                      ENDIF
                    XOldCubed= XOld*XOld*XOld 
                    dtNew    = (XOldCubed*C3New + VarA*DSQRT(Y)) /
     &                          DSQRT(mu)

                    ! --------  Readjust upper and lower bounds -------
                    IF ( dtNew .lt. Dtsec ) THEN
                        Lower= PsiOld 
                      ENDIF
                    IF ( dtNew .gt. Dtsec ) THEN
                        Upper= PsiOld 
                      ENDIF
                    PsiNew= (Upper+Lower) * 0.5D0 

                    ! ------------- Find c2 and c3 functions ----------
                    CALL FINDC2C3( PsiNew, C2New,C3New )
                    PsiOld = PsiNew
                    Loops = Loops + 1

                    ! --- Make sure the first guess isn't too close ---
                    IF ( (DABS(dtNew - Dtsec) .lt. Small)
     &                        .and. (Loops .eq. 1) ) THEN
                        dtNew= Dtsec-1.0D0
                      ENDIF
                  ENDIF  ! If  YNegKtr .lt. 10
c               write(20,'(4f14.6)') y,xold,dtnew,psinew
              ENDDO ! Do While Loop

            IF ( (Loops .ge. NumIter) .or. (YNegKtr .ge. 10) ) THEN
                Error= 'GNotConv'
                IF ( YNegKtr .ge. 10 ) THEN
                    Error= 'Y negative' 
                  ENDIF
              ELSE
                ! --- Use F and G series to find Velocity Vectors -----
                F   = 1.0D0 - Y/magro
                GDot= 1.0D0 - Y/magr
                G   = 1.0D0 / (VarA*DSQRT( Y/mu ))  ! 1 over G
                DO i= 1 , 3
                    vo(i)= ( r(i) - F*ro(i) )*G
                    v(i) = ( GDot*r(i) - ro(i) )*G
                  ENDDO
              ENDIF   ! If  the answer has converged
          ELSE
            Error= 'impos 180 deg' 
          ENDIF  ! If  Var A .gt. 0.0D0

      RETURN
      END
*
* --------- Two recursion algorithms needed by the LambertBattin routine
*
      REAL*8 FUNCTION SEE    ( v )
        IMPLICIT NONE
        REAL*8 v
* -----------------------------  Locals  ------------------------------
        REAL*8 c(0:20),term, termold, del, delold, sum1, eta, SQRTopv
        INTEGER i

        ! --------------------  Implementation   ----------------------
          c(0) =    0.2D0
          c(1) =    9.0D0 /  35.0D0
          c(2) =   16.0D0 /  63.0D0
          c(3) =   25.0D0 /  99.0D0
          c(4) =   36.0D0 / 143.0D0
          c(5) =   49.0D0 / 195.0D0
          c(6) =   64.0D0 / 255.0D0
          c(7) =   81.0D0 / 323.0D0
          c(8) =  100.0D0 / 399.0D0
          c(9) =  121.0D0 / 483.0D0
          c(10)=  144.0D0 / 575.0D0
          c(11)=  169.0D0 / 675.0D0
          c(12)=  196.0D0 / 783.0D0
          c(13)=  225.0D0 / 899.0D0
          c(14)=  256.0D0 /1023.0D0
          c(15)=  289.0D0 /1155.0D0
          c(16)=  324.0D0 /1295.0D0
          c(17)=  361.0D0 /1443.0D0
          c(18)=  400.0D0 /1599.0D0
          c(19)=  441.0D0 /1763.0D0
          c(20)=  484.0D0 /1935.0D0
          SQRTOpv= DSQRT(1.0D0 + v) 
          eta    = v / ( ( 1.0D0+SQRTOpv )**2 )

          ! ------------------- Process Forwards ----------------------
          delold = 1.0D0
          termold= c(0)   ! * eta}
          sum1   = termold 
          i= 1 
          DO WHILE ((i.le.20) .and. (DABS(Termold) .gt. 0.000001D0 ))
              del  = 1.0D0 / ( 1.0D0 + c(i)*eta*delold )
              term = termold * (del - 1.0D0) 
              sum1 = sum1 + term 
              i    = i + 1
              delold = del
              termold= term 
            ENDDO

c          See= (1.0D0 / (8.0D0*(1.0D0+SQRTOpv))) *
c     &         ( 3.0D0 + Sum1 / ( 1.0D0+eta*sum1 ) )
          See= 1.0D0 / ((1.0D0/(8.0D0*(1.0D0+sqrtopv))) *
     &         ( 3.0D0 + sum1 / ( 1.0D0+eta*sum1 ) ) )
      RETURN
      END


      REAL*8 FUNCTION k      ( v )
        IMPLICIT NONE
        REAL*8 v
* -----------------------------  Locals  ------------------------------
        REAL*8 d(0:20), del,delold,term,termold, sum1
        INTEGER i

        ! --------------------  Implementation   ----------------------
          d(0) =     1.0D0 /    3.0D0
          d(1) =     4.0D0 /   27.0D0
          d(2) =     8.0D0 /   27.0D0
          d(3) =     2.0D0 /    9.0D0
          d(4) =    22.0D0 /   81.0D0
          d(5) =   208.0D0 /  891.0D0
          d(6) =   340.0D0 / 1287.0D0
          d(7) =   418.0D0 / 1755.0D0
          d(8) =   598.0D0 / 2295.0D0
          d(9) =   700.0D0 / 2907.0D0
          d(10)=   928.0D0 / 3591.0D0
          d(11)=  1054.0D0 / 4347.0D0
          d(12)=  1330.0D0 / 5175.0D0
          d(13)=  1480.0D0 / 6075.0D0
          d(14)=  1804.0D0 / 7047.0D0
          d(15)=  1978.0D0 / 8091.0D0
          d(16)=  2350.0D0 / 9207.0D0
          d(17)=  2548.0D0 /10395.0D0
          d(18)=  2968.0D0 /11655.0D0
          d(19)=  3190.0D0 /12987.0D0
          d(20)=  3658.0D0 /14391.0D0

          ! ----------------- Process Forwards ------------------------
          sum1   = d(0)
          delold = 1.0D0 
          termold= d(0)
          i      = 1
          DO WHILE ((i.le.20) .and. (DABS(Termold) .gt. 0.000001D0 ))
              del  = 1.0D0 / ( 1.0D0 + d(i)*v*delold )
              term = termold * ( del - 1.0D0 )
              sum1 = sum1 + term 
              i    = i + 1
              delold = del 
              termold= term 
            ENDDO

          k= Sum1
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LAMBERBATTIN
*
*  this subroutine solves Lambert's problem using Battins method. The method is
*    developed in Battin (1987). It uses contiNued fractions to speed the
*    solution and has several parameters that are defined differently than
*    the traditional Gaussian technique.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ro          - IJK Position vector 1          km
*    R           - IJK Position vector 2          km
*    DM          - direction of motion            'L','S'
*    Dtsec        - Time between R1 and R2         s
*
*  OutPuts       :
*    Vo          - IJK Velocity vector            km / s
*    V           - IJK Velocity vector            km / s
*    Error       - Error flag                     'ok',...
*
*  Locals        :
*    i           - Index
*    Loops       -
*    u           -
*    b           -
*    Sinv        -
*    Cosv        -
*    rp          -
*    x           -
*    xn          -
*    y           -
*    l           -
*    m           -
*    CosDeltaNu  -
*    SinDeltaNu  -
*    DNu         -
*    a           -
*    Tan2w       -
*    RoR         -
*    h1          -
*    h2          -
*    Tempx       -
*    eps         -
*    denom       -
*    chord       -
*    k2          -
*    s           -
*    f           -
*    g           -
*    gDot        -
*    am          -
*    ae          -
*    be          -
*    tm          -
*    arg1        -
*    arg2        -
*    tc          -
*    AlpE        -
*    BetE        -
*    AlpH        -
*    BetH        -
*    DE          -
*    DH          -
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    ASINH     - Inverse hyperbolic sine
*    ARCCOSH     - Inverse hyperbolic cosine
*    SINH        - Hyperbolic sine
*
*  References    :
*    Vallado       2001, 464-467, Ex 7-5
*
* ------------------------------------------------------------------------------  
*
      SUBROUTINE LAMBERTBATTIN ( ro,r, dm,OverRev, Dtsec, mu,
     &                          vo,v, Error )
        IMPLICIT NONE
        REAL*8 Ro(3), r(3), Dtsec, Vo(3),v(3),mu
        CHARACTER Dm, OverRev
        CHARACTER*12 Error
        EXTERNAL Dot, See, K, ASINH, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i, Loops
        REAL*8 RCrossR(3),Dot, See, K, ASINH, MAG,
*         lambda,bigt,     testamt,
     &   u, b, Sinv,Cosv, rp, x, xn, y, l, m, CosDeltaNu, SinDeltaNu,
     &   DNu, a, tan2w, RoR, h1, h2, tempx, eps, denom, chord, k2, s,
     &   f, g, am, ae, be, tm, gDot, arg1, arg2, AlpE, BetE,
     &   AlpH, BetH, DE, DH, magr, magro, magrcrossr,y1, lim1

!         INCLUDE 'astmath.cmn'
!         INCLUDE 'astconst.cmn'
        REAL*8 Small,    Rad2Deg,  Deg2Rad,  HalfPi,
     &         Pi,       TwoPi,    Infinite, Undefined

        REAL*8     rekm,     omegaearth, flat,     EESqrd, auer

        DATA rekm       /6378.137D0/
!         DATA mu         /398600.4418D0/
!         DATA mu         /1.32712440018D11/       ! Sun
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

        ! --------------------  Implementation   ----------------------
        Error = 'ok'
        magr = MAG(r)
        magro = MAG(ro)
        CosDeltaNu= DOT(ro,r)/(magro*magr)
        CALL CROSS( ro,r, RCrossR )
        magrcrossr = MAG(rcrossr)
        SinDeltaNu= magRCrossr/(magro*magr)
        DNu   = DATAN2( SinDeltaNu,CosDeltaNu )

        RoR   = magr/magro
        eps   = RoR - 1.0D0
        tan2w = 0.25D0*eps*eps / ( DSQRT( RoR ) + RoR *
     &                             ( 2.0D0 + DSQRT( RoR ) ) )
        rp    = DSQRT( magro*magr )*( (DCOS(DNu*0.25D0))**2 + tan2w )

        IF ( DNu .lt. Pi ) THEN
            L = ( (DSIN(DNu*0.25D0))**2 + tan2w ) /
     &          ( (DSIN(DNu*0.25D0))**2 + tan2w + DCOS( DNu*0.5D0 ) )
          ELSE
            L = ( (DCOS(DNu*0.25D0))**2 + tan2w - DCOS( DNu*0.5D0 ) ) /
     &            ( (DCOS(DNu*0.25D0))**2 + tan2w )
          ENDIF

        m    = mu*Dtsec*Dtsec / ( 8.0D0*rp*rp*rp )
        x    = 10.0D0
        xn   = L  !0.0D0   !L    ! 0 for par and hyp
        chord= DSQRT( magro*magro + magr*magr - 2.0D0*magro*magr*
     &         DCOS( DNu ) )
        s    = ( magro + magr + chord )*0.5D0
        lim1 = dsqrt(m/l)

        Loops= 1
        DO WHILE ((DABS(xn-x) .ge. Small) .and. (Loops .le. 30))
            x    = xn 
            Tempx= See(x) 
            Denom= 1.0D0 / ( (1.0D0+2.0D0*x+L) * (4.0D0*x +
     &             tempx*(3.0D0+x) ) )
            h1   = ( L+x )**2 * ( 1.0D0+ 3.0D0*x + Tempx )*Denom
            h2   = m*( x - L + Tempx )*Denom

            ! ----------------------- Evaluate CUBIC ------------------
            b = 0.25D0*27.0D0*h2 / ((1.0D0+h1)**3 )
            if (b .lt. -1.0D0) THEN ! reset the initial condition
                xn = 1.0D0 - 2.0D0*l
             else
                if (y1 .gt. lim1) THEN
                    xn = xn * (lim1/y1)
                else
                  u = 0.5D0*b / ( 1.0D0 + DSQRT( 1.0D0 + b ) ) 
                  K2= K(u) 

                   y = ( ( 1.0D0+h1 ) / 3.0D0 )*
     &                 ( 2.0D0 + DSQRT( 1.0D0+b ) /
     &                 ( 1.0D0+2.0D0*u*k2*k2 ) )
                   xn= DSQRT( ( (1.0D0-L)*0.5D0 )**2 + m/(y*y) ) -
     &                 ( 1.0D0+L )*0.5D0
                endif
             endif

            Loops = Loops + 1
        y1=dsqrt(m/((L+x)*(1.0D0+x)) )
        write(*,'(i3,6f11.7)') loops,y,x,k2,b,u,y1
          ENDDO

        a=  mu*Dtsec*Dtsec / (16.0D0*Rp*rp*xn*y*y )
*
        ! ------------------ Find Eccentric anomalies -----------------
        ! ------------------------ Hyperbolic -------------------------
        IF ( a .lt. -Small ) THEN
            arg1 = DSQRT( s / ( -2.0D0*a ) )
            arg2 = DSQRT( ( s-chord ) / ( -2.0D0*a ) )
            ! ------- Evaluate f and g functions --------
            AlpH = 2.0D0 * ASINH( arg1 )
            BetH = 2.0D0 * ASINH( arg2 )
            DH   = AlpH - BetH
            F    = 1.0D0 - (a/magro)*(1.0D0 - COSH(DH) )
            GDot = 1.0D0 - (a/magr) *(1.0D0 - COSH(DH) )
            G    = Dtsec - DSQRT(-a*a*a/mu)*(SINH(DH)-DH)
          ELSE
            ! ------------------------ Elliptical ---------------------
            IF ( a .gt. small ) THEN
                arg1 = DSQRT( s / ( 2.0D0*a ) )
                arg2 = DSQRT( ( s-chord ) / ( 2.0D0*a ) )
                Sinv = arg2
                Cosv = DSQRT( 1.0D0 - (magro+magr-chord)/(4.0D0*a) )
                BetE = 2.0D0*DACOS(Cosv)
                BetE = 2.0D0*DASIN(Sinv)
                IF ( DNu .gt. Pi ) THEN
                    BetE= -BetE
                  ENDIF

                Cosv= DSQRT( 1.0D0 - s/(2.0D0*a) ) 
                Sinv= arg1 

                am  = s*0.5D0 
                ae  = Pi 
                be  = 2.0D0*DASIN( DSQRT( (s-chord)/s ) ) 
                tm  = DSQRT(am*am*am/mu)*(ae - (be-DSIN(be)))
                IF ( Dtsec .gt. tm ) THEN
                    AlpE= 2.0D0*pi-2.0D0*DASIN( Sinv )
                  ELSE
                    AlpE= 2.0D0*DASIN( Sinv )
                  ENDIF
                DE  = AlpE - BetE 
                F   = 1.0D0 - (a/magro)*(1.0D0 - DCOS(DE) )
                GDot= 1.0D0 - (a/magr)* (1.0D0 - DCOS(DE) )
                G   = Dtsec - DSQRT(a*a*a/mu)*(DE - DSIN(DE))
              ELSE
                ! --------------------- Parabolic ---------------------
                arg1 = 0.0D0
                arg2 = 0.0D0 
                Error= 'a = 0 '
                Write(10,*) ' a parabolic orbit '
              ENDIF
          ENDIF

        DO i= 1 , 3
            vo(i)= ( r(i) - F*ro(i) )/G
            v(i) = ( GDot*r(i) - ro(i) )/G
          ENDDO

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE TARGET
*
*  this subroutine accomplishes the targeting problem using KEPLER/PKEPLER .and.
*    Lambert.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInt        - Initial Position vector of Int km
*    VInt        - Initial Velocity vector of Int km/s
*    RTgt        - Initial Position vector of Tgt km
*    VTgt        - Initial Velocity vector of Tgt km/s
*    dm          - Direction of Motion for Gauss  'L','S'
*    Kind        - Type of propagator             'K','P'
*    Dtsec        - Time of flight to the int      s
*    mu          - Grav. constant                 km3/s2
*
*  Outputs       :
*    V1t         - Initial Transfer Velocity vec  km/s
*    V2t         - Final Transfer Velocity vec    km/s
*    DV1         - Initial Change Velocity vec    km/s
*    DV2         - Final Change Velocity vec      km/s
*    Error       - Error flag from Gauss          'ok', ...
*
*  Locals        :
*    TransNormal - CROSS product of trans orbit   km
*    IntNormal   - CROSS product of int orbit     km
*    R1Tgt       - Position vector after Dt, Tgt  km
*    V1Tgt       - Velocity vector after Dt, Tgt  km/s
*    RIRT        - RInt(4) * R1Tgt(4)
*    CosDeltaNu  - Cosine of DeltaNu              rad
*    SinDeltaNu  - Sine of DeltaNu                rad
*    DeltaNu     - ANGLE between position vectors rad
*    i           - Index
*
*  Coupling      :
*    KEPLER      - Find R and V at future time
*    LAMBERTUNIV - Find velocity vectors at each ENDIF of transfer
*    LNCOM2      - Linear combination of two vectors and constants
*
*  References    :
*    Vallado       2001, 468-474, Alg 58
*
* ------------------------------------------------------------------------------

      SUBROUTINE TARGET      ( RInt,VInt,RTgt,VTgt, Dm,Kind, Dtsec,mu,
     &                         V1t,V2t,DV1,DV2, Error  )
        IMPLICIT NONE
        REAL*8 RInt(3),VInt(3),RTgt(3),VTgt(3),Dtsec,mu,V1t(3),V2t(3),
     &         DV1(3),DV2(3)
        CHARACTER Kind, Dm
        CHARACTER*12 Error
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 R1Tgt(3), V1Tgt(3)

        ! --------------------  Implementation   ----------------------
        ! ----------- Propogate TARGET forward by time ----------------
        IF (Kind.eq.'K') THEN
            CALL KEPLER ( RTgt,VTgt,Dtsec,  R1Tgt,V1Tgt,Error )
          ENDIF
*        IF (Kind.eq.'P') THEN
*            CALL PKEPLER( RTgt,VTgt,Dtsec,  R1Tgt,V1Tgt )
*          ENDIF

        ! ----------- Calculate transfer orbit between r's ------------
        IF ( Error .eq. 'ok' ) THEN
            CALL LAMBERTUNIV( RInt,R1Tgt,dm,'N',Dtsec,mu,V1t,V2t,Error )

            IF ( Error .eq. 'ok' ) THEN
                CALL LNCOM2( -1.0D0, 1.0D0,VInt, V1t,  DV1 )
                CALL LNCOM2(  1.0D0,-1.0D0,V1Tgt,V2t,  DV2 )
              ELSE
                DO i= 1 , 3
                    V1t(i)= 0.0D0
                    V2t(i)= 0.0D0
                    DV1(i)= 0.0D0
                    dV2(i)= 0.0D0
                  ENDDO
              ENDIF 
          ENDIF 
      RETURN
      END
*
