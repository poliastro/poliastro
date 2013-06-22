*   -------------------------------------------------------------------
*
*                              ASTTIME.FOR
*
*   This file contains fundamental Astrodynamic procedures and functions
*   relating to the time functions. These routines are discussed in Ch 3
*   and Ch 5.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                   2007
*                             by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*               7 sep 07  david vallado
*                           fix some headers
*    changes :
*              15 mar 07  david vallado
*                           3rd edition baseline
*              21 jul 05  david vallado
*                           2nd printing baseline
*              28 Jan 04  David Vallado
*                           Update headers                                    
*              14 Mar 03  David Vallado
*                           Fixes to IDINT
*              28 Feb 03  David Vallado
*                           New baseline
*              14 May 01  David Vallado
*                           2nd edition baseline
*              23 Nov 87  David Vallado
*                           Original baseline
*
*   -------------------------------------------------------------------
*
*     Uses object files:
*         Astmath
*     Uses common files:
*         Astmath.cmn
*
*
*
*      SUBROUTINE INITTIME
*
*      INTEGER FUNCTION GETINTMON  ( MonStr )
*
*      INTEGER FUNCTION GETINTDAY  ( DayStr )
*
*      INTEGER FUNCTION DAYOFWEEK  ( JD )
*
*      SUBROUTINE DAYLIGHTST  ( Year, StartDay, StopDay, JDStartDST,
*     &                         JDStopDST )
*
*      SUBROUTINE JDay        ( Year,Mon,Day,Hr,minute, Sec, JD )
*
*      SUBROUTINE JDayALL     ( Year,Mon,Day,Hr,minute,Sec, WhichType, JD )
*
*      SUBROUTINE DAY2SMDHMS  ( Year,Days,  Mon,Day,Hr,minute,Sec )
*
*      SUBROUTINE INVJDay     ( JD, Year,Mon,Day,Hr,minute, Sec )
*
*      SUBROUTINE FINDDAYS    ( Year,Month,Day,Hr,minute, Sec,  Days )
*
*      SUBROUTINE LSTIME      ( Lon, JD, LST, GST )
*
*      SUBROUTINE SUNRISESET  ( JD,Latgd,Lon, WhichKind, UTSunRise,
*     &                         UTSunSet, Error )
*
*      SUBROUTINE MOONRISESET ( JD,Latgd,Lon, UTMoonRise, UTMoonSet,
*     &                         MoonPhaseAng, Error )
*
*      SUBROUTINE HMS_SEC     ( Hr,minute, Sec, Direction, UTSec )
*
*      SUBROUTINE HMS_UT      ( Hr,minute, Sec, Direction, UT )
*
*      SUBROUTINE HMS_RAD     ( Hr,minute, Sec, Direction, HMS )
*
*      SUBROUTINE DMS_RAD     ( Deg,minute, Sec, Direction, DMS )
*
*      SUBROUTINE jd2sse      ( jd,Direction, sse )
*
*      SUBROUTINE CONVTIME    ( Year, Mon, Day, Hr, minute, SEC,
*     &                         TimeZone, TypeUTIn, DUT1, DAT, xp, yp,
*     &                         UT1, TUT1, JDUT1, UTC, TAI, TT, TTT,
*     &                         JDTT, TDB, TTDB, JDTDB, DDPSi, DDEps,
*     &                         LOD, Error )
*
      SUBROUTINE INITTIME
        IMPLICIT NONE
* ----------------------------  Locals  -------------------------------
        INTEGER i, LMonth(12)
        CHARACTER*3 MonthTitle(12), DayTitle(12)

        ! --------------------  Implementation   ----------------------
        DO i = 1,12
            LMonth(i) = 31
          ENDDO
        LMonth( 2) = 28
        LMonth( 4) = 30
        LMonth( 6) = 30
        LMonth( 9) = 30
        LMonth(11) = 30

        MonthTitle( 1)= 'Jan'
        MonthTitle( 2)= 'Feb'
        MonthTitle( 3)= 'Mar'
        MonthTitle( 4)= 'Apr'
        MonthTitle( 5)= 'May'
        MonthTitle( 6)= 'Jun'
        MonthTitle( 7)= 'Jul'
        MonthTitle( 8)= 'Aug'
        MonthTitle( 9)= 'Sep'
        MonthTitle(10)= 'Oct'
        MonthTitle(11)= 'Nov'
        MonthTitle(12)= 'Dec'

        DayTitle(1)= 'Sun'
        DayTitle(2)= 'Mon'
        DayTitle(3)= 'Tue'
        DayTitle(4)= 'Wed'
        DayTitle(5)= 'Thr'
        DayTitle(6)= 'Fri'
        DayTitle(7)= 'Sat'

      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           FUNCTION GETINTMON
*
*  this function finds the INTEGER equivalent of the 3 character string
*    representation of month.
*
*  Inputs          Description                    Range / Units
*    MonStr      - Month name                     'Jan','Feb' ...
*
*  OutPuts       :
*    GETINTMON   - INTEGER Month equivalent       1 .. 12
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
* -----------------------------------------------------------------------------

      INTEGER FUNCTION GETINTMON       ( MonStr )
        CHARACTER*3 MonStr
* ----------------------------  Locals  -------------------------------
        INTEGER i
        CHARACTER*3 MonthTitle(12),MonthlTitle(12)

        ! --------------------  Implementation   ----------------------
        MonthTitle( 1)= 'JAN'
        MonthTitle( 2)= 'FEB'
        MonthTitle( 3)= 'MAR'
        MonthTitle( 4)= 'APR'
        MonthTitle( 5)= 'MAY'
        MonthTitle( 6)= 'JUN'
        MonthTitle( 7)= 'JUL'
        MonthTitle( 8)= 'AUG'
        MonthTitle( 9)= 'SEP'
        MonthTitle(10)= 'OCT'
        MonthTitle(11)= 'NOV'
        MonthTitle(12)= 'DEC'

        MonthlTitle( 1)= 'Jan'
        MonthlTitle( 2)= 'Feb'
        MonthlTitle( 3)= 'Mar'
        MonthlTitle( 4)= 'Apr'
        MonthlTitle( 5)= 'May'
        MonthlTitle( 6)= 'Jun'
        MonthlTitle( 7)= 'Jul'
        MonthlTitle( 8)= 'Aug'
        MonthlTitle( 9)= 'Sep'
        MonthlTitle(10)= 'Oct'
        MonthlTitle(11)= 'Nov'
        MonthlTitle(12)= 'Dec'

        i = 1
        DO WHILE (( MonthTitle(i).ne.MonStr ).and.( i .lt. 12 ))
            i= i+1
          ENDDO

        IF ( (i.eq.12) .and. (MonthTitle(i).ne.MonStr) ) THEN
            i = 1
            DO WHILE (( MonthlTitle(i).ne.MonStr ).and.( i .lt. 12 ))
                i= i+1
              ENDDO
          ENDIF

        GETINTMON= i
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           FUNCTION GETINTDAY
*
*  this function finds the INTEGER equivalent of the 3 character string
*    representation of the day of the week.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    DayStr      - Day name string                'Sun','Mon' ...
*
*  OutPuts       :
*    GETINTDAY   - INTEGER Day equivalent         1 .. 7
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
* -----------------------------------------------------------------------------

      INTEGER FUNCTION GETINTDAY       ( DayStr )
        CHARACTER*3 DayStr
* ----------------------------  Locals  -------------------------------
        INTEGER i
        CHARACTER*3 DayTitle(12)

        ! --------------------  Implementation   ----------------------
        DayTitle(1)= 'Sun'
        DayTitle(2)= 'Mon'
        DayTitle(3)= 'Tue'
        DayTitle(4)= 'Wed'
        DayTitle(5)= 'Thr'
        DayTitle(6)= 'Fri'
        DayTitle(7)= 'Sat'

        i= 1
        DO WHILE (( (DayTitle(i)) .ne. (DayStr) ) .and.
     &           ( i .lt. 12 ))
            i= i+1
          ENDDO
        GETINTDAY= i
      RETURN
      END

* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DAYOFWEEK
*
*  this function finds the day of the week. Integers are used for the days,
*    1 = 'SUN', 2 = 'Mon', ... 7 = 'Sat'.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian date of interest        days from 4713 BC
*
*  OutPuts       :
*    DAYOFWEEK   - answer                         1 to 7
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
*  References    :
*    vallado       2007, 188, Eq 3-39
*
* -----------------------------------------------------------------------------

      INTEGER FUNCTION DAYOFWEEK       ( JD )
        REAL*8 JD

        ! --------------------  Implementation   ----------------------
        ! ------- Be sure JD is at 0.0D0 h on the day of interest -----
        JD       = DINT(JD+0.5D0)
        DAYOFWEEK= IDINT( JD - 7 * IDINT( (JD+1)/7 ) + 2 )
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DAYLIGHTST
*
*  this subroutine finds the dates for switiching to daylight savings time in
*    a given year. The date is set as the 2nd sunday in march and the first
*    sunday in november. The DST dates are adjusted -10hr to get 0200, -Zone to get
*    the local time zone, and -1.0hr to process the stop because the local time is
*    on DST before a stop.
*
*  Author        : David Vallado                  719-573-2600   17 Mar 2007
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1900 .. 2100
*    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
*
*  Outputs       :
*    StartDay    - Day in April when DST begins   1 .. 28,29,30,31
*    StopDay     - Day in October when DST ends   1 .. 28,29,30,31
*
*  Locals        :
*    DW          - Day of the week                1 .. 7
*    JDStartDST  - Julian date of start           Days from 4713 BC
*    JDStopDST   - Julian date of stop            Days from 4713 BC
*    Zone        - Time zone of site. Default
*                  of 0.0 gives Greenwich         hrs
*
*  Coupling      :
*    JDAY   - Find the Julian Date
*
*  References    :
*    Vallado       2007, 188
*
* -----------------------------------------------------------------------------

      SUBROUTINE DAYLIGHTST ( Year, Lon, StartDay, StopDay, JDStartDST,
     &                        JDStopDST )
        IMPLICIT NONE
        INTEGER Year, StartDay, StopDay
        REAL*8  Lon, JDStartDST, JDStopDST
* ----------------------------  Locals  -------------------------------
        REAL*8 Zone
        INTEGER DW
        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        ! ------- Find time zone information to adjust to a site ------ 
        Zone= DINT( Lon*Rad2Deg/15.0D0 )
        IF (Zone .gt. 0.0D0) THEN
            Zone= Zone - 24.0D0
          ENDIF

        IF (year .lt. 2007) THEN
             StartDay= 0
            Dw = 0
            DO WHILE ((DW .ne. 1) .and. (StartDay .ne. 7)) ! 1 is Sunday
                StartDay= StartDay + 1
                CALL JDay( Year,4,StartDay,12,0,0.0D0, JDStartDST )
                DW= IDINT( JDStartDST - 7* IDINT( (JDStartDST+1)/7 )+2)
              ENDDO
            JDStartDST= JDStartDST - 10.0D0/24.0D0   ! set to 0200 UTC

            StopDay= 32
            DO WHILE ((DW .ne. 1) .and. (StopDay .ne. 25)) ! 1 is Sunday
                StopDay= StopDay - 1
                CALL JDay( Year,10,StopDay,12,0,0.0D0, JDStopDST )
                DW= IDINT( JDStopDST - 7* IDINT( (JDStopDST+1)/7 ) + 2 )
              ENDDO
            JDStopDST= JDStopDST - 10.0D0/24.0D0   ! set to 0200 UTC
          ELSE
           StartDay = 7
            Dw = 0
            DO WHILE ((DW .ne. 1) .and. (StartDay .ne. 15)) ! 1 is Sunday
                StartDay= StartDay + 1
                CALL JDay( Year,3,StartDay,12,0,0.0D0, JDStartDST )
                DW= IDINT( JDStartDST - 7* IDINT( (JDStartDST+1)/7 )+2)
              ENDDO
            JDStartDST= JDStartDST - 10.0D0/24.0D0   ! set to 0200 UTC

            StopDay = 0
            DO WHILE ((DW .ne. 1) .and. (StopDay .ne. 8)) ! 1 is Sunday
                StopDay= StopDay + 1
                CALL JDay( Year,11,StopDay,12,0,0.0D0, JDStopDST )
                DW= IDINT( JDStopDST - 7* IDINT( (JDStopDST+1)/7 ) + 2 )
              ENDDO
            JDStopDST= JDStopDST - 10.0D0/24.0D0   ! set to 0200 UTC
          ENDIF

      RETURN
      END
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE JDay
*
*  this subroutine finds the Julian date given the Year, Month, Day, and Time.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Universal Time Hour            0 .. 23
*    minute         - Universal Time minute             0 .. 59
*    Sec         - Universal Time Sec             0.0D0 .. 59.999D0
*    WhichType   - Julian .or. Gregorian calender   'J' .or. 'G'
*
*  Outputs       :
*    JD          - Julian Date                    days from 4713 BC
*
*  Locals        :
*    B           - Var to aid Gregorian dates
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 189, Alg 14, Ex 3-14
* -----------------------------------------------------------------------------

      SUBROUTINE JDay        ( Year,Mon,Day,Hr,minute, Sec, JD )
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, minute
        REAL*8  Sec, JD

        ! --------------------  Implementation   ----------------------
        JD= 367.0D0 * Year
     &        - INT( (7* (Year+INT ( (Mon+9)/12) ) ) * 0.25D0 )
     &        + INT( 275*Mon / 9 )
     &        + Day + 1721013.5D0
     &        + ( (Sec/60.0D0 + minute ) / 60.0D0 + Hr ) / 24.0D0
*     &      - 0.5D0*DSIGN(1.0D0, 100.0D0*Year + Mon - 190002.5D0) + 0.5D0
      RETURN
      END

*
*  References    :
*    Vallado       2007, 190

      SUBROUTINE JDayALL   ( Year,Mon,Day,Hr,minute,Sec, WhichType, JD )
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, minute
        CHARACTER WhichType
        REAL*8  Sec, JD
* ----------------------------  Locals  -------------------------------
        REAL*8 B

        ! --------------------  Implementation   ----------------------
        IF ( Mon .le. 2 ) THEN
            Year= Year - 1
            Mon = Mon + 12
          ENDIF
        IF ( WhichType .eq. 'J' ) THEN
            ! -------- Use for Julian calender, every 4 years ---------
            B= 0.0D0
          ELSE
            ! --------------------- Use for Gregorian -----------------
            B= 2 - DINT(Year*0.01D0) + DINT(DINT(Year*0.01D0)*0.25D0)
          ENDIF
        JD= DINT( 365.25D0*(Year + 4716) )
     &       + DINT( 30.6001D0*(Mon+1) )
     &       + Day + B - 1524.5D0
     &       + ( (Sec/60.0D0 + minute ) / 60.0D0 + Hr ) / 24.0D0
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DAYS2MDHMS
*
*  this subroutine converts the day of the year, days, to the equivalent month
*    day, hour, Minute and second.
*
*  Algorithm     : Set up array for the Number of days per month
*                  Find Leap Year - be sure to account for the 400 years
*                  Loop through a Temp value for WHILE the value is .lt. the days
*                  Perform INTEGER conversions to the correct day and month
*                  Convert remainder into H M S using type conversions
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                          +1900 .. 2100+
*    Days        - Julian Day of the year         0.0D0  .. 366.0D0
*
*  OutPuts       :
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*
*  Locals        :
*    DayofYr     - Day of year
*    Temp        - Temporary REAL*8 values
*    IntTemp     - Temporary INTEGER value
*    i           - Index
*    LMonth[12]  - INTEGER Array containing the Number of days per month
*
*  Coupling      :
*    None.
* -----------------------------------------------------------------------------

      SUBROUTINE DAYS2MDHMS  ( Year,Days,  Mon,Day,Hr,minute,Sec )
        IMPLICIT NONE
        REAL*8 Days,Sec
        INTEGER Year, Mon, Day, Hr, minute
* ----------------------------  Locals  -------------------------------
        INTEGER IntTemp,i,DayofYr, LMonth(12)
        REAL*8 Temp

        ! --------------------  Implementation   ----------------------
        ! -------------- Set up array of days in month  ---------------
        DO i = 1,12
            LMonth(i) = 31
          ENDDO
        LMonth( 2) = 28
        LMonth( 4) = 30
        LMonth( 6) = 30
        LMonth( 9) = 30
        LMonth(11) = 30

        DayofYr= IDINT(Days )

        ! ---------------- Find month and Day of month ----------------
        IF (MOD(Year,4).eq.0) THEN
            LMonth(2)= 29
          ENDIF
        i= 1
        IntTemp= 0
        DO WHILE ( (DayofYr.gt.IntTemp + LMonth(i) ) .and. ( i.lt.12 ))
            IntTemp= IntTemp + LMonth(i)
            i= i+1
          ENDDO
        Mon= i
        Day= DayofYr - IntTemp

        ! ---------------- Find hours Minutes and seconds -------------
        Temp= (Days - DayofYr )*24.0D0
        Hr  = IDINT( Temp )
        Temp= (Temp-Hr) * 60.0D0
        minute = IDINT( Temp )
        Sec = (Temp-minute) * 60.0D0

        ! ---- Check for roundoff errors
c        IF (Sec .ge. 59.9999D0) THEN
c            Sec = 0.0D0
c            minute = minute + 1
c            IF (minute .gt. 59) THEN
c                minute = 0
c                Hr = Hr + 1
c                IF (Hr .gt. 23) THEN
c                    Hr = 0
c                    Day = Day + 1
c                  ENDIF
c              ENDIF
c          ENDIF
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE INVJDay
*
*  this subroutine finds the Year, month, day, hour, Minute and second
*  given the Julian date. TU can be UT1, TDT, TDB, etc.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*
*  Locals        :
*    Days        - Day of year plus fractional
*                  portion of a day               days
*    Tu          - Julian Centuries from 0 h
*                  Jan 0, 1900
*    Temp        - Temporary real values
*    LeapYrs     - Number of Leap years from 1900
*
*  Coupling      :
*    DAYS2MDHMS  - Finds MD HMS given Days and Year
*
*  References    :
*    Vallado       2007, 208, Alg 22, Ex 3-13
* -----------------------------------------------------------------------------

      SUBROUTINE INVJDay     ( JD, Year,Mon,Day,Hr,minute, Sec )
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, minute
        REAL*8  Sec, JD
* ----------------------------  Locals  -------------------------------
        INTEGER LeapYrs
        REAL*8  Days, Tu, Temp

        ! --------------------  Implementation   ----------------------
        ! ---------------- Find Year and Days of the year -------------
        Temp   = JD-2415019.5D0
        Tu     = Temp / 365.25D0
        Year   = 1900 + IDINT( Tu )
        LeapYrs= IDINT( ( Year-1901 )*0.25D0 )
        Days   = Temp - ((Year-1900)*365.0D0 + LeapYrs )

        ! -------------- Check for case of beginning of a year --------
        IF ( Days .lt. 1.0D0 ) THEN
            Year   = Year - 1
            LeapYrs= IDINT( ( Year-1901 )*0.25D0 )
            Days   = Temp - ((Year-1900)*365.0D0 + LeapYrs )
          ENDIF

        ! ------------------ Find remaing data  -----------------------
        CALL DAYS2MDHMS( Year,Days, Mon,Day,Hr,minute,Sec )

      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE FINDDAYS
*
*  this subroutine finds the fractional days through a year given the year,
*    month, day, hour, Minute and second.
*
*  Algorithm     : Set up array for the Number of days per month
*                  Find Leap Year - be sure to account for the 400 years 
*                  Check for a leap year
*                  Loop to find the elapsed days in the year
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*
*  OutPuts       :
*    Days        - Day of year plus fraction of a
*                    day                          days
*
*  Locals        :
*    LMonth      - Length of months of year
*    i           - Index
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 207, Ex 3-12
*
* -----------------------------------------------------------------------------

      SUBROUTINE FINDDAYS    ( Year,Month,Day,Hr,minute, Sec,  Days )
        IMPLICIT NONE
        INTEGER Year, Month, Day, Hr, minute
        REAL*8 Sec, Days
* ----------------------------  Locals  -------------------------------
        INTEGER i, LMonth(12)

        ! --------------------  Implementation   ----------------------
        DO i = 1,12
            LMonth(i) = 31
          ENDDO
        LMonth( 2) = 28
        LMonth( 4) = 30
        LMonth( 6) = 30
        LMonth( 9) = 30
        LMonth(11) = 30
        IF (MOD(Year,4).eq.0) THEN
            LMonth(2)= 29
          ENDIF

        i   = 1
        Days= 0.0D0
        DO WHILE ((i .lt. Month) .and. ( i .lt. 12 ))
            Days= Days + LMonth(i)
            i= i + 1
          ENDDO

        Days= Days + Day + Hr/24.0D0 + minute/1440.0D0 + Sec/86400.0D0
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           FUNCTION GSTIME
*
*  this function finds the Greenwich sidereal time (iau-82). 
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    GSTIME      - Greenwich SIDEREAL Time        0 to 2Pi rad
*
*  Locals        :
*    Temp        - Temporary variable for reals   rad
*    TUT1        - Julian Centuries from the
*                  Jan 1, 2000 12 h epoch (UT1)
*
*  Coupling      :
*
*  References    :
*    vallado       2007, 193, Eq 3-43
* -----------------------------------------------------------------------------

      REAL*8 FUNCTION GSTIME ( JD )
        IMPLICIT NONE
        REAL*8 JD
* ----------------------------  Locals  -------------------------------
        REAL*8 Temp, TUT1

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------

        TUT1= ( JD - 2451545.0D0 ) / 36525.0D0
        Temp= - 6.2D-6*TUT1*TUT1*TUT1
     &        + 0.093104D0*TUT1*TUT1
     &        + (876600.0D0*3600.0D0 + 8640184.812866D0)*TUT1
     &        + 67310.54841D0
        Temp= DMOD( Temp*Deg2Rad/240.0D0,TwoPi ) ! 360/86400 = 1/240, to deg, to rad

        ! ------------------------ Check quadrants --------------------
        IF ( Temp .lt. 0.0D0 ) THEN
            Temp= Temp + TwoPi
          ENDIF

        GSTIME= Temp

      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           FUNCTION GSTIM0
*
*  this function finds the Greenwich sidereal time at the beginning of a year.
*    This formula is derived from the Astronomical Almanac and is good only
*    0 hr UT1, Jan 1 of a year.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1998, 1999, etc.
*
*  OutPuts       :
*    GSTIM0      - Greenwich sidereal Time        0 to 2Pi rad
*
*  Locals        :
*    JD          - Julian Date                    days from 4713 BC
*    Temp        - Temporary variable for Reals   rad
*    TUT1        - Julian Centuries from the
*                  Jan 1, 2000 12 h epoch (UT1)
*
*  Coupling      :
*
*  References    :
*    Vallado       2007, 195, Eq 3-46
*
* -----------------------------------------------------------------------------

      REAL*8 FUNCTION GSTIM0 ( Year )
        IMPLICIT NONE
        INTEGER Year
* ----------------------------  Locals  -------------------------------
        REAL*8 JD, Temp, TUT1

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        JD  = 367.0D0 * Year - ( INT((7*(Year+INT(10/12)))*0.25D0) )+
     &                           ( INT(275/9) ) + 1721014.5D0
        TUT1= ( JD - 2451545.0D0 ) / 36525.0D0

        Temp= - 6.2D-6*TUT1*TUT1*TUT1
     &        + 0.093104D0*TUT1*TUT1
     &        + (876600.0D0*3600.0D0 + 8640184.812866D0)*TUT1
     &        + 67310.54841D0

        ! ------------------------ Check quadrants --------------------
        Temp= DMOD( Temp,TwoPi )
        IF ( Temp .lt. 0.0D0 ) THEN
            Temp= Temp + TwoPi
          ENDIF

        GSTIM0= Temp

      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE LSTIME
*
*  this subroutine finds the Local sidereal time at a given location.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    LST         - Local sidereal Time            0.0D0 to 2Pi rad
*    GST         - Greenwich sidereal Time        0.0D0 to 2Pi rad
*
*  Locals        :
*    None.
*
*  Coupling      :
*    GSTIME        Finds the Greenwich sidereal Time
*
*  References    :
*    Vallado       2007, 194, alg 15, ex 3-5
*
* -----------------------------------------------------------------------------

      SUBROUTINE LSTIME      ( Lon, JD, LST, GST )
        IMPLICIT NONE
        REAL*8 Lon, JD, LST, GST

* ----------------------------  Locals  -------------------------------
        EXTERNAL GSTIME
        REAL*8 GSTIME

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        GST = GSTIME( JD )
        LST = Lon + GST

        ! ----------------------- Check quadrants ---------------------
        LST = DMOD( LST,TwoPi )
        IF ( LST .lt. 0.0D0 ) THEN
            LST= LST + TwoPi
          ENDIF

      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE SUNRISESET
*
*  this subroutine finds the Universal time for Sunrise and Sunset given the
*    day and SITE location.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*    Latgd       - SITE latitude (SOUTH -)        -65 to 65 rad (deg?)
*    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
*    WhichKind   - Character for which rise/set   'S' 'C' 'N' 'A'
*
*  OutPuts       :
*    UTSunRise   - Universal time of sunrise      hrs
*    UTSunSet    - Universal time of sunset       hrs
*    Error       - Error Parameter
*
*  Locals        :
*    SunAngle    - ANGLE between the SUN vector
*                  and a point on the Earth     rad
*    JDTemp      - Julian date for sunrise/set    days from 4713 BC
*    UTTemp      - Temporary UT time              days
*    TUT1        - Julian Centuries from the
*                  Jan 1, 2000 12 h epoch (UT1)
*    Ra          - Right ascension                rad
*    Decl        - Declination                    rad
*    MeanLonSun  -                                rad
*    MeanAnomalySun                               rad
*    LonEcliptic - Longitude of the ecliptic      rad
*    Obliquity   - Obliquity of the ecliptic      rad
*    GST         - for 0 h UTC of each day        rad
*    LHA         - Local hour ANGLE               rad
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*    Opt         - Idx to for rise and set calc    1,2
*
*  Coupling      :
*    INVJDay- Finds the Year day mon hr minute Sec from the Julian Date
*    JDay   - Finds the Julian date given Year, mon day, hr, minute, Sec
*
*  References    :
*    Vallado       2007, 283, Alg 30, Ex 5-2
*
* -----------------------------------------------------------------------------

      SUBROUTINE SUNRISESET  ( JD,Latgd,Lon, WhichKind, UTSunRise,
     &                         UTSunSet, Error )
        IMPLICIT NONE
        REAL*8 JD, Latgd, Lon, UTSunRise, UTSunSet
        CHARACTER WhichKind
        CHARACTER*12 Error
* ----------------------------  Locals  -------------------------------
        INTEGER Opt, Year, Month, Day, Hr, minute
        REAL*8 MeanAnomalySun,
     &         JDTemp, UTTemp, SunAngle, TUT1, Ra, Sec, MeanLonSun,
     &         LonEcliptic, Decl, Obliquity, GST, LHA

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Error= 'ok'

        ! -------------- Make sure lon is within +- 180 deg -----------
        IF ( Lon .gt. Pi ) THEN
            Lon= Lon - 2.0D0*Pi
          ENDIF
        IF ( Lon .lt. -Pi ) THEN
            Lon= Lon + 2.0D0*Pi
          ENDIF
*
        IF (WhichKind .eq. 'S') THEN
            SunAngle= (90.0D0+50.0D0/60.0D0 )*Deg2Rad
          ENDIF
        IF (WhichKind .eq. 'C') THEN
            SunAngle=  96.0D0 *Deg2Rad
          ENDIF
        IF (WhichKind .eq. 'N') THEN
            SunAngle= 102.0D0 *Deg2Rad
          ENDIF
        IF (WhichKind .eq. 'A') THEN
            SunAngle= 108.0D0 *Deg2Rad
          ENDIF
        CALL INVJDay( JD,  Year,Month,Day,Hr,minute,Sec )
        DO Opt= 1 , 2
            IF ( Opt .eq. 1 ) THEN
                CALL JDay( Year,Month,Day, 6,0,0.0D0, JDTemp )
              ELSE
                CALL JDay( Year,Month,Day,18,0,0.0D0, JDTemp )
              ENDIF
            JDTemp= JDTemp - Lon*Rad2Deg/15.0D0/24.0D0

            TUT1 = (JDTemp - 2451545.0D0)/36525.0D0 
            MeanLonSun = 280.4606184D0 + 36000.77005361D0*TUT1
            MeanAnomalySun= 357.5277233D0+35999.05034D0*TUT1
            MeanAnomalySun= DMOD( MeanAnomalySun*Deg2Rad,TwoPi )
            IF ( MeanAnomalySun .lt. 0.0D0 ) THEN
                MeanAnomalySun= MeanAnomalySun + TwoPi
              ENDIF
            LonEcliptic= MeanLonSun + 1.914666471D0*DSIN(MeanAnomalySun)
     &                     + 0.019994643D0*DSIN(2.0D0*MeanAnomalySun)
            LonEcliptic= DMOD( LonEcliptic*Deg2Rad,TwoPi )
            IF ( LonEcliptic .lt. 0.0D0 ) THEN
                LonEcliptic= LonEcliptic + TwoPi
              ENDIF
            Obliquity= 23.439291D0 - 0.0130042D0*TUT1
            Obliquity= Obliquity *Deg2Rad
            Ra  = DATAN( DCOS(Obliquity) * DTAN(LonEcliptic) )
            Decl= DASIN( DSIN(Obliquity) * DSIN(LonEcliptic) )
            IF ( Ra .lt. 0.0D0 ) THEN
                Ra= Ra + TwoPi
              ENDIF
            IF ( (LonEcliptic .gt. Pi) .and. (Ra .lt. Pi) ) THEN
                Ra= Ra + Pi
              ENDIF
            IF ( (LonEcliptic .lt. Pi) .and. (Ra .gt. Pi) ) THEN
                Ra= Ra - Pi
              ENDIF
            LHA= (DCOS(SunAngle) - DSIN(Decl)*DSIN(Latgd)) /
     &           (DCOS(Decl)*DCOS(Latgd) )
            IF ( DABS(LHA) .le. 1.0D0 ) THEN
                LHA= DACOS( LHA )
              ELSE
                Error= 'Not ok'
              ENDIF
            IF ( Error .eq. 'ok' ) THEN
                IF ( Opt .eq. 1 ) THEN
                    LHA= TwoPi - LHA
                  ENDIF
                GST= 1.75336855923327D0 + 628.331970688841D0*TUT1
     &                 + 6.77071394490334D-06*TUT1*TUT1
     &                 - 4.50876723431868D-10*TUT1*TUT1*TUT1
                GST= DMOD( GST,TwoPi )
                IF ( GST .lt. 0.0D0 ) THEN
                    GST= GST + TwoPi
                  ENDIF
                UTTemp= LHA + Ra  - GST
                UTTemp= UTTemp * Rad2Deg/15.0D0
                UTTemp= DMOD( UTTemp,24.0D0 )
                UTTemp= UTTemp - Lon*Rad2Deg/15.0D0
                IF ( UTTemp .lt. 0.0D0 ) THEN
                    UTTemp= UTTemp + 24.0D0
                    Error= 'Day before'
                  ENDIF
                IF ( UTTemp .gt. 24.0D0 ) THEN
                    UTTemp= UTTemp - 24.0D0
                    Error= 'Day After'
                  ENDIF
              ELSE
                UTTemp= 99.99D0
              ENDIF
            IF ( Opt .eq. 1 ) THEN
                UTSunRise= UTTemp
              ELSE
                UTSunSet = UTTemp
              ENDIF
          ENDDO
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE MOONRISESET
*
*  this subroutine finds the Universal time for Moonrise and Moonset given the
*    day and SITE location.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*    Latgd       - SITE latitude (SOUTH -)        -65 to 65 rad (deg?)
*    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
*
*  OutPuts       :
*    UTMoonRise  - Universal time of Moonrise     hrs
*    UTMoonSet   - Universal time of Moonset      hrs
*    MoonPhaseAng- Phase angle of the Moon        deg
*    Error       - Error Parameter
*
*  Locals        :
*    MoonAngle   - ANGLE between the Moon vector
*                  and a point on the Earth       rad
*    JDTemp      - Julian date for Moonrise/set   days from 4713 BC
*    UTTemp      - Temporary UT time              days
*    TUT1        - Julian Centuries from the
*                  Jan 1, 2000 12 h epoch (UT1)
*    RtAsc       - Right ascension                rad
*    Decl        - Declination                    rad
*    MeanLonMoon -                                rad
*    MeanAnomaly -                                rad
*    EclpLong    - Longitude of the ecliptic      rad
*    Obliquity   - Obliquity of the ecliptic      rad
*    RMoon
*    RMoonRS
*    RV
*    RhoSat
*    Try
*    l, m, n     - Direction cosines
*    EclpLat
*    MoonGHA, MoonGHAn
*    DGHA, DGHAn
*    LHAn
*    LST
*    DeltaUT, DeltaUTn
*    t, tn
*    HzParal
*    LonEclSun
*    LonEclMoon
*    TTDB
*    GST         - for 0 h UTC of each day        rad
*    LHA         - Local hour ANGLE               rad
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*    Opt         - Idx to for rise and set calc    1,2
*
*  Coupling      :
*    INVJDay- Finds the Year day mon hr minute Sec from the Julian Date
*    JDay   - Finds the Julian date given Year, mon day, hr, minute, Sec
*
*  References    :
*    Vallado       2007, 292, Alg 32, Ex 5-4
*
* -----------------------------------------------------------------------------

      SUBROUTINE MOONRISESET ( JD,Latgd,Lon, UTMoonRise, UTMoonSet,
     &                         MoonPhaseAng, Error )
        IMPLICIT NONE
        REAL*8 JD, Latgd, Lon, UTMoonRise, UTMoonSet, MoonPhaseAng
        CHARACTER*12 Error
* ----------------------------  Locals  -------------------------------
        INTEGER Opt, i, Year, Month, Day, Hr, minute, Try
        REAL*8 DeltaUT, tn, GST, t,
     &    l,m,n, EclpLong, EclpLat, Obliquity, MoonGHA, DGHA, LHAn,
     &    MoonGHAn, LHA, LST,
     &    LonEclSun, LonEclMoon, MeanAnomaly, MeanLong, ttdb,
     &    Sec, JDTemp, UTTemp, RtAsc, Decl

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Error = 'ok'

        ! ---------- for once for MoonRise (1), ) THEN set (2) --------
        ! -------------- Make sure lon is within +- 180 deg -----------
        IF ( Lon .gt. Pi ) THEN
            Lon= Lon - 2.0D0*Pi
          ENDIF
        IF ( Lon .lt. -Pi ) THEN
            Lon= Lon + 2.0D0*Pi
          ENDIF

        Try= 1
        Opt= 1
        DO WHILE (Opt .le. 2)
            CALL INVJDay( JD,  Year,Month,Day,Hr,minute,Sec )
            CALL JDay( Year,Month,Day,0,0,0.0D0, JDTemp )
            UTTemp= 0.5D0

            IF ( Try .eq. 2 ) THEN
                IF ( Opt .eq. 1 ) THEN
                    UTTemp= 0.25D0
                  ELSE
                    UTTemp= 0.75D0
                  ENDIF
              ENDIF

            i = 0
            tn= UTTemp
            t = tn + 10.0D0
            JDTemp= JDTemp + UTTemp

            DO WHILE ( (DABS(tn-t).ge.0.008D0) .and. (i .le. 5) )
                TTDB = ( JDTemp - 2451545.0D0 ) / 36525.0D0
                EclpLong= 218.32D0 + 481267.883D0*TTDB
     &              + 6.29D0*DSIN( (134.9D0+477198.85D0*TTDB)*Deg2Rad )
     &              - 1.27D0*DSIN( (259.2D0-413335.38D0*TTDB)*Deg2Rad )
     &              + 0.66D0*DSIN( (235.7D0+890534.23D0*TTDB)*Deg2Rad )
     &              + 0.21D0*DSIN( (269.9D0+954397.70D0*TTDB)*Deg2Rad )
     &              - 0.19D0*DSIN( (357.5D0+ 35999.05D0*TTDB)*Deg2Rad )
     &              - 0.11D0*DSIN( (186.6D0+966404.05D0*TTDB)*Deg2Rad )
            EclpLat = 5.13D0*DSIN( ( 93.3D0+483202.03D0*TTDB)*Deg2Rad )
     &              + 0.28D0*DSIN( (228.2D0+960400.87D0*TTDB)*Deg2Rad )
     &              - 0.28D0*DSIN( (318.3D0+  6003.18D0*TTDB)*Deg2Rad )
     &              - 0.17D0*DSIN( (217.6D0-407332.20D0*TTDB)*Deg2Rad )
                EclpLong = DMOD( EclpLong*Deg2Rad, TwoPi )
                EclpLat  = DMOD( EclpLat*Deg2Rad, TwoPi )
                Obliquity= 23.439291D0 - 0.0130042D0*TTDB
                Obliquity= Obliquity *Deg2Rad
                ! ------- Find the geocentric direction cosines -------
                l= DCOS( EclpLat ) * DCOS( EclpLong )
                m= DCOS(Obliquity)*DCOS(EclpLat)*DSIN(EclpLong)
     &               - DSIN(Obliquity)*DSIN(EclpLat)
                n= DSIN(Obliquity)*DCOS(EclpLat)*DSIN(EclpLong)
     &               + DCOS(Obliquity)*DSIN(EclpLat)
                RtAsc= DATAN2( m,l )
                ! - Check that RtAsc is in the same quadrant as EclpLong
                IF ( EclpLong .lt. 0.0D0 ) THEN
                    EclpLong= EclpLong + TwoPi
                  ENDIF   
                IF ( DABS( EclpLong - RtAsc ) .gt. Pi*0.5D0 ) THEN
                    RtAsc= RtAsc + 0.5D0*Pi*
     &                     DINT( 0.5D0 + (EclpLong-RtAsc) / (0.5D0*Pi) )
                  ENDIF
                Decl = DASIN( n )
                CALL LSTIME( Lon,JDTemp,LST,GST )
                MoonGHAn= LST - Lon - RtAsc
                IF ( i .eq. 0 ) THEN
                    LHA = MoonGHAn + Lon
                    DGHA= 347.8D0 * Deg2Rad
                  ELSE
                    DGHA= (MoonGHAn - MoonGHA) / DeltaUT 
                  ENDIF
                IF ( DGHA .lt. 0.0D0 ) THEN
                    DGHA= DGHA + TwoPi/DABS(DeltaUT)
                  ENDIF
                LHAn= 0.00233D0 - (DSIN(Latgd)*DSIN(Decl)) /
     &                            (DCOS(Latgd)*DCOS(Decl))
                IF ( LHAn .gt. 1.0D0 ) THEN
                    LHAn= 0.0D0
                  ENDIF
                IF ( LHAn .lt. -1.0D0 ) THEN
                    LHAn= -1.0D0
                  ENDIF
                LHAn= DACOS( LHAn )
                IF ( Opt .eq. 1 ) THEN
                    LHAn= TwoPi - LHAn 
                  ENDIF 
                IF ( DABS( DGHA ) .gt. 0.0001D0 ) THEN
                    DeltaUT= (LHAn - LHA ) / DGHA
                  ELSE
                    DeltaUT= (LHAn - LHA )
                    DeltaUT= 1.0D0
                    Write( *,*)  'Fileout,x'
                  ENDIF
                t= tn 
                IF ( DABS( DeltaUT ) .gt. 0.5D0 ) THEN
                    IF ( DABS( DGHA ) .gt. 0.001D0 ) THEN
                        IF ( DeltaUT .lt. 0.0D0 ) THEN
                            DeltaUT= DeltaUT + TwoPi/DGHA
                            IF ( DABS( DeltaUT ) .gt. 0.51D0 ) THEN
                                i= 6 
                              ENDIF
                          ELSE
                            DeltaUT= DeltaUT - TwoPi/DGHA
                            IF ( DABS( DeltaUT ) .gt. 0.51D0 ) THEN
                                i= 6 
                              ENDIF
                          ENDIF
                      ELSE
                        DeltaUT= DeltaUT
                        Write(*,*) 'Fileout,y'
                      ENDIF 
                  ENDIF
                tn     = UTTemp + DeltaUT
                JDTemp = JDTemp - UTTemp + tn
                i = i + 1
                MoonGHA= MoonGHAn 

              ENDDO

            UTTemp= tn*24.0D0
            IF ( i .gt. 5 ) THEN
                UTTemp= 9999.99D0 
              ENDIF
            IF ( UTTemp .lt. 9999.0D0 ) THEN
                UTTemp= DMOD( UTTemp,24.0D0)
              ENDIF   
            IF ( UTTemp .lt. 0.0D0 ) THEN
                UTTemp= UTTEmp + 24.0D0
              ENDIF   
            IF ( UTTemp .gt. 900 ) THEN
                UTTemp= 24.0D0
              ENDIF

            IF (Opt .eq. 1 ) THEN
                UTMoonRise= UTTemp
              ENDIF
            IF (Opt .eq. 2 ) THEN
                UTMoonSet = UTTemp
              ENDIF

            Try= Try + 1 
            IF ( (i .gt. 5) .and. (Try .lt. 3) ) THEN
                Write(*,*) 'try #2 ',opt
              ELSE
                IF ( (i.gt.5) .and. (Try.gt.2) ) THEN
                    IF (Opt .eq. 1 ) THEN
                        Error = 'No Rise'
                      ENDIF
                    IF (Opt .eq. 2 ) THEN
                        Error = 'No Set'
                      ENDIF
                  ENDIF
                Opt= Opt + 1
                Try= 1
              ENDIF

          ENDDO

        ! ------------- determine phase ANGLE of the MOON --------------
        MeanLong= 280.4606184D0 + 36000.77005361D0*TTDB
        MeanLong= DMOD( MeanLong,360.0D0 )

        MeanAnomaly= 357.5277233D0 + 35999.05034D0*TTDB
        MeanAnomaly= DMOD( MeanAnomaly*Deg2Rad,TwoPi )
        IF ( MeanAnomaly .lt. 0.0D0 ) THEN
            MeanAnomaly= TwoPi + MeanAnomaly
          ENDIF

        LonEclSun= MeanLong + 1.914666471D0*DSIN(MeanAnomaly)
     &              + 0.019994643D0*DSIN(2.0D0*MeanAnomaly)

        LonEclMoon=   218.32D0 + 481267.883D0*TTDB
     &              + 6.29D0*DSIN( (134.9D0+477198.85D0*TTDB)*Deg2Rad )
     &              - 1.27D0*DSIN( (259.2D0-413335.38D0*TTDB)*Deg2Rad )
     &              + 0.66D0*DSIN( (235.7D0+890534.23D0*TTDB)*Deg2Rad )
     &              + 0.21D0*DSIN( (269.9D0+954397.70D0*TTDB)*Deg2Rad )
     &              - 0.19D0*DSIN( (357.5D0+ 35999.05D0*TTDB)*Deg2Rad )
     &              - 0.11D0*DSIN( (186.6D0+966404.05D0*TTDB)*Deg2Rad )
        LonEclMoon= DMOD( LonEclMoon, 360.0D0 )

        MoonPhaseAng= LonEclMoon - LonEclSun

        IF ( MoonPhaseAng .lt. 0.0D0 ) THEN
            MoonPhaseAng= 360.0D0 + MoonPhaseAng
          ENDIF

      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE HMS_SEC
*
*  this subroutine converts Hours, Minutes and Seconds into seconds from the
*    beginning of the day.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Hr          - Hours                          0 .. 24
*    minute         - Minutes                        0 .. 59
*    Sec         - Seconds                        0.0D0 .. 59.99D0
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    Sec         - Seconds                        0.0D0 .. 86400.0D0
*
*  Locals        :
*    Temp        - Temporary variable
*
*  Coupling      :
*    None.
*
* -----------------------------------------------------------------------------

      SUBROUTINE HMS_SEC     ( Hr,minute, Sec, Direction, UTSec )
        IMPLICIT NONE
        INTEGER Hr, minute
        REAL*8  Sec, UTSec
        CHARACTER*4 Direction
* ----------------------------  Locals  -------------------------------
        Real*8 Temp
        ! --------------------  Implementation   ----------------------
        IF ( Direction .eq. 'FROM' ) THEN
            Temp= UTSec / 3600.0D0
            Hr  = IDINT( Temp )
            minute = IDINT( (Temp - Hr)*60.0D0 )
            Sec = (Temp - Hr - minute/60.0D0 ) * 3600.0D0
          ELSE
            UTSec= Hr*3600.0D0 + minute*60.0D0 + Sec
          ENDIF
      RETURN
      END

* -----------------------------------------------------------------------------
*
*                           SUBROUTINE HMS_UT
*
*  this subroutine converts Hours, Minutes and Seconds into Universal Time.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Hr          - Hours                          0 .. 24
*    minute         - Minutes                        0 .. 59
*    Sec         - Seconds                        0.0D0 .. 59.99D0
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    UT          - Universal Time                 HrMin.Sec
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 205, alg 21, ex 3-10
*
* -----------------------------------------------------------------------------

      SUBROUTINE HMS_UT      ( Hr,minute, Sec, Direction, UT )
        IMPLICIT NONE
        INTEGER Hr, minute
        REAL*8  Sec, UT
        CHARACTER*4 Direction

        ! --------------------  Implementation   ----------------------
        IF ( Direction .eq. 'FROM' ) THEN
            Hr = IDINT( UT*0.01D0 )
            minute= IDINT( UT - Hr*100.0D0 )
            Sec= ( UT-DINT(UT) ) * 100.0D0
          ELSE
            UT= Hr*100.0D0 + minute + Sec*0.01D0
          ENDIF
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE HMS_RAD
*
*  this subroutine converts Hours, Minutes and seconds into radians.  Notice
*    the conversion 0.2617D0 is simply the radian equivalent of 15 degrees.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Hr          - Hours                          0 .. 24
*    minute         - Minutes                        0 .. 59
*    Sec         - Seconds                        0.0D0 .. 59.99D0
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    HMS         - Result                         rad
*
*  Locals        :
*    Temp        - Conversion from hours to rad   0.261799D0
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 204, alg 19 alg 20, ex 3-9
*
* -----------------------------------------------------------------------------

      SUBROUTINE HMS_RAD     ( Hr,minute, Sec, Direction, HMS )
        IMPLICIT NONE
        INTEGER Hr, minute
        REAL*8  Sec, HMS
        CHARACTER*4 Direction
* ----------------------------  Locals  ------------------------------
        REAL*8  Temp

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Temp= 15.0D0*Pi/180.0D0
        IF ( Direction .eq. 'FROM' ) THEN
            Temp= HMS   / Temp
            Hr  = IDINT( Temp   )
            minute = IDINT( (Temp - Hr)*60.0D0 )
            Sec = (Temp - Hr - minute/60.0D0 ) * 3600.0D0
          ELSE
            HMS= ( Hr + minute/60.0D0 + Sec/3600.0D0 )*Temp
          ENDIF
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DMS_RAD
*
*  this subroutine converts Degrees, Minutes and seconds into radians.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Deg         - Degrees                        0 .. 360
*    minute         - Minutes                        0 .. 59
*    Sec         - Seconds                        0.0D0 .. 59.99D0
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    DMS         - Result                         rad
*
*  Locals        :
*    Temp        - Temporary variable
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 203, alg 17 alg 18, ex 3-8
*
* -----------------------------------------------------------------------------

      SUBROUTINE DMS_RAD     ( Deg,minute, Sec, Direction, DMS )
        IMPLICIT NONE
        INTEGER Deg, minute
        REAL*8  Sec, DMS
        CHARACTER*4 Direction
* ----------------------------  Locals  ------------------------------
        REAL*8  Temp

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        IF ( Direction.eq.'FROM' ) THEN
            Temp= DMS * Rad2Deg
            Deg = IDINT( Temp )
            minute = IDINT( (Temp - Deg)*60.0D0 )
            Sec = (Temp - Deg - minute/60.0D0 ) * 3600.0D0
          ELSE
            DMS= ( Deg + minute/60.0D0 + Sec/3600.0D0 ) * Deg2Rad
          ENDIF
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE jd2sse.m
*
*  this subroutine finds the seconds since epoch (1 Jan 2000) given the julian date
*
*  author        : david vallado                  719-573-2600   12 dec 2002
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    jd          - julian date                    days from 4713 bc
*    Direction   - Which set of vars to output    FROM  TOO
*
*  outputs       :
*    sse         - seconds since epoch 1 jan 2000
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
*
*  references    :
*    none.
*
* -----------------------------------------------------------------------------

      SUBROUTINE jd2sse      ( jd,Direction, sse )
        IMPLICIT NONE
        REAL*8  sse, jd
        CHARACTER*4 Direction

        ! --------------------  Implementation   ----------------------
        IF ( Direction.eq.'FROM' ) THEN
            jd = 2451544.5D0 + sse/86400.0D0
          ELSE
            sse = (jd - 2451544.5D0) * 86400.0D0
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE CONVTIME
*
*  this subroutine finds the time parameters and Julian century values for inputs
*    of UTC or UT1. Numerous outputs are found as shown in the local variables.
*    Because calucations are in UTC, you must include TimeZone IF ( you enter a
*    local time, otherwise it should be zero.
*
*  Algorithm     : A file of record contains the timing data
*                  Seeks are performed to obtain the data
*                    Data starts Jan 1, 1980, thus JD = 2444238.5D0 in the code
*                  Calculate the answer depending on initial time type
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Universal Time Hour            0 .. 23
*    minute         - Universal Time minute             0 .. 59
*    SEC         - Universal Time SEC             0.0D0 .. 59.999D0
*    TimeZone    - Offset to UTC from local SITE  0 .. 23 hr
*    TypeUTIn    - Type of input UT               1 (UT1), else UTC
*
*  Outputs       :
*    DUT1        - Delta of UTC - UT1             SEC
*    DAT         - Delta of UTC - TAI             SEC
*    xp          - Polar motion coefficient       arcsec
*    yp          - Polar motion coefficient       arcsec
*    UT1         - Universal time                 SEC
*    TUT1        - Julian centuries of UT1
*    JDUT1       - Julian Date of UT1             days from 4713 BC
*    UTC         - Coordinated Universal Time     SEC
*    TAI         - Atomic time                    SEC
*    TDT         - Terrestrial Dynamical time     SEC
*    TTDT        - Julian centuries of TDT
*    JDTDT       - Julian Date of TDT             days from 4713 BC
*    TDB         - Terrestrial Barycentric time   SEC
*    TTDB        - Julian centuries of TDB
*    JDTDB       - Julian Date of TDB             days from 4713 BC
*    Error       - Error flag DO SUBROUTINE       'ok'
*
*  Locals        :
*    HrTemp      - Temporary hours                hr
*    MinTemp     - Temporary miNutes              minute
*    SecTemp     - Temporary seconds              SEC
*    LocalHr     - Difference to local time       hr
*    JD          - Julian Date of request         days from 4713 BC
*    ME          - Mean Anomaly of the Earth      rad
*    TimeFile    - File of record with time data
*    CurrTimeRec - Current Time record
*
*  Coupling      :
*    HMS_SEC     - Conversion between hr-minute-SEC .and. seconds
*    jday   - Find the Julian date
*
*  References    :
*    vallado       2007, 201, alg 16, ex 3-7
*
* ------------------------------------------------------------------------------

      SUBROUTINE CONVTIME    ( Year, Mon, Day, Hr, minute, SEC,
     &                         TimeZone, TypeUTIn, DUT1, DAT, xp, yp,
     &                         UT1, TUT1, JDUT1, UTC, TAI, TT, TTT,
     &                         JDTT, TDB, TTDB, JDTDB, DDPSi, DDEps,
     &                         LOD, Error )
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, minute, TimeZone
        REAL*8 DUT1, DAT, xp, yp, UT1, TUT1, JDUT1, UTC, TAI, TT, TTT,
     &  Sec, JDTT, TDB, TTDB, JDTDB, DDEPs, DDPsi, LOD, MFME, MJD
        CHARACTER TypeUTIn
        CHARACTER*12  Error
* ----------------------------  Locals  -------------------------------
        INTEGER hrTemp, minTemp, LocalHr
        REAL*8 secTemp, ME, JD

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Error= 'ok'

        ! ------------------ Start IF ( UT1 is known ------------------
        LocalHr= TimeZone + Hr
        IF ( TypeUTIn .eq. '1' ) THEN
            CALL HMS_SEC( LocalHr,minute,SEC, 'TOO',  UT1 )
            CALL jday( Year,Mon,Day, LocalHr, minute, SEC, JDUT1 )
            TUT1= (JDUT1 - 2451545.0D0 )/ 36525.0D0

            UTC= UT1 - dUT1
*            CALL HMS_SEC( HrTemp,MinTemp,SecTemp,'FROM', UTC )
*            CALL jday( Year,Mon,Day, HrTemp, MinTemp, SecTemp,
*     &                      JDUTC )
          ELSE
            ! ---------------- Start IF ( UTC is known ----------------
            CALL HMS_SEC( LocalHr,minute,SEC,'TOO ',  UTC )
*            CALL jday( Year,Mon,Day, LocalHr, minute, SEC, JDUTC )

            UT1= UTC + DUT1
            CALL HMS_SEC( HrTemp,MinTemp,SecTemp,'FROM', UT1 )
            CALL jday( Year,Mon,Day, HrTemp, MinTemp, SecTemp,
     &                      JDUT1 )
            TUT1= (JDUT1 - 2451545.0D0 )/ 36525.0D0
          ENDIF

        TAI= UTC + DAT
*        CALL HMS_SEC( HrTemp,MinTemp,SecTemp,'FROM', TAI )
*        CALL jday( Year,Mon,Day, HrTemp, MinTemp, SecTemp, JDTAI)

        TT= TAI + 32.184D0  ! SEC
        CALL HMS_SEC( HrTemp,MinTemp,SecTemp,'FROM', TT )
        CALL jday( Year,Mon,Day, HrTemp, MinTemp, SecTemp, JDTT)
        TTT= (JDTT - 2451545.0D0 )/ 36525.0D0

        ME= 357.5277233D0 + 35999.05034D0*TTT   ! approx - should do with TTDB
        ME= DMOD( ME,360.0D0 )
        ME= ME * Deg2Rad
        TDB= TT + 0.001658D0 * DSIN(ME) + 0.00001385D0*DSIN(2.0D0*ME)
        CALL HMS_SEC( HrTemp,MinTemp,SecTemp,'FROM', TDB )
        CALL jday( Year,Mon,Day, HrTemp, MinTemp, SecTemp, JDTDB )
        TTDB= (JDTDB - 2451545.0D0 )/ 36525.0D0

      RETURN
      END   ! SUBROUTINE CONVTIME
