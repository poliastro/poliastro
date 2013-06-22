        implicit none

        real*8 k, ro(3), vo(3), r(3), v(3), tof
        character*12 error

        external KEPLER

        ro = (/3.20712964848534613848D07,  1.02947987170683801174D08,
     &           -4.57704247282716271002D05 /)
        vo(1) = -33.55207844069929024045D0
        vo(2) = 10.24876605422645070576D0
        vo(3) = 2.07652696079603060753D0
        tof = 13305600.0D0
        k = 1.32712440018D11

        call KEPLER(ro, vo, tof, k, r, v, error)

        print *, r
        print *, v
        print *, error

        end
