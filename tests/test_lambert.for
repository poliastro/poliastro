        implicit none

*         Exactly the same results (success!)
*         vo = 2.05899   2.91597  -0.00000
*         v = -3.45153   0.91031  -0.00000

        real*8 ro(3),r(3), Dtsec, k_Earth,vo(3),v(3)
        character dm, OverRev
        character*12 Error

        external LAMBERTBATTIN, LAMBERTUNIV

        ro = (/15945.15, 0.0, 0.0/)
        r = (/12214.83899, 10249.46731, 0.0/)
        dm = 'S'
        OverRev = 'N'
        Dtsec = 76.0 * 60
        k_Earth = 398600.4418

        call LAMBERTBATTIN( ro,r, dm,OverRev, Dtsec, k_Earth,
     &                     vo,v, Error )

        print *, vo
        print *, v

        call LAMBERTUNIV  ( ro,r, dm,OverRev, Dtsec, k_Earth,
     &                     vo,v, Error )

        print *, vo
        print *, v

        end
