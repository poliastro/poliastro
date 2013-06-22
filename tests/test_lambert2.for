        implicit none

!         ?
!         v0 = [-14.21190144,  18.80586495,   1.59734319]
!         v  = [ 22.36280009,  27.66651827,   0.31027719]

        real*8 ro(3),r(3), Dtsec, k_Sun,vo(3),v(3)
        character dm, OverRev
        character*12 Error

        external LAMBERTBATTIN, LAMBERTUNIV

        ro = (/-1.09038587D08,  -1.04221731D08,   0.00000000D00/)
        r = (/-1.02959046D08,   3.05523274D07,   6.35994599D06/)
        dm = 'S'
        OverRev = 'N'
        Dtsec = 5184000.0D0
        k_Sun = 132712440018.0D0

!         call LAMBERTBATTIN( ro,r, 'S',OverRev, Dtsec, k_Sun,
!      &                     vo,v, Error )
! 
!         print *, vo
!         print *, v
! 
!         call LAMBERTUNIV  ( ro,r, 'S',OverRev, Dtsec, k_Sun,
!      &                     vo,v, Error )
! 
!         print *, vo
!         print *, v

        call LAMBERTUNIV  ( ro,r, 'L',OverRev, Dtsec, k_Sun,
     &                     vo,v, Error )

        print *, vo
        print *, v

        end
