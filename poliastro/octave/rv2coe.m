%
% ------------------------------------------------------------------------------
%
%                           function rv2coe
%
%  this function finds the classical orbital elements given the geocentric
%    equatorial position and velocity vectors.
%
%  author        : david vallado                  719-573-2600   21 jun 2002
%
%  revisions
%    vallado     - fix special cases                              5 sep 2002
%    vallado     - delete extra check in inclination code        16 oct 2002
%    vallado     - add constant file use                         29 jun 2003
%    vallado     - add mu                                         2 apr 2007
%
%  inputs          description                    range / units
%    r           - ijk position vector            km
%    v           - ijk velocity vector            km / s
%    mu          - gravitational parameter        km3 / s2
%
%  outputs       :
%    p           - semilatus rectum               km
%    a           - semimajor axis                 km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omega       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    m           - mean anomaly                   0.0  to 2pi rad
%    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
%    truelon     - true longitude            (ce) 0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
%
%  locals        :
%    hbar        - angular momentum h vector      km2 / s
%    ebar        - eccentricity     e vector
%    nbar        - line of nodes    n vector
%    c1          - v**2 - u/r
%    rdotv       - r dot v
%    hk          - hk unit vector
%    sme         - specfic mechanical energy      km2 / s2
%    i           - index
%    e           - eccentric, parabolic,
%                  hyperbolic anomaly             rad
%    temp        - temporary variable
%    typeorbit   - type of orbit                  ee, ei, ce, ci
%
%  coupling      :
%    mag         - magnitude of a vector
%    angl        - find the angl between two vectors
%    newtonnu    - find the mean anomaly
%
%  references    :
%    vallado       2007, 121, alg 9, ex 2-5
%
% [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r,v);
% ------------------------------------------------------------------------------

function [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r,v, mu);

        constmath;
%          constastro;  % don't overwrite mu

        % -------------------------  implementation   -----------------
        magr= mag( r );
        magv= mag( v );
        % ------------------  find h n and e vectors   ----------------
        [hbar] = cross( r,v );
        magh= mag( hbar );
        if ( magh > small )
            nbar(1)= -hbar(2);
            nbar(2)=  hbar(1);
            nbar(3)=   0.0;
            magn = mag( nbar );
            c1 = magv*magv - mu /magr;
            rdotv= dot( r,v );
            for i= 1 : 3
                ebar(i)= (c1*r(i) - rdotv*v(i))/mu;
              end
            ecc = mag( ebar );

            % ------------  find a e and semi-latus rectum   ----------
            sme= ( magv*magv*0.5  ) - ( mu /magr );
            if ( abs( sme ) > small )
                a= -mu  / (2.0 *sme);
              else
                a= infinite;
              end
            p = magh*magh/mu;

            % -----------------  find inclination   -------------------
            hk= hbar(3)/magh;
            incl= acos( hk );

            % --------  determine type of orbit for later use  --------
            % ------ elliptical, parabolic, hyperbolic inclined -------
            typeorbit= 'ei';
            if ( ecc < small )
                % ----------------  circular equatorial ---------------
                if  (incl<small) | (abs(incl-pi)<small)
                    typeorbit= 'ce';
                  else
                    % --------------  circular inclined ---------------
                    typeorbit= 'ci';
                  end
              else
                % - elliptical, parabolic, hyperbolic equatorial --
                if  (incl<small) | (abs(incl-pi)<small)
                    typeorbit= 'ee';
                  end
              end

            % ----------  find longitude of ascending node ------------
            if ( magn > small )
                temp= nbar(1) / magn;
                if ( abs(temp) > 1.0  )
                    temp= sign(temp);
                  end
                omega= acos( temp );
                if ( nbar(2) < 0.0  )
                    omega= twopi - omega;
                  end
              else
                omega= undefined;
              end

            % ---------------- find argument of perigee ---------------
            if ( typeorbit == 'ei' )
                argp = angl( nbar,ebar);
                if ( ebar(3) < 0.0  )
                    argp= twopi - argp;
                  end
              else
                argp= undefined;
              end

            % ------------  find true anomaly at epoch    -------------
            if ( typeorbit(1:1) == 'e' )
                nu =  angl( ebar,r);
                if ( rdotv < 0.0  )
                    nu= twopi - nu;
                  end
              else
                nu= undefined;
              end

            % ----  find argument of latitude - circular inclined -----
            if ( typeorbit == 'ci' )
                arglat = angl( nbar,r );
                if ( r(3) < 0.0  )
                    arglat= twopi - arglat;
                  end
                m = arglat;
              else
                arglat= undefined;
              end

            % -- find longitude of perigee - elliptical equatorial ----
            if  ( ecc>small ) & (typeorbit=='ee')
                temp= ebar(1)/ecc;
                if ( abs(temp) > 1.0  )
                    temp= sign(temp);
                  end
                lonper= acos( temp );
                if ( ebar(2) < 0.0  )
                    lonper= twopi - lonper;
                  end
                if ( incl > halfpi )
                    lonper= twopi - lonper;
                  end
              else
                lonper= undefined;
              end

            % -------- find true longitude - circular equatorial ------
            if  ( magr>small ) & ( typeorbit=='ce' )
                temp= r(1)/magr;
                if ( abs(temp) > 1.0  )
                    temp= sign(temp);
                  end
                truelon= acos( temp );
                if ( r(2) < 0.0  )
                    truelon= twopi - truelon;
                  end
                if ( incl > halfpi )
                    truelon= twopi - truelon;
                  end
                m = truelon;
              else
                truelon= undefined;
              end

            % ------------ find mean anomaly for all orbits -----------
            if ( typeorbit(1:1) == 'e' )
                [e,m] = newtonnu(ecc,nu );
              end

         else
           p    = undefined;
           a    = undefined;
           ecc  = undefined;
           incl = undefined;
           omega= undefined;
           argp = undefined;
           nu   = undefined;
           m    = undefined;
           arglat = undefined;
           truelon= undefined;
           lonper = undefined;
         end

