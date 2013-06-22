% ------------------------------------------------------------------------------
%
%                           function newtonm
%
%  this function performs the newton rhapson iteration to find the
%    eccentric anomaly given the mean anomaly.  the true anomaly is also
%    calculated.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    ecc         - eccentricity                   0.0  to
%    m           - mean anomaly                   -2pi to 2pi rad
%
%  outputs       :
%    e0          - eccentric anomaly              0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%
%  locals        :
%    e1          - eccentric anomaly, next value  rad
%    sinv        - sine of nu
%    cosv        - cosine of nu
%    ktr         - index
%    r1r         - cubic roots - 1 to 3
%    r1i         - imaginary component
%    r2r         -
%    r2i         -
%    r3r         -
%    r3i         -
%    s           - variables for parabolic solution
%    w           - variables for parabolic solution
%
%  coupling      :
%    cubic       - solves a cubic polynomial
%
%  references    :
%    vallado       2001, 72-75, alg 2, ex 2-1
%
% [e0,nu] = newtonm ( ecc,m );
% ------------------------------------------------------------------------------

function [e0,nu] = newtonm ( ecc,m );

        % -------------------------  implementation   -----------------
        numiter =    50;
        small   =     0.00000001;
        halfpi  = pi * 0.5;
        
        % -------------------------- hyperbolic  ----------------------
        if ( (ecc-1.0 ) > small )
           % -------------------  initial guess -----------------------
            if ( ecc < 1.6  )
                if ( ((m<0.0 ) & (m>-pi)) | (m>pi) )
                    e0= m - ecc;
                  else
                    e0= m + ecc;
                  end
              else
                if ( (ecc < 3.6 ) & (abs(m) > pi) )
                    e0= m - sign(m)*ecc;
                  else
                    e0= m/(ecc-1.0 );
                  end
              end
            ktr= 1;
            e1 = e0 + ( (m-ecc*sinh(e0)+e0) / (ecc*cosh(e0) - 1.0 ) );
            while ((abs(e1-e0)>small ) & ( ktr<=numiter ))
                e0= e1;
                e1= e0 + ( ( m - ecc*sinh(e0) + e0 ) / ( ecc*cosh(e0) - 1.0  ) );
                ktr = ktr + 1;
              end
            % ----------------  find true anomaly  --------------------
            sinv= -( sqrt( ecc*ecc-1.0  ) * sinh(e1) ) / ( 1.0  - ecc*cosh(e1) );
            cosv= ( cosh(e1) - ecc ) / ( 1.0  - ecc*cosh(e1) );
            nu  = atan2( sinv,cosv );
          else
            % --------------------- parabolic -------------------------
            if ( abs( ecc-1.0  ) < small )
%                c = [ 1.0/3.0; 0.0; 1.0; -m];
%                [r1r] = roots (c);
%                e0= r1r;
                 s = 0.5  * (halfpi - atan( 1.5 *m ) );
                 w = atan( tan( s )^(1.0 /3.0 ) );
                 e0= 2.0 *cot(2.0 *w);
                ktr= 1;
                nu = 2.0  * atan(e0);
              else
                % -------------------- elliptical ----------------------
                if ( ecc > small )
                    % -----------  initial guess -------------
                    if ( ((m < 0.0 ) & (m > -pi)) | (m > pi) )
                        e0= m - ecc;
                      else
                        e0= m + ecc;
                      end
                    ktr= 1;
                    e1 = e0 + ( m - e0 + ecc*sin(e0) ) / ( 1.0  - ecc*cos(e0) );
                    while (( abs(e1-e0) > small ) & ( ktr <= numiter ))
                        ktr = ktr + 1;
                        e0= e1;
                        e1= e0 + ( m - e0 + ecc*sin(e0) ) / ( 1.0  - ecc*cos(e0) );
                      end
                    % -------------  find true anomaly  ---------------
                    sinv= ( sqrt( 1.0 -ecc*ecc ) * sin(e1) ) / ( 1.0 -ecc*cos(e1) );
                    cosv= ( cos(e1)-ecc ) / ( 1.0  - ecc*cos(e1) );
                    nu  = atan2( sinv,cosv );
                  else
                    % -------------------- circular -------------------
                    ktr= 0;
                    nu= m;
                    e0= m;
                  end
              end
          end

