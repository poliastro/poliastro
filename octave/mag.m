% ------------------------------------------------------------------------------
%
%                            function mag
%
%  this function finds the magnitude of a vector.  the tolerance is set to
%    0.000001, thus the 1.0e-12 for the squared test of underflows.
%
%  author        : david vallado                  719-573-2600   30 may 2002
%
%  revisions
%    vallado     - fix tolerance to match coe, eq, etc            3 sep 2002
%
%  inputs          description                    range / units
%    vec         - vector
%
%  outputs       :
%    mag         - magnitude
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
% mag = ( vec );
% ----------------------------------------------------------------------------- }

function mag = mag ( vec );

        temp= vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3);

        if abs( temp ) >= 1.0e-16
            mag= sqrt( temp );
          else
            mag= 0.0;
          end

