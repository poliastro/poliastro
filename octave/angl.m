% ------------------------------------------------------------------------------
%
%                            function angl
%
%  this function calculates the angle between two vectors.  the output is
%    set to 999999.1 to indicate an undefined value.  be sure to check for
%    this at the output phase.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%    vallado     - fix tolerances                                 5 sep 2002
%
%  inputs          description                    range / units
%    vec1        - vector number 1
%    vec2        - vector number 2
%
%  outputs       :
%    theta       - angle between the two vectors  -pi to pi
%
%  locals        :
%    temp        - temporary real variable
%
%  coupling      :
%
% [theta] = angl ( vec1,vec2 );
% ----------------------------------------------------------------------------- }

function [theta] = angl ( vec1,vec2 );

        small     = 0.00000001;
        undefined = 999999.1;

        magv1 = mag(vec1);
        magv2 = mag(vec2);

        if magv1*magv2 > small^2
            temp= dot(vec1,vec2) / (magv1*magv2);
            if abs( temp ) > 1.0
                temp= sign(temp) * 1.0;
              end
            theta= acos( temp );
          else
            theta= undefined;
          end

