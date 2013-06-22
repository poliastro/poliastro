% Compares if two floats are almost equal within a tolerance
%
% Handles corner cases, such as 0 / 0
% Translated from http://floating-point-gui.de/errors/comparison/
%
% Author: Juan Luis Cano
function [res] = almost_equal(a, b, tol);

diff = abs(a - b);
if a == b
	res = true;
	return
elseif (a == 0) || (b == 0) || (diff < eps(1.0))
	res = diff < tol * eps(1.0);
	return
else
	res = diff / (abs(a) + abs(b)) < tol;
	return
end
