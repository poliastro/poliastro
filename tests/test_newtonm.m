% Test `newtonm`
%
% Schlesinger & Udick
% Author: Juan Luis Cano

addpath('../')

% Constants
tol = 1.0e-4;
rad = 180.0 / pi;  % deg / rad

% Data: ecc, M (deg), nu (deg)
data = [
	0.0, 0.0, 0.0;
	0.05, 10.0, 11.06;
	0.06, 30.0, 33.67;
	0.04, 120.0, 123.87;
	0.14, 65.0, 80.50;
	0.19, 21.0, 30.94;
	0.35, 65.0, 105.71;
	0.48, 180.0, 180.0;
	0.75, 125.0, 167.57;
];

for ii = 1:size(data, 1)
	ecc0 = data(ii, 1);
	M0 = data(ii, 2) / rad;
	nu0 = data(ii, 3) / rad;
	[_, nu] = newtonm(ecc0, M0);
	assert(almost_equal(nu, nu0, tol));
end

disp('All tests passed')
