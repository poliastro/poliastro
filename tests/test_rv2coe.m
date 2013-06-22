% Test `rv2coe`
%
% Vallado, 3rd ed., 2007, ex. 2.5
% Author: Juan Luis Cano

addpath('../')

% Constants
tol = 1.0e-4;
rad = 180.0 / pi;  % deg / rad

% Initial data
k_Earth = 3.986e5;  % km3 / s2
r = [6524.834, 6862.875, 6448.296];  % km
v = [4.901327, 5.533756, -1.976341];  % km / s

% Solution
p0 = 11067.790;  % km
a0 = 36127.343;  % km
ecc0 = 0.832953;
inc0 = 87.870 / rad;  % rad
omega0 = 227.89 / rad;  % rad
argp0 = 53.38 / rad;  % rad
nu0 = 92.335 / rad;  % rad

[p, a, ecc, inc, omega, argp, nu] = rv2coe(r, v, k_Earth);

% Assertions
assert(almost_equal(p, p0, tol));
assert(almost_equal(a, a0, tol));
assert(almost_equal(ecc, ecc0, tol));
assert(almost_equal(inc, inc0, tol));
assert(almost_equal(omega, omega0, tol));
assert(almost_equal(argp, argp0, tol));
assert(almost_equal(nu, nu0, tol));

disp('All tests passed')
