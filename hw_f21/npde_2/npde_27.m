% A script that runs a semiimplicit and an implicit method on the Burgers
% equation.
h = 1/100;
T = 1;
k = 1/200;
c = 1;
x = transpose(0:h:2*pi);
f = sin(x);
g = cos(x);

% u = leapfrog(f, g, T, k, c);
u = newmark(f, g, T, k, c);