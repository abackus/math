% A script that runs a semiimplicit and an implicit method on the Burgers
% equation.
h = 1/100;
T = 50;
k = 1/300;
eta = .005;
x = transpose(0:h:2*pi);
f = x - pi;
D = center_diff(length(f));
Laplace = right_diff(length(f)) * left_diff(length(f));

tic
% Comment the next line out when using the implicit method 
% inverse = inv(eye(length(f)) - eta * k * Laplace);

u = f;
for i=1:round(T/k)
    % Implicit
    u = (eye(length(f)) - eta * k * Laplace - k * diag(u) * D) \ u;
    % Semi-implicit
    % u = inverse * (eye(length(f)) + k * diag(u) * D) * u;
end
toc