h = 2^(-6);
x = 0:h:2 * pi;
f = zeros(length(x),1);
u = zeros(length(x),1); % Represents u at time 1.

for j = 1:length(x)
    if x(j) <= pi 
        f(j) = x(j);
        u(mod(j - 1/h, length(x)) + 1) = x(j);
    else 
        f(j) = 2 * pi - x(j);
        u(mod(j - 1/h, length(x)) + 1) = 2 * pi - x(j);
    end
end

% v = upwind_scheme(f, 0.5, 1); % Upwind scheme
% v = naive_scheme(f, 0.5, 1); % Naive scheme
% v = lax_schemes(f, 0.5, "Fred", 1); % Lax-Friedrichs scheme
% v = lax_schemes(f, 0.5, "Wend", 1); % Lax-Wendroff scheme
% v = theta_scheme(f, 0.5, 0.5, 1); % Crank-Nicholson scheme
v = theta_scheme(f, 1, 0.5, 1); % Backwards Euler scheme

uapprox = v(end,:);

err = sqrt(h) * norm(u - transpose(uapprox));

plot(x, u, x, uapprox);
legend({'Analytic solution', 'Numeric solution'});
disp(err);