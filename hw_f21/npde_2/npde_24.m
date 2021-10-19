h = 1/100;
T = 2;

x = transpose(0:h:2*pi);
x_strong = transpose(0:1/500:2*pi);

f = x - pi;
f_strong = x_strong - pi;

u_strong = crank_nicholson_heat(f_strong, 1/3 * (1/500)^2, T);
% u_strong = crank_nicholson_heat(f_strong, 1/3 * (1/500), T);
plot(x_strong, u_strong);

k = 1/3 * h^2;

% u = dufort_frankel_heat(f, k, 2);
% plot(x, u);


axis([0 2*pi -4 4]);
