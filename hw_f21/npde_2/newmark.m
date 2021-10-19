% A wave equation Newmark scheme.
function u = newmark(f, g, T, k, c)
    P = k * (c^2)/4 * right_diff(length(f)) * left_diff(length(f));
    u = f;
    w = g;
    Q = k * inv(eye(length(f)) - k * P);
    R = (eye(length(f)) - k * P) \ (eye(length(f)) + k * P);
    for i=1:(round(T/k) - 1)
        u_next = Q * w + R * u;
        w = P * (u_next + u) + (1/k) * (u_next - u);
        u = u_next;
    end
end