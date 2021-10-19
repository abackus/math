% A wave equation leapfrog scheme.
function u = leapfrog(f, g, T, k, c)
    P = 2 * eye(length(f)) + k^2 * c^2 * right_diff(length(f)) * left_diff(length(f));
    prev = f;
    curr = f + k * g;
    for i=1:(round(T/k) - 1)
        next = P * curr - prev;
        prev = curr;
        curr = next;
    end
    u = curr;
end