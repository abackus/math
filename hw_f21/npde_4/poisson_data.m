% Create the data used to solve Poisson's problem
% Trace is top, left, bottom, right

function F = poisson_data(RHS, trace, N)
    F = RHS;
    for i = 1:N
        F(1, i) = F(1, i) + N^2 * trace(1, i);
    end
    for i = 1:N
        F(i, 1) = F(i, 1) + N^2 * trace(2, i);
    end
    for i = 1:N
        F(N, i) = F(N, i) + N^2 * trace(3, i);
    end
    for i = 1:N
        F(i, N) = F(i, N) + N^2 * trace(4, i);
    end
    F = F .* N^(-2);
end