% Generate the Poisson data for Problem 3

function F = pr3_poisson_data(N)
    x = [1/(N+2):1/(N+1):1-1/(N+2)];
    row = 4 .* sin(x) .* cos(x) - 2;
    RHS = zeros(N);
    for i=1:N
        for j=1:N
            RHS(j,i) = row(i);
        end
    end
    trace = zeros([4 N]);
    trig = sin(x) .* cos(x);
    square = (1 - x) .^ 2;
    for i=1:N
        trace(1, i) = trig(i) + 1;
        trace(2, i) = square(i);
        trace(3, i) = trig(i);
        trace(4, i) = sin(1) * cos(1) + square(i);
    end
    F = poisson_data(RHS, trace, N);
end