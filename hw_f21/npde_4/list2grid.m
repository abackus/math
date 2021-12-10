% Converts a list to a function on 2-space

function u = list2grid(v, N)
    u = zeros(N);
    for i=1:N
        for j = 1:N
            u(i, j) = v((i-1)*N + j);
        end
    end
end