% Converts a function on 2-space into a list

function v = grid2list(u, N)
    v = zeros([1 N^2]);
    for i = 1:N
        for j = 1:N
            v((i - 1) * N + j) = u(i, j);
        end
    end
end