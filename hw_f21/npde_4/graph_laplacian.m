% Returns a graph Laplacian.
% Our convention is that this Laplacian is positive-definite!

function Laplace = graph_laplacian(N)
    Laplace = 4 .* eye(N^2);
    for i = 1:N^2
        for j = (i+1):N^2
            if is_adjacent(i, j, N)
                Laplace(i, j) = -1;
                Laplace(j, i) = -1;
            end
        end
    end
end