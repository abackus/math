% Returns true iff two indices refer to adjacent vertices in a NxN grid.

function bool = is_adjacent(i, j, N)
    bool = (i == j);
    if (j < i) % Enforce that i < j
        jtmp = i;
        i = j;
        j = jtmp;
    end
    bool = bool || (j == i + N) || ((j == i + 1) && ~(mod(i, N) == 0));
end