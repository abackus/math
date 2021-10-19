function u = dufort_frankel_heat(f, k, T)
    % Runs a Dufort-Frankel method on the heat equation.
    n = length(f);
    Laplace = right_diff(n) * left_diff(n);
    
    prev = f;
    curr = prev + k * Laplace * prev;
    for i=1:(round(T/k) - 1)
        next = prev + 2 * k * Laplace * curr;
        prev = curr;
        curr = next;
    end 
    u = curr;
end