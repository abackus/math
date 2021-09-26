% An implementation of an upwind scheme for the transport equation 
% on a circle.
% f is the initial data, lambda is the ratio
% of spatial to temporal grid sizes, T is the final time.

function v = upwind_scheme(f, lambda, T)
    h = 2 * pi/length(f); % Spatial grid size
    k = h * lambda; % Temporal grid size
    timesteps = round(T / k);
    
    v = zeros(timesteps, length(f));
    
    % Initialize v
    for j = 1:length(f)
        v(1, j) = f(j);
    end
    
    for n = 2:timesteps
        for j = 1:length(f)
            next = obob_mod(j + 1, length(f));
            v(n, j) = v(n - 1, j) + lambda * (v(n - 1, next) - v(n - 1, j));
        end
    end
end
