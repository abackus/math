% An implementation of an implicit \theta-scheme for the transport equation 
% on a circle.
% f is the initial data, theta is the \theta-parameter, lambda is the ratio
% of spatial to temporal grid sizes, T is the final time.
% Need \theta \in [0.5, 1] to get stability.

function v = theta_scheme(f, theta,lambda, T)
    h = 2 * pi/length(f); % Spatial grid size
    k = h * lambda; % Temporal grid size
    timesteps = round(T / k);
    
    v = zeros(timesteps, length(f));
    
    % Initialize v
    for j = 1:length(f)
        v(1, j) = f(j);
    end
    
    % Create the difference operator
    x = zeros(length(f));
    x(1) = 1;
    x(2) = -theta * lambda/2;
    x(length(f)) = theta * lambda/2;
    A = circulant(x);
    x(2) = (1 - theta) * lambda/2;
    x(length(f)) = - (1 - theta) * lambda/2;
    B = circulant(x);
    Q = A\B;
    
    for n = 2:timesteps
        x = zeros(length(f));
        for j = 1:length(f)
            x(j) = v(n - 1, j);
        end
        x = Q * x;
        for j = 1:length(f)
            v(n, j) = x(j);
        end
    end
end