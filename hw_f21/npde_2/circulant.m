% A function that takes in a vector x and returns the circulant matrix
% whose first row is x.

function A = circulant(x)
    A = zeros(length(x), length(x));
    for n = 1:length(x)
        for m = 1:length(x)
            A(n, m) = x(m);
        end
        x = circshift(x, 1);
    end
end