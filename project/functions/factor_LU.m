function [L, U] = factor_LU(A)

    rows = height(A);

    [U, ~, L] = gauss(A, zeros(rows, 1));

end