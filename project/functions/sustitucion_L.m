function x = sustitucion_L(A, b)

    rows = height(A);
    cols = length(A);

    x = zeros(rows, 1);
    
    for j = 1:(cols - 1)
        x(j) = b(j)/A(j, j);
        b(j + 1) = b(j + 1) - A(j + 1, j)*x(j);
    end
    x(cols) = b(cols)/A(cols, cols);

end