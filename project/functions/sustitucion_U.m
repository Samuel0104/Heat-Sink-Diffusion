function x = sustitucion_U(A, b)

    rows = height(A);
    cols = length(A);

    x = zeros(rows, 1);
    
    for j = cols:-1:2
        x(j) = b(j)/A(j, j);
        b(j - 1) = b(j - 1) - A(j - 1, j)*x(j);
    end
    x(1) = b(1)/A(1, 1);

end