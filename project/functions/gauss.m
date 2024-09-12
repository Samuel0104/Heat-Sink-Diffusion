function [A, b, L] = gauss(A, b)

    rows = height(A);
    cols = length(A);

    a = A(1, 2);
    d = A(1, 1);

    L = eye(rows, cols);
    
    for i = 2:(rows - 1)
        pivot = a/A(i - 1, i - 1);
        A(i, i - 1) = 0;
        A(i, i) = d - a*pivot;
        b(i) = b(i) - pivot*b(i - 1);
        L(i, i - 1) = pivot;
    end

    pivot = 2*a/A(rows - 1, cols - 1);
    A(rows, cols - 1) = 0;
    A(rows, cols) = d - a*pivot;
    b(rows) = b(rows) - pivot*b(rows - 1);
    L(rows, cols - 1) = pivot;

end