function x = resolver(A, b)

    [U, b, ~] = gauss(A, b);
    x = sustitucion_U(U, b);

end