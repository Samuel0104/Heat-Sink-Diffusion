function x = solve_LU(L, U, b)

    y = sustitucion_L(L, b);

    x = sustitucion_U(U, y);

end