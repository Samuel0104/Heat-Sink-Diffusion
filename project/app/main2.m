function main2(axes1, axes2, L, t_max, n_x, n_t, alpha, T_1j, T_i1)

    d_x = L/(n_x - 1);
    d_t = t_max/(n_t - 1);
    
    
    
    
    
    %% Estado transitorio

    a = d_t;
    b = -(2*d_t + d_x^2*d_t*alpha^2 + d_x^2);
    c = d_t;
    d = -d_x^2;

    A_1 = [zeros(1, n_x - 1); a*eye(n_x - 2, n_x - 1)];
    A_2 = b*eye(n_x - 1);
    A_3 = [zeros(n_x - 1, 1) c*eye(n_x - 1, n_x - 2)];
    A = A_1 + A_2 + A_3;
    A(n_x - 1, n_x - 2) = A(n_x - 1, n_x - 2) + c;
    
    temp = T_i1*ones(n_x - 1, 1);
    T_A = [T_1j; temp];
    
    p = plot(axes1, 0:d_x:L, T_A);
    p.LineWidth = 2;
    hold(axes1, "on")
    
    [L_A, U_A] = factor_LU(A);

    sol_A = d*T_i1*ones(n_x - 1, 1);
    sol_A(1) = sol_A(1) - a*T_1j;
    
    for i = 2:n_t
        temp = solve_LU(L_A, U_A, sol_A);
        T_A = [T_1j; temp];
    
        p = plot(axes1, 0:d_x:L, T_A);
        p.LineWidth = 2;
    
        sol_A = d*T_A;
        sol_A(1) = sol_A(1) - a*T_1j;
    end
    hold(axes1, "off")
    
    
    
    
    
    %% Estado estable
    
    a = -(2 + d_x^2*alpha^2);
    
    B_1 = [zeros(1, n_x - 1); eye(n_x - 2, n_x - 1)];
    B_2 = a*eye(n_x - 1);
    B_3 = [zeros(n_x - 1, 1) eye(n_x - 1, n_x - 2)];
    B = B_1 + B_2 + B_3;
    B(n_x - 1, n_x - 2) = 2;
    
    sol_B = zeros(n_x - 1, 1);
    sol_B(1) = -T_1j;
    
    T_B = resolver(B, sol_B);
    T_B = [T_1j; T_B];
    
    p = plot(axes2, 0:d_x:L, T_B);
    p.Color = "#4DBEEE";
    p.LineWidth = 2;

end