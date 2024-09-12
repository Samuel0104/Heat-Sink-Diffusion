L = input("\nLongitud de la superficie de difusión: ");
t_max = input("\nTiempo total de la simulación: ");
n_x = input("\nNúmero de nodos en el espacio: ");
n_t = input("\nNúmero de nodos en el tiempo: ");
alpha = input("\nConstante del material: ");
T_1j = input("\nTemperatura del contorno: ");
T_i1 = input("\nTemperatura inicial: ");

d_x = L/(n_x - 1); % Se calcula delta de x
d_t = t_max/(n_t - 1); % Se calcula delta de t





%% Estado transitorio

fprintf("\n\n########## Estado transitorio ##########\n\n");

% Coeficientes para la matriz A y el vector b
a = d_t;
b = -(2*d_t + d_x^2*d_t*alpha^2 + d_x^2);
c = d_t;
d = -d_x^2;

% Se crea la matriz A de la forma:
% b  c  0  0  0
% a  b  c  0  0
% 0  a  b  c  0
% 0  0  a  b  c
% 0  0  0 a+c b
A_1 = [zeros(1, n_x - 1); a*eye(n_x - 2, n_x - 1)];
A_2 = b*eye(n_x - 1);
A_3 = [zeros(n_x - 1, 1) c*eye(n_x - 1, n_x - 2)];
A = A_1 + A_2 + A_3;
A(n_x - 1, n_x - 2) = A(n_x - 1, n_x - 2) + c;

% Se crea la solución inicial (t = 0) para graficar
temp = T_i1*ones(n_x - 1, 1); % Temperatura inicial
T_A = [T_1j; temp]; % Temperatura de contorno

figure(1)
p = plot(0:d_x:L, T_A);
p.LineWidth = 2;
hold on

tic
[L_A, U_A] = factor_LU(A); % Se factoriza la matriz A

% Se crea el vector b de la forma:
% d*T_(2,j-1) - a*T_(1,j)
%       d*T_(3,j-1)
%       d*T_(4,j-1)
%       d*T_(5,j-1)
%       d*T_(6,j-1)
sol_A = d*T_i1*ones(n_x - 1, 1);
sol_A(1) = sol_A(1) - a*T_1j;

for i = 2:n_t
    % Se halla la nueva solución
    temp = solve_LU(L_A, U_A, sol_A);
    T_A = [T_1j; temp];

    p = plot(0:d_x:L, T_A);
    p.LineWidth = 2;

    % Se actualiza el vector b para usar en la siguiente iteración
    sol_A = d*T_A;
    sol_A(1) = sol_A(1) - a*T_1j;
end
hold off
toc




%% Estado estable

fprintf("\n\n########## Estado estable ##########\n\n");

% Coeficientes para la matriz B
a = -(2 + d_x^2*alpha^2);

% Se crea la matriz B de la forma:
% a  1  0  0  0
% 1  a  1  0  0
% 0  1  a  1  0
% 0  0  1  a  1
% 0  0  0  2  a
B_1 = [zeros(1, n_x - 1); eye(n_x - 2, n_x - 1)];
B_2 = a*eye(n_x - 1);
B_3 = [zeros(n_x - 1, 1) eye(n_x - 1, n_x - 2)];
B = B_1 + B_2 + B_3;
B(n_x - 1, n_x - 2) = 2;

% Se crea el vector b de la forma:
% -T_(1,j)
%    0
%    0
%    0
%    0
sol_B = zeros(n_x - 1, 1);
sol_B(1) = -T_1j;
tic
% Se halla la solución para graficar
T_B = resolver(B, sol_B);
T_B = [T_1j; T_B];

figure(2)
p = plot(0:d_x:L, T_B);
p.Color = "#4DBEEE";
p.LineWidth = 2;
toc