% Parâmetros do feixe
delta_rho = 60e-6; % Spot inicial (60 micrômetros)
lambda = 1064e-9; % Comprimento de onda (metros)
c = 3e8; % Velocidade da luz (m/s)
ct = 0.18; % Constante c * t (metros)

% Constantes derivadas
k0 = 2 * pi / lambda; % Número de onda

%% a) Calculando s
w0 = sqrt(sqrt(3)*k0*(delta_rho^2)/2);
s_a = 1/(k0*w0)

%% b) Calculando s para lambda = 0.5um e w0 = 50um
lambda = 0.5e-6;
w0 = 50e-6;
k0 = 2 * pi / lambda; % Número de onda
s_b = 1/(k0*w0)

%% c) Calculando s para w0 = lambda, lambda = 532nm e lambda = 1064nm
lambda1 = 523e-9;
lambda2 = 1064e-9;

w0_1 = lambda1;
w0_2 = lambda2;

k0_1 = 2 * pi / lambda1;
k0_2 = 2 * pi / lambda2;

s_c1 = 1/(k0_1*w0_1)
s_c2 = 1/(k0_2*w0_2)

