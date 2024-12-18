% Parâmetros do feixe
delta_rho = 60e-6; % Spot inicial (60 micrômetros)
lambda = 1064e-9; % Comprimento de onda (metros)
c = 3e8; % Velocidade da luz (m/s)
ct = 0.18; % Constante c * t (metros)

% Constantes derivadas
a = delta_rho / 2; % Amplitude inicial relacionada ao spot
k0 = 2 * pi / lambda; % Número de onda
b = ct / (2 * c); % Cálculo da constante b

% Coordenadas
rho = linspace(-3e-3, 3e-3, 500); % Coordenada radial limitada de -0.4mm a 0.4mm
z = linspace(0, 0.4, 500); % Coordenada longitudinal de 0 até 1.8cm
[rho1, z1] = meshgrid(rho, z);

% Cálculo da função Psi
denominator = a^2 + 1i * z1 / (2 * k0); % Denominador
spatial_gaussian = exp(-rho1.^2 ./ (4 * denominator)); % Termo espacial
temporal_gaussian = exp(-((z1 - ct).^2) ./ (4 * c^2 * b^2)); % Termo temporal
psi = (a^2 ./ denominator) .* spatial_gaussian .* temporal_gaussian;

% Normaliza Psi
psi_normalized = abs(psi) / max(abs(psi(:)));

%% Figura 1: Gráfico 3D do Pulso Escalar Gaussiano
figure;
surf(rho, z, psi_normalized, 'EdgeColor', 'none');
colormap jet;
colorbar;
xlabel('\rho (mm)', 'FontSize', 12);
ylabel('z (cm)', 'FontSize', 12);
zlabel('|Psi Normalizado|', 'FontSize', 12);
title('Pulso Escalar Gaussiano - 3D', 'FontSize', 14);
shading interp;

%% Figura 2: Gráficos 2D para z = 0, z = 0.9 cm, z = 1.8 cm
z_values = [0, 0.9e-2, 1.8e-2]; % Valores específicos de z (em metros)
%rho = linspace(-0.4e-3, 0.4e-3, 500);
figure;
xlim([-0.4 0.4]);
hold on;
for i = 1:length(z_values)
    % Seleciona o índice de z mais próximo do valor desejado
    [~, idx] = min(abs(z - z_values(i)));
    z_fixed = z(idx);
    
    % Extração da linha correspondente
    psi_line = psi_normalized(idx, :);
    
    % Plot da linha
    plot(rho * 1e3, psi_line, 'LineWidth', 1.5, 'DisplayName', sprintf('z = %.1f cm', z_fixed * 100));
end
hold off;

% Melhorar visualização
xlabel('\rho (mm)', 'FontSize', 12);

ylabel('|Psi Normalizado|', 'FontSize', 12);
title('Perfil do feixe gaussiano psi(rho)', 'FontSize', 14);
legend show;
grid on;
