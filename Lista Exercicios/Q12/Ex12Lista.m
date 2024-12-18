% Parâmetros do feixe
delta_rho = 60e-6; % Spot inicial (60 micrômetros)
lambda = 1064e-9; % Comprimento de onda (metros)
c = 3e8; % Velocidade da luz (m/s)
a = delta_rho / 2; % Amplitude inicial relacionada ao spot
k0 = 2 * pi / lambda; % Número de onda

% Constantes "ct" para os gráficos
ct_values = [0.00001, 0.09, 0.18, 0.27]; % Valores de c*t (metros)

% Coordenadas
rho = linspace(-0.4e-2, 0.4e-2, 1000); % Coordenada radial limitada de -0.4mm a 0.4mm
z = linspace(0, 0.4, 1000); % Coordenada longitudinal de 0 até 4cm
[rho1, z1] = meshgrid(rho, z);

% Loop para gerar os gráficos

for j = 1:length(ct_values)
    figure;
    ct = ct_values(j); % Valor atual de ct
    b = delta_rho*10e2/(2*c); % Constante b dependente de ct, evita valores muito pequenos
    
    % Cálculo da função Psi
    denominator = a^2 + (1i * z1 / (2 * k0)); % Denominador
    spatial_gaussian = exp(-rho1.^2 ./ (4 * denominator)); % Termo espacial
    temporal_gaussian = exp(-((z1 - ct).^2) ./ (4 * c^2 * b^2)); % Termo temporal
    psi = (a^2 ./ denominator) .* spatial_gaussian .* temporal_gaussian;

    % Normaliza Psi
    psi_normalized = real(abs(psi) / max(abs(psi(:))));

    % Plot 3D - Superfície
    surf(rho, z, psi_normalized, 'EdgeColor', 'none'); % rho em mm, z em cm
    colormap jet;
    shading interp;
    xlabel('\rho (mm)', 'FontSize', 12);
    ylabel('z (cm)', 'FontSize', 12);
    zlabel('|Psi Normalizado|', 'FontSize', 12);
    title(sprintf('Pulso Gaussiano - ct = %.2f m', ct), 'FontSize', 12);
    view(45, 30);
end