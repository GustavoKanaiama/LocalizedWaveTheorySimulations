import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0  # Função Bessel de ordem zero
from tqdm import tqdm
from PIL import Image

# Parâmetros principais
wavelength = 632e-9         # Comprimento de onda (632 nm)
k = 2 * np.pi / wavelength  # Número de onda
L = 0.06                    # Comprimento longitudinal (6 cm)
N = 15                      # Número de Bessel Beams (BBs) por onda congelada
Q = 0.9998 * k              # Grau de paraxialidade
P = 276                     # Número de ondas congeladas (igual ao número de colunas da imagem)
delta_x = 48.38e-6          # Espaçamento entre ondas congeladas (48.38 μm)

# Carregar a imagem binarizada
imagem_bin = Image.open("imagem_binarizada_276x75.png").convert('L')
imagem_array = np.array(imagem_bin) / 255.0  # Normalizar para [0, 1]

# Inverter a imagem: onde for 0 vira 1, e onde for 1 vira 0
imagem_invertida = 1 - imagem_array  # Isso inverte a imagem binarizada

# Converter de volta para imagem (se necessário)
imagem_invertida_bin = (imagem_invertida * 255).astype(np.uint8)  # Converter para 0-255
imagem_invertida_pil = Image.fromarray(imagem_invertida_bin)

# Salvar ou mostrar a imagem invertida
imagem_invertida_pil.save("imagem_invertida.png")
imagem_invertida_pil.show()

imagem_bin = Image.open("imagem_invertida.png").convert('L')
imagem_array = np.array(imagem_bin) / 255.0  # Normalizar para [0, 1]

# Coordenadas da grade
x_points, z_points = imagem_array.shape
x_vals = np.linspace(0, delta_x * P, x_points)  # Ajuste para começar de 0
z_vals = np.linspace(0, L, imagem_array.shape[0])  # Garantir correspondência com as linhas da imagem

# Função auxiliar para calcular kρq e βq
def k_rho_beta(q, Q, k, L):
    """ Calcula kρq e βq para o índice q. """
    beta_q = Q + 2 * np.pi * q / L
    k_rho_q = np.sqrt(k**2 - beta_q**2)
    return k_rho_q, beta_q

# Calcular os coeficientes Aqpp
def calculate_Aqpp(image_column, q, z_vals, L):
    """ Calcula os coeficientes Aqpp com base na coluna da imagem binarizada. """
    phase = np.exp(1j * 2 * np.pi * q * z_vals / L)
    integral = np.trapz(image_column * phase, z_vals)  # Integração numérica
    return integral / L

# Construção de ΨSFW
def calculate_Psi_SFW(x_vals, z_vals, image_array, delta_x, N, Q, k, L):
    """ Calcula ΨSFW(x, z) com base nas equações (4) e (6). """
    P = image_array.shape[1]
    psi_sfw = np.zeros((len(x_vals), len(z_vals)), dtype=complex)

    for p, rho_0p in tqdm(enumerate(np.linspace(0, delta_x * P, P))):  # Para cada onda congelada
        image_column = image_array[:, p]

        for q in range(-N, N + 1):  # Para cada Bessel Beam
            k_rho_q, beta_q = k_rho_beta(q, Q, k, L)
            Aqpp = calculate_Aqpp(image_column, q, z_vals, L)

            for i, rho in enumerate(x_vals):
                for j, z in enumerate(z_vals):
                    arg = np.sqrt(rho**2 + rho_0p**2 - 2 * rho * rho_0p)
                    bessel = j0(k_rho_q * arg) if arg > 0 else 1
                    psi_sfw[i, j] += Aqpp * bessel * np.exp(-1j * beta_q * z)

    return np.real(psi_sfw)**2 + np.imag(psi_sfw)**2 # Intensidade (magnitude ao quadrado)

# Calcular ΨSFW
psi_sfw = calculate_Psi_SFW(x_vals, z_vals, imagem_array, delta_x, N, Q, k, L)

# Criar um colormap do preto ao verde
import matplotlib.colors as mcolors
cmap = mcolors.LinearSegmentedColormap.from_list("blackgreen", ["black", "green"])

# Plotar o resultado
plt.figure(figsize=(12, 4))
plt.imshow(psi_sfw.T, extent=[0, L * 1e2, 0, delta_x * P * 1e2], cmap=cmap, aspect='auto')
plt.title('Density plot of $|\Psi_{SFW}(x, z)|^2$')
plt.xlabel('z (cm)')
plt.ylabel('x (cm)')
plt.colorbar(label='Intensidade Normalizada')
plt.show()