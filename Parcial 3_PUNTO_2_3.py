# PUNTO 2 - Ajuste de un modelo biológico
# Modelo: X(D) = Yxs * (Sin - D*Ks / (mu_max - D))
# Parámetros biológicos a ajustar: mu_max, Ks, Yxs  (3 parámetros)

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")

# 2.1 Datos "experimentales" simulados a partir de la Figura 2a del articulo

# Tasas de dilución (h^-1) aproximadas visualmente del artículo
D_exp = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14])

# Biomasa X (g/L) aproximada de la curva de la Figura 2a
X_exp = np.array([1.45, 1.38, 1.33, 1.30, 1.22, 1.05, 0.55])

df_biomasa = pd.DataFrame({"Dilucion_h-1": D_exp, "Biomasa_gL": X_exp})
print("\n PUNTO 2: Datos simulados de biomasa vs dilución ")
print(df_biomasa)

#
# 2.2 Modelo biológico X(D) basado en Monod + chemostato
#     X(D) = Yxs * (Sin - D*Ks / (mu_max - D))

Sin = 10.0  # g/L, concentración de glucosa en la alimentación (tomada del medio del artículo)

def modelo_X(D, mu_max, Ks, Yxs):

    # evitar división por cero
    return Yxs * (Sin - (D * Ks) / (mu_max - D))

# 2.3 Ajuste no lineal de los parámetros mu_max, Ks, Yxs

# Valores iniciales razonables para el ajuste
p0 = [0.15, 0.05, 0.5]   # [mu_max, Ks, Yxs]

param_opt, param_cov = curve_fit(modelo_X, D_exp, X_exp, p0=p0, maxfev=10000)
mu_max_fit, Ks_fit, Yxs_fit = param_opt

print("\n=== Parámetros ajustados (Punto 2) ===")
print(f"mu_max ajustado = {mu_max_fit:.4f} h^-1")
print(f"Ks ajustado     = {Ks_fit:.4f} g/L")
print(f"Yxs ajustado    = {Yxs_fit:.4f} gX/gS")

# 2.4 Gráfico: datos vs modelo ajustado

D_fino = np.linspace(0.019, 0.145, 300)
X_modelo = modelo_X(D_fino, mu_max_fit, Ks_fit, Yxs_fit)

plt.figure(figsize=(8,5))
plt.scatter(D_exp, X_exp, color="black", label="Datos simulados (Figura 2a)")
plt.plot(D_fino, X_modelo, color="red", linewidth=2,
         label="Modelo ajustado X(D)")

plt.xlabel("Tasa de dilución D [h$^{-1}$]")
plt.ylabel("Biomasa X [g/L]")
plt.title("PUNTO 2 - Ajuste de biomasa vs dilución\nModelo biológico basado en Monod")
plt.legend()
plt.tight_layout()
plt.show()

# PUNTO 3 - Dinámica de un sistema biológico con 3 componentes

from scipy.integrate import solve_ivp

# Usamos los parámetros ajustados del Punto 2:
mu_max = mu_max_fit  # h^-1
Ks = Ks_fit  # g/L
Yxs = Yxs_fit  # gX/gS
Sin = Sin  # ya definido (10 g/L)

# Definimos un parámetro de producción de riboflavina
qp = 0.6  # mg de P producidos / (g de X · h)

# Dos tasas de dilución para comparar escenarios
D_baja = 0.06  # h^-1  (cultivo estable, buena biomasa)
D_alta = 0.13  # h^-1  (cerca del washout)

# 3.1 Sistema de EDOs

def modelo_quimiostato(t, y, D):

    X, S, P = y

    # cinética de Monod
    mu = mu_max * S / (Ks + S)

    dXdt = (mu - D) * X
    dSdt = D * (Sin - S) - (1.0 / Yxs) * mu * X
    dPdt = qp * X - D * P

    return [dXdt, dSdt, dPdt]

# 3.2 Condiciones iniciales y tiempo

X0 = 0.1  # g/L
S0 = Sin  # g/L
P0 = 0.0  # mg/L

y0 = [X0, S0, P0]

t_span = (0, 200)  # horas
t_eval = np.linspace(t_span[0], t_span[1], 500)

# 3.3 Simulaciones para D_baja y D_alta


sol_baja = solve_ivp(modelo_quimiostato, t_span, y0, t_eval=t_eval, args=(D_baja,))
sol_alta = solve_ivp(modelo_quimiostato, t_span, y0, t_eval=t_eval, args=(D_alta,))

# 3.4 Gráficas de la dinámica

def plot_dinamica(sol, D, titulo_extra=""):
    t = sol.t
    X, S, P = sol.y

    plt.figure(figsize=(9, 6))
    plt.plot(t, X, label="Biomasa X [g/L]")
    plt.plot(t, S, label="Sustrato S [g/L]")
    plt.plot(t, P, label="Riboflavina P [mg/L]")
    plt.xlabel("Tiempo [h]")
    plt.ylabel("Concentración")
    plt.title(f"PUNTO 3 - Dinámica del quimiostato (D = {D:.3f} h$^{{-1}}$){titulo_extra}")
    plt.legend()
    plt.tight_layout()
    plt.show()


plot_dinamica(sol_baja, D_baja, " - Dilución baja")
plot_dinamica(sol_alta, D_alta, " - Dilución alta")

# Diagrama de fase X vs S para ver el estado estacionario
plt.figure(figsize=(7, 5))
plt.plot(sol_baja.y[1], sol_baja.y[0], label=f"D baja = {D_baja:.3f}")
plt.plot(sol_alta.y[1], sol_alta.y[0], label=f"D alta = {D_alta:.3f}")
plt.xlabel("Sustrato S [g/L]")
plt.ylabel("Biomasa X [g/L]")
plt.title("PUNTO 3 - Diagrama de fase X vs S")
plt.legend()
plt.tight_layout()
plt.show()

## =D ##