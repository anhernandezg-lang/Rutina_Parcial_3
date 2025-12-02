# PARCIAL 3
# PUNTO 1 - Reproducir un análisis estadístico de un artículo
# Dataset: "Pathogen detection Salmonella enterica" (Kaggle)
# Prueba estadística: Kruskal–Wallis Articulo

import pandas as pd
from scipy.stats import kruskal
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Cargar el dataset

ruta = "PARCIAL 3/archive/Pathogen detection Salmonella enterica.csv"
df = pd.read_csv(ruta)

print("Columnas del archivo:")
print(df.columns.tolist())
print()

# 2. Limpieza para tener solo las columnas necesarias para el analisis

columnas_necesarias = ["Serovar", "Min-diff", "Isolation type"]
df = df[columnas_necesarias].copy()

# Quitar filas con datos faltantes en las columnas clave
df = df.dropna(subset=columnas_necesarias)

# Asegurar que Min-diff es numérico
df["Min-diff"] = pd.to_numeric(df["Min-diff"], errors="coerce")
df = df.dropna(subset=["Min-diff"])

print(f"Número de filas después de limpieza: {len(df)}")
print()

# 3 ANÁLISIS ESTADISTICO: KRUSKAL–WALLIS POR TIPO DE AISLAMIENTO
#   Se utiliza Kruskal–Wallis para comparar un conteo entre varios grupos.
groups = []
labels = []

for iso_type, subdf in df.groupby("Isolation type"):
    if len(subdf) >= 5:  # mínimo de datos por grupo
        groups.append(subdf["Min-diff"].values)
        labels.append(iso_type)

print("Tipos de aislamiento incluidos en Kruskal:")
for lab, g in zip(labels, groups):
    print(f"  {lab}: n = {len(g)}")
print()

# 4. Prueba de Kruskal–Wallis
stat, p_value = kruskal(*groups)

print("Kruskal–Wallis test for SNP differences by isolation type")
print(f"H = {stat:.3f}")
print(f"p-value = {p_value:.6e}")
print()

# 5. Boxplot general por tipo de aislamiento

plt.figure(figsize=(8, 5))
sns.boxplot(data=df, x="Isolation type", y="Min-diff")
plt.xticks(rotation=45, ha="right")
plt.title("Min-diff (SNP differences) por tipo de aislamiento")
plt.tight_layout()
plt.show()


# VISUALIZACIÓN POR SEROVAR

conteos_serovar = df["Serovar"].value_counts()

print("los 20 serovares más frecuentes:")
print(conteos_serovar.head(20))
print()

# Umbral mínimo de tamaño de muestra para considerar un serovar
n_min = 500

serovares_candidatos = conteos_serovar[conteos_serovar >= n_min].index

# De esos candidatos, se seleccionan los 8 más frecuentes
if len(serovares_candidatos) > 0:
    serovares_seleccionados = conteos_serovar[serovares_candidatos].head(8).index
else:
    # Si ningún serovar supera n_min, tomamos simplemente los 8 más frecuentes
    serovares_seleccionados = conteos_serovar.head(8).index
    print(f"No hubo serovares con >= {n_min} aislamientos, se usan los 8 más frecuentes.")

print("\nSerovares seleccionados para las figuras:")
for s in serovares_seleccionados:
    print(f"  {s}: n = {conteos_serovar[s]}")
print()

# DataFrame solo con esos serovares para las gráficas
df_plot = df[df["Serovar"].isin(serovares_seleccionados)].copy()

# Ordenar los serovares por la mediana de Min-diff
orden_serovar = (
    df_plot.groupby("Serovar")["Min-diff"]
    .median()
    .sort_values()
    .index
)

# 7. Boxplot por serovar (serovares frecuentes seleccionados)

plt.figure(figsize=(10, 5))
sns.boxplot(
    data=df_plot,
    x="Serovar",
    y="Min-diff",
    order=orden_serovar
)
plt.xticks(rotation=45, ha="right")
plt.title("Min-diff por serovar (serovares frecuentes)")
plt.tight_layout()
plt.show()

# 8. Violinplot por serovar
#
plt.figure(figsize=(10, 5))
sns.violinplot(
    data=df_plot,
    x="Serovar",
    y="Min-diff",
    order=orden_serovar,
    inner="box"
)
plt.xticks(rotation=45, ha="right")
plt.title("Distribución de Min-diff por serovar (violin plot, serovares frecuentes)")
plt.tight_layout()
plt.show()
