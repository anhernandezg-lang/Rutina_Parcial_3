📘 README

Este repositorio contiene el desarrollo de tres ejercicios independientes para el examen final de modelado y análisis de datos biológicos:

1️⃣ Análisis estadístico (Kruskal–Wallis)

Se reproduce el enfoque estadístico utilizado en un artículo sobre resistencia antimicrobiana en Salmonella.
Se aplica la prueba Kruskal–Wallis a un dataset real de Kaggle para evaluar si las diferencias SNP (Min-diff) cambian según el tipo de aislamiento.
Incluye filtrado, justificación metodológica y visualizaciones (boxplot y violin plot).

2️⃣ Ajuste de un modelo biológico

Se implementa el ajuste de un modelo tipo Monod con tres parámetros (µ_max, K_s, Y_X/S) usando datos simulados basados en estudios reales de producción de riboflavina en bioprocesos.
Se utiliza curve_fit para estimar parámetros y comparar datos experimentales vs modelo.

3️⃣ Dinámica de un sistema biológico

Se simula un sistema de tres componentes: biomasa X(t), sustrato S(t) y producto P(t) en un quimiostato, resolviendo un sistema de EDOs con solve_ivp.
Incluye comparación de escenarios (D bajo vs D alto), curvas dinámicas y análisis en estado estacionario.

📦 Requisitos

Python 3.9+

numpy

pandas

scipy

matplotlib

seaborn
