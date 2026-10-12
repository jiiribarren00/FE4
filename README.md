# PyODR

*Creado por Juan I Iribarren para los estudiantes de Física Experimental IV. Se hizo integramente con IA, no es una excusa, pero ante cualquier error avisar al autor con un email a juan.iribarren@fisica.unlp.edu.ar*

**PyODR** es una implementación nativa en Python moderna, ligera y altamente optimizada del algoritmo de **Regresión por Distancia Ortogonal (ODR)**. Está inspirada en el paquete clásico ODRPACK y diseñada con una interfaz 100% compatible con `scipy.odr`.

Permite realizar ajustes de modelos no lineales y lineales cuando **tanto la variable independiente ($X$) como la dependiente ($Y$) poseen incertezas experimentales**.

---

## 🚀 Características Principales

* **Compatibilidad Total con SciPy:** Comparte la misma estructura de clases e interfaz que `scipy.odr` (`Data`, `RealData`, `Model`, `ODR`, `Output` y la función rápida `odr`), facilitando la migración de código existente.
* **Alto Rendimiento:** Bucle interno optimizado con preasignación de memoria estática (*in-place updates*), reduciendo drásticamente el *overhead* de recolección de basura durante iteraciones masivas de Levenberg-Marquardt.
* **Robustez Numérica:** Implementación de inversiones matriciales seguras mediante descomposición SVD y pseudoinversas (`scipy.linalg`), previniendo errores por matrices singulares (`LinAlgError`).
* **Tipado Estricto (*Type Hinting*):** Anotaciones de tipo completas en todo el código para potenciar el autocompletado y análisis estático en IDEs como VSCode o PyCharm.

---

## 📦 Requisitos e Instalación

PyODR requiere únicamente dependencias estándar de cómputo científico en Python:
* Python 3.8 o superior.
* `numpy >= 1.20`
* `scipy >= 1.6`

Puedes integrar el archivo `pyodr.py` directamente en la raíz de tu proyecto.

---

## 📖 Guía de Inicio Rápido

A continuación se muestra un ejemplo completo de cómo ajustar un modelo lineal simple con incertezas en ambos ejes utilizando `pyodr`:

```python
import numpy as np
from pyodr import RealData, Model, ODR
```
### 1. Definir la función del modelo: f(beta, x) -> y
beta[0]: Ordenada al origen | beta[1]: Pendiente
```
def modelo_lineal(beta, x):
    return beta[0] + beta[1] * x
```

### 2. Datos experimentales y sus desviaciones estándar (sx, sy)
```
x_data = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
y_data = np.array([2.1, 3.9, 6.2, 7.8, 10.1])
s_x = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
s_y = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
```

### 3. Empaquetar los datos usando RealData (acepta desviaciones estándar)
```
datos = RealData(x_data, y_data, sx=s_x, sy=s_y)
```

### 4. Crear la instancia del modelo
```
modelo = Model(modelo_lineal)
```

### 5. Configurar el objeto ODR con estimaciones iniciales para beta
```
problema_odr = ODR(datos, modelo, beta0=[0.0, 2.0])
```

### 6. Ejecutar la optimización
```
resultado = problema_odr.run()
```

### 7. Imprimir informe resumido
```
resultado.pprint()
```
---

## 🧩 Referencia de la API

1. Data
Contenedor base para los datos crudos y sus pesos de precisión estadística.
 * Parámetros:
   * `x` (array_like): Coordenadas independientes. Forma (M, N) o (N,).
   * `y` (array_like, optional): Coordenadas dependientes. Forma (NQ, N) o (N,).
   * `we` (array_like, optional): Pesos estadísticos para Y (inverso de la varianza).
   * `wd` (array_like, optional): Pesos estadísticos para X (inverso de la varianza).
   * `fix` (array_like, optional): Arreglo de enteros para fijar observaciones específicas.
   * `meta` (dict, optional): Diccionario con metadatos personalizados.

2. RealData
Hereda de Data. Permite ingresar directamente las desviaciones estándar o incertidumbres experimentales habituales en laboratorios de física e ingeniería.
 * Parámetros:
   * `x`, `y` (array_like): Datos centrales de las mediciones.
   * `sx` (array_like, optional): Desviación estándar en X (se convierte automáticamente a wd = 1/sx^2).
   * `sy` (array_like, optional): Desviación estándar en Y (se convierte automáticamente a we = 1/sy^2).
   * `covx`, `covy` (array_like, optional): Covarianzas alternativas en lugar de varianzas simples.

3. Model
Define la función matemática que modela el fenómeno físico o estadístico a ajustar.
 * Parámetros:
   * `fcn` (Callable): Función con la firma fcn(beta, x).
   * `fjacb` (Callable, optional): Jacobiano analítico respecto a los parámetros \beta.
   * `fjacd` (Callable, optional): Jacobiano analítico respecto a las variables X.
   * `estimate` (Callable, optional): Función auxiliar para estimar parámetros iniciales a partir de los datos.

4. ODR
Motor principal del algoritmo de optimización y control de convergencia.
 * Método clave: `run()` -> `Output`
   Ejecuta el proceso iterativo y devuelve un objeto contenedor con los resultados detallados del ajuste.
   
5. Output
Estructura de datos que almacena el resultado final del ajuste estadístico.
 * Atributos Principales:
   * `beta` (np.ndarray): Valores óptimos de los parámetros ajustados.
   * `sd_beta` (np.ndarray): Errores estándar (incertezas 1\sigma) de cada parámetro.
   * `cov_beta` (np.ndarray): Matriz de covarianza de los parámetros.
   * `res_var` (float): Varianza residual (\chi^2 reducido).
   * `xplus` (np.ndarray): Puntos óptimos corregidos en el eje X sobre la curva.
   * `ynew` (np.ndarray): Valores evaluados del modelo en los puntos corregidos.

## 📐 Fundamentación Matemática del Algoritmo ODR

A diferencia de los Mínimos Cuadrados Ordinarios (OLS), que asumen que la variable independiente x está libre de error y solo minimizan los residuos verticales en y, la Regresión por Distancia Ortogonal (ODR) considera que tanto x como y contienen incertidumbres experimentales.

1. Formulación del Problema

Para un conjunto de n observaciones independientes, denotadas por coordenadas $(x_i, y_i)$ con $i = 1, \dots, n$, se busca ajustar un modelo no lineal dado por:
Donde:
 * $\beta$ es el vector de parámetros globales del modelo a determinar de tamaño $p$ .
 * $\delta_i$ representa el vector de ajuste u optimización (error implícito) en la variable independiente $x_i$.
 * $\varepsilon_i$ representa el residuo en la variable dependiente $y_i$.

2. La Función Objetivo Ponderada

El algoritmo ODR minimiza la suma ponderada de los cuadrados de las distancias ortogonales (o generalizadas) tanto en el eje x como en el eje y. La función de costo $S(\beta, \delta)$ a minimizar es:
Donde:
 * $w_{e,i} = \sigma_{y,i}^{-2}$ son los pesos asociados a las observaciones de y.
 * $w_{d,i} = \sigma_{x,i}^{-2}$ son los pesos asociados a las observaciones de x.

3. Resolución Numérica (Levenberg-Marquardt Generalizado)

Este problema de mínimos cuadrados no lineales combinados se resuelve reestructurándolo como un sistema extendido donde se optimizan simultáneamente los parámetros libres $\beta$ y los incrementos de posición $\delta_i$. PyODR implementa esto mapeando el problema hacia el algoritmo de Levenberg-Marquardt a través de la función de costo vectorizada:

5. Estimación de Errores y Matriz de Covarianza

Una vez alcanzada la convergencia en el punto óptimo, se calcula la Varianza Residual (`res_var` o $\chi^2_{reducido}$) para cuantificar la calidad global del ajuste frente a los grados de libertad del sistema ($N_{obs} - p$):
La Matriz de Covarianza de los parámetros \beta se obtiene a partir del Jacobiano J de los residuos evaluado respecto a los parámetros libres en el óptimo:
Finalmente, la incertidumbre estándar (desviación estándar de cada parámetro, reportada como `sd_beta`) se extrae directamente de la raíz cuadrada de la diagonal principal de esta matriz:

## 📚 Referencias Bibliográficas
 * Boggs, P. T., Byrd, R. H., & Schnabel, R. B. (1987). A stable and efficient algorithm for nonlinear orthogonal distance regression. SIAM Journal on Scientific and Statistical Computing, 8(6), 1052-1078.
 * Boggs, P. T., Donaldson, J. R., Byrd, R. H., & Schnabel, R. B. (1989). Algorithm 676: ODRPACK: Software for orthogonal distance regression. ACM Transactions on Mathematical Software (TOMS), 15(4), 348-364.
 * National Institute of Standards and Technology (NIST). Guide to the Expression of Uncertainty in Measurement (GUM) / ODRPACK Documentation and Reference Manuals.

