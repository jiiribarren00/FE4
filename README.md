Aquí tienes una documentación completa, profesional y lista para ser utilizada como archivo README.md (o cualquier documento de especificación técnica) para tu librería PyODR.
PyODR: Documentación Técnica de la Librería
PyODR es una implementación nativa en Python moderna, ligera y altamente optimizada del algoritmo de Regresión por Distancia Ortogonal (ODR), inspirada en el paquete clásico ODRPACK y diseñada con una interfaz 100% compatible con scipy.odr.
Permite realizar ajustes de modelos no lineales y lineales cuando tanto la variable independiente (X) como la dependiente (Y) poseen incertezas experimentales.
🚀 Características Principales
 * Compatibilidad Total con SciPy: Comparte exactamente la misma estructura de clases e interfaz que scipy.odr (Data, RealData, Model, ODR, Output, y la función rápida odr), facilitando la migración de código existente.
 * Alto Rendimiento (Optimización de Memoria): Bucle interno optimizado con preasignación de memoria estática (in-place updates), reduciendo drásticamente el overhead de recolección de basura (garbage collection) durante iteraciones masivas de Levenberg-Marquardt.
 * Robustez Numérica: Implementación de inversiones matriciales seguras mediante descomposición SVD y pseudoinversas (scipy.linalg), previniendo errores por matrices singulares (LinAlgError).
 * Tipado Estricto (Type Hinting): Anotaciones de tipo completas en todo el código para potenciar el autocompletado y análisis estático en IDEs como VSCode o PyCharm.
📦 Requisitos e Instalación
PyODR requiere únicamente dependencias estándar de cómputo científico en Python:
 * Python 3.8 o superior.
 * numpy >= 1.20
 * scipy >= 1.6
📖 Guía de Inicio Rápido
A continuación se muestra un ejemplo completo de cómo ajustar un modelo lineal simple con incertezas en ambos ejes utilizando pyodr:
import numpy as np
from pyodr import RealData, Model, ODR

# 1. Definir la función del modelo: f(beta, x) -> y
# beta[0]: Ordenada al origen | beta[1]: Pendiente
def modelo_lineal(beta, x):
    return beta[0] + beta[1] * x

# 2. Datos experimentales y sus desviaciones estándar (sx, sy)
x_data = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
y_data = np.array([2.1, 3.9, 6.2, 7.8, 10.1])
s_x = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
s_y = np.array([0.2, 0.2, 0.2, 0.2, 0.2])

# 3. Empaquetar los datos usando RealData (acepta desviaciones estándar)
datos = RealData(x_data, y_data, sx=s_x, sy=s_y)

# 4. Crear la instancia del modelo
modelo = Model(modelo_lineal)

# 5. Configurar el objeto ODR con estimaciones iniciales para beta ([ordenada, pendiente])
problema_odr = ODR(datos, modelo, beta50=[0.0, 2.0] if 'beta50' in globals() else [0.0, 2.0]) # o beta0=[0.0, 2.0]

# 6. Ejecutar la optimización
resultado = problema_odr.run()

# 7. Imprimir informe resumido
resultado.pprint()

🧩 Referencia de la API
1. Data
Contenedor base para los datos crudos y sus pesos de precisión estadística.
 * Parámetros:
   * x (array_like): Coordenadas independientes. Forma (M, N) o (N,).
   * y (array_like, optional): Coordenadas dependientes. Forma (NQ, N) o (N,).
   * we (array_like, optional): Pesos estadísticos para Y (inverso de la varianza).
   * wd (array_like, optional): Pesos estadísticos para X (inverso de la varianza).
   * fix (array_like, optional): Arreglo de enteros para fijar observaciones específicas.
2. RealData
Hereda de Data. Permite ingresar directamente las desviaciones estándar o incertidumbres experimentales comunes en laboratorios de física e ingeniería.
 * Parámetros:
   * x, y (array_like): Datos centrales.
   * sx (array_like, optional): Desviación estándar en X (se convierte automáticamente a wd = 1/sx^2).
   * sy (array_like, optional): Desviación estándar en Y (se convierte automáticamente a we = 1/sy^2).
3. Model
Define la función matemática que modela el fenómeno físico o estadístico.
 * Parámetros:
   * fcn (Callable): Función con la firma fcn(beta, x).
   * fjacb (Callable, optional): Jacobiano analítico respecto a los parámetros \beta.
   * fjacd (Callable, optional): Jacobiano analítico respecto a las variables X.
   * estimate (Callable, optional): Función auxiliar para estimar parámetros iniciales a partir de los datos.
4. ODR
Motor principal del algoritmo de optimización.
 * Método clave: run() -> Output
   Ejecuta el proceso iterativo y devuelve un objeto contenedor con los resultados detallados.
5. Output
Estructura de datos que almacena el resultado del ajuste.
 * Atributos Principales:
   * beta (np.ndarray): Valores óptimos de los parámetros ajustados.
   * sd_beta (np.ndarray): Errores estándar (incertezas 1\sigma) de cada parámetro.
   * cov_beta (np.ndarray): Matriz de covarianza de los parámetros.
   * res_var (float): Varianza residual (\chi^2 reducido).
   * xplus (np.ndarray): Puntos óptimos corregidos en el eje X sobre la curva.
   * ynew (np.ndarray): Valores evaluados del modelo en los puntos corregidos.
💡 Notas Técnicas y Buenas Prácticas
 * Relación entre Incertezas y Matriz de Covarianza:
   El error estándar reportado en sd_beta[i] combina la curvatura geométrica del ajuste (cov_beta) escalada por la varianza residual (res_var):
   
 * Control de Grados de Libertad:
   Asegúrate de que el número de observaciones independientes supere estrictamente al número de parámetros a ajustar (N > p), de lo contrario la estimación de la varianza residual carecerá de grados de libertad válidos.
