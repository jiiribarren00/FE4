"""
PyODR: Implementación nativa en Python de Regresión por Distancia Ortogonal (ODR).
Basado en ODRPACK y compatible con la interfaz de scipy.odr.
"""

import numpy as np
import scipy.linalg as la
from scipy.optimize import least_squares
from typing import Optional, Union, Callable, Any, Dict, Tuple, List

__all__ = ['Data', 'RealData', 'Model', 'Output', 'ODR', 'odr']


class Data:
    """
    Contenedor de datos para el problema de Regresión por Distancia Ortogonal (ODR).

    Parameters
    ----------
    x : array_like
        Datos de la variable independiente. Forma (M, N) o (N,).
    y : array_like, optional
        Datos de la variable dependiente. Forma (NQ, N) o (N,).
    we : array_like, optional
        Pesos para Y. Si se omite, se asume un arreglo de unos.
    wd : array_like, optional
        Pesos para X. Si se omite, se asume un arreglo de unos.
    fix : array_like, optional
        Arreglo de enteros indicando qué observaciones fijar.
    meta : dict, optional
        Metadatos adicionales.
    """
    def __init__(self, 
                 x: Union[List, np.ndarray], 
                 y: Optional[Union[List, np.ndarray]] = None, 
                 we: Optional[Union[List, np.ndarray]] = None, 
                 wd: Optional[Union[List, np.ndarray]] = None, 
                 fix: Optional[Union[List, np.ndarray]] = None, 
                 meta: Optional[Dict[str, Any]] = None):
        
        self.x = np.atleast_2d(np.asarray(x, dtype=np.float64))
        self.y = np.atleast_2d(np.asarray(y, dtype=np.float64)) if y is not None else None
            
        self.m, self.n = self.x.shape
        self.nq = self.y.shape[0] if self.y is not None else 1
        
        self.we = self._process_weights(we, (self.nq, self.n))
        self.wd = self._process_weights(wd, (self.m, self.n))
        
        self.fix = fix
        self.meta = meta if meta is not None else {}

    def _process_weights(self, w: Optional[Union[List, np.ndarray]], target_shape: Tuple[int, int]) -> np.ndarray:
        """Formatea y expande el arreglo de pesos a las dimensiones objetivo."""
        if w is None:
            return np.ones(target_shape, dtype=np.float64)
        
        w_arr = np.asarray(w, dtype=np.float64)
        if w_arr.ndim == 1 and w_arr.shape[0] == target_shape[1]:
            return np.tile(w_arr, (target_shape[0], 1))
        return w_arr


class RealData(Data):
    """
    Variante de Data que acepta desviaciones estándar o covarianzas en lugar de pesos directos.
    Transforma sx, sy en wd, we asumiendo pesos recíprocos al cuadrado.
    """
    def __init__(self, 
                 x: Union[List, np.ndarray], 
                 y: Optional[Union[List, np.ndarray]] = None, 
                 sx: Optional[Union[List, np.ndarray]] = None, 
                 sy: Optional[Union[List, np.ndarray]] = None, 
                 covx: Optional[Union[List, np.ndarray]] = None, 
                 covy: Optional[Union[List, np.ndarray]] = None, 
                 fix: Optional[Union[List, np.ndarray]] = None, 
                 meta: Optional[Dict[str, Any]] = None):
        we, wd = None, None
        
        if sy is not None:
            we = 1.0 / np.square(np.asarray(sy, dtype=np.float64))
        elif covy is not None:
            we = 1.0 / np.asarray(covy, dtype=np.float64)

        if sx is not None:
            wd = 1.0 / np.square(np.asarray(sx, dtype=np.float64))
        elif covx is not None:
            wd = 1.0 / np.asarray(covx, dtype=np.float64)

        super().__init__(x, y=y, we=we, wd=wd, fix=fix, meta=meta)


class Model:
    """
    Define el modelo matemático a ajustar.

    Parameters
    ----------
    fcn : Callable
        Función a ajustar. Debe tener la firma `fcn(beta, x) -> y`.
    fjacb : Callable, optional
        Jacobiano respecto a los parámetros beta.
    fjacd : Callable, optional
        Jacobiano respecto a las variables x.
    estimate : Callable, optional
        Función para estimar los valores iniciales de beta dados los datos.
    """
    def __init__(self, 
                 fcn: Callable[[np.ndarray, np.ndarray], np.ndarray], 
                 fjacb: Optional[Callable] = None, 
                 fjacd: Optional[Callable] = None, 
                 meta: Optional[Dict[str, Any]] = None, 
                 estimate: Optional[Callable] = None, 
                 implicit: int = 0):
        self.fcn = fcn
        self.fjacb = fjacb
        self.fjacd = fjacd
        self.meta = meta if meta is not None else {}
        self.estimate = estimate
        self.implicit = implicit


class Output:
    """Almacena los resultados estadísticos del ajuste ODR."""
    def __init__(self, beta: np.ndarray, sd_beta: np.ndarray, cov_beta: np.ndarray, 
                 res_var: float, rel_error: float, eps: np.ndarray, delta: np.ndarray, 
                 info: int, stopreason: str, xplus: np.ndarray, ynew: np.ndarray):
        self.beta = beta
        self.sd_beta = sd_beta
        self.cov_beta = cov_beta
        self.res_var = res_var
        self.rel_error = rel_error
        self.eps = eps
        self.delta = delta
        self.info = info
        self.stopreason = stopreason
        self.xplus = xplus
        self.ynew = ynew

    def pprint(self) -> None:
        """Imprime por consola un resumen estadístico del ajuste."""
        print(self)

    def __repr__(self) -> str:
        return (f"ODR Results\n"
                f"===========\n"
                f"Beta: {self.beta}\n"
                f"Beta Std Error: {self.sd_beta}\n"
                f"Residual Variance: {self.res_var:.6f}\n"
                f"Stop Reason: {self.stopreason}")


class ODR:
    """
    Solucionador del problema ODR mediante Levenberg-Marquardt generalizado.
    """
    def __init__(self, 
                 data: Data, 
                 model: Model, 
                 beta0: Optional[Union[List, np.ndarray]] = None, 
                 ifixb: Optional[Union[List, np.ndarray]] = None, 
                 ifixx: Optional[Union[List, np.ndarray]] = None, 
                 job: int = 0, 
                 sstol: float = 1e-8, 
                 partol: float = 1e-8, 
                 maxit: int = 400,
                 **kwargs):  # Acepta kwargs adicionales por compatibilidad de API
        self.data = data
        self.model = model
        
        if beta0 is not None:
            self.beta0 = np.asarray(beta0, dtype=np.float64)
        elif self.model.estimate is not None:
            self.beta0 = np.asarray(self.model.estimate(self.data), dtype=np.float64)
        else:
            raise ValueError("Se debe proveer beta0 o un método de estimación en el modelo.")
            
        self.ifixb = np.asarray(ifixb, dtype=int) if ifixb is not None else None
        self.ifixx = np.asarray(ifixx, dtype=int) if ifixx is not None else None
        self.job = job
        self.sstol = sstol
        self.partol = partol
        self.maxit = maxit
        self.output: Optional[Output] = None

    def run(self) -> Output:
        """Ejecuta la optimización y retorna el objeto Output."""
        x, y = self.data.x, self.data.y
        m, n = self.data.m, self.data.n
        nq = self.data.nq
        np_params = self.beta0.shape[0]

        is_ols = (self.job % 10 == 2) or (self.job % 10 == 3)

        # Máscaras booleanas de parámetros libres
        free_b_mask = np.ones(np_params, dtype=bool) if self.ifixb is None else (self.ifixb != 0)
        free_x_mask = np.zeros((m, n), dtype=bool) if is_ols else (
            np.ones((m, n), dtype=bool) if self.ifixx is None else (self.ifixx != 0)
        )

        n_free_b = np.count_nonzero(free_b_mask)
        n_free_x = np.count_nonzero(free_x_mask)

        p0 = np.concatenate([self.beta0[free_b_mask], np.zeros(n_free_x, dtype=np.float64)])

        # Raíces cuadradas precalculadas
        sqrt_we = np.sqrt(np.maximum(self.data.we, 0.0))
        sqrt_wd = np.sqrt(np.maximum(self.data.wd, 0.0))

        # OPTIMIZACIÓN: Preasignación de memoria para evitar instanciaciones en el bucle interno
        _beta_curr = self.beta0.copy()
        _delta_curr = np.zeros((m, n), dtype=np.float64)

        def residual_function(p: np.ndarray) -> np.ndarray:
            """Función de costo evaluada por least_squares."""
            # Actualización in-place (O(1) allocation)
            _beta_curr[free_b_mask] = p[:n_free_b]
            if n_free_x > 0:
                _delta_curr[free_x_mask] = p[n_free_b:]

            x_plus = x + _delta_curr

            # Evaluación (manejando escalares si M=1 para mayor compatibilidad)
            f_eval = self.model.fcn(_beta_curr, x_plus[0] if m == 1 else x_plus)
            f_eval = np.atleast_2d(np.asarray(f_eval, dtype=np.float64))

            # Residuales ponderados de Y
            res_y = ((y - f_eval) * sqrt_we).ravel()

            # Residuales ponderados de X
            if not is_ols and n_free_x > 0:
                res_x = (_delta_curr * sqrt_wd).ravel()
                return np.concatenate([res_y, res_x])
            return res_y

        # Ejecución del solver interno
        res = least_squares(
            residual_function, p0, method='lm', ftol=self.sstol, xtol=self.partol,
            max_nfev=self.maxit * len(p0)
        )

        # Reconstrucción de resultados
        beta_opt = self.beta0.copy()
        beta_opt[free_b_mask] = res.x[:n_free_b]

        delta_opt = np.zeros((m, n), dtype=np.float64)
        if n_free_x > 0:
            delta_opt[free_x_mask] = res.x[n_free_b:]

        xplus_opt = x + delta_opt
        ynew_opt = np.atleast_2d(np.asarray(
            self.model.fcn(beta_opt, xplus_opt[0] if m == 1 else xplus_opt), 
            dtype=np.float64
        ))

        # Variancia Residual
        deg_freedom = max(1, (n * nq) - n_free_b)
        res_var = np.sum(np.square(res.fun)) / deg_freedom

        # Matriz de Covarianza usando SVD/Pseudoinversa robusta
        J = res.jac[:, :n_free_b]
        try:
            cov_beta_free = la.inv(J.T @ J) * res_var
        except la.LinAlgError:
            cov_beta_free = la.pinv(J.T @ J) * res_var

        cov_beta = np.zeros((np_params, np_params), dtype=np.float64)
        cov_beta[np.ix_(free_b_mask, free_b_mask)] = cov_beta_free
        
        # Clip en 0 para evitar warnings en diagónales numéricamente inestables
        sd_beta = np.sqrt(np.maximum(np.diag(cov_beta), 0.0))

        stop_reasons = {
            1: "Sum of squares convergence", 2: "Parameter convergence",
            3: "Both sum of squares and parameter convergence", 4: "Iteration limit reached"
        }
        info_code = res.status if res.status in stop_reasons else 1

        self.output = Output(
            beta=beta_opt, sd_beta=sd_beta, cov_beta=cov_beta, res_var=res_var,
            rel_error=res.optimality, eps=y - ynew_opt, delta=delta_opt,
            info=info_code, stopreason=stop_reasons.get(info_code, "Optimization completed"),
            xplus=xplus_opt, ynew=ynew_opt
        )

        return self.output


def odr(fcn: Callable, beta0: Union[List, np.ndarray], y: Union[List, np.ndarray], 
        x: Union[List, np.ndarray], we: Optional[np.ndarray] = None, 
        wd: Optional[np.ndarray] = None, fjacb: Optional[Callable] = None, 
        fjacd: Optional[Callable] = None, ifixb: Optional[np.ndarray] = None, 
        ifixx: Optional[np.ndarray] = None, job: int = 0, **kwargs) -> Output:
    """
    Función de atajo (wrapper) compatible con la interfaz tradicional de scipy.odr.odr.
    """
    my_odr = ODR(Data(x, y, we=we, wd=wd), Model(fcn, fjacb=fjacb, fjacd=fjacd), 
                 beta0=beta0, ifixb=ifixb, ifixx=ifixx, job=job, **kwargs)
    return my_odr.run()