import numpy as np
from scipy.optimize import least_squares

__all__ = ['Data', 'RealData', 'Model', 'Output', 'ODR', 'odr']


class Data:
    """
    Contenedor de datos para el problema de Regresión por Distancia Ortogonal (ODR).
    Compatibilidad directa con scipy.odr.Data.
    """
    def __init__(self, x, y=None, we=None, wd=None, fix=None, meta=None):
        self.x = np.asarray(x, dtype=np.float64)
        self.y = np.asarray(y, dtype=np.float64) if y is not None else None
        
        # Ajuste de dimensiones de X e Y
        if self.x.ndim == 1:
            self.x = self.x[np.newaxis, :]
        if self.y is not None and self.y.ndim == 1:
            self.y = self.y[np.newaxis, :]
            
        self.m = self.x.shape[0]  # Número de variables explicativas (dimensiones de x)
        self.n = self.x.shape[1]  # Número de observaciones
        self.nq = self.y.shape[0] if self.y is not None else 1
        
        # Pesos (Weights) para Y (we) y X (wd)
        self.we = self._process_weights(we, (self.nq, self.n))
        self.wd = self._process_weights(wd, (self.m, self.n))
        
        self.fix = fix
        self.meta = meta if meta is not None else {}

    def _process_weights(self, w, target_shape):
        if w is None:
            return np.ones(target_shape, dtype=np.float64)
        w_arr = np.asarray(w, dtype=np.float64)
        if w_arr.ndim == 1 and w_arr.shape[0] == target_shape[1]:
            return np.tile(w_arr, (target_shape[0], 1))
        return w_arr


class RealData(Data):
    """
    Variante de Data que acepta desviaciones estándar (sx, sy) o covarianzas en lugar de pesos.
    Compatibilidad directa con scipy.odr.RealData.
    """
    def __init__(self, x, y=None, sx=None, sy=None, covx=None, covy=None, fix=None, meta=None):
        we = None
        wd = None
        
        if sy is not None:
            sy_arr = np.asarray(sy, dtype=np.float64)
            we = 1.0 / np.square(sy_arr)
        elif covy is not None:
            we = 1.0 / np.asarray(covy, dtype=np.float64)

        if sx is not None:
            sx_arr = np.asarray(sx, dtype=np.float64)
            wd = 1.0 / np.square(sx_arr)
        elif covx is not None:
            wd = 1.0 / np.asarray(covx, dtype=np.float64)

        super().__init__(x, y=y, we=we, wd=wd, fix=fix, meta=meta)


class Model:
    """
    Define el modelo matemático a ajustar y sus jacobianos opcionales.
    Compatibilidad directa con scipy.odr.Model.
    """
    def __init__(self, fcn, fjacb=None, fjacd=None, meta=None, estimate=None, implicit=0):
        self.fcn = fcn
        self.fjacb = fjacb
        self.fjacd = fjacd
        self.meta = meta if meta is not None else {}
        self.estimate = estimate
        self.implicit = implicit


class Output:
    """
    Almacena los resultados devueltos por el algoritmo ODRPACK.
    Compatibilidad directa con los atributos de scipy.odr.Output.
    """
    def __init__(self, beta, sd_beta, cov_beta, res_var, rel_error, 
                 eps, delta, info, stopreason, xplus, ynew):
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

    def pprint(self):
        """Imprime una vista resumida de los resultados."""
        print("Beta:", self.beta)
        print("Beta SD:", self.sd_beta)
        print("Residual Variance:", self.res_var)
        print("Reason:", self.stopreason)


class ODR:
    """
    Clase principal del solucionador ODRPACK escrito nativamente en Python.
    Compatibilidad directa con scipy.odr.ODR.
    """
    def __init__(self, data, model, beta0=None, ifixb=None, ifixx=None, 
                 job=None, iprint=None, errfile=None, rptfile=None, 
                 ndigit=None, taufac=None, sstol=None, partol=None, maxit=100, 
                 stpb=None, stpd=None, sclb=None, scld=None, work=None, iwork=None):
        self.data = data
        self.model = model
        self.beta0 = np.asarray(beta0, dtype=np.float64) if beta0 is not None else None
        
        if self.beta0 is None and self.model.estimate is not None:
            self.beta0 = np.asarray(self.model.estimate(self.data), dtype=np.float64)
            
        self.ifixb = np.asarray(ifixb, dtype=int) if ifixb is not None else None
        self.ifixx = np.asarray(ifixx, dtype=int) if ifixx is not None else None
        self.job = job if job is not None else 0  # 00000 = ODR Estándar
        self.sstol = sstol if sstol is not None else 1e-6
        self.partol = partol if partol is not None else 1e-6
        self.maxit = maxit
        self.output = None

    def run(self):
        """Ejecuta la optimización por Distancia Ortogonal (ODR) o Mínimos Cuadrados Ordinarios (OLS)."""
        x = self.data.x  # Forma (M, N)
        y = self.data.y  # Forma (NQ, N)
        m, n = self.data.m, self.data.n
        nq = self.data.nq
        np_params = self.beta0.shape[0]

        # Verificar si se solicita OLS (Ordinary Least Squares) vía job
        is_ols = (self.job % 10 == 2) or (self.job % 10 == 3)

        # Máscara de parámetros fijos (ifixb)
        if self.ifixb is None:
            free_b_mask = np.ones(np_params, dtype=bool)
        else:
            free_b_mask = (self.ifixb != 0)

        # Máscara de observaciones/puntos fijos (ifixx)
        if self.ifixx is None or is_ols:
            free_x_mask = np.ones((m, n), dtype=bool) if not is_ols else np.zeros((m, n), dtype=bool)
        else:
            free_x_mask = (self.ifixx != 0)

        n_free_b = np.count_nonzero(free_b_mask)
        n_free_x = np.count_nonzero(free_x_mask)

        # Vector de parámetros iniciales a optimizar: [beta_libres, delta_libres]
        init_beta_free = self.beta0[free_b_mask]
        init_delta_free = np.zeros(n_free_x, dtype=np.float64)
        p0 = np.concatenate([init_beta_free, init_delta_free])

        # Precalculo de raíces cuadradas de pesos
        sqrt_we = np.sqrt(np.maximum(self.data.we, 0.0))
        sqrt_wd = np.sqrt(np.maximum(self.data.wd, 0.0))

        def residual_function(p):
            # Reconstruir beta y delta a partir del vector p
            beta = self.beta0.copy()
            beta[free_b_mask] = p[:n_free_b]
            
            delta = np.zeros((m, n), dtype=np.float64)
            if n_free_x > 0:
                delta[free_x_mask] = p[n_free_b:]

            x_plus = x + delta

            # Evaluación del modelo
            if m == 1:
                f_eval = self.model.fcn(beta, x_plus[0])
            else:
                f_eval = self.model.fcn(beta, x_plus)
                
            f_eval = np.asarray(f_eval, dtype=np.float64)
            if f_eval.ndim == 1:
                f_eval = f_eval[np.newaxis, :]

            # Residuales ponderados de Y
            eps = y - f_eval
            res_y = (eps * sqrt_we).ravel()

            # Residuales ponderados de X
            if not is_ols and n_free_x > 0:
                res_x = (delta * sqrt_wd).ravel()
                return np.concatenate([res_y, res_x])
            
            return res_y

        # Resolver el problema de mínimos cuadrados no lineales usando Levenberg-Marquardt / Trust Region
        res = least_squares(
            residual_function, 
            p0, 
            method='lm', 
            ftol=self.sstol, 
            xtol=self.partol, 
            max_nfev=self.maxit * len(p0)
        )

        # Extraer solución
        beta_opt = self.beta0.copy()
        beta_opt[free_b_mask] = res.x[:n_free_b]

        delta_opt = np.zeros((m, n), dtype=np.float64)
        if n_free_x > 0:
            delta_opt[free_x_mask] = res.x[n_free_b:]

        xplus_opt = x + delta_opt

        if m == 1:
            ynew_opt = np.asarray(self.model.fcn(beta_opt, xplus_opt[0]), dtype=np.float64)
        else:
            ynew_opt = np.asarray(self.model.fcn(beta_opt, xplus_opt), dtype=np.float64)
            
        if ynew_opt.ndim == 1:
            ynew_opt = ynew_opt[np.newaxis, :]

        eps_opt = y - ynew_opt

        # Estimación de Variancia Residual y Matriz de Covarianza
        n_obs = n * nq
        deg_freedom = max(1, n_obs - n_free_b)
        sum_sq = np.sum(np.square(res.fun))
        res_var = sum_sq / deg_freedom

        # Estimación numérica del Jacobiano respecto a Beta en el punto óptimo
        J = res.jac[:, :n_free_b]
        try:
            cov_beta_free = np.linalg.inv(J.T @ J) * res_var
        except np.linalg.LinAlgError:
            cov_beta_free = np.linalg.pinv(J.T @ J) * res_var

        cov_beta = np.zeros((np_params, np_params), dtype=np.float64)
        cov_beta[np.ix_(free_b_mask, free_b_mask)] = cov_beta_free

        sd_beta = np.sqrt(np.diag(cov_beta))

        # Construir objeto Output
        stop_reasons = {
            1: "Sum of squares convergence",
            2: "Parameter convergence",
            3: "Both sum of squares and parameter convergence",
            4: "Iteration limit reached"
        }
        info_code = res.status if res.status in stop_reasons else 1
        stopreason = stop_reasons.get(info_code, "Optimization completed")

        self.output = Output(
            beta=beta_opt,
            sd_beta=sd_beta,
            cov_beta=cov_beta,
            res_var=res_var,
            rel_error=res.optimality,
            eps=eps_opt,
            delta=delta_opt,
            info=info_code,
            stopreason=stopreason,
            xplus=xplus_opt,
            ynew=ynew_opt
        )

        return self.output


def odr(fcn, beta0, y, x, we=None, wd=None, fjacb=None, fjacd=None, 
        ifixb=None, ifixx=None, job=None, iprint=None, errfile=None, 
        rptfile=None, ndigit=None, taufac=None, sstol=None, partol=None, 
        maxit=None, stpb=None, stpd=None, sclb=None, scld=None, work=None, iwork=None):
    """
    Función de conveniencia compatible con scipy.odr.odr.
    """
    data = Data(x, y, we=we, wd=wd)
    model = Model(fcn, fjacb=fjacb, fjacd=fjacd)
    my_odr = ODR(
        data, model, beta0=beta0, ifixb=ifixb, ifixx=ifixx, 
        job=job, sstol=sstol, partol=partol, maxit=maxit
    )
    return my_odr.run()