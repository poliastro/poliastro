from scipy.integrate import OdeSolver, DenseOutput
from copy import copy
import numpy as np
from poliastro import integrator_params
from scipy.integrate._ivp.common import (warn_extraneous, validate_max_step, validate_tol, select_initial_step)
from poliastro.jit import jit


def validate_max_nsteps(max_nsteps):
    if max_nsteps <= 0:
        raise ValueError("`max_nsteps` must be positive.")
    return max_nsteps


def validate_safety_factor(safety_factor):
    if safety_factor >= 1.0 or safety_factor <= 1e-4:
        raise ValueError("`safety_factor` must lie within 1e-4 and 1.0.")
    return safety_factor


def validate_beta_stabilizer(beta_stabilizer):
    if beta_stabilizer < 0 or beta_stabilizer > 0.2:
        raise ValueError("`beta_stabilizer` must lie within 0 and 0.2.")
    return beta_stabilizer


class DOP835(OdeSolver):
    A = integrator_params.A
    C = integrator_params.C
    B = integrator_params.B
    E = integrator_params.E
    BHH = integrator_params.BHH
    D = integrator_params.D

    def __init__(self, fun, t0, y0, t_bound, max_step=np.inf,
                 rtol=1e-7, atol=1e-12, safety_factor=0.9,
                 min_step_change=0.333, max_step_change=6.0,
                 beta_stabilizer=0.00, max_nsteps=100000,
                 vectorized=False, **extraneous):
        warn_extraneous(extraneous)
        super(DOP835, self).__init__(fun, t0, y0, t_bound, vectorized,
                                     support_complex=True)
        self.y_old = None
        self.max_step = validate_max_step(max_step)
        self.beta_stabilizer = validate_beta_stabilizer(beta_stabilizer)
        self.max_nsteps = validate_max_nsteps(max_nsteps)
        self.safety_factor = validate_safety_factor(safety_factor)
        self.rtol, self.atol = validate_tol(rtol, atol, self.n)
        self.min_step_change = min_step_change
        self.max_step_change = max_step_change
        self.order = 7

        self.f = self.fun(self.t, self.y)
        self.h_abs = select_initial_step(self.fun, self.t, self.y, self.f, self.direction,
                                         self.order, self.rtol, self.atol)
        self.nfev += 2

        self.n_steps = 0
        self.n_accepted = 0
        self.n_rejected = 0
        self.factor_old = 1e-4  # Lund-stabilization factor
        self.K = np.zeros((16, self.n))
        self.interpolation = np.zeros((8, self.n))

    def _step_impl(self):
        t = self.t
        y = self.y
        f = self.f
        K = self.K

        max_step = self.max_step
        rtol = self.rtol
        atol = self.atol

        min_step = 10 * np.abs(np.nextafter(t, self.direction * np.inf) - t)

        if self.h_abs > max_step:
            h_abs = max_step
        elif self.h_abs < min_step:
            h_abs = min_step
        else:
            h_abs = self.h_abs

        order = self.order
        step_accepted = False

        while not step_accepted:
            if self.h_abs < min_step:
                return False, self.TOO_SMALL_STEP

            h = self.h_abs * self.direction
            t_new = t + h

            if self.direction * (t_new - self.t_bound) > 0:
                t_new = self.t_bound

            h = t_new - t
            h_abs = np.abs(h)

            K[0] = f
            for s in range(1, 12):
                a, c = self.A[s], self.C[s]
                dy = np.dot(K[:s].T, a) * h
                K[s] = self.fun(t + c * h, y + dy)
            self.nfev += 11

            f_B = np.dot(K[:12].T, self.B)
            y_final = y + h * f_B

            scale = atol + np.maximum(np.abs(y), np.abs(y_final)) * rtol
            err_BHH = f_B - self.BHH[0] * K[0] - self.BHH[1] * K[8] - self.BHH[2] * K[11]
            err_BHH = np.sum((err_BHH / scale) ** 2)

            err_E = np.dot(K[:12].T, self.E)
            err_E = np.sum((err_E / scale) ** 2)

            denominator = err_E + 1e-2 * err_BHH
            err_E = h_abs * err_E / np.sqrt(self.n * denominator)

            err_exp = err_E ** (0.125 - self.beta_stabilizer * 0.2)
            dh_factor = err_exp / (self.factor_old ** self.beta_stabilizer)
            dh_factor = np.max([1.0 / self.max_step_change,
                                np.min([1.0 / self.min_step_change,
                                        dh_factor / self.safety_factor])])
            h_new_abs = h_abs / dh_factor

            if err_E < 1.0:
                step_accepted = True
                self.factor_old = np.max([err_E, 1e-4])
                self.n_accepted += 1
                K[12] = self.fun(t + h, y_final)

                for s in range(13, 16):
                    a, c = self.A[s], self.C[s]
                    dy = np.dot(K[:s].T, a) * h
                    K[s] = self.fun(t + c * h, y + dy)
                self.nfev += 4

                # prepare the dense output
                self.interpolation[0] = y
                self.interpolation[1] = y_final - y
                self.interpolation[2] = h * K[0] - self.interpolation[1]
                self.interpolation[3] = self.interpolation[1] - h * K[12] - self.interpolation[2]
                for n in range(4, 8):
                    self.interpolation[n] = h * np.dot(K[:16].T, self.D[n])

                self.y_old = y
                self.t = t_new
                self.y = y_final
                self.f = K[12]
                self.h_abs = h_new_abs

                return True, None
            else:
                self.n_rejected += 1
                self.h_abs /= np.min([1.0 / self.min_step_change,
                                      err_exp / self.safety_factor])

    def _dense_output_impl(self):
        return DOP835DenseOutput(self.t_old, self.t, self.y_old, self.interpolation)


def get_coeffs(s):
    coeffs = np.zeros((8))
    s_back = 1.0 - s
    coeffs[0] = 1.0
    for i in range(7):
        if i % 2 == 0:
            coeffs[i + 1] = coeffs[i] * s
        else:
            coeffs[i + 1] = coeffs[i] * s_back
    return np.array(coeffs)


class DOP835DenseOutput(DenseOutput):
    def __init__(self, t_old, t_new, y_old, interpolation):
        super(DOP835DenseOutput, self).__init__(t_old, t_new)
        self.h = t_new - t_old
        self.interpolation = copy(interpolation)
        self.y_old = y_old

    def _call_impl(self, t_eval):
        s = (t_eval - self.t_old) / self.h
        coeffs = get_coeffs(s)
        return np.dot(self.interpolation.T, coeffs)
