import numpy as np
import Parameter

class Environment(object):
    def __init__(self, param: Parameter):
        super().__init__()
        self.param = param
        self.init()

    def init(self):
        param = self.param
        param.D_L      = np.zeros(param.N, dtype=float)
        param.D_O      = np.zeros(param.N, dtype=float)
        param.Q        = np.zeros(param.N, dtype=float)
        param.A        = np.random.uniform(param.A_lb, param.A_ub)
        param.B_remain = np.zeros(param.N, dtype=float)
        param.B_minus  = param.B_max
        #param.B_remain = 0.5 * param.B_max
        #param.B_minus  = param.B_max - param.B_remain

    def step(self):
        param = self.param
        param.E_H = np.zeros(param.N, dtype=float)
        param.D_O = np.zeros(param.N, dtype=float)
        param.A = np.random.uniform(param.A_lb, param.A_ub)

        # energy related
        param.E_H = param.mu * (param.h_D @ (param.a_T * param.P_T * param.tau_T))
        param.E_L = param.kappa * param.f**3 * param.tau_L
        param.E_O = param.P_O * param.tau_O
        param.B_remain = np.minimum(
            param.B_max,
            param.B_remain - param.E_L - param.E_O + param.E_H
        )
        assert np.all(param.B_remain >= 0 - param.ABS_TOL)
        """
        Due to computation accuracy, the value of B_remain[i] may be a very small number (e.g., +/- 1e-22) when it is supposed to zero
        If this small number is negative, it will yield a negative f_ub and cause the checkConstraints fail
        To avoid this problem, we manually set all negative small values in B_remain to zero
        """
        param.B_remain[param.B_remain < 0] = 0
        param.B_minus = param.B_max - param.B_remain

        # queue related
        param.D_L = np.minimum(param.Q_act, param.f * param.tau_L / param.phi)
        param.D_O = np.sum(
            param.a_O * param.B * param.tau_O[:, None] / param.v[:, None] *
            np.log2(1 + param.P_O[:, None] * param.h_U / param.sigma**2),
            axis=1
        )
        """
        The calculation of E_C depends on D_O.
        The actual D_O should not exceed Q - D_L.
        """
        param.D_O = np.minimum(param.D_O, param.Q_act - param.D_L)
        param.Q_act = np.maximum(0, param.Q_act - param.D_L - param.D_O) + param.A

        # update place-holder
        if param.placeholder:
            param.q_hat = (1-param.alpha) * param.q_hat + param.alpha * param.Q
            param.q = np.maximum(0, param.q_hat - param.r * np.log(param.V)**2)
            param.Q = param.Q_act + param.q
        else:
            param.Q = param.Q_act

        param.calChannelGain()
