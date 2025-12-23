import numpy as np
import scipy as sp
from scipy.optimize import linear_sum_assignment, root, fsolve, least_squares
import math
import Parameter

'''
execute: 执行算法
'''

class Algorithm(object):
    def __init__(self, param: Parameter):
        super().__init__()
        self.param = param

    def WPT(self):
        param = self.param
        param.c_T = np.array([0. for _ in range(param.M)])
        param.a_T = np.array([0 for _ in range(param.M)])
        param.P_T = np.array([0. for _ in range(param.M)])
        param.tau_T = np.array([0. for _ in range(param.M)])

        param.c_T = param.V - (param.B_minus * param.B_scalar * param.mu) @ param.h_D

        j_star = np.argmin(param.c_T * param.P_T_max)
        cp_min = (param.c_T * param.P_T_max)[j_star]

        if cp_min < 0:
            param.a_T[j_star] = 1
            param.P_T[j_star] = param.P_T_max[j_star]
            param.tau_T[j_star] = param.delta_t

    def LocalComputation(self):
        param = self.param
        param.f = np.zeros(param.N, dtype=float)

        f_ub = np.minimum(
            np.cbrt(param.B_remain / (param.kappa * param.delta_t)),
            param.f_max
        )

        stationary_point = np.sqrt(
            param.Q * param.Q_scalar /
            ( 3 * param.kappa * param.phi * np.where(param.B_minus == 0, np.inf, param.B_minus) * param.B_scalar)
        )

        param.f = np.where(
            param.B_minus == 0,
            f_ub,
            np.minimum(f_ub, stationary_point)
        )

        assert np.all(param.f >= 0)
        assert np.all(param.f <= param.f_max)

    def ComputationOffloading(self):
        param       = self.param
        param.a_O   = np.zeros((param.N, param.M), dtype=int)
        param.P_O   = np.zeros(param.N, dtype=float)
        param.tau_O = np.zeros(param.N, dtype=float)
        param.c_O   = np.zeros((param.N, param.M), dtype=float)
        P_O_ub      = np.minimum(param.P_O_max, param.B_remain / param.delta_t)
        w           = np.zeros((param.N, param.M), dtype=float)
        P_O_record  = np.zeros((param.N, param.M), dtype=float)
        f_record    = np.zeros((param.N, param.M), dtype=float)
        f_ub = np.minimum(
            np.cbrt(param.B_remain / (param.kappa * param.delta_t)),
            param.f_max
        )

        for j in range(param.M):
            for i in range(param.N):
                # the objective of this loop is to calculate c_O, and we should not update the real P_O
                # hence, we introduce a local variable P_O_i
                # we add param.ABS_TOL to avoid dividing by zero when B_minus is zero
                P_hat_O_i = ( (param.Q[i] * param.Q_scalar - param.V * param.phi[i] * param.eta) * param.B
                    / (param.B_minus[i] * param.B_scalar * param.v[i] * np.log(2) + param.ABS_TOL ) ) \
                    - ( param.sigma[j]**2 / param.h_U[i][j] )
                P_O_i = np.maximum( 0, np.minimum(P_O_ub[i], P_hat_O_i))
                f_star_i = param.f[i]
                E_L_i = param.kappa[i] * param.f[i]**3 * param.delta_t
                E_O_i = P_O_i * param.delta_t
                if E_L_i + E_O_i > param.B_remain[i]:
                    # obtain f_star_i
                    roots = np.roots([
                        param.v[i] * np.log(2) * param.kappa[i] / param.B,
                        3 * param.kappa[i] * param.phi[i],
                        0,
                        - param.v[i] * np.log(2) / param.B * (param.sigma[j]**2 / param.h_U[i][j] + param.B_remain[i] / param.delta_t) \
                        - param.eta * param.phi[i]
                    ])
                    real_roots = roots[np.isreal(roots)].real
                    positive_real_roots = real_roots[real_roots > 0]
                    # only roots in [0, f_max_i] is meaningful
                    lower_bound, upper_bound = 0, f_ub[i]
                    actual_roots = np.clip(positive_real_roots, lower_bound, upper_bound)
                    if param.debug:
                        assert len(actual_roots) == 1
                    f_star_i = actual_roots[0]

                    # calculate P_O_i
                    P_O_i = (param.B_remain[i] / param.delta_t) - (param.kappa[i] * f_star_i**3)
                    assert P_O_i >= 0 - param.ABS_TOL

                # record the value of P_O_i for this particular (i,j), so that we do not need to re-calculate it afterwards
                P_O_record[i][j] = P_O_i
                f_record[i][j] = f_star_i

                # calculate c_O_ij
                param.c_O[i][j] = ( (param.V * param.eta * param.phi[i] - param.Q[i] * param.Q_scalar) * param.B / param.v[i] * \
                                np.log2(1 + P_O_i * param.h_U[i][j] / param.sigma[j]**2) ) \
                            - param.B_minus[i] * param.B_scalar * P_O_i

        for i in range(param.N):
            for j in range(param.M):
                if param.c_O[i][j] < 0:
                    w[i][j] = param.c_O[i][j] * param.delta_t

        row_ind, col_ind = linear_sum_assignment(w)
        param.a_O[row_ind, col_ind] = 1
        
        for j in range(param.M):
            pos = np.where(col_ind == j)[0][0]
            i_star = row_ind[pos]
            if param.debug:
                assert param.a_O[i_star][j] == 1
            if param.c_O[i_star][j] < 0:
                """
                The implementation here is slightly different with the pseudocode in the paper.
                Specifically, we use P_O_record instead of y.
                But the logic is the same.
                """
                param.tau_O[i_star] = param.delta_t
                param.P_O[i_star] = P_O_record[i_star][j]
                param.f[i_star] = f_record[i_star][j]


    def AdjustVariableValues(self):
        param = self.param
        if not param.a_T.any():
            # no wpt
            return
        else:
            j_star = np.where(param.a_T == 1)[0][0]
            i_star = np.where(param.a_O[:, j_star] == 1)[0][0]
            if param.debug:
                assert param.a_O[i_star, j_star] == 1
            if param.tau_T[j_star] + param.tau_O[i_star] > param.delta_t:
                if param.c_T[j_star] * param.P_T_max[j_star] < param.c_O[i_star, j_star]:
                    param.tau_T[j_star] = param.delta_t
                    param.tau_O[i_star] = 0
                else:
                    param.tau_T[j_star] = 0
                    param.tau_O[i_star] = param.delta_t


    def CheckConstraints(self):
        param = self.param
        assert np.sum(param.a_T) <= 1
        for i in range(param.N):
            E_L_i = param.kappa[i] * param.f[i]**3 * param.tau_L[i]
            E_O_i = param.P_O[i] * param.tau_O[i]
            assert E_L_i + E_O_i <= param.B_remain[i] + param.ABS_TOL
        for i in range(param.N):
            assert np.sum(param.a_O[i]) <= 1
        for j in range(param.M):
            assert param.a_T[j] * param.tau_T[j] + np.sum(param.a_O[:,j] * param.tau_O) <= param.delta_t + param.ABS_TOL
        for j in range(param.M):
            assert 0 <= param.tau_T[j] <= param.delta_t
            assert 0 <= param.P_T[j] <= param.P_T_max[j]
        for i in range(param.N):
            assert 0 - param.ABS_TOL <= param.P_O[i] <= param.P_O_max[i]
            assert 0 - param.ABS_TOL <= param.tau_O[i] <= param.delta_t
            assert 0 - param.ABS_TOL <= param.f[i] <= param.f_max[i]
    
    def Show(self):
        print("B_remain:", end=''); print(self.param.B_remain);
        print("q_T:", end=''); print(self.param.a_T * self.param.P_T * self.param.tau_T);
        print(self.param.a_T); print(self.param.tau_T);
        print("Q:", end=''); print(self.param.Q);
        print("P:", end=''); print(self.param.P_O);
        print("f:", end=''); print(self.param.f);
        
    '''proposed algorithm'''
    def executeProp(self):
        self.WPT()
        self.LocalComputation()
        self.ComputationOffloading()
        self.AdjustVariableValues()
        if self.param.debug:
            self.CheckConstraints()

    '''benchmark-1-LCO'''
    def executeLCO(self):
        self.WPT()
        self.LocalComputation()
        # guarantee that the actually processed data do not exceed the remaining data in the queue
        self.param.f = np.minimum(self.param.f, self.param.Q * self.param.phi / self.param.delta_t)
        if self.param.debug:
            self.CheckConstraints()

    '''benchmark-2-FO'''
    def executeFO(self):
        self.param.f_max = np.zeros(self.param.N, dtype=float)
        self.WPT()
        #self.LocalComputation()
        self.ComputationOffloading()
        """
        Because the FO has limited data-processing capability, the queue length Q grows rapidly when the arrival rate is high. 
        This triggers the AdjustVariableValues procedure to allocate time to data offloading. 
        However, doing so prevents the WD from harvesting energy, which further accelerates the growth of Q, creating a vicious cycle. 
        There are two possible solutions to this problem. One is to modify system parameters, such as increasing B_max. 
        The other is to guarantee a minimum amount of WPT time. Here, we adopt the second approach.
        """
        param = self.param
        if not param.a_T.any():
            # no wpt
            return
        else:
            j_star = np.where(param.a_T == 1)[0][0]
            i_star = np.where(param.a_O[:, j_star] == 1)[0][0]
            if param.tau_T[j_star] + param.tau_O[i_star] > param.delta_t:
                param.tau_T[j_star] = param.delta_t
                param.tau_O[i_star] = 0
        if self.param.debug:
            self.CheckConstraints()

    '''benchmark-3-Greey'''
    def executeGreedy(self):
        param = self.param
        param.c_T = np.array([0. for _ in range(param.M)])
        param.a_T = np.array([0 for _ in range(param.M)])
        param.P_T = np.array([0. for _ in range(param.M)])
        param.tau_T = np.array([0. for _ in range(param.M)])
        param.f = np.zeros(param.N, dtype=float)
        param.a_O   = np.zeros((param.N, param.M), dtype=int)
        param.P_O   = np.zeros(param.N, dtype=float)
        param.tau_O = np.zeros(param.N, dtype=float)
        param.c_O   = np.zeros((param.N, param.M), dtype=float)
        
        self.WPT()
        """
        Each WD is associated with the AP that has the best channel gain.
        The offloading time is solely allocated to the WD with the largest Q_i(t).
        Each WD tries to process as much data as possible using the remaining energy.
        """
        # calculate AP-WD association
        association = np.zeros((param.N, param.M), dtype=int)
        best_ap = np.argmax(param.h_U, axis=1)
        association[np.arange(param.N), best_ap] = 1

        # allocate offloading time to the WD with the largest Q_i(t)
        for j in range(param.M):
            if param.a_T[j] == 1:
                continue
            associated_wd = np.where(association[:,j] == 1)[0]
            if associated_wd.size == 0:
                continue
            Q_values = param.Q[associated_wd]
            Q_max_index = np.argmax(Q_values)
            wd_index = associated_wd[Q_max_index]
            param.a_O[wd_index, j] = 1
            param.tau_O[wd_index] = param.delta_t

        # greedily determine f_i(t) and P^O_i(t)
        for i in range(param.N):
            f_ub = min(
                np.cbrt(param.B_remain[i] / (param.kappa[i] * param.delta_t)),
                param.Q[i] * param.phi[i] / param.delta_t,
                param.f_max[i]
            )
            if np.sum(param.a_O[i, :]) == 0:
                # only local computing
                param.f[i] = f_ub
            else:
                # both local computing and data offloading
                """
                Greedy seeks to maximize the processed bits so we removed the energy consumption term at the AP
                """
                roots = np.roots([
                    param.v[i] * np.log(2) * param.kappa[i] / param.B,
                    3 * param.kappa[i] * param.phi[i],
                    0,
                    - param.v[i] * np.log(2) / param.B * (param.sigma[j]**2 / param.h_U[i][j] + param.B_remain[i] / param.delta_t)
                ])
                real_roots = roots[np.isreal(roots)].real
                positive_real_roots = real_roots[real_roots > 0]
                lower_bound, upper_bound = 0, f_ub
                actual_roots = np.clip(positive_real_roots, lower_bound, upper_bound)
                param.f[i] = actual_roots[0]
                param.P_O[i] = np.minimum(
                    (param.B_remain[i] / param.delta_t) - (param.kappa[i] * param.f[i]**3),
                    param.P_O_max[i]
                )
                assert param.P_O[i] >= 0 - param.ABS_TOL

        if self.param.debug:
            self.CheckConstraints()
