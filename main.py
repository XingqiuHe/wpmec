from Parameter import Parameter
from Environment import Environment
from Algorithm import Algorithm
from scipy.io import savemat
import numpy as np
import sys
import matplotlib.pyplot as plt

def simulation(alg='Prop', V=10000, N=30, M=5, seed=0, placeholder=False, arrivalRateDiscount=1.0):
    if alg not in ['Prop', 'FO', 'LCO', 'Greedy']:
        print("Algorithm not implemented!")
        return 0
    if alg == 'Prop':
        placeholder = True
    param = Parameter(V=V, N=N, M=M, seed=seed, placeholder=placeholder)
    param.A_lb = param.A_lb * arrivalRateDiscount
    param.A_ub = param.A_ub * arrivalRateDiscount
    env = Environment(param=param)
    algorithm = Algorithm(param=param)

    E_history = np.zeros(param.T, dtype=float)
    Q_history = np.zeros(param.T, dtype=np.int64)
    A_history = np.zeros(param.T, dtype=np.int64)
    
    for t in range(param.T):
        param.t = t
        if alg == 'Prop':
            algorithm.executeProp()
        elif alg == 'LCO':
            algorithm.executeLCO()
        elif alg == 'FO':
            algorithm.executeFO()
        else:
            algorithm.executeGreedy()
        env.step()

        E_T = np.sum(param.a_T * param.P_T * param.tau_T)
        E_C = param.eta * np.sum(param.a_O * (param.phi[:, None] * param.D_O[:, None]))
        E_history[t] = E_T + E_C
        Q_history[t] = np.sum(param.Q_act)
        A_history[t] = np.sum(param.A)

    E_av = np.mean(E_history)
    L_av = np.sum(Q_history) / np.sum(A_history) * param.delta_t * 1000

    return E_av, L_av


mode = {1}

if 1 in mode:
    print('Queue Dynamics')
    data = {'Q_PH':[], 'Q_act_PH':[], 'q_PH':[], 'B_remain_PH':[], 'Q':[], 'B_remain':[]}
    param = Parameter(V=10000, placeholder=True)
    env = Environment(param=param)
    algorithm = Algorithm(param=param)
    for t in range(param.T):
        param.t = t
        data['Q_PH'].append(param.Q)
        data['Q_act_PH'].append(param.Q_act)
        data['q_PH'].append(param.q)
        data['B_remain_PH'].append(param.B_remain)
        algorithm.executeProp()
        env.step()

    param = Parameter(V=10000, placeholder=False)
    env = Environment(param=param)
    algorithm = Algorithm(param=param)
    for t in range(param.T):
        param.t = t
        data['Q'].append(param.Q)
        data['B_remain'].append(param.B_remain)
        algorithm.executeProp()
        env.step()

    savemat('data/queue_dynamics.mat', data, False)

if 2 in mode:
    print("Different V")
    arrivalRateDiscount = 0.75
    V_list = list(range(5000, 15001, 1000))
    alg_list = ['Prop', 'LCO', 'FO', 'Greedy']
    seed_list = list(range(10))
    E_data = dict.fromkeys(alg_list)
    E_data = {k: [] for k in E_data}
    L_data = dict.fromkeys(alg_list)
    L_data = {k: [] for k in L_data}
    for V in V_list:
        for alg in alg_list:
            for seed in seed_list:
                E_av, L_av = simulation(alg=alg, V=V, seed=seed, arrivalRateDiscount=arrivalRateDiscount)
                E_data[alg].append(E_av)
                L_data[alg].append(L_av)

    data = {
        'seed_list': seed_list,
        'V_list': V_list,
        'E_data': E_data,
        'L_data': L_data
    }
    savemat('data/different_V.mat', data, False)

if 3 in mode:
    print("Arrival Rate")
    discount_list = [round(i * 0.1, 1) for i in range(1, 11)]
    alg_list = ['Prop', 'LCO', 'FO', 'Greedy']
    seed_list = list(range(10))
    E_data = dict.fromkeys(alg_list)
    E_data = {k: [] for k in E_data}
    L_data = dict.fromkeys(alg_list)
    L_data = {k: [] for k in L_data}
    for discount in discount_list:
        for alg in alg_list:
            for seed in seed_list:
                E_av, L_av = simulation(alg=alg, seed=seed, arrivalRateDiscount=discount)
                E_data[alg].append(E_av)
                L_data[alg].append(L_av)

    data = {
        'seed_list': seed_list,
        'discount_list': discount_list,
        'E_data': E_data,
        'L_data': L_data
    }
    savemat('data/arrival_rate.mat', data, False)

if 4 in mode:
    print("Different M")
    M_list = [i for i in range(1, 11)]
    alg_list = ['Prop', 'LCO', 'FO', 'Greedy']
    rate_discount = 0.75
    seed_list = list(range(30))
    E_data = dict.fromkeys(alg_list)
    E_data = {k: [] for k in E_data}
    L_data = dict.fromkeys(alg_list)
    L_data = {k: [] for k in L_data}
    for M in M_list:
        for alg in alg_list:
            for seed in seed_list:
                E_av, L_av = simulation(alg=alg, M=M, seed=seed, arrivalRateDiscount=rate_discount)
                E_data[alg].append(E_av)
                L_data[alg].append(L_av)

    data = {
        'seed_list': seed_list,
        'M_list': M_list,
        'E_data': E_data,
        'L_data': L_data
    }
    savemat('data/different_M.mat', data, False)

if 5 in mode:
    print("Different N")
    N_list = [10+i*4 for i in range(0, 11)]
    alg_list = ['Prop', 'LCO', 'FO', 'Greedy']
    rate_discount = 0.75
    seed_list = list(range(30))
    E_data = dict.fromkeys(alg_list)
    E_data = {k: [] for k in E_data}
    L_data = dict.fromkeys(alg_list)
    L_data = {k: [] for k in L_data}
    for N in N_list:
        for alg in alg_list:
            for seed in seed_list:
                E_av, L_av = simulation(alg=alg, N=N, seed=seed, arrivalRateDiscount=rate_discount)
                E_data[alg].append(E_av)
                L_data[alg].append(L_av)

    data = {
        'seed_list': seed_list,
        'N_list': N_list,
        'E_data': E_data,
        'L_data': L_data
    }
    savemat('data/different_N.mat', data, False)
