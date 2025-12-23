import numpy as np
import os
from scipy.io import savemat

'''
Parameter List:
    M              = 5                                    ; number of APs
    N              = 30                                   ; number of WDs
    T              = 10000                                ; maximum time slots
    delta_t        = 10 ms                                ; duration of each time slot
    B              = 0.1 MHz                              ; wireless bandwidth
    A[i]           = random~[1, 2] kb                     ; arrived computation data

    mu[i]          = 0.51                                 ; energy conversion efficiency
    k[i]           = 10^(-28)                             ; energy efficiency coefficient of cpu
    phi[i]         = 1000 cycles/bit                      ; number of cpu cycles required to process one bit of data
    sigma[j]^2     = 10^(-9) W                            ; noise power at APs
    v[i]           = 1.1                                  ; wireless communication overhead
    eta            = 10^(-9) J                            ; energy consumption per cycle in MEC server

    f_max[i]       = 0.5 GHz                              ; maximum cpu frequency of WDs
    P_T_max[j]     = 3 W                                  ; maximum power of WPT
    P_O_max[i]     = 0.1 W                                ; maximum power of computation offloading
    B_max[i]       = 2e-3 J                               ; battery capacity

    Randomly generate the locations of WDs and APs in a 10m*10m square
    theta_U        = 5e-4                                 ; uplink gain coefficient
    theta_D        = 1e-3                                 ; downlink gain coefficient
    dist[i][k]     = random                               ; distance between WD i and AP j
    h_random[i][j] = CN(0,1)                              ; complex normal distribution
    h_U[i][j]      = theta_U * dist[i][k]^(-2) * h_random ; uplink channel gain
    h_D[i][j]      = theta_D * dist[i][k]^(-2) * h_random ; downlink channel gain
'''


class Parameter(object):
    def __init__(self, V=1000, N=30, M=5, seed=0, placeholder=False):
        super().__init__()

        '''parameters'''

        self.mapsize      = 10
        self.debug        = 0
        self.placeholder  = placeholder
        self.ABS_TOL      = 1e-10
        self.t            = 0
        self.seed         = seed
        
        # scalars are used to unifies the magnitudes of queue lengths
        self.B_scalar = 1e10
        self.Q_scalar = 3e-7

        self.delta_t = 0.01
        self.B       = 1e5
        self.N       = N
        self.M       = M
        self.T       = 10000
        self.alpha   = 0.0003
        self.r       = 50

        self.V       = V
        self.mu      = np.full(self.N, 0.51)
        self.kappa   = np.full(self.N, 1e-28)
        self.phi     = np.full(self.N, 1e3)
        self.v       = np.full(self.N, 1.1)
        self.sigma   = np.full(self.M, 10**(-4.5))
        self.eta     = 1e-9
        self.theta_U = 5e-4
        self.theta_D = 1e-3

        self.f_max   = np.full(self.N, 5e8)
        self.P_O_max = np.full(self.N, 0.1)
        self.P_T_max = np.full(self.M, 3.0)
        self.B_max   = np.full(self.N, 2e-3)
        self.tau_L   = np.full(self.N, self.delta_t)
        self.A_lb    = np.full(self.N, 1000)
        self.A_ub    = np.full(self.N, 2000)

        '''network condition'''

        self.locWD = np.zeros((self.N, 2), dtype=float)
        self.locAP = np.zeros((self.M, 2), dtype=float)
        self.h_D   = np.zeros((self.N, self.M), dtype=float)
        self.h_U   = np.zeros((self.N, self.M), dtype=float)

        '''variables'''

        # updated by the environment
        self.D_L   = np.zeros(self.N, dtype=np.int64)
        self.D_O   = np.zeros(self.N, dtype=np.int64)
        self.Q     = np.zeros(self.N, dtype=np.int64)
        self.A     = np.zeros(self.N, dtype=np.int64)
        self.Q_act = np.zeros(self.N, dtype=np.int64)
        self.q_hat = np.zeros(self.N, dtype=np.int64)
        self.q     = np.zeros(self.N, dtype=np.int64)

        # B_remain is the remaining energy in the battery
        # it is denoted by B_i(t) in the paper
        # we name it B_remain to distinguish it from the spectrum bandwidth
        self.B_remain = np.zeros(self.N, dtype=float)
        self.B_minus  = np.zeros(self.N, dtype=float)

        self.E_H = np.zeros(self.N, dtype=float)
        self.E_L = np.zeros(self.N, dtype=float)
        self.E_O = np.zeros(self.N, dtype=float)

        # updated by the algorithm
        self.c_T   = np.zeros(self.M, dtype=float)
        self.a_T   = np.zeros(self.M, dtype=int)
        self.P_T   = np.zeros(self.M, dtype=float)
        self.tau_T = np.zeros(self.M, dtype=float)

        self.c_O   = np.zeros((self.N, self.M), dtype=float)
        self.f     = np.zeros(self.N, dtype=float)
        self.tau_O = np.zeros(self.N, dtype=float)
        self.a_O   = np.zeros((self.N, self.M), dtype=int)
        self.P_O   = np.zeros(self.N, dtype=float)

        # generate map and calculate channel gain
        np.random.seed(self.seed)
        self.newMap()


    def newMap(self):
        filename = 'data/MAP_' + str(self.M) + '_' + str(self.N) + '_' + str(self.seed) + '.npz'
        if os.path.exists(filename):
            self.loadMap(filename)
        else:
            self.generateMap()
            self.saveMap(filename)
        self.calChannelGain()

    def generateMap(self):
        # Generate APs
        """
        Finding APs' coordinates to minimize the maximum distance from any given location is the Circle Covering Problem,
        which is actually hard to solve and admits no analytical solutions.
        Hence, I just hardcore the optimal coordinates for some small numbers of APs.
        For larger numbers of APs, one can obtain approximate solutions via Lloyd's Algorithm.
        """
        OPTIMAL_COORDS = {
            1: [
                (0.500, 0.500)
            ],
            2: [
                (0.250, 0.500), (0.750, 0.500)
            ],
            3: [
                (0.500, 0.634),
                (0.146, 0.146),
                (0.854, 0.146)
            ],
            4: [
                (0.250, 0.250), (0.750, 0.250),
                (0.250, 0.750), (0.750, 0.750)
            ],
            5: [
                (0.500, 0.500),
                (0.195, 0.195), (0.805, 0.195),
                (0.195, 0.805), (0.805, 0.805)
            ],
            6: [
                (0.250, 0.167), (0.750, 0.167),
                (0.250, 0.500), (0.750, 0.500),
                (0.250, 0.833), (0.750, 0.833)
            ],
            7: [
                (0.500, 0.500),
                (0.150, 0.250), (0.850, 0.250),
                (0.150, 0.750), (0.850, 0.750),
                (0.500, 0.900), (0.500, 0.100)
            ],
            8: [
                (0.190, 0.190), (0.810, 0.190),
                (0.190, 0.810), (0.810, 0.810),
                (0.500, 0.350), (0.500, 0.650),
                (0.150, 0.500), (0.850, 0.500)
            ],
            9: [
                (0.167, 0.167), (0.500, 0.167), (0.833, 0.167),
                (0.167, 0.500), (0.500, 0.500), (0.833, 0.500),
                (0.167, 0.833), (0.500, 0.833), (0.833, 0.833)
            ],
            10: [
                (0.160, 0.200), (0.500, 0.200), (0.840, 0.200),
                (0.250, 0.500), (0.750, 0.500),
                (0.120, 0.500), (0.880, 0.500),
                (0.160, 0.800), (0.500, 0.800), (0.840, 0.800)
            ]
        }
        if self.M not in OPTIMAL_COORDS:
            raise ValueError(f"Only support 1 <= M <= 10.")
        else:
            self.locAP = np.array(OPTIMAL_COORDS[self.M]) * self.mapsize

        # Generate WDs
        for i in range(self.N):
            while 1:
                xi = np.random.randint(0, self.mapsize*10) / 10.
                yi = np.random.randint(0, self.mapsize*10) / 10.
                isExist = False
                for [x, y] in self.locWD:
                    if x == xi and y == yi:
                        isExist = True
                        break
                for [x, y] in self.locAP:
                    if x == xi and y == yi:
                        isExist = True
                        break
                if not isExist:
                    self.locWD[i] = [xi, yi]
                    break

    def loadMap(self, filename: str):
        npzfile = np.load(filename)
        self.locAP = npzfile['arr_0']
        self.locWD = npzfile['arr_1']

    def saveMap(self, filename: str):
        np.savez(filename, self.locAP, self.locWD)
        savemat(filename[:-3]+'mat', {
            'locAP': self.locAP,
            'locWD': self.locWD
            })

    def calChannelGain(self):
        xn = self.locWD[:, 0].reshape(-1, 1)
        yn = self.locWD[:, 1].reshape(-1, 1)
        
        xm = self.locAP[:, 0].reshape(1, -1)
        ym = self.locAP[:, 1].reshape(1, -1)
        
        dist = np.sqrt( (xn - xm)**2. + (yn - ym)**2. )
        path_loss = dist**(-2)
        
        h_bar_square_U = np.random.exponential(1.0, size=(self.N, self.M))
        h_bar_square_D = np.random.exponential(1.0, size=(self.N, self.M))
        
        self.h_U = self.theta_U * path_loss * h_bar_square_U
        self.h_D = self.theta_D * path_loss * h_bar_square_D
