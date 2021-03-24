# coding=utf-8
import numpy as np
import scipy
from scipy.sparse import csr_matrix
import pdb

def min_ydot_least_sq_init(neq,eps_min,yinit,block_list,args,dt,rho,eps_factor=5.0):
    # System : min (over y) ||Fy+C||^2 + eps||Ay-yinit||^2
    # Inversion equation:
    # y <-- inv(F'F+eps*D) (-F'C+eps*D*yinit)
    # yinit : Desired set of initial conditions
    # D : diag([(i in yinit) for i in range(neq) ] )
    # If yinit == zeros, set D = I

    # ydot is solved in an approximate sense: ydot <-- np.linalg.lstsq(E,-Fy-C)

    # Solve as a sequence of problems : eps ---> eps_min
    eps = 10.0
    iit = 0
    args['Time'] = 0
    # y0 = np.zeros(neq)

    y0 = yinit

    E,F,C,dE,dF,dC = initialize_solution_matrices(neq)

    if np.linalg.norm(yinit) == 0.:
        D = np.eye(neq)
    else :
        D = np.diag([_!=0 for _ in yinit])

    print("Approximate consistent initialization : \n\n")

    while eps > eps_min :
        iit +=1
        args['Solution'] = y0
        assemble_structures(E,F,C,dE,dF,dC,args,block_list)

        M = np.dot(F.transpose(),F)+eps*D
        v = -np.dot(F.transpose(),C) + eps*np.dot(D,yinit)
        y0,_,_,_ = np.linalg.lstsq(M,v)

        ydot0,_,_,_ = np.linalg.lstsq(E,-np.dot(F,y0)-C)

        print("Iteration ",iit,", Initializing residual: ",np.linalg.norm(form_rhs_NR(E,F,C,y0,ydot0)))
        eps = eps/eps_factor


    return y0,ydot0

def min_ydot_cons_least_sq_init(neq,eps_min,yinit,block_list,args,dt,rho,eps_factor=5.0):
    # System : min (over y) ||Fy+C||^2 + sum_j(lambda_j(y_j-yinit_j))
    # Inversion equation:
    # [2(F'F)  A' ][  y0  ] = [ -2F'C ]
    # [  A     0  ][lambda] = [ yinit ]

    # ydot is solved in an approximate sense: ydot <-- np.linalg.lstsq(E,-Fy-C)

    # Solve as a sequence of problems : eps ---> eps_min
    eps = 1.0
    iit = 0
    args['Time'] = 0

    y0 = yinit
    ydot0 = np.zeros(neq)

    indx = np.nonzero(yinit)

    E,F,C,dE,dF,dC = initialize_solution_matrices(neq)

    print("Approximate consistent initialization : \n\n")

    isize = indx[0].size

    if isize == 0:
        while eps > eps_min :
            iit +=1
            args['Solution'] = y0
            assemble_structures(E,F,C,dE,dF,dC,args,block_list)
            M =  np.dot(F.transpose(),F)
            v = -np.dot(F.transpose(),C+np.dot(E,ydot0))
            y0,_,_,_ = np.linalg.lstsq(M,v)
            ydot0,_,_,_ = np.linalg.lstsq(E,-np.dot(F,y0)-C)
            print("Iteration ",iit,", Initializing residual: ",np.linalg.norm(form_rhs_NR(E,F,C,y0,ydot0)))
            print("Iteration ",iit,", Initializing time derivative size: ",np.linalg.norm(ydot0))
            eps = eps/eps_factor

    else :
        A = np.zeros((isize,neq))
        i = 0
        print(indx[0])
        for j in indx[0]:
            A[i,j] = 1
            i+=1

        while eps > eps_min :
            iit +=1
            args['Solution'] = y0
            assemble_structures(E,F,C,dE,dF,dC,args,block_list)

            T = 2*np.dot(F.transpose(),F)

            M = np.block([
                [T,A.transpose()],
                [A,np.zeros((isize,isize))]
                ])

            v = np.block([-2*np.dot(F.transpose(),C+np.dot(E,ydot0)),yinit[indx]]).transpose()

            yM,_,_,_ = np.linalg.lstsq(M,v)

            y0 = yM[:neq]
            ydot0,_,_,_ = np.linalg.lstsq(E,-np.dot(F,y0)-C)

            print("Iteration ",iit,", Initializing residual: ",np.linalg.norm(form_rhs_NR(E,F,C,y0,ydot0)))
            print("Iteration ",iit,", Initializing time derivative size: ",np.linalg.norm(ydot0))
            eps = eps/eps_factor


    return y0,ydot0


# Equation: E*ydot + F*y + C = 0
class GenAlpha():
    def __init__(self, rho, y):
        # Constants for generalized alpha
        self.alpha_m = 0.5*(3.0-rho)/(1.0+rho)
        self.alpha_f = 1.0/(1.0+rho)
        self.gamma = 0.5 + self.alpha_m - self.alpha_f
        self.n = y.shape[0]

        self.damping_step = 1.5

        self.mat = {}
        self.M = []
        self.res = []
        self.res0 = []

        self.initialize_solution_matrices()

    def initialize_solution_matrices(self):
        mats = ['E', 'F', 'dE', 'dF', 'dC']
        vecs = ['C']

        # create empty dense matrices
        for m in mats:
            self.mat[m] = np.zeros((self.n, self.n))
        for v in vecs:
            self.mat[v] = np.zeros(self.n)

    def assemble_structures(self, block_list):
        # assemble local into global matrices
        trg = ['E', 'F', 'C', 'dE', 'dF', 'dC']

        for bl in block_list:
            src = [bl.emxcoe, bl.fmxcoe, bl.cveccoe, bl.demxcoe, bl.dfmxcoe, bl.dcmxcoe]
            for S, T in zip(src, trg):
                # vectors
                if len(self.mat[T].shape) == 1:
                    for i in range(len(S)):
                        self.mat[T][bl.global_row_id[i]] = S[i]
                # matrices
                else:
                    for i in range(len(S)):
                        for j in range(len(S[i])):
                            self.mat[T][bl.global_row_id[i], bl.global_col_id[j]] = S[i][j]

    def form_matrix_NR(self, dt):
        self.M = (self.mat['F'] + (self.mat['dE'] + self.mat['dF'] + self.mat['dC'] + self.mat['E'] * self.alpha_m / (self.alpha_f * self.gamma * dt)))

    def form_rhs_NR(self, y, ydot):
        self.res = - np.dot(self.mat['E'], ydot) - np.dot(self.mat['F'], y) - self.mat['C']
        # return - csr_matrix(E).dot(ydot) - csr_matrix(F).dot(y) - C

    def step(self, y, ydot, t, block_list, args, dt, nit=30):
        # Initial guess for n+1-th step -- explicit euler type guess, half step
        curr_y = y+0.5*dt*ydot
        curr_ydot = np.copy(ydot) * ((self.gamma - 0.5) / self.gamma)

        # Substep level quantities
        yaf = y + self.alpha_f * (curr_y - y)
        ydotam = ydot + self.alpha_m * (curr_ydot - ydot)

        # print t
        args['Time'] = t + self.alpha_f * dt

        iit = 0

        args['Solution'] = yaf

        # initialize blocks
        for b in block_list:
            b.update_constant()
            b.update_time(args)
            b.update_solution(args)

        self.assemble_structures(block_list)

        self.form_rhs_NR(yaf, ydotam)
        self.res0 = self.res

        # print "time = ", t, " , Max residual (outside while loop) = ", max(abs(res0))
        while max(abs(self.res0)) > 5e-4 and iit < nit:

            damping = 1.

            if iit > 0:
                # update solution-dependent blocks
                for b in block_list:
                    b.update_solution(args)

                # update newton
                self.assemble_structures(block_list)
                self.form_rhs_NR(yaf, ydotam)

            self.form_matrix_NR(dt)
            dy = scipy.sparse.linalg.spsolve(csr_matrix(self.M), self.res)
            while np.linalg.norm(self.res) >= np.linalg.norm(self.res0) and damping > 1e-5:
                yaf2 = yaf + damping * dy
                ydotam2 = ydotam + damping * self.alpha_m * dy / (self.alpha_f * self.gamma * dt)
                damping /= self.damping_step
                self.form_rhs_NR(yaf2, ydotam2)

            yaf = yaf + self.damping_step * damping * dy
            ydotam = ydotam + self.damping_step * damping * self.alpha_m * dy / (self.alpha_f * self.gamma * dt)

            self.res0 = self.res
            if np.any(np.isnan(self.res0)):
                raise RuntimeError('Solution nan')
            # print "time = ", t, " , Max residual (in while loop) = ", max(abs(res0))

            # Check this equation up
            # ydotam = (1-alpha_m/gamma)*ydot + (alpha_m/(gamma*dt*alpha_f))*(yaf-y)

            args['Solution'] = yaf
            iit += 1

        if iit >= nit:
            print("Max NR iterations reached at time: ", t, " , max error: ", max(abs(res0))) # NOTE: "max error" = max residual here
            # print "Condition number of F ", np.linalg.cond(F)
            # print "Condition number of NR matrix: ",np.linalg.cond(M)
            # print M

        curr_y = y + (yaf - y) / self.alpha_f
        curr_ydot = ydot + (ydotam - ydot) / self.alpha_m

        args['Time'] = t+dt

        return curr_y, curr_ydot
