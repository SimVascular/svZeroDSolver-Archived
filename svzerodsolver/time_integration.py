# coding=utf-8

# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import scipy
import scipy.sparse.linalg
from scipy.sparse import csr_matrix
import copy

try:
    import matplotlib.pyplot as plt
except:
    pass

class GenAlpha:
    """
    Solves system E*ydot + F*y + C = 0 with generalized alpha and Newton-Raphson for non-linear residual
    """
    def __init__(self, rho, y, sparse=False):
        # Constants for generalized alpha
        self.alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho)
        self.alpha_f = 1.0 / (1.0 + rho)
        self.gamma = 0.5 + self.alpha_m - self.alpha_f

        # problem dimension
        self.n = y.shape[0]

        # stores matrices E, F, vector C, and tangent matrices dE, dF, dC
        self.mat = {}

        # jacobian matrix
        self.M = np.zeros((self.n, self.n))
        self.sparse = sparse
        self.is_sparse = False

        # residual vector
        self.res = np.zeros(self.n)

        self.mats = ['E', 'F', 'dE', 'dF', 'dC']
        self.vecs = ['C']
        for m in self.mats:
            self.mat[m] = np.zeros((self.n, self.n))
        for v in self.vecs:
            self.mat[v] = np.zeros(self.n)

    def setup_sparse(self, block_list):
        """Setup sparsity pattern."""

        for m in self.mats:
            self.mat[m] *= 0.0
        for v in self.vecs:
            self.mat[v] *= 0.0
        for bl in block_list:
            for n in bl.mat.keys():
                if (self.mat[n].ndim == 1):
                    self.mat[n][bl.global_row_id] = 1.0
                else:
                    for i in range(len(bl.global_row_id)):
                        self.mat[n][bl.global_row_id[i], bl.global_col_id] = 1.0

        self.M = csr_matrix(self.M) * 0.0
        for m in self.mats:
            self.mat[m] = csr_matrix(self.mat[m]) * 0.0

    def assemble_structures(self, block_list):
        """
        Assemble block matrices into global matrices
        """
        if self.sparse and not self.is_sparse:
            self.setup_sparse(block_list)
            self.is_sparse = True
        for bl in block_list:
            for n, emat in bl.mat.items():
                # vectors
                if (self.mat[n].ndim == 1):
                    self.mat[n][bl.global_row_id] = emat
                # matrices
                else:
                    for i in range(len(bl.global_row_id)):
                        self.mat[n][bl.global_row_id[i], bl.global_col_id] = emat[i]

    def form_matrix_NR(self, dt):
        """
        Create Jacobian matrix
        """
        self.M = (self.mat['F'] + (self.mat['dE'] + self.mat['dF'] + self.mat['dC'] + self.mat['E'] * self.alpha_m / (
                    self.alpha_f * self.gamma * dt)))

    def form_rhs_NR(self, y, ydot):
        """
        Create residual vector
        """
        self.res = - np.dot(self.mat['E'], ydot) - np.dot(self.mat['F'], y) - self.mat['C']

    def form_matrix_NR_numerical(self, res_i, ydotam, args, block_list, epsilon):
        """
        Numerically compute the Jacobian by computing the partial derivatives of the residual using forward finite differences
        """
        # save original values for restoration later
        yaf_original = copy.deepcopy(args['Solution']) # yaf_i

        # compute numerical Jacobian
        J_numerical = np.zeros((self.n, self.n))
        for jj in range(self.n):

            yaf_step_size = np.zeros(self.n)
            yaf_step_size[jj] = np.abs(yaf_original[jj])  * epsilon

            # get solution at the i+1 step
            args['Solution'] = yaf_original  + yaf_step_size # yaf_ip1

            for b in block_list:
                b.update_solution(args)
            self.assemble_structures(block_list)
            self.form_rhs_NR(args['Solution'], ydotam)

            # use forward finite differences (multiply by -1 b/c form_rhs_NR creates the negative residual)
            J_numerical[:, jj] = (self.res - res_i) / yaf_step_size[jj] * -1

        # restore original quantities
        args['Solution'] = yaf_original

        for b in block_list:
            b.update_solution(args)
        self.assemble_structures(block_list)
        self.form_rhs_NR(args['Solution'], ydotam)

        return J_numerical

    def check_jacobian(self, res_i, ydotam, args, block_list):
        """
        Check if the analytical Jacobian (computed from form_matrix_NR) matches the numerical Jacobian
        """

        epsilon_list = np.power(10, np.linspace(-6, 4, 25))

        fig, axs = plt.subplots(self.n, self.n, figsize = (20, 20))
        for epsilon in epsilon_list:
            J_numerical = self.form_matrix_NR_numerical(res_i, ydotam, args, block_list, epsilon)
            error = np.abs(self.M - J_numerical)
            for ii in range(self.n):
                for jj in range(self.n):
                    axs[ii, jj].loglog(epsilon, error[ii, jj], 'k*-')

        for ax in axs.flat:
            ax.set(xlabel='epsilon', ylabel='error')

        fig.suptitle('absolute error vs epsilon')
        plt.show()

    def step(self, y, ydot, t, block_list, args, dt, nit=30):
        """
        Perform one time step
        """
        # initial guess for time step
        curr_y = y.copy() + 0.5 * dt * ydot
        curr_ydot = ydot.copy() * ((self.gamma - 0.5) / self.gamma)

        # Substep level quantities
        yaf = y + self.alpha_f * (curr_y - y)
        ydotam = ydot + self.alpha_m * (curr_ydot - ydot)

        # initialize solution
        args['Time'] = t + self.alpha_f * dt
        args['Solution'] = yaf

        # initialize blocks
        for b in block_list:
            b.update_constant()
            b.update_time(args)

        iit = 0
        while (np.abs(self.res).max() > 5e-4 or iit == 0) and iit < nit:
            # update solution-dependent blocks
            for b in block_list:
                b.update_solution(args)

            # update residual and jacobian
            self.assemble_structures(block_list)
            self.form_rhs_NR(yaf, ydotam)
            self.form_matrix_NR(dt)

            # perform finite-difference check of jacobian if requested
            if args['check_jacobian']:
                if args['Time'] > dt:
                    self.check_jacobian(copy.deepcopy(self.res), ydotam, args, block_list)

            # solve for Newton increment
            if self.sparse:
                dy = scipy.sparse.linalg.spsolve(self.M, self.res)
            else:
                dy = scipy.linalg.solve(self.M, self.res)

            # update solution
            yaf += dy
            ydotam += self.alpha_m * dy / (self.alpha_f * self.gamma * dt)

            if np.any(np.isnan(self.res)):
                raise RuntimeError('Solution nan')

            args['Solution'] = yaf
            iit += 1

        if iit >= nit:
            print("Max NR iterations reached at time: ", t, " , max error: ", max(abs(self.res)))

        # update time step
        curr_y = y + (yaf - y) / self.alpha_f
        curr_ydot = ydot + (ydotam - ydot) / self.alpha_m

        args['Time'] = t + dt

        return curr_y, curr_ydot
