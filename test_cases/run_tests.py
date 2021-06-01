#!/usr/bin/env python3
# coding=utf-8

import os
import pdb
import subprocess
import numpy as np

# last here - then run this code to see if it works for the below single test case.
# Then add more test cases if the below R_R test works

def get_tests():
    """
    Add new test cases here
    How to define a new test case <NAME>:
      - create new file with name <NAME>.in with the following naming conventions:
        MODEL results_<NAME>_
        ELEMENT branch0_seg0 ...
        ELEMENT branch1_seg0 ...
      - add test case to the dictionary below
      - all dictionary values are lists, so add as many results as you like
      - see the documentation of the Test class below on how to set up a test case
    """
    tests = {}

    # field, branch, branch_node, time_step, res, tol, fun -- last here: delete this line when I finish making all test cases
    tests['steadyFlow_R_R'] = [ Test('pressure', 0,  0, -1, 1000.0, 1.0e-10, 'point'),
                                Test('pressure', 0, -1, -1,  500.0, 1.0e-10, 'point'),
                                Test('flow',     0,  0, -1,    5.0, 1.0e-16, 'point'),
                                Test('flow',     0, -1, -1,    5.0, 1.0e-16, 'point')   ]

    return tests


class Test:
    """
    GOOD
    Class to define (and check) test cases
    """

    def __init__(self, field, branch, branch_node, time_step, res, tol, fun):
        """
        GOOD
        Args:
            field: field to check ('flow', 'pressure')
            branch: branch to check
            branch_node: branch node to check (usually 0 or -1)
            time_step: individual time step to check (usually -1) or list of time steps (e.g. np.arange(5,10))
            res: result to check (float)
            tol: relative tolerance for result check
            fun: type of result comparison to perform in time ('point', 'mean', 'max', 'min')
                 specify an interval for 'time_step' when chosing 'mean', 'max', or 'min'
        """
        # sanity checks
        if field not in ['flow', 'pressure']:
            raise ValueError('Field ' + field + ' unknown. Please select from flow, pressure')
        if fun not in ['point', 'mean', 'max', 'min']:
            raise ValueError('Function ' + fun + ' unknown. Please select from point, mean, max, min')
        if not np.isscalar(time_step) and fun == 'point':
            raise ValueError('Specify a single time point when selecting result type ' + fun)
        if np.isscalar(time_step) and (fun == 'mean' or fun == 'max' or fun == 'min'):
            raise ValueError('Specify a time interval when selecting result type ' + fun)

        self.field = field
        self.branch = branch
        self.branch_node = branch_node
        self.time_step = time_step
        self.res = res
        self.tol = tol
        self.fun = fun

    def check(self, results):
        """
        GOOD
        Perform the actual result check
        """
        # read result from svZeroDSolver
        res = self.read_result(results)

        # calculate relative difference
        diff = np.abs((res - self.res) / self.res)

        # check if difference in results is larger than given tolerance
        if diff > self.tol:
            return self.print_err(res, diff)

        return False

    def read_result(self, results):
        """
        GOOD
        Read results and select function
        """
        # extract result
        res = results[self.field][self.branch][self.branch_node, self.time_step]

        # select result type
        if self.fun == 'point':
            return res
        elif self.fun == 'mean':
            return np.mean(res)
        elif self.fun == 'max':
            return np.max(res)
        elif self.fun == 'min':
            return np.min(res)
        else:
            raise ValueError('Unknown result type ' + self.fun)

    def print_err(self, res, diff):
        """
        GOOD
        Create error string for user
        """
        err = 'Test failed. ' + self.field + ' in branch ' + str(self.branch)
        err += ', branch_node ' + str(self.branch_node) + ', time_step ' + str(self.time_step) + '. expected: ' + str(self.res)
        err += '. got: ' + str(res) + '. abs rel diff: ' + str(diff) + ' > ' + str(self.tol)
        return err


def read_results_0d(zero_d_simulation_results_file_path):
    """
    GOOD
    Read results from svZeroDSolver
    Args:
        str zero_d_simulation_results_file_path: path to 0d simulation results (saved from svZeroDSolver using branching structure)
    Returns:
        dictionary res[result field][branch id][branch_node, time step]
    """
    res = np.load(zero_d_simulation_results_file_path, allow_pickle = True).item()
    os.remove(zero_d_simulation_results_file_path)
    os.remove(zero_d_simulation_results_file_path[:-18] + "all_results.npy")
    return res


def run_check(results, result_checks):
    """
    GOOD
    Check the results of a test
    """
    # loop all results
    for test in result_checks:
        err = test.check(results)
        if err:
            # test failed
            return err

    # all tests passed
    return False


def run_test(build_dir, test_dir, name, check):
    """
    GOOD
    Run a test case and check the results
    """
    # name of svZeroDSolver executable
    exe = os.path.join(build_dir, 'svZeroDSolver.py')

    # name of input file
    inp = os.path.join(test_dir, name + '.in')

    # run simulation
    try:
        subprocess.check_output([exe, inp])
    except subprocess.CalledProcessError as err:
        return 'Test failed. svZeroDSolver returned error:\n' + err.output.decode("utf-8")

    # extract results
    try:
        zero_d_simulation_results_file_path = inp[:-3] + "_branch_results.npy"
        res = read_results_0d(zero_d_simulation_results_file_path)
    except Exception as err:
        return 'Test failed. Result extraction failed:\n' + str(err)

    # compare to stored results
    try:
        res = run_check(res, check)
    except Exception as err:
        return 'Test failed. Result check failed:\n' + str(err)

    return res


def main():
    """
    GOOD
    Loop over all test cases and check if all results match
    """
    # set paths
    if 'BUILD_DIR' not in os.environ and 'TEST_DIR' not in os.environ:
        # run locally
        fpath = os.path.dirname(os.path.realpath(__file__))
        build_dir = os.path.join(fpath, '..')
        test_dir = fpath
    else:
        # run on Travis
        build_dir = os.environ['BUILD_DIR']
        test_dir = os.environ['TEST_DIR']

    # get test cases
    try:
        test_cases = get_tests()
    except Exception as err:
        print(err)
        return True

    # loop all test cases
    for name, check in test_cases.items():
        print('Running test ' + name)
        err = run_test(build_dir, test_dir, name, check)

        # check if errors occured
        if err:
            print(err)
            return True
        else:
            print('Test passed')

    # no tests failed
    return False


if __name__ == '__main__':
    # tests fail
    if main():
        exit(1)

    # tests passs
    else:
        exit(0)
