from svzerodsolver.solver import set_up_and_run_0d_simulation
import os

import shutil
import numpy as np

RTOL_PRES = 1.0e-7
RTOL_FLOW = 1.0e-8


def run_test_case_by_name(name, testdir):
    """Run a test case by its case name.
    
    Args:
        name: Name of the test case.
        testdir: Directory for performing the simulation.
    """
    testfile = os.path.join(os.path.dirname(__file__), "cases", name + ".json")
    shutil.copyfile(testfile, os.path.join(testdir, name + ".json"))
    set_up_and_run_0d_simulation(os.path.join(testdir, name + ".json"))
    result_file = os.path.join(testdir, name + "_branch_results.npy")
    return np.load(result_file, allow_pickle=True).item()

def get_result(result_array, field, branch, branch_node, time_step):
    """"Get results at specific field, branch, branch_node and time step."""
    # extract result
    return result_array[field][branch][branch_node, time_step]


def test_steady_flow_R_R(tmpdir):
    results = run_test_case_by_name("steadyFlow_R_R", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 1100.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 600.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # outlet flow


def test_steady_flow_r_coronary(tmpdir):
    results = run_test_case_by_name("steadyFlow_R_coronary", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 2000.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 1500.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # outlet flow


def test_steady_flow_rlc_r(tmpdir):
    results = run_test_case_by_name("steadyFlow_RLC_R", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 1100.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 600.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # outlet flow


def test_steady_flow_rc_r(tmpdir):
    results = run_test_case_by_name("steadyFlow_RC_R", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 1100.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 600.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # outlet flow


def test_steady_flow_rl_r(tmpdir):
    results = run_test_case_by_name("steadyFlow_RL_R", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 1100.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 600.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # outlet flow


def test_steady_flow_r_rcr(tmpdir):
    results = run_test_case_by_name("steadyFlow_R_RCR", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 10500.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 10000.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # outlet flow


def test_steady_flow_r_steady_pressure(tmpdir):
    results = run_test_case_by_name("steadyFlow_R_steadyPressure", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 1500.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 1000.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # outlet flow


def test_steady_flow_stenosis_r(tmpdir):
    results = run_test_case_by_name("steadyFlow_stenosis_R", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 3600.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 600.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # outlet flow


def test_steady_flow_bifurcationr_r1(tmpdir):
    results = run_test_case_by_name("steadyFlow_bifurcationR_R1", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 1100.0, rtol=RTOL_PRES
    )  # parent inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 600.0, rtol=RTOL_PRES
    )  # parent outlet pressure
    assert np.isclose(
        get_result(results, "pressure", 1, 0, -1), 600.0, rtol=RTOL_PRES
    )  # daughter1 inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 1, -1, -1), 350.0, rtol=RTOL_PRES
    )  # daughter1 outlet pressure
    assert np.isclose(
        get_result(results, "pressure", 2, 0, -1), 600.0, rtol=RTOL_PRES
    )  # daughter2 inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 2, -1, -1), 350.0, rtol=RTOL_PRES
    )  # daughter2 outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # parent inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # parent outlet flow
    assert np.isclose(
        get_result(results, "flow", 1, 0, -1), 2.5, rtol=RTOL_FLOW
    )  # daughter1 inlet flow
    assert np.isclose(
        get_result(results, "flow", 1, -1, -1), 2.5, rtol=RTOL_FLOW
    )  # daughter1 outlet flow
    assert np.isclose(
        get_result(results, "flow", 2, 0, -1), 2.5, rtol=RTOL_FLOW
    )  # daughter2 inlet flow
    assert np.isclose(
        get_result(results, "flow", 2, -1, -1), 2.5, rtol=RTOL_FLOW
    )  # daughter2 outlet flow


def test_steady_flow_bifurcationr_r2(tmpdir):
    results = run_test_case_by_name("steadyFlow_bifurcationR_R2", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 3462.5, rtol=RTOL_PRES
    )  # parent inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 1962.5, rtol=RTOL_PRES
    )  # parent outlet pressure
    assert np.isclose(
        get_result(results, "pressure", 1, 0, -1), 1962.5, rtol=RTOL_PRES
    )  # daughter1 inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 1, -1, -1), 432.5, rtol=RTOL_PRES
    )  # daughter1 outlet pressure
    assert np.isclose(
        get_result(results, "pressure", 2, 0, -1), 1962.5, rtol=RTOL_PRES
    )  # daughter2 inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 2, -1, -1), 1375.0, rtol=RTOL_PRES
    )  # daughter2 outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # parent inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # parent outlet flow
    assert np.isclose(
        get_result(results, "flow", 1, 0, -1), 3.825, rtol=RTOL_FLOW
    )  # daughter1 inlet flow
    assert np.isclose(
        get_result(results, "flow", 1, -1, -1), 3.825, rtol=RTOL_FLOW
    )  # daughter1 outlet flow
    assert np.isclose(
        get_result(results, "flow", 2, 0, -1), 1.175, rtol=RTOL_FLOW
    )  # daughter2 inlet flow
    assert np.isclose(
        get_result(results, "flow", 2, -1, -1), 1.175, rtol=RTOL_FLOW
    )  # daughter2 outlet flow


def test_pulsatile_flow_r_rcr(tmpdir):
    results = run_test_case_by_name("pulsatileFlow_R_RCR", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, 0), 4620.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, 0), 4400.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, 0), 2.2, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, 0), 2.2, rtol=RTOL_FLOW
    )  # outlet flow


def test_pulsatile_flow_r_coronary(tmpdir):
    results = run_test_case_by_name("pulsatileFlow_R_coronary", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, 0), 880.0, rtol=RTOL_PRES
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, 0), 660.0, rtol=RTOL_PRES
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, 0), 2.2, rtol=RTOL_FLOW
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, 0), 2.2, rtol=RTOL_FLOW
    )  # outlet flow


def test_pusatile_flow_cstenosis_steady_pressure(tmpdir):
    results = run_test_case_by_name("pusatileFlow_CStenosis_steadyPressure", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -439),
        0.5937169800360568,
        rtol=1.0e-5,
    )  # inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -439), 0.1, rtol=1.0e-5
    )  # outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -439),
        0.7026499697830042,
        rtol=1.0e-5,
    )  # inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -439),
        0.7026499697830042,
        rtol=1.0e-5,
    )  # outlet flow


def test_steady_flow_confluencer_r(tmpdir):
    results = run_test_case_by_name("steadyFlow_confluenceR_R", tmpdir)
    assert np.isclose(
        get_result(results, "pressure", 0, 0, -1), 6600.0, rtol=RTOL_PRES
    )  # parent inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 0, -1, -1), 6100.0, rtol=RTOL_PRES
    )  # parent outlet pressure
    assert np.isclose(
        get_result(results, "pressure", 1, 0, -1), 8100.0, rtol=RTOL_PRES
    )  # daughter1 inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 1, -1, -1), 6100.0, rtol=RTOL_PRES
    )  # daughter1 outlet pressure
    assert np.isclose(
        get_result(results, "pressure", 2, 0, -1), 6100.0, rtol=RTOL_PRES
    )  # daughter2 inlet pressure
    assert np.isclose(
        get_result(results, "pressure", 2, -1, -1), 1600.0, rtol=RTOL_PRES
    )  # daughter2 outlet pressure
    assert np.isclose(
        get_result(results, "flow", 0, 0, -1), 5.0, rtol=RTOL_FLOW
    )  # parent inlet flow
    assert np.isclose(
        get_result(results, "flow", 0, -1, -1), 5.0, rtol=RTOL_FLOW
    )  # parent outlet flow
    assert np.isclose(
        get_result(results, "flow", 1, 0, -1), 10.0, rtol=RTOL_FLOW
    )  # daughter1 inlet flow
    assert np.isclose(
        get_result(results, "flow", 1, -1, -1), 10.0, rtol=RTOL_FLOW
    )  # daughter1 outlet flow
    assert np.isclose(
        get_result(results, "flow", 2, 0, -1), 15.0, rtol=RTOL_FLOW
    )  # daughter2 inlet flow
    assert np.isclose(
        get_result(results, "flow", 2, -1, -1), 15.0, rtol=RTOL_FLOW
    )  # daughter2 outlet flow
