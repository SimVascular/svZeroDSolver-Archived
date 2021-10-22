## About svZeroDSolver

svZeroDSolver is a Python code that simulates the hemodynamics in zero-dimensional (0D) lumped parameter models of 
vascular networks. These 0D models are governed by differential algebraic equations (DAEs).

The solver uses a highly modular framework to model the vascular anatomy, using individual 0D elements to represent different parts of 
the vascular anatomy (and boundary conditions). The individual 0D elements and their associated governing equations defined in `blocks.py`. 
In `solver.py`, the blocks are assembled and simulated using the generalized-alpha time-stepping method defined in `time_integration.py`.

The svZeroDSolver Python files are in the `svzerodsolver` Python package directory. 

<!-- add link to the 0D solver and theory documentation on SimVascular website when it is available -->

svZeroDSolver currently supports the following vascular 0D modeling options and boundary conditions:

#### Vascular 0D elements:
- Resistor
- Resistor-capacitor
- Resistor-inductor
- Resistor-capacitor-inductor

#### Boundary conditions:
- Pressure
- Resistor
- RCR
- Coronary
- Flow


### Prerequisites

The following software is required:

- Python 3

The required Python packages are listet in `requirements.txt`

### Execution

The solver can be executed four ways.

1) Execute from the top level Git repository using

~~~bash
python -m svzerodsolver.solver SOLVER_INPUT_FILE 
~~~

2) Setting the `PYTHONPATH` environent variable

The solver can be executed from any directory by setting the `PYTHONPATH` environent variable 
to the top level Git repository 

~~~bash
export PYTHONPATH=$PYTHONPATH:/$HOME/svZeroDSolver/
~~~

3) Installing the `svzerodsolver` Python package 

The `svzerodsolver` Python package is installed using the `setup.py` script

~~~bash
python setup.py install
~~~

4) Importing the `svzerodsolver` Python package

~~~python
import svzerodsolver
svzerodsolver.solver.set_up_and_run_0d_simulation('input.json')
~~~

This variant enables running svZeroDSolver within a user-defined Python code (e.g. parameter optimization, uncertainty quantification)

### Solver options 

The options supported by the solver are listed using

~~~
python -m svzerodsolver.solver --help
~~~



