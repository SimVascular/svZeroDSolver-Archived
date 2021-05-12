import oneD_to_zeroD_convertor
import run_0d_solver

sys.path.append("/home/jonathanpham/Documents/Marsden_Lab/3DModelConversion")
import use_mean_bcs # todo: need to import this code to svZeroDSolver and update the code to make sure that it works with the new 0d input file format (e.g. inclusion of the SOLVEROPTIONS_0D) and updated oneD_to_zeroD_convertor.py
sys.path.pop()

def create_zero_d_input_file(one_d_solver_input_file_path):
    zero_d_solver_input_file_path = oneD_to_zeroD_convertor.convert_1d_to_0d(one_d_solver_input_file_path, destination_folder = ".", default_element_type = "R")
    return zero_d_solver_input_file_path

def run_zero_d_simulation_with_zero_ICs(zero_d_solver_input_file_path, model_number):
    save_y_ydot_to_npy = True
    y_ydot_file_path = "./" + str(model_number) + "_ics.npy"

    run_0d_solver.set_up_and_run_0d_simulation(zero_d_solver_input_file_path, draw_directed_graph, last_cycle, save_results_all, save_results_branch, save_y_ydot_to_npy = save_y_ydot_to_npy, y_ydot_file_path = y_ydot_file_path)

def run_zero_d_simulation_with_mean_ICs(zero_d_solver_input_file_path, model_number):
    use_ICs_from_npy_file = True
    ICs_npy_file_path = "./" + str(model_number) + "_ics.npy"

    run_0d_solver.set_up_and_run_0d_simulation(zero_d_solver_input_file_path, draw_directed_graph, last_cycle, save_results_all, save_results_branch, use_ICs_from_npy_file = use_ICs_from_npy_file, ICs_npy_file_path = ICs_npy_file_path)

if __name__ == "__main__":

    zero_d_solver_input_file_path = "/path/to/0d/input/file"

    model_number = int(os.path.basename(one_d_solver_input_file_path).split("_")[0])

    rather than creating a new 0d input file and saving it to a file, maybe I should just update the 0d parameters to use the mean? but the problem with this is that, it would be harder to check if the 0d model with mean BCs is correct or not... (like visually check, since i wouldnt have any input file to visually check, but i would have to write actual code to check "parameters" (from oneD_to_zeroD_convertor.extract_info_from_solver_input_file), but actually, thats okay. It will be better if I can prevent the code from creating so many unnecessary auxiliary files that i wont ever use) but note that if I decide to go down this path of not creating a new 0d input file, i will have to modify the run_0d_solver code instead, because i will have to manually update "parameters" in which case, this file is irrelevant

    last here - delete this file bc it is irrelevant as mentioned in the above comment

    zero_d_solver_input_file_path = use_mean_bcs.use_mean_values_for_bcs(zero_d_solver_input_file_path)

    run_zero_d_simulation_with_zero_ICs(zero_d_solver_input_file_path, model_number)

    zero_d_solver_input_file_path = create_zero_d_input_file(one_d_solver_input_file_path)

    run_zero_d_simulation_with_mean_ICs(zero_d_solver_input_file_path, model_number)
