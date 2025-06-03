# api/simulation/run_model.py
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import traceback # For detailed error logging

from .params import get_model_constants
from .initial_values import get_initial_values
from .simulation_setup import get_simulation_times, base_feed_intake_pattern, scale_feed_intake_pattern
from .ode_function import ode_function_py

DMI_OF_BASE_PATTERN = 22.50015 # Sum of (FI_rate * duration) for the base_feed_intake_pattern

def execute_single_simulation(sim_params, diet_params):
    """
    Executes a single run of the rumen simulation model.

    Args:
        sim_params (dict): Parameters for the simulation run itself (e.g., hours, points_per_hour).
        diet_params (dict): Diet-specific parameters (e.g., DMI, NDF, St, WSC in g/kg DM).

    Returns:
        dict: A dictionary containing simulation results or an error message.
              Format: {"success": True/False, "results": [...], "error": "message"}
    """
    try:
        # 1. Prepare Model Constants based on diet
        #    `diet_params` from frontend should contain NDF, St, etc. in g/kg DM
        #    `DMI` from `diet_params` is kg/day for scaling feed intake.
        model_nutrient_overrides = {
            key: diet_params[key] for key in diet_params if key != "DMI" # Exclude DMI for now
        }
        # The 'DMI' key in _base_model_constants is a default; it might or might not be used
        # if the model's equations directly use a DMI parameter beyond feed intake rate.
        # For now, we assume nutrient composition overrides are the main dynamic part.
        current_model_constants = get_model_constants(model_nutrient_overrides)

        # 2. Prepare Initial Values
        #    Initial values might depend on some base parameters, so get_initial_values can take them.
        initial_values_data = get_initial_values() # Using default constants for initial values
        y0 = np.array(initial_values_data["initialValuesVec"])
        state_variable_names = initial_values_data["stateVariableNames"]

        # 3. Prepare Simulation Time Points
        hours = sim_params.get("hours", 24)
        points_per_hour = sim_params.get("points_per_hour", 1000)
        time_data = get_simulation_times(hours, points_per_hour)
        times_eval = time_data["times"]
        
        if len(times_eval) < 2: # Should be handled by get_simulation_times but defensive check
            return {"success": False, "error": "Invalid simulation time setup, too few points."}
        t_span = (times_eval[0], times_eval[-1])

        # 4. Prepare Feed Intake Data (Scaled)
        target_dmi_for_scaling = diet_params.get("DMI") # This is kg/day from frontend
        if target_dmi_for_scaling is None or target_dmi_for_scaling <= 0:
            return {"success": False, "error": f"Invalid DMI provided for scaling: {target_dmi_for_scaling}"}
        
        scaled_feed_intake = scale_feed_intake_pattern(
            base_feed_intake_pattern,
            target_dmi_for_scaling,
            DMI_OF_BASE_PATTERN
        )
        
        # delt_ext (stpsz_value) - not directly used by solve_ivp if ode_function_py doesn't need it
        stpsz_value = times_eval[1] - times_eval[0] if len(times_eval) > 1 else 0


        # 5. Prepare Arguments for ODE Function
        #    The order must match ode_function_py: (t, state_vec, parms_dict, feedintake_data_ext, state_variable_names_arg, delt_ext)
        ode_args = (current_model_constants, scaled_feed_intake, state_variable_names, stpsz_value)

        # 6. Run Simulation
        print(f"Running simulation with DMI: {target_dmi_for_scaling}, NDF: {current_model_constants.get('NDF')}, St: {current_model_constants.get('St')}")
        solution = solve_ivp(
            fun=ode_function_py,
            t_span=t_span,
            y0=y0,
            t_eval=times_eval,
            args=ode_args,
            method='RK45', # RK45 LSODA Good for potentially stiff problems
            rtol=1e-6,
            atol=1e-8
        )

        if not solution.success:
            print(f"ODE solver failed: {solution.message}")
            return {"success": False, "error": f"ODE solver failed: {solution.message}"}

        # 7. Process and Return Results
        results_array = solution.y.T
        result_df = pd.DataFrame(results_array, columns=state_variable_names)
        result_df.insert(0, 'time', solution.t)
        
        # Convert DataFrame to list of dictionaries for JSON response
        simulation_output = result_df.to_dict(orient='records')

        return {"success": True, "results": simulation_output}

    except Exception as e:
        print(f"Error during simulation: {str(e)}")
        print(traceback.format_exc()) # Log full traceback to server console
        return {"success": False, "error": f"An unexpected error occurred in simulation: {str(e)}"}

if __name__ == '__main__':
    print("--- Testing execute_single_simulation ---")
    test_sim_params = {"hours": 24, "points_per_hour": 100} # Lower res for quick test
    test_diet_params = {
        "DMI": 20.0, # kg/day
        "NDF": 300.0, # g/kg DM
        "St": 250.0,  # g/kg DM
        "WSC": 100.0, # g/kg DM
        "Acin": 10.0, # g/kg DM
        "Prin": 1.5,  # g/kg DM
        "Buin": 1.5,  # g/kg DM
        "Lain": 15.0  # g/kg DM
    }
    result = execute_single_simulation(test_sim_params, test_diet_params)
    if result["success"]:
        print(f"Simulation successful. Number of result points: {len(result['results'])}")
        print("First result point:", result['results'][0])
        print("Last result point:", result['results'][-1])
    else:
        print(f"Simulation failed: {result['error']}")