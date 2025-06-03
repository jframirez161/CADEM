# api/simulation/simulation_setup.py
import numpy as np
import copy # For deepcopy

# --- Base Feed Intake Pattern (as per JavaScript version) ---
base_feed_intake_pattern = {
    "feedintake": [
        {"T_begin": 0.0, "T_eind": 0.5, "FI_rate": 6.1899},
        {"T_begin": 0.5, "T_eind": 1.0, "FI_rate": 3.5438},
        {"T_begin": 1.0, "T_eind": 1.5, "FI_rate": 2.8811},
        {"T_begin": 1.5, "T_eind": 2.0, "FI_rate": 1.8753},
        {"T_begin": 2.0, "T_eind": 2.5, "FI_rate": 1.0058},
        {"T_begin": 2.5, "T_eind": 3.0, "FI_rate": 0.7567},
        {"T_begin": 3.0, "T_eind": 4.0, "FI_rate": 0.3384},
        {"T_begin": 4.0, "T_eind": 5.0, "FI_rate": 0.2538},
        {"T_begin": 5.0, "T_eind": 6.0, "FI_rate": 0.1692},
        {"T_begin": 6.0, "T_eind": 7.0, "FI_rate": 0.0611},
        {"T_begin": 7.0, "T_eind": 12.0, "FI_rate": 0.0},
        {"T_begin": 12.0, "T_eind": 12.5, "FI_rate": 6.1899},
        {"T_begin": 12.5, "T_eind": 13.0, "FI_rate": 3.5438},
        {"T_begin": 13.0, "T_eind": 13.5, "FI_rate": 2.8811},
        {"T_begin": 13.5, "T_eind": 14.0, "FI_rate": 1.8753},
        {"T_begin": 14.0, "T_eind": 14.5, "FI_rate": 1.0058},
        {"T_begin": 14.5, "T_eind": 15.0, "FI_rate": 0.7567},
        {"T_begin": 15.0, "T_eind": 16.0, "FI_rate": 0.3384},
        {"T_begin": 16.0, "T_eind": 17.0, "FI_rate": 0.2538},
        {"T_begin": 17.0, "T_eind": 18.0, "FI_rate": 0.1692},
        {"T_begin": 18.0, "T_eind": 19.0, "FI_rate": 0.0611},
        {"T_begin": 19.0, "T_eind": 24.0, "FI_rate": 0.0}
    ]
}

def get_simulation_times(hours=24, num_points_per_hour=1000):
    """
    Generates an array of time points for the simulation.
    Args:
        hours (int): Total simulation time in hours.
        num_points_per_hour (int): Number of output points per hour.
    Returns:
        dict: Containing {"times": numpy_array_of_time_points}.
    """
    # Ensure num_points_per_hour results in at least 2 points for linspace if hours > 0
    # The total number of points will be hours * num_points_per_hour.
    # linspace's `num` argument is the total number of samples.
    if hours <= 0:
        return {"times": np.array([0.0])}
    if num_points_per_hour <=0: # Should not happen with Pydantic validation, but good check
        num_points_per_hour = 1

    total_num_points = int(hours * num_points_per_hour)
    if total_num_points < 2: # linspace needs at least 2 points if start != stop
        total_num_points = 2
        
    # endpoint=True means the stop value (hours) is included.
    t_eval = np.linspace(start=0, stop=hours, num=total_num_points, endpoint=True)
    return {"times": t_eval}

def scale_feed_intake_pattern(base_pattern, target_daily_dmi, original_pattern_dmi):
    """
    Scales the FI_rate in a feed intake pattern to match a target daily DMI.
    Args:
        base_pattern (dict): The base feed intake pattern structure.
        target_daily_dmi (float): The desired total DMI over 24 hours (kg/day).
        original_pattern_dmi (float): The total DMI represented by the base_pattern over 24 hours (kg/day).
    Returns:
        dict: A new feed intake pattern with scaled FI_rate values.
    """
    if not (isinstance(target_daily_dmi, (int, float)) and target_daily_dmi > 0):
        print(f"Warning: Invalid target_daily_dmi: {target_daily_dmi}. Using base pattern without scaling.")
        return copy.deepcopy(base_pattern)
    if not (isinstance(original_pattern_dmi, (int, float)) and original_pattern_dmi > 0):
        print(f"Warning: Original pattern DMI is zero or invalid ({original_pattern_dmi}). Cannot scale. Using base pattern.")
        return copy.deepcopy(base_pattern)

    scaling_factor = target_daily_dmi / original_pattern_dmi
    
    # Deepcopy to avoid modifying the original base_pattern if it's reused
    scaled_pattern = copy.deepcopy(base_pattern) 
    
    for entry in scaled_pattern["feedintake"]:
        entry["FI_rate"] = entry["FI_rate"] * scaling_factor
            
    return scaled_pattern

if __name__ == '__main__':
    print("--- Testing get_simulation_times ---")
    time_data_default = get_simulation_times()
    print(f"Default times (first 5): {time_data_default['times'][:5]}, last: {time_data_default['times'][-1]}, len: {len(time_data_default['times'])}")
    time_data_custom = get_simulation_times(hours=12, num_points_per_hour=2)
    print(f"Custom times (12h, 2pts/h): {time_data_custom['times']}, len: {len(time_data_custom['times'])}")
    time_data_1pt = get_simulation_times(hours=1, num_points_per_hour=1) # will default to 2 points
    print(f"Custom times (1h, 1pt/h -> forced 2): {time_data_1pt['times']}, len: {len(time_data_1pt['times'])}")


    print("\n--- Testing scale_feed_intake_pattern ---")
    DMI_OF_BASE_PATTERN_TEST = 22.50015 
    scaled_half = scale_feed_intake_pattern(base_feed_intake_pattern, DMI_OF_BASE_PATTERN_TEST / 2, DMI_OF_BASE_PATTERN_TEST)
    scaled_double = scale_feed_intake_pattern(base_feed_intake_pattern, DMI_OF_BASE_PATTERN_TEST * 2, DMI_OF_BASE_PATTERN_TEST)
    
    original_sum_fi_rate_duration = sum([(e["T_eind"] - e["T_begin"]) * e["FI_rate"] for e in base_feed_intake_pattern["feedintake"]])
    scaled_half_sum_fi_rate_duration = sum([(e["T_eind"] - e["T_begin"]) * e["FI_rate"] for e in scaled_half["feedintake"]])
    
    print(f"Original pattern DMI (calculated from pattern): {original_sum_fi_rate_duration:.4f}")
    print(f"Test DMI_OF_BASE_PATTERN_TEST: {DMI_OF_BASE_PATTERN_TEST}")
    print(f"First FI_rate in base: {base_feed_intake_pattern['feedintake'][0]['FI_rate']}")
    print(f"First FI_rate in scaled_half: {scaled_half['feedintake'][0]['FI_rate']:.4f} (expected: {base_feed_intake_pattern['feedintake'][0]['FI_rate']/2:.4f})")
    print(f"Calculated DMI for scaled_half pattern: {scaled_half_sum_fi_rate_duration:.4f} (expected: {DMI_OF_BASE_PATTERN_TEST/2:.4f})")