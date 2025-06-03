# api/simulation/initial_values.py
from .params import get_model_constants

# Ordered list of state variable names
state_variable_names = [
    "Q_Fd", "QSd", "QWr", "QNa", "QHe", "QMi", "Q_Ac", "Q_Pr", "Q_Bu",
    "QHy", "Q_Me", "fNADH", "QMt"
]

def get_initial_values(dynamic_params_for_constants=None):
    """
    Calculates initial values for state variables.
    Some initial values might depend on model parameters (e.g., fNADH).
    Args:
        dynamic_params_for_constants (dict, optional): Parameters to override base model constants
                                                     if initial values depend on them.
    Returns:
        dict: Containing initialValuesDict, initialValuesVec, and stateVariableNames.
    """
    # Get model constants, possibly influenced by dynamic parameters
    parms = get_model_constants(dynamic_params_for_constants)

    # Base initial values for state variables
    # QMi_initial is fixed here for the fNADH calculation as per original setup
    qmi_initial_fixed = 1.680e+03

    _initial_values_dict_base = {
        "Q_Fd": 8.617267e+02,  # Amount of degradable Fiber (g)
        "QSd": 2.052867e+02,   # Amount of degradable Starch (g)
        "QWr": 3.663752e-57,  # Amount of Water Soluble Carbohydrates (g) - effectively zero
        "QNa": 1.000000e-05,  # Amount of Nitrate (mol)
        "QHe": 4.221166e-02,  # Amount of Hexose (mol)
        "QMi": qmi_initial_fixed,     # Amount of general Microbes (g)
        "Q_Ac": 3.030543e+00,  # Amount of Acetate (mol)
        "Q_Pr": 6.449737e-01,  # Amount of Propionate (mol)
        "Q_Bu": 5.346567e-01,  # Amount of Butyrate (mol)
        "QHy": 4.195278e-04,  # Amount of dissolved Hydrogen (mol)
        "Q_Me": 2.253747e+01,  # Amount of Methanogens (g)
        # fNADH state is Q_NADH (amount of NADH in mol)
        # Initial Q_NADH = initial_fNADH_fraction * c_NAD_param * QMi_initial
        "fNADH": 0.2709195 * parms["c_NAD_param"] * qmi_initial_fixed,
        "QMt": 0.0             # Amount of Methane (mol)
    }

    initial_values_vec_ordered = [_initial_values_dict_base[name] for name in state_variable_names]

    return {
        "initialValuesDict": _initial_values_dict_base,
        "initialValuesVec": initial_values_vec_ordered,
        "stateVariableNames": state_variable_names
    }

# Global instances for convenience or legacy imports, using default constants
default_initial_values_data = get_initial_values()
initial_values_dict = default_initial_values_data["initialValuesDict"]
initial_values_vec = default_initial_values_data["initialValuesVec"]

if __name__ == '__main__':
    data_default = get_initial_values()
    print("Default Initial values dictionary:")
    for key, value in data_default["initialValuesDict"].items():
        print(f"  {key}: {value:.4e}")
    
    print("\nInitial values vector (for ODE solver):")
    for val in data_default["initialValuesVec"]:
        print(f"  {val:.4e}")
    
    print(f"\nUsing overridden c_NAD_param for initial fNADH calculation example:")
    data_dynamic_params = get_initial_values(dynamic_params_for_constants={"c_NAD_param": 8e-06})
    print(f"  fNADH with c_NAD_param=7e-06 (default): {data_default['initialValuesDict']['fNADH']:.4e}")
    print(f"  fNADH with c_NAD_param=8e-06 (dynamic): {data_dynamic_params['initialValuesDict']['fNADH']:.4e}")