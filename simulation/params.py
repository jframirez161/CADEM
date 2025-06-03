# api/simulation/params.py

# Base model constants
_base_model_constants = {
    "DMI": 17.9, "NDF": 287.123, "St": 158.6, "WSC": 125.26, "Acin": 10.68, "Prin": 1.78, "Buin": 1.78, "Lain": 21.36,
    "k_SdHe": 0.078218, "k_FdHe": 0.042862, "k_PdPs": 0.053613, "k_WrHe": 14.0, "c_For": 0.7,
    "c_HNO3_param": 0.16,
    "V_headspace": 40.0, "V_mol": 25.0, "C_starMi": 12.5, "k_H": 1382.0, "f_b": 703.0,
    "M_He_PsMi": 0.0224, "M_He_AmMi": 0.00135, "J_Am_HeVf": 0.00861, "J_Ps_HeVf": 0.01465,
    "v_AcAb": 0.07271622, "v_PrAb": 0.226362, "v_BuAb": 1.012668, "k_NPAb": 0.3, "M_Ac_AcAb": 0.0791,
    "M_Pr_PrAb": 0.112, "M_Bu_BuAb": 0.4934, "M_He_HeMi": 0.02, "M_He_HeVa": 0.055,
    "q_GM": 3.0, "M_NAD_Ac": 9.0, "J_NAD_AP": 1.0, "k_NADH_Fe": 202.0, "k_NaAb": 0.3,
    "k_HyEm": 7.490259, "v_HeAc": 0.05277139, "v_HeAP": 0.007528171, "v_HeBu": 0.002435758, "v_HyMe": 0.2935437,
    "M_Hy_HyMe": 6.289578e-07, "v_VfAb": 0.8384857, "V_fluid": 79.34609999999999, "k_FlEx": 0.1281074, "k_SoEx": 0.04528851,
    "k_LiEx": 0.10615711252653928, "k_NaAm_param": 6.987,
    "Y_Ac_HeAc": 2.0, "Y_Bu_HeBu": 1.0, "Y_Ac_HeAP": 0.6666667, "Y_Pr_HeAP": 1.3333333,
    "Y_AcMi_HeMi": 84.25, "Y_BuMi_HeMi": 68.95, "Y_APMi_HeMi": 77.86,
    "Y_He_LaHe": 0.00278, "Y_Hy_PsMi": 0.00058, "X_Hy_AmMi": 0.00041, "Y_Hy_HeAc": 4.0,
    "Y_Hy_HeBu": 2.0, "Y_Hy_NaAm_param": 4.0,
    "Y_Hy_NaNi": 1.0, "Y_Hy_NiAm": 3.0,
    "c_PdMi": 0.55, "c_NAD_param": 7e-06,
    "Y_Mt_HyMe": 0.25, "Y_Me_HyMt": 2.0,
    "Y_Mt_FmLi": 0.002, "f_He_HeAc": 0.645, "f_He_HeBu": 0.706, "f_He_HeAP": 0.67,
    "f_SdSi": 0.1, "f_PdSi": 0.25
    # Add any other base parameters your model uses
}

def get_model_constants(dynamic_params=None):
    """
    Returns the model constants, allowing for overrides with dynamic_params.
    These dynamic_params typically come from diet composition (NDF, St, WSC, etc. in g/kg DM).
    The DMI parameter within dynamic_params is usually handled separately for scaling feed intake,
    but if the model itself has a DMI parameter that influences other constants directly, it could be included.
    """
    constants = _base_model_constants.copy() # Start with a copy of the base constants
    if dynamic_params:
        constants.update(dynamic_params) # Override with any dynamic parameters passed
    return constants

# This global instance is for convenience if some modules absolutely need to import
# a default set, but it's better to call get_model_constants() for specific simulation runs.
model_constants = get_model_constants()

if __name__ == '__main__':
    print(f"Number of base parameters defined: {len(_base_model_constants)}")
    default_set = get_model_constants()
    print(f"V_fluid (base): {default_set['V_fluid']} L")
    
    test_overrides = {"NDF": 300.0, "St": 200.0, "NewParamTest": 123}
    test_constants = get_model_constants(test_overrides)
    print(f"NDF (overridden): {test_constants['NDF']}")
    print(f"St (overridden): {test_constants['St']}")
    print(f"WSC (base, not overridden): {test_constants['WSC']}")
    print(f"NewParamTest (dynamic): {test_constants['NewParamTest']}")
    print(f"Original DMI still present if not overridden: {test_constants['DMI']}")