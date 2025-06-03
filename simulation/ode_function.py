# api/simulation/ode_function.py
import numpy as np
# Removed direct imports of model_constants and initial_values_dict as they are passed in.
# from .params import model_constants
# from .initial_values import initial_values_dict, state_variable_names

def ode_function_py(t, state_vec, parms_dict, feedintake_data_ext, state_variable_names_arg, delt_ext=None):
    """
    Defines the system of ordinary differential equations for the rumen model.

    Args:
        t (float): Current time.
        state_vec (np.ndarray): Array of current state variable values.
        parms_dict (dict): Dictionary of model parameters.
        feedintake_data_ext (dict): Feed intake schedule data.
        state_variable_names_arg (list): Ordered list of state variable names.
        delt_ext (float, optional): Time step (not directly used by solve_ivp's adaptive steps).

    Returns:
        np.ndarray: Array of derivatives (d(state)/dt).
    """

    state = dict(zip(state_variable_names_arg, state_vec))

    Q_Fd = state["Q_Fd"]
    QSd = state["QSd"]
    QWr = state["QWr"]
    QNa = state["QNa"]
    QHe = state["QHe"]
    QMi = state["QMi"]
    Q_Ac = state["Q_Ac"]
    Q_Pr = state["Q_Pr"]
    Q_Bu = state["Q_Bu"]
    QHy = state["QHy"]
    Q_Me = state["Q_Me"]
    fNADH = state["fNADH"] # Represents Q_NADH (mol)
    QMt = state["QMt"]

    # --- Concentrations ---
    C_Na = QNa / parms_dict["V_fluid"] if parms_dict["V_fluid"] > 1e-9 else 0
    C_He = QHe / parms_dict["V_fluid"] if parms_dict["V_fluid"] > 1e-9 else 0
    C_Mi = QMi / parms_dict["V_fluid"] if parms_dict["V_fluid"] > 1e-9 else 0
    C_Ps = 3e-3
    C_Am = 5e-3
    C_Ac = Q_Ac / parms_dict["V_fluid"] if parms_dict["V_fluid"] > 1e-9 else 0
    C_Pr = Q_Pr / parms_dict["V_fluid"] if parms_dict["V_fluid"] > 1e-9 else 0
    C_Bu = Q_Bu / parms_dict["V_fluid"] if parms_dict["V_fluid"] > 1e-9 else 0
    
    C_VFA_sum = C_Ac + C_Pr + C_Bu
    
    P_Hy = (QHy * (1e3 * 8.3145 * 312 / 1.01325e5)) / parms_dict["V_headspace"] if parms_dict["V_headspace"] > 1e-9 else 0
    C_Hy = P_Hy / parms_dict["k_H"] if parms_dict["k_H"] > 1e-9 else 0
    
    pH = 7.73 - (0.014 * 1e3 * C_VFA_sum)
    pH_cell = 6.43 + 3.62e-8 * np.exp(2.4 * pH)
    
    Q_NAD_total = parms_dict["c_NAD_param"] * QMi
    Q_NAD_plus = max(0, Q_NAD_total - fNADH) # fNADH is Q_NADH

    if fNADH > 1e-9:
        r_NAD = Q_NAD_plus / fNADH
    else:
        r_NAD = float('inf')

    # --- DMI Rate ---
    D_DMI = 0.0
    current_t_mod_24 = t % 24
    for fi_entry in feedintake_data_ext["feedintake"]:
        if fi_entry["T_begin"] <= current_t_mod_24 < fi_entry["T_eind"]:
            D_DMI = fi_entry["FI_rate"]
            break

    # --- Process Rates (P, U variables) ---
    P_Fd_InFd = D_DMI * parms_dict["NDF"]
    U_Fd_FdHe = parms_dict["k_FdHe"] * (C_Mi / parms_dict["C_starMi"] if parms_dict["C_starMi"] > 1e-9 else 0) * Q_Fd
    U_Fd_FdEx = parms_dict["k_SoEx"] * Q_Fd

    P_Sd_InSd = D_DMI * parms_dict["St"]
    U_Sd_SdHe = parms_dict["k_SdHe"] * (C_Mi / parms_dict["C_starMi"] if parms_dict["C_starMi"] > 1e-9 else 0) * QSd
    U_Sd_SdEx = parms_dict["k_SoEx"] * QSd

    P_Wr_InWr = D_DMI * parms_dict["WSC"]
    U_Wr_WrHe = parms_dict["k_WrHe"] * (C_Mi / parms_dict["C_starMi"] if parms_dict["C_starMi"] > 1e-9 else 0) * QWr
    U_Wr_WrEx = parms_dict["k_FlEx"] * QWr

    c_HNO3_local = 0.0
    P_Na_InNa = D_DMI * c_HNO3_local
    U_Na_NaAb = parms_dict["k_NaAb"] * QNa
    U_Na_NaEx = parms_dict["k_FlEx"] * QNa
    U_Na_NaAm = parms_dict["k_NaAm_param"] * QMi * QNa * QHy

    P_He_LaHe = parms_dict["Y_He_LaHe"] * D_DMI * parms_dict["Lain"] / 89.0
    P_He_XdHe = (U_Fd_FdHe + U_Sd_SdHe + U_Wr_WrHe) / 162.0

    den_AcHe_HeMi = (1 + parms_dict["M_He_HeMi"] / (C_He if C_He > 1e-9 else 1e-9)) * \
                    (1 + parms_dict["M_NAD_Ac"] / (r_NAD if r_NAD > 1e-9 else 1e-9))
    den_BuHe_HeMi = (1 + parms_dict["M_He_HeMi"] / (C_He if C_He > 1e-9 else 1e-9))
    den_APHe_HeMi = (1 + parms_dict["M_He_HeMi"] / (C_He if C_He > 1e-9 else 1e-9)) * \
                    (1 + r_NAD / (parms_dict["J_NAD_AP"] if parms_dict["J_NAD_AP"] > 1e-9 else 1e-9))
    
    U_AcHe_HeMi = (parms_dict["v_HeAc"] / (parms_dict["q_GM"] if parms_dict["q_GM"] > 1e-9 else 1e-9)) * QMi / (den_AcHe_HeMi if den_AcHe_HeMi > 1e-9 else 1e-9)
    U_BuHe_HeMi = (parms_dict["v_HeBu"] / (parms_dict["q_GM"] if parms_dict["q_GM"] > 1e-9 else 1e-9)) * QMi / (den_BuHe_HeMi if den_BuHe_HeMi > 1e-9 else 1e-9)
    U_APHe_HeMi = (parms_dict["v_HeAP"] / (parms_dict["q_GM"] if parms_dict["q_GM"] > 1e-9 else 1e-9)) * QMi / (den_APHe_HeMi if den_APHe_HeMi > 1e-9 else 1e-9)

    den_He_HeVf_common = (1 + parms_dict["M_He_HeVa"] / (C_He if C_He > 1e-9 else 1e-9)) * \
                         (1 + C_Am / (parms_dict["J_Am_HeVf"] if parms_dict["J_Am_HeVf"] > 1e-9 else 1e-9)) * \
                         (1 + C_Ps / (parms_dict["J_Ps_HeVf"] if parms_dict["J_Ps_HeVf"] > 1e-9 else 1e-9))

    U_He_HeAc = parms_dict["v_HeAc"] * QMi / ((den_He_HeVf_common if den_He_HeVf_common > 1e-9 else 1e-9) * \
                (1 + parms_dict["M_NAD_Ac"] / (r_NAD if r_NAD > 1e-9 else 1e-9)))
    U_He_HeBu = parms_dict["v_HeBu"] * QMi / (den_He_HeVf_common if den_He_HeVf_common > 1e-9 else 1e-9)
    U_He_HeAP = parms_dict["v_HeAP"] * QMi / ((den_He_HeVf_common if den_He_HeVf_common > 1e-9 else 1e-9) * \
                (1 + r_NAD / (parms_dict["J_NAD_AP"] if parms_dict["J_NAD_AP"] > 1e-9 else 1e-9)))
    U_He_HeEx = parms_dict["k_FlEx"] * QHe

    P_Mi_HeMi = parms_dict["Y_AcMi_HeMi"] * U_AcHe_HeMi + \
                parms_dict["Y_BuMi_HeMi"] * U_BuHe_HeMi + \
                parms_dict["Y_APMi_HeMi"] * U_APHe_HeMi
    U_Mi_MiEx = (0.15 * parms_dict["k_FlEx"] + 0.65 * parms_dict["k_SoEx"]) * QMi

    P_Ac_in = D_DMI * parms_dict["Acin"] / 59.0
    P_Ac_HeAc = parms_dict["Y_Ac_HeAc"] * (parms_dict["f_He_HeAc"] * U_AcHe_HeMi + U_He_HeAc) + \
                parms_dict["Y_Ac_HeAP"] * (parms_dict["f_He_HeAP"] * U_APHe_HeMi + U_He_HeAP)
    den_Ac_Ab = (1 + (parms_dict["M_Ac_AcAb"] / (C_Ac if C_Ac > 1e-9 else 1e-9))**1.17) * \
                (1 + (pH / (6.02 if 6.02 > 1e-9 else 1e-9))**3.91)
    U_Ac_AcAb = (parms_dict["v_VfAb"] * parms_dict["v_AcAb"] * Q_Ac * (parms_dict["V_fluid"]**0.75)) / \
                (den_Ac_Ab if den_Ac_Ab > 1e-9 else 1e-9)
    U_Ac_AcEx = parms_dict["k_FlEx"] * Q_Ac

    P_Pr_in = D_DMI * parms_dict["Prin"] / 73.0
    P_Pr_HePr = parms_dict["Y_Pr_HeAP"] * (parms_dict["f_He_HeAP"] * U_APHe_HeMi + U_He_HeAP)
    den_Pr_Ab = (1 + (parms_dict["M_Pr_PrAb"] / (C_Pr if C_Pr > 1e-9 else 1e-9))**0.95) * \
                (1 + (pH / (6.02 if 6.02 > 1e-9 else 1e-9))**4.61)
    U_Pr_PrAb = (parms_dict["v_VfAb"] * parms_dict["v_PrAb"] * Q_Pr * (parms_dict["V_fluid"]**0.75)) / \
                (den_Pr_Ab if den_Pr_Ab > 1e-9 else 1e-9)
    U_Pr_PrEx = parms_dict["k_FlEx"] * Q_Pr

    P_Bu_in = D_DMI * parms_dict["Buin"] / 87.0
    P_Bu_HeBu = parms_dict["Y_Bu_HeBu"] * (parms_dict["f_He_HeBu"] * U_BuHe_HeMi + U_He_HeBu)
    den_Bu_Ab = (1 + (parms_dict["M_Bu_BuAb"] / (C_Bu if C_Bu > 1e-9 else 1e-9))**0.99) * \
                (1 + (pH / (6.02 if 6.02 > 1e-9 else 1e-9))**5.13)
    U_Bu_BuAb = (parms_dict["v_VfAb"] * parms_dict["v_BuAb"] * Q_Bu * (parms_dict["V_fluid"]**0.75)) / \
                (den_Bu_Ab if den_Bu_Ab > 1e-9 else 1e-9)
    U_Bu_BuEx = parms_dict["k_FlEx"] * Q_Bu

    P_Hy_HeAc = parms_dict["Y_Hy_HeAc"] * (parms_dict["f_He_HeAc"] * U_AcHe_HeMi + U_He_HeAc)
    P_Hy_HeBu = parms_dict["Y_Hy_HeBu"] * (parms_dict["f_He_HeBu"] * U_BuHe_HeMi + U_He_HeBu)
    U_Hy_HyMe = (parms_dict["v_HyMe"] * Q_Me) / (1 + parms_dict["M_Hy_HyMe"] / (C_Hy if C_Hy > 1e-9 else 1e-9))
    U_Hy_NaAm = parms_dict["Y_Hy_NaAm_param"] * U_Na_NaAm

    term_HyEm_factor = 0
    if parms_dict["V_headspace"] > 1e-9 and parms_dict["k_H"] > 1e-9:
        term_HyEm_factor = parms_dict["f_b"] * parms_dict["V_mol"] / (parms_dict["V_headspace"] * parms_dict["k_H"])
    U_Hy_HyEm = (parms_dict["k_HyEm"] + term_HyEm_factor) * QHy
    
    term_HyEx_factor = 0
    if parms_dict["V_headspace"] > 1e-9 and parms_dict["k_H"] > 1e-9:
         term_HyEx_factor = parms_dict["V_fluid"] * parms_dict["k_FlEx"] * parms_dict["V_mol"] / (parms_dict["V_headspace"] * parms_dict["k_H"])
    U_Hy_HyEx = term_HyEx_factor * QHy

    P_Me_HyMt = parms_dict["Y_Me_HyMt"] * parms_dict["Y_Mt_HyMe"] * U_Hy_HyMe
    U_Me_MeEx = (0.4 * parms_dict["k_SoEx"] + 0.4 * parms_dict["k_FlEx"]) * Q_Me

    den_FdLi = parms_dict["k_FdHe"] + parms_dict["k_LiEx"]
    den_SdLi = parms_dict["k_SdHe"] + parms_dict["k_LiEx"]
    den_PdLi = parms_dict["k_PdPs"] + parms_dict["k_LiEx"]
    term1_FmLi = U_Fd_FdEx * (parms_dict["k_FdHe"] / (den_FdLi if den_FdLi > 1e-9 else 1e-9))
    term2_FmLi = U_Sd_SdEx * parms_dict["f_SdSi"] * (parms_dict["k_SdHe"] / (den_SdLi if den_SdLi > 1e-9 else 1e-9))
    term3_FmLi = 2 * (U_Mi_MiEx + U_Me_MeEx) * parms_dict["c_PdMi"] * parms_dict["f_PdSi"] * \
                 (parms_dict["k_PdPs"] / (den_PdLi if den_PdLi > 1e-9 else 1e-9))
    P_FmLi = term1_FmLi + term2_FmLi + term3_FmLi

    P_Mt_HyMe = parms_dict["Y_Mt_HyMe"] * U_Hy_HyMe
    P_Mt_MtLi = P_FmLi * parms_dict["Y_Mt_FmLi"]

    P_NADH_Ac = 2 * (parms_dict["f_He_HeAc"] * U_AcHe_HeMi + U_He_HeAc)
    U_NADH_AP = 0.67 * (parms_dict["f_He_HeAP"] * U_APHe_HeMi + U_He_HeAP)
    term_F_T_Fd_exp = np.exp((-102e3) / (2 * 8.3145 * 312))
    val_for_sqrt = 0
    pH_cell_term_cubed = (10**(-1 * pH_cell))**3
    if pH_cell_term_cubed > 1e-12:
         val_for_sqrt = (9 * P_Hy**2) / (9 * pH_cell_term_cubed)
         if val_for_sqrt < 0: val_for_sqrt = 0
    F_T_Fd = 1 - (np.sqrt(val_for_sqrt) * term_F_T_Fd_exp)
    U_NADH_Fe = 0.0
    if F_T_Fd > 0:
        U_NADH_Fe = parms_dict["k_NADH_Fe"] * fNADH * F_T_Fd # fNADH is Q_NADH (mol)

    # --- Derivatives ---
    dFd_dt = P_Fd_InFd - U_Fd_FdHe - U_Fd_FdEx
    dSd_dt = P_Sd_InSd - U_Sd_SdHe - U_Sd_SdEx
    dWr_dt = P_Wr_InWr - U_Wr_WrHe - U_Wr_WrEx
    dNa_dt = D_DMI #P_Na_InNa - U_Na_NaAb - U_Na_NaEx - U_Na_NaAm 
    dHe_dt = P_He_LaHe + P_He_XdHe - U_AcHe_HeMi - U_APHe_HeMi - U_BuHe_HeMi - \
             U_He_HeAc - U_He_HeAP - U_He_HeBu - U_He_HeEx
    dMi_dt = P_Mi_HeMi - U_Mi_MiEx
    dAc_dt = P_Ac_HeAc + P_Ac_in - U_Ac_AcAb - U_Ac_AcEx
    dPr_dt = P_Pr_HePr + P_Pr_in - U_Pr_PrAb - U_Pr_PrEx
    dBu_dt = P_Bu_HeBu + P_Bu_in - U_Bu_BuAb - U_Bu_BuEx
    dfNADH_dt = P_NADH_Ac - U_NADH_AP - U_NADH_Fe # This is dQ_NADH/dt
    dHy_dt = P_Hy_HeAc + P_Hy_HeBu - U_Hy_HyMe - U_Hy_NaAm - U_Hy_HyEm - U_Hy_HyEx
    dMe_dt = P_Me_HyMt - U_Me_MeEx
    dMt_dt = P_Mt_HyMe + P_Mt_MtLi

    derivatives = np.array([
        dFd_dt, dSd_dt, dWr_dt, dNa_dt, dHe_dt, dMi_dt, dAc_dt, dPr_dt, dBu_dt,
        dHy_dt, dMe_dt, dfNADH_dt, dMt_dt
    ])
    return derivatives

if __name__ == '__main__':
    print("Testing ode_function_py with initial values and t=0...")
    # To run this test block standalone, you need params.py, initial_values.py, 
    # and simulation_setup.py in the same directory or adjust Python path.
    # For this example, assume they are in a "simulation" subfolder and we adjust imports.
    # This test block is more illustrative when part of a package structure.
    try:
        # Need to import from the current package for the test to run standalone
        from .params import get_model_constants
        from .initial_values import get_initial_values
        from .simulation_setup import base_feed_intake_pattern, get_simulation_times

        test_sim_params = {"hours": 24, "points_per_hour": 10} # Low res for test
        test_diet_params = {"NDF": 287, "St": 150, "WSC":120, "Acin":10, "Prin":2, "Buin":2, "Lain":20} # Example diet

        parms_for_test = get_model_constants(test_diet_params)
        initial_values_data_test = get_initial_values() # Uses default params for initial values
        y0_test = np.array(initial_values_data_test["initialValuesVec"])
        state_names_test = initial_values_data_test["stateVariableNames"]
        
        time_data_test = get_simulation_times(test_sim_params["hours"], test_sim_params["points_per_hour"])
        stpsz_value_test = time_data_test["times"][1] - time_data_test["times"][0] if len(time_data_test["times"]) > 1 else 0


        derivatives_at_t0 = ode_function_py(0.0, y0_test, parms_for_test, base_feed_intake_pattern, state_names_test, stpsz_value_test)
        print("\nCalculated derivatives at t=0:")
        for name, deriv in zip(state_names_test, derivatives_at_t0):
            print(f"  d{name}/dt: {deriv:.4e}")

    except ImportError as e:
        print(f"ImportError during test: {e}. Ensure files are in correct package structure (e.g. 'simulation' folder with __init__.py).")
    except Exception as e:
        import traceback
        print(f"An error occurred during the test run: {e}")
        traceback.print_exc()