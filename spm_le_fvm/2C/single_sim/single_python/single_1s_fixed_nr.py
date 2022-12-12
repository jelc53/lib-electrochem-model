
import numpy as np
import pybamm
import matplotlib.pyplot as plt
model = pybamm.lithium_ion.SPM()
parameter_values = pybamm.ParameterValues("Chen2020")
import pandas as pd
Nr=10;
use_par=pd.read_csv('use_var.csv')
use_par=use_par.values[:,0]

parameter_values['Upper voltage cut-off [V]'] = 5.5
parameter_values['Lower voltage cut-off [V]'] = 1.5
parameter_values['Positive particle radius [m]'] = use_par[0]
parameter_values['Negative particle radius [m]'] = use_par[1]


parameter_values[ 'Maximum concentration in positive electrode [mol.m-3]'] = use_par[2]
parameter_values[ 'Maximum concentration in negative electrode [mol.m-3]'] = use_par[3]

parameter_values['Initial concentration in positive electrode [mol.m-3]'] = use_par[4]
parameter_values['Initial concentration in negative electrode [mol.m-3]'] = use_par[5]

parameter_values['Positive electrode thickness [m]'] = use_par[6]
parameter_values['Separator thickness [m]'] = use_par[7]
parameter_values['Negative electrode thickness [m]'] = use_par[8]



parameter_values['Positive electrode porosity'] = use_par[9]
parameter_values['Separator porosity'] = use_par[10]
parameter_values['Negative electrode porosity'] = use_par[11]

parameter_values['Positive electrode active material volume fraction'] = use_par[12]
parameter_values['Negative electrode active material volume fraction'] = use_par[13]

parameter_values['Positive electrode diffusivity [m2.s-1]'] = use_par[14]
parameter_values['Negative electrode diffusivity [m2.s-1]'] = use_par[15]





def my_kp(c_e, c_s_surf, c_s_max, T):
    # m_ref = 3.42e-6  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
    m_ref=use_par[16]*96487
    E_r = 17800
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))
    return (m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5)



def my_kn(c_e, c_s_surf, c_s_max, T):
    # m_ref = 6.48e-7  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
    m_ref = use_par[17]*96487
    E_r = 35000
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))
    return (m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5)

parameter_values['Negative electrode exchange-current density [A.m-2]'] = my_kn
parameter_values['Positive electrode exchange-current density [A.m-2]'] = my_kp


parameter_values["Electrode height [m]"] = 1
parameter_values["Electrode width [m]"] = use_par[18]


 
parameter_values['Current function [A]'] = -9.74
tot_sim_time=1200




# experiment = pybamm.Experiment(
#     [
#         ("Charge at 19.5 A for 50 seconds (1 second period)"),
#     ] 
# )


geometry = model.default_geometry

# load parameter values and process model and geometry
param = parameter_values
param.process_model(model)
param.process_geometry(geometry)

var = pybamm.standard_spatial_vars
var_pts = {var.x_n: 15, var.x_s: 10, var.x_p: 15, var.r_n: Nr, var.r_p: Nr}
mesh = pybamm.Mesh(geometry, model.default_submesh_types, var_pts)

# discretise model
spatial_methods = {
    "macroscale": pybamm.FiniteVolume(),
    "negative particle": pybamm.FiniteVolume(),
    "positive particle": pybamm.FiniteVolume(),
    "current collector": pybamm.ZeroDimensionalSpatialMethod()
}
disc = pybamm.Discretisation(mesh, spatial_methods)


disc.process_model(model)

sim = pybamm.Simulation(model, parameter_values=parameter_values)



t_eval = np.linspace(0,tot_sim_time,tot_sim_time+1)
sim.solve(t_eval=t_eval)


solution = sim.solution
t = solution["Time [s]"]
V = solution["Terminal voltage [V]"]
V=V.entries
t=t.entries






#==============E
U_p = solution["X-averaged positive electrode open circuit potential [V]"]
U_p=U_p.entries
U_n = solution["X-averaged negative electrode open circuit potential [V]"]
U_n=U_n.entries


E=U_p-U_n
phisp=solution["X-averaged positive electrode reaction overpotential [V]"]
phisp=phisp.entries

phisn=solution["X-averaged negative electrode reaction overpotential [V]"]
phisn=phisn.entries

plt.plot(V)






#=======All concentraion
csp=solution["X-averaged positive particle concentration [mol.m-3]"]
csp=csp.entries

csn=solution["X-averaged negative particle concentration [mol.m-3]"]
csn=csn.entries

#=======Surface concentraion
cspsur=solution["X-averaged positive particle surface concentration [mol.m-3]"]
cspsur=cspsur.entries

csnsur=solution["X-averaged negative particle surface concentration [mol.m-3]"]
csnsur=csnsur.entries




#=======Volume average concentraion
avecsp=solution["Positive electrode volume-averaged concentration [mol.m-3]"]
avecsp=avecsp.entries

avecsn=solution["Negative electrode volume-averaged concentration [mol.m-3]"]
avecsn=avecsn.entries


#==========t_spanfor pybamm
tspan=t



np.savetxt('pybamm_t_single_nr'+str(Nr)+'.txt',tspan)
np.savetxt('pybamm_V_single_nr'+str(Nr)+'.txt',V)
np.savetxt('pybamm_cspsur_single_nr'+str(Nr)+'.txt',cspsur)
np.savetxt('pybamm_csnsur_single_nr'+str(Nr)+'.txt',csnsur)
np.savetxt('pybamm_cspave_single_nr'+str(Nr)+'.txt',avecsp)
np.savetxt('pybamm_csnave_single_nr'+str(Nr)+'.txt',avecsn)

