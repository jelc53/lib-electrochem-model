function eps_el_p = lam_effect(param,aina_p,af_p)

% Modified electrolyte porosity
eps_el_p = param.eps_el_p + (aina_p-af_p)*param.Rs_p/3;

end