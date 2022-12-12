function i_lpl = side_current_plating_aftersolve(R_sei,T_core,phi_sn,phi_e,param,input_crt)
            
for i = 1:param.Nc
    i_lpl(i,1) = -param.i0_li*exp(-param.alpha_li*param.F./(param.Rg*T_core(i)).*(phi_sn(i)-phi_e(i)-R_sei(i)*input_crt));
end

end