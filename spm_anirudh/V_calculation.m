function [V_cell] = V_calculation(ocp_p,ocp_n,eta_p,eta_n,I,param)


V_cell = ocp_p - ocp_n + eta_p - eta_n  - param.R_l*I;

for i = 1:param.Nc
    if V_cell(i)<2.49
        V_cell(i) = 2.49;
    else
        V_cell(i) = V_cell(i);
    end
    if isreal(V_cell(i))
        V_cell(i) = V_cell(i);
    else
        V_cell(i) = 2.49;
    end      
end
end
