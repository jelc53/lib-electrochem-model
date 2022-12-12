function [Asei_mat, Bsei_mat] = matrices_sei(param)

% ODE matrix formulation of SEI layer (one-state) for all cells in module %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% ODE: State-space formulation
Asei_mat = zeros(param.Nc);
b_sei = [-param.M_sei/(2*param.F*param.rho_sei) ...
         -param.beta_li2sei*param.M_sei/(2*param.F*param.rho_sei)];
Bsei_mat = repmat(b_sei,param.Nc,1); %Create Nc x 1 vector w/ b_sei

end
