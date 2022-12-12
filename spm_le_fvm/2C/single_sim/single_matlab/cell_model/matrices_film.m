function [Afilm_mat, Bfilm_mat] = matrices_film(param)

% ODE matrix formulation of fils layer (one-state) for all cells in module %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% ODE: State-space formulation
Afilm_mat = zeros(param.Nc);
b_film = [-param.M_sei/(2*param.F*param.rho_sei) ...
         -param.beta_li2sei*param.M_sei/(2*param.F*param.rho_sei)-(1-param.beta_li2sei)*param.M_li/(2*param.F*param.rho_li)];
Bfilm_mat = repmat(b_film,param.Nc,1); %Create Nc x 1 vector w/ b_sei

end
