function [A_en_d,A_es_d,A_ep_d,B_e_d] = matrices_elecphase(param)

%%% ODE matrix formulation for Electrolyte Phase of all cells in module %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Create State Space Coefficient Matrices corresponding to D(i-1), D(i), and D(i+1)
%1) D(i-1/2) <--> A1
%2) D(i+1/2) <--> A2


%Negative Electrolyte region
A1_mat_n = zeros(param.Nx_n);
A2_mat_n = zeros(param.Nx_n);
Be_mat_n = zeros(param.Nx_n,1);
 
for k = 1:(param.Nx_n-1) %Final state will have all zero coefficients, so leave it out
    if k == 1
        A2_mat_n(k,k) = -1;
        A2_mat_n(k,k+1) = 1;    
        Be_mat_n(k,1) = 1;
    else
        A1_mat_n(k,k-1) = 1;
        A1_mat_n(k,k) = -1;
        
        A2_mat_n(k,k) = -1;
        A2_mat_n(k,k+1) = 1;

        Be_mat_n(k,1) = 1;
    end
end


%Separator Region
A1_mat_s = zeros(param.Nx_s);
A2_mat_s = zeros(param.Nx_s);
Be_mat_s = zeros(param.Nx_s,1);

for k = 2:(param.Nx_s-1) %Final state will have all zero coefficients, so leave it out
        A1_mat_s(k,k-1) = 1;
        A1_mat_s(k,k) = -1;
        
        A2_mat_s(k,k) = -1;
        A2_mat_s(k,k+1) = 1;
end

%Positive Electrolyte Region
A1_mat_p = zeros(param.Nx_p);
A2_mat_p = zeros(param.Nx_p);
Be_mat_p = zeros(param.Nx_p,1);

for k = 2:(param.Nx_p) %Final state will have all zero coefficients, so leave it out
    if k == param.Nx_p
        A1_mat_p(k,k) = -1;
        A1_mat_p(k,k-1) = 1;    
        Be_mat_p(k,1) = 1;
    else
        A1_mat_p(k,k-1) = 1;
        A1_mat_p(k,k) = -1;
        
        A2_mat_p(k,k) = -1;
        A2_mat_p(k,k+1) = 1;
        
        Be_mat_p(k,1) = 1;
    end
end

% A_mat_n = [A1_mat_n; A2_mat_n; A3_mat_n];
% A_mat_s = [A1_mat_s; A2_mat_s; A3_mat_s];
% A_mat_p = [A1_mat_p; A2_mat_p; A3_mat_p];
% Be_mat = [Be_mat_n; Be_mat_s; Be_mat_p];

% Creating a block-diagonal matrix with A matrices for all cells
%kron replaces each non-zero element of Nc x Nc identity matrix w/ (Nr - 1)x(Nr - 1) A_mat
A1_en_d = kron(eye(param.Nc),A1_mat_n);
A2_en_d = kron(eye(param.Nc),A2_mat_n);
B_en_d = repmat(Be_mat_n,param.Nc,1);

A1_es_d = kron(eye(param.Nc),A1_mat_s);
A2_es_d = kron(eye(param.Nc),A2_mat_s);
B_es_d = repmat(Be_mat_s,param.Nc,1);

A1_ep_d = kron(eye(param.Nc),A1_mat_p);
A2_ep_d = kron(eye(param.Nc),A2_mat_p);
B_ep_d = repmat(Be_mat_p,param.Nc,1);

%CONSOLIDATE BEFORE SENDING TO WORKSPACE
A_en_d = sparse([A1_en_d; A2_en_d]);
A_es_d = sparse([A1_es_d; A2_es_d]);
A_ep_d = sparse([A1_ep_d; A2_ep_d]);
B_e_d = sparse([B_en_d; B_es_d; B_ep_d]);

end