function [At_mat, Bt_mat] = matrices_thermal(param)

%%%% ODE matrix formulation for two-state thermal model of all cells %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%NOTE: EDITED FROM ORIGINAL
% ODE: State-space formulation
At_mat = zeros(param.T_states); %Ex. 12 x 12 matrix
Bt_mat = zeros(param.T_states,2); %Ex. 12 x 2 matrix
for k = 1:(param.T_states)
    %mod = remainder of k / 2
    if mod(k,2) == 0 % even (Surface temperature ODEs)
        At_mat(k,k-1) = 1/(param.Rc*param.Cs);
        Bt_mat(k,2) = 1/(param.Ru*param.Cs); 
        
        if param.T_states == 2
            %Applies for Nc = 1 (Single cell, no interconnection)
            At_mat(k,k) = -(1/param.Ru + 1/param.Rc)*(1/param.Cs);
        else
            if k == 2 %Skin temp. of 1st cell
                At_mat(k,k+2) = 1/(param.Rm*param.Cs);
                At_mat(k,k) = -(1/param.Ru + 1/param.Rc + 1/param.Rm)*(1/param.Cs);
            elseif k == param.T_states %Skin temp. of last cell
                At_mat(k,k-2) = 1/(param.Rm*param.Cs);
                At_mat(k,k) = -(1/param.Ru + 1/param.Rc + 1/param.Rm)*(1/param.Cs);
            else %Skin temp. of internal cells (2 adjacent neighbors)
                At_mat(k,k-2) = 1/(param.Rm*param.Cs);            
                At_mat(k,k+2) = 1/(param.Rm*param.Cs);
%                 At_mat(k,k) = -(1/param.Ru + 1/param.Rc + 2/param.Rm)*(1/param.Cs);
                At_mat(k,k) = -((4/5)*1/param.Ru + 1/param.Rc + 2/param.Rm)*(1/param.Cs);
                Bt_mat(k,2) = 4/5*(1/(param.Ru*param.Cs));
                %NOTICE FACTOR OF 4/5 - come up w/ a reasonable value for
                %this scaling
            end
        end
    else % odd (Core temperature ODEs)
        At_mat(k,k) = -1/(param.Rc*param.Cc);
        At_mat(k,k+1) = 1/(param.Rc*param.Cc);
        Bt_mat(k,1) = 1/(param.Cc);
    end
end
At_mat = sparse(At_mat);
Bt_mat = sparse(Bt_mat);

end