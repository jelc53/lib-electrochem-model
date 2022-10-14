function [value,isterminal,direction] = physical_event_function_v2(t,x,param)
% when value is equal to zero, an event is triggered.
% set isterminal to 1 to stop the solver at the first event, or 0 to
% get all the events.
%  direction=0 if all zeros are to be computed (the default), +1 if
%  only zeros where the event function is increasing, and -1 if only
%  zeros where the event function is decreasing.

value1 = min(x(1:param.Nc*(param.Nr-1))) - 60; %Cs_n < 10, end integration
% value1 = min(x(1:param.Nc*(param.Nr-1))) - param.theta0_n*param.c_n_max; %Cs_n < 10, end integration
value2 = max(x(1:param.Nc*(param.Nr-1))) - (param.theta100_n*param.c_n_max-10);
value3 = min(x(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1))) - 18000;
% value3 = min(x(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1))) - param.theta100_p*param.c_p_max;
value4 = max(x(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1))) - (param.theta0_p*param.c_p_max-10); 
value5 = isreal(x);

value = [value1; value2; value3; value4; value5];  % when value = 0, an event is triggered
isterminal = [1 ; 1; 1; 1; 1]; % terminate after the first event
direction = [-1; 1; -1; 1; 0];  % get all the zeros