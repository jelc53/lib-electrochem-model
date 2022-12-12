function [value,isterminal,direction] = physical_event_function_tk(t,x,param)
% when value is equal to zero, an event is triggered.
% set isterminal to 1 to stop the solver at the first event, or 0 to
% get all the events.
% direction=0 if all zeros are to be computed (the default), +1 if
% only zeros where the event function is increasing, and -1 if only
% zeros where the event function is decreasing.
TIMENEW = datetime('now','Format','HHmmss');

value1 = min(x(1:param.Nc*(param.Nr-1))) - (param.theta0_n*param.c_n_max); %Cs_n < 10, end integration
value1 = 1;
value2 = max(x(1:param.Nc*(param.Nr-1))) - (param.theta100_n*param.c_n_max-10);
value3 = min(x(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1))) - (param.theta100_p*param.c_p_max+10);
value4 = max(x(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1))) - (param.theta0_p*param.c_p_max); 
value4 = 1;
value5 = isreal(x);
% value6 = seconds(TIMENEW-param.TIMEOLD) < 10; 
value6 = seconds(TIMENEW-param.TIMEOLD) < 360;
%
value = [value1; value2; value3; value4; value5; value6];  % when value = 0, an event is triggered
isterminal = [1 ; 1; 1; 1; 1; 1]; % terminate after the first event
direction = [-1; 1; -1; 1; 0; 0]; % get all the zeros

% New stop conditions - KM June 6 2022
value = [value5; value6];  % when value = 0, an event is triggered
isterminal = [1; 1]; % terminate after the first event
direction = [0; 0]; % get all the zeros