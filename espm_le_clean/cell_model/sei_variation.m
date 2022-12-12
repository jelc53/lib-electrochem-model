function [Lsei] = sei_variation(param)   

% Randomize the variation in SEI between cells in a series module
%Using rand instead of randn because it's bounded between 0 and 1
Lsei = param.Lsei_ref*ones(param.Nc,1) + (param.SEI_deviation)*rand(param.Nc,1);

end