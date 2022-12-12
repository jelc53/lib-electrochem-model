function y = U_n(x)

% Analytical expression to calculate Open Circuit Potential (OCP) of
% anode. Expression borrowed from Reference:
% Ji, Y., Zhang, Y., & Wang, C. Y. (2013). 
% Li-ion cell operation at low temperatures. 
% Journal of The Electrochemical Society, 160(4), A636-A649.

y = 0.1493 + 0.8493*exp(-61.79*x)+ ...  
    0.3824*exp(-665.8*x)-exp(39.42*x-41.92)-...
    0.03131*atan(25.59*x-4.099)-...
    0.009434*atan(32.49*x- 15.74);
end