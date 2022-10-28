function y = U_n(x)

% [chen] Development of Experimental Techniques for Parameterization of Multiscale Lithium-ion Battery Models

y = 1.9793*exp(-39.3631*x) + 0.2482 -0.0909*tanh(29.8538*(x-0.1234))...
    -0.04478*tanh(14.9159*(x-0.2769))-0.0205*tanh(30.4444*(x-0.6103));
end