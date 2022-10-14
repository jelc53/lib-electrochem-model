function y = U_p(x)

% Polynomial expression to calculate Open Circuit Potential (OCP) of
% cathode. Expression has been derived from experimental data.

y_fit = [10188.54, -66535.82, 189316.65, -307780.79, 314825.52, ...
    -209988.26, 91295.30, -24944.01, 3884.75, -258.27];

y = polyval(y_fit,x);
end