function [SOC_cells] = soc_bound(x)

% Function to bound the values of SOC to its physical limits of 0 and 1
h = length(x);
SOC_cells = zeros(h,1);

for i = 1:h
    if x(i )> 1
        SOC_cells(i) = 1;
    elseif x(i) < 0
        SOC_cells(i) = 0;
    else
        SOC_cells(i) = x(i);
    end
end

end
