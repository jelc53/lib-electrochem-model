function [soc_bulk_n,soc_bulk_p] = soc_calculation_fvm(cs_n,cs_p,param)




for j = 1:size(cs_n,2) 
for i = 1:param.Nc

dr_n=param.Rs_n/(param.Nr-1);
r_n=linspace(dr_n,param.Rs_n,param.Nr-1);
r_n=[0,r_n];
for kk=2:1:length(r_n)
avecsni(kk-1)=cs_n(kk-1,j)*4/3*pi*(r_n(kk)^3-r_n(kk-1)^3);
end
xbulk_n(i,j)=sum(avecsni)/(4/3*pi*(param.Rs_n^3))/param.c_n_max;   

dr_p=param.Rs_p/(param.Nr-1);
r_p=linspace(dr_p,param.Rs_p,param.Nr-1);
r_p=[0,r_p];
for kk=2:1:length(r_p)
avecspi(kk-1)=cs_p(kk-1,j)*4/3*pi*(r_p(kk)^3-r_p(kk-1)^3);
end
xbulk_p(i,j)=sum(avecspi)/(4/3*pi*(param.Rs_p^3))/param.c_p_max;



soc_bulk_n(i,j) = (xbulk_n(i,j) - param.theta0_n)/(param.theta100_n - param.theta0_n);
soc_bulk_p(i,j) = (param.theta0_p - xbulk_p(i,j))/(param.theta0_p - param.theta100_p);

if soc_bulk_n(i,j)>1
    soc_bulk_n(i,j) = 1;
end
if soc_bulk_p(i,j)<0
    soc_bulk_p(i,j) = 0;
end
    if soc_bulk_n(i,j)>1
        soc_bulk_n(i,j) = 1;
    elseif soc_bulk_n(i,j)<0
        soc_bulk_n(i,j) = 0;
    else
        soc_bulk_n(i,j) = soc_bulk_n(i,j);           
    end
    if soc_bulk_p(i,j)>1
        soc_bulk_p(i,j) = 1;
    elseif soc_bulk_p(i,j)<0
        soc_bulk_p(i,j) = 0;
    else
        soc_bulk_p(i,j) = soc_bulk_p(i,j);           
    end

end
end

end
