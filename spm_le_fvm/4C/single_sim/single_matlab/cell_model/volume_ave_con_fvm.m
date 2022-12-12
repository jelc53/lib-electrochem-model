function csi_ave = volume_ave_con_fvm(param,cs_i,flag)
%============The flag indicate cathode (0) or anode(1)


if flag==0 %===========Cathode
cs_p=cs_i;    
dr_p=param.Rs_p/(param.Nr-1);
r_p=linspace(dr_p,param.Rs_p,param.Nr-1);
r_p=[0,r_p];
for kk=2:1:length(r_p)
avecspi(kk-1,:)=cs_p(kk-1,:)*4/3*pi*(r_p(kk)^3-r_p(kk-1)^3);
end
csi_ave=sum(avecspi)/(4/3*pi*(param.Rs_p^3));

elseif flag==1  %===========Anode
cs_n=cs_i;   
dr_n=param.Rs_n/(param.Nr-1);
r_n=linspace(dr_n,param.Rs_n,param.Nr-1);
r_n=[0,r_n];
for kk=2:1:length(r_n)
avecsni(kk-1,:)=cs_n(kk-1,:)*4/3*pi*(r_n(kk)^3-r_n(kk-1)^3);
end
csi_ave=sum(avecsni)/(4/3*pi*(param.Rs_n^3));       
end

end

