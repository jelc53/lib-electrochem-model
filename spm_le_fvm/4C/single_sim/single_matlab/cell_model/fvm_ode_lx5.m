function dudt = fvm_ode_lx5(t_in, u_in, param)
us_n=u_in(1:param.Nr-1);
us_p=u_in(param.Nr:end);
tuse=ceil(t_in);
input_crt=param.I_data(tuse+1);
% %=========Anode dr
dr_n=param.Rs_n/(param.Nr-1);
% r_n=linspace(dr_n,param.Rs_n,param.Nr-1);
% % us_n=cs_n.*r_n';
% %=========Cathode dr
dr_p=param.Rs_p/(param.Nr-1);
% r_p=linspace(dr_p,param.Rs_p,param.Nr-1);
% us_p=cs_p.*r_p';
ap=param.a_sp;
an=param.a_sn;
F=param.F;
Lp=param.Lp;
Ln=param.Ln;
jp=-input_crt/(F*ap*Lp*param.A);
jn=input_crt/(F*an*Ln*param.A);
% jp=-input_crt/(F*ap*Lp);
% jn=input_crt/(F*an*Ln);


r_n=linspace(dr_n/2,param.Rs_n-dr_n/2,param.Nr-1);
r_nplus=r_n+dr_n/2;
r_nminus=r_n-dr_n/2;

%=========Cathode dr

r_p=linspace(dr_p/2,param.Rs_p-dr_p/2,param.Nr-1);
r_pplus=r_p+dr_p/2;
r_pminus=r_p-dr_p/2;
usurn=(us_n(end)-jn*param.Rs_n*dr_n/(2*param.Dsn_ref))/(1-dr_n/param.Rs_n/2);
% dudt_un(1)=(us_n(2)/(dr_n)^2-3*us_n(1)/(dr_n)^2)*param.Dsn_ref;
% dudt_un(end)=((usurn-us_n(end))/(dr_n/2)-(us_n(end)-us_n(end-1))/(dr_n))*param.Dsn_ref/dr_n;

usurp=(us_p(end)-jp*param.Rs_p*dr_p/(2*param.Dsp_ref))/(1-dr_p/param.Rs_p/2);
% dudt_up(1)=(us_p(2)/(dr_p)^2-3*us_p(1)/(dr_p)^2)*param.Dsp_ref;
% dudt_up(end)=((usurp-us_p(end))/(dr_p/2)-(us_p(end)-us_p(end-1))/(dr_p))*param.Dsp_ref/dr_p;
%====Anode ode
for kk=1:1:size(us_n,1)
    if kk==1
%   dudt_un(kk)=(us_n(2)/(dr_n)^2-3*us_n(1)/(dr_n)^2)*param.Dsn_ref;
dudt_un(kk)=((us_n(2)-us_n(1))/dr_n*4*pi*r_n(kk)^2-(us_n(1))/(dr_n/2)*4*pi*r_n(kk)^2)*param.Dsn_ref/(4*pi*(r_nplus(kk)^3-r_nminus(kk)^3)/3);

    elseif kk~=size(us_n,1)&&kk~=1
%       dudt_un(kk)=(us_n(kk+1)/(dr_n)^2-2*us_n(kk)/(dr_n)^2+us_n(kk-1)/(dr_n)^2)*param.Dsn_ref;
dudt_un(kk)=((us_n(kk+1)-us_n(kk))/dr_n*4*pi*(r_n(kk)^2)-...
    (us_n(kk)-us_n(kk-1))/dr_n*4*pi*(r_n(kk)^2))*param.Dsn_ref/(4*pi*(r_nplus(kk)^3-r_nminus(kk)^3)/3);
    else
%        dudt_un(kk)=((usurn-us_n(end))/(dr_n/2)-(us_n(end)-us_n(end-1))/(dr_n))*param.Dsn_ref/dr_n;
dudt_un(kk)=((usurn-us_n(end))/(dr_n/2)*4*pi*(r_n(kk)^2)-(us_n(end)-us_n(end-1))/(dr_n)*4*pi*(r_n(kk)^2))*param.Dsn_ref/(4*pi*(r_nplus(kk)^3-r_nminus(kk)^3)/3);
    end    
end

%====Cathode ode
for kk=1:1:size(us_p,1)
    if kk==1
%   dudt_up(kk)=(us_p(2)/(dr_p)^2-3*us_p(1)/(dr_p)^2)*param.Dsp_ref;
dudt_up(kk)=((us_p(2)-us_p(1))/dr_p*4*pi*r_p(kk)^2-(us_p(1))/(dr_p/2)*4*pi*r_p(kk)^2)*param.Dsp_ref/(4*pi*(r_pplus(kk)^3-r_pminus(kk)^3)/3);

    elseif kk~=size(us_p,1)&&kk~=1
%      dudt_up(kk)=(us_p(kk+1)/(dr_p)^2-2*us_p(kk)/(dr_p)^2+us_p(kk-1)/(dr_p)^2)*param.Dsp_ref;
dudt_up(kk)=((us_p(kk+1)-us_p(kk))/dr_p*4*pi*(r_p(kk)^2)-...
    (us_p(kk)-us_p(kk-1))/dr_p*4*pi*(r_p(kk)^2))*param.Dsp_ref/(4*pi*(r_pplus(kk)^3-r_pminus(kk)^3)/3);
    else
%        dudt_up(kk)=((usurp-us_p(end))/(dr_p/2)-(us_p(end)-us_p(end-1))/(dr_p))*param.Dsp_ref/dr_p;
dudt_up(kk)=((usurp-us_p(end))/(dr_p/2)*4*pi*(r_p(kk)^2)-(us_p(end)-us_p(end-1))/(dr_p)*4*pi*(r_p(kk)^2))*param.Dsp_ref/(4*pi*(r_pplus(kk)^3-r_pminus(kk)^3)/3);    
    end    
end
dudt=[dudt_un dudt_up]';

end

