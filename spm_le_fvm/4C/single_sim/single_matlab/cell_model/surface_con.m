function csi_surf = surface_con(param,cs_i,flag)
%============The flag indicate cathode (0) or anode(1)


if flag==0 %===========Cathode
cs_p=cs_i;    
x1p=(cs_p(end,:)+cs_p(end-1,:))/2;
xmp=cs_p(end,:);
x2p=2*xmp-x1p;
csi_surf=x2p;

elseif flag==1  %===========Anode
cs_n=cs_i;   
x1n=(cs_n(end,:)+cs_n(end-1,:))/2;
xmn=cs_n(end,:);
x2n=2*xmn-x1n;
csi_surf=x2n; 
end

end

