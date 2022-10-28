function Dsolv0 = Solv_Diff_update(input_crt)
input_crt = abs(input_crt);

p1 =   1.634e-21;
p2 =  -2.171e-20; 
p3 =    1.76e-19; 
p4 =   9.041e-20;

p = [p1 p2 p3 p4];

Dsolv0 = polyval(p,input_crt);
end