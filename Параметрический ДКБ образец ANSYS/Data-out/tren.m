clc
clear all
close all


fid = fopen('coeff.txt', 'r'); 
a= fscanf(fid, '%f'); 
fclose(fid);

E=3000;
mu=0.33;

j=1;

syms r t
sx=0;
sy=0;
sxy=0;

for n=1:length(a)

sx=sx+n/2*r^(n/2-1)*a(n)*((2+n/2+(-1)^n)*cos((n/2-1)*t)-(n/2-1)*cos((n/2-3)*t));
sy=sy+n/2*r^(n/2-1)*a(n)*((2-n/2-(-1)^n)*cos((n/2-1)*t)+(n/2-1)*cos((n/2-3)*t));
sxy=sxy+n/2*r^(n/2-1)*a(n)*((n/2-1)*sin((n/2-3)*t)-(n/2+(-1)^n)*sin((n/2-1)*t));
end

U=1/(2*E)*(sx^2+sy^2-2*mu*sx*sy+2*(1+mu)*sxy^2);

dUt=diff(U, t);
d2Ut=diff(dUt, t);
dUr=diff(U, r);
d2Ur=diff(dUr, r);




r=0.5*j;
dUt_=eval(dUt);
t_mas=solve(dUt_==0,t);

for k=1:length(t_mas)
if (isreal(eval(t_mas(k)))==1)
    t=eval(t_mas(k));
    if (eval(d2Ut)>0)
    tc=t
    
    
 
    end
end
end

