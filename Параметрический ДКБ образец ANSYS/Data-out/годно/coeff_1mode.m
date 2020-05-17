clc
clear all
close all


fid = fopen('coeff_1mode.txt', 'r'); 
a= fscanf(fid, '%f'); 
fclose(fid);

E=3000;
mu=0.3;

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
d2Ut=diff(dUt,t);

s=r
s=s+r
kk=0;
for j=1:100
r=0.1*j;
t=0:pi/18000:pi/2;
U_=eval(U);
dUt_=eval(dUt);
d2Ut_=eval(d2Ut);
check=0;
for k=1:(length(t)-1)
    if (dUt_(k)<=0)&&(dUt_(k+1)>0)&&(check==0)
        kk=kk+1;
        tc(kk)=t(k)/pi*180;
        rc(kk)=r;
        check=1;
    end
end
end
tc
plot(t,dUt_)
figure
plot(t,U_)
figure
plot(rc,tc)




