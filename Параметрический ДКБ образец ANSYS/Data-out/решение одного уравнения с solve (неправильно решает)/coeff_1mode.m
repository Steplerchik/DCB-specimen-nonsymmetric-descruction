clc
clear all
close all


fid = fopen('coeff.txt', 'r'); 
a= fscanf(fid, '%f'); 
fclose(fid);

E=3000;
mu=0.33;

for j=1:10

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

kk=0;



r=1.5*j;
dUt_=eval(dUt);
t_mas=solve(dUt_==0,t);

for k=1:length(t_mas)
if (isreal(eval(t_mas(k)))==1)
    t=eval(t_mas(k));
    d2Ut_=eval(d2Ut);
    if (d2Ut_>0)
    kk=kk+1
    tt(j,kk)=t/pi*180;
    rr(j,kk)=r;
    else
    kk=kk+1
    tt(j,kk)=-1;
    rr(j,kk)=-1;    
    end
end
end

end




