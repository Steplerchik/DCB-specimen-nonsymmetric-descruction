clc
clear all
close all


fid = fopen('coeff_04.txt', 'r'); 
a_mas(:,1)= fscanf(fid, '%f'); 
fclose(fid);
fid = fopen('coeff_06.txt', 'r'); 
a_mas(:,2)= fscanf(fid, '%f'); 
fclose(fid);
fid = fopen('coeff_08.txt', 'r'); 
a_mas(:,3)= fscanf(fid, '%f'); 
fclose(fid);
fid = fopen('coeff_10.txt', 'r'); 
a_mas(:,4)= fscanf(fid, '%f'); 
fclose(fid);
fid = fopen('coeff_12.txt', 'r'); 
a_mas(:,5)= fscanf(fid, '%f'); 
fclose(fid);
fid = fopen('coeff_14.txt', 'r'); 
a_mas(:,6)= fscanf(fid, '%f'); 
fclose(fid);





for gg=1:6
gstr=0;
gstb=1;
    for g=1:length(a_mas(:,gg))
    gstr=gstr+1;
    if (gstr==11)
        gstr=1;
        gstb=gstb+1;
    end
    a_(gstr,gstb,gg)=a_mas(g,gg);
    end
end
    
E=3000;
mu=0.3;

gg=2;
g=8;
    
a=a_(:,g,gg);
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

dUr=diff(U, r);
d2Ur=diff(dUr, r);

s=r
s=s+r
kk=0;
sto=100;
if (g==1)
    sto=50;
end
if (g==19)
    sto=50;
end

R=0.1*5:0.1:0.1*sto;
T=0:pi/1800:pi/2;
[r,t]=meshgrid(R,T);
U_=eval(U);

x=r.*cos(t);
y=r.*sin(t);

mesh(x,y,U_)
xlabel('x');
ylabel('y');

% max=U_(1);
% for j=1:length(U_)
%     if (U_(j)>=max)
%         rc_(gg,g)=r(j);
%         tc_(gg,g)=t(j);
%         max=U_(j);
%     end
% end



% rc_
% tc_








