clc
clear all
close all


fid = fopen('C_0.4_matl.txt', 'r'); 
C_mas(:,1)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_0.6_matl.txt', 'r'); 
C_mas(:,2)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_0.8_matl.txt', 'r'); 
C_mas(:,3)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.0_matl.txt', 'r'); 
C_mas(:,4)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.2_matl.txt', 'r'); 
C_mas(:,5)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.4_matl.txt', 'r'); 
C_mas(:,6)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.6_matl.txt', 'r'); 
C_mas(:,7)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.8_matl.txt', 'r'); 
C_mas(:,8)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_2.0_matl.txt', 'r'); 
C_mas(:,9)= fscanf(fid, '%f'); 
fclose(fid);
schetchik=0;
for p=1:9


for j=1:length(C_mas(:,1))
    
C=C_mas(j,p);
E=1;
mu=0.33;
K1=1;

syms r t

f1x=cos(t/2)*(1-sin(t/2)*sin(3/2*t));
f1y=cos(t/2)*(1+sin(t/2)*sin(3/2*t));
f1xy=cos(t/2)*sin(t/2)*cos(3/2*t);

sx=K1*(f1x/sqrt(2*pi*r)+C);
sy=K1*f1y/sqrt(2*pi*r);
sxy=K1*f1xy/sqrt(2*pi*r);

U=1/(2*E)*(sx^2+sy^2-2*mu*sx*sy+2*(1+mu)*sxy^2);

dU=diff(U, t);
d2U=diff(dU, t);

r_sim=solve(dU==0,r);

dr=diff(r_sim, t);
d2r=diff(dr, t);

t_mas=solve(dr==0,t);

kk=0; k=0;
for n=1:length(t_mas)
    if (isreal(eval(t_mas(n)))==1)
        k=k+1;
        t=eval(t_mas(n));
        r=eval(r_sim);
        d2U_=eval(d2U);
        d2r_=eval(d2r);
        if (d2U_>0)&&(d2r_>0)
            kk=kk+1;
            tt(j,p,kk)=t/pi*180;
            rr(j,p,kk)=r;
        end
    end
end

kk
schetchik=schetchik+1
end
end
tt
rr
    



