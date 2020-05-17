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

for gg=1:6
for g=1:19
rc=[];
tc=[];
a=a_(:,g,gg);
syms r t
sx=0;
sy=0;
sxy=0;

for n=1:2

sx=sx+n/2*r^(n/2-1)*a(n)*((2+n/2+(-1)^n)*cos((n/2-1)*t)-(n/2-1)*cos((n/2-3)*t));
sy=sy+n/2*r^(n/2-1)*a(n)*((2-n/2-(-1)^n)*cos((n/2-1)*t)+(n/2-1)*cos((n/2-3)*t));
sxy=sxy+n/2*r^(n/2-1)*a(n)*((n/2-1)*sin((n/2-3)*t)-(n/2+(-1)^n)*sin((n/2-1)*t));
end

U=1/(2*E)*(sx^2+sy^2-2*mu*sx*sy+2*(1+mu)*sxy^2);
dUt=diff(U, t);
d2Ut=diff(dUt, t);

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
for j=10:sto*2
r=0.05*j
t=0:pi/180000:pi/2;
U_=eval(U);
%d2Ut_=eval(d2Ut);
dUt_=eval(dUt);
check=0;
for k=1:(length(t)-1)
    if (dUt_(k)<=0)&&(dUt_(k+1)>=0)&&(check==0)
% for k=1:(length(t))
%     if (dUt_(k)==0)&&(d2Ut_(k)>0)&&(check==0)
        kk=kk+1;
        tc(kk)=t(k)/pi*180;
        rc(kk)=r;
        check=1;
    end
end
end
tc;
%plot(t,dUt_)
%figure
%plot(t,U_)
%figure
%plot(rc,tc)
r=[];
t=[];
for kk=1:length(rc)
    r(kk)=rc(kk);
    t(kk)=tc(kk);
end

%r=rc(1):(rc(length(rc))-rc(1))/200:rc(length(rc));
%t=spline(rc,tc,r);

%dUr_=eval(dUr);
%d2Ur_=eval(d2Ur);
U_=eval(U);


%figure
%plot(r,dUr_)
%figure
%plot(r,d2Ur_)
%figure
%plot(r,U_)

max=U_(1);
for j=1:length(U_)
    if (U_(j)>=max)
        rc_(gg,g)=r(j);
        tc_(gg,g)=t(j);
        max=U_(j);
    end
end

end
end
rc_
tc_

fid = fopen('tc_.txt', 'w'); 
if fid == -1 
    error('File is not opened'); 
end 
  

fprintf(fid, '%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\r\n', tc_'); 
fclose(fid);






