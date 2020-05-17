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
s=0;

for gg=1:6
for g=1:19
rc=[];
tc=[];
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
d2Ut=diff(dUt,t);

s=s+1

kk=0;

for j=1:11
   syms r t 
r=0.375+j*0.0375-0.0375;

dUt_=eval(dUt);
t_mas=solve(dUt_==0,t);

for k=1:length(t_mas)
if (isreal(eval(t_mas(k)))==1)
    t=eval(t_mas(k));
    if (eval(d2Ut)>0)
    kk=kk+1;
        tc(kk)=t;
        rc(kk)=r;
    end
end
end
end

if (length(tc)~=0)
max=tc(1);
for j=1:length(tc)
    if (tc(j)>=max)
        rc_(gg,g)=rc(j);
        tc_(gg,g)=tc(j);
        max=tc(j);
    end
end
else
    rc_(gg,g)=0;
    tc_(gg,g)=-1;
end

end
end

fid = fopen('tc_.txt', 'w'); 
if fid == -1 
    error('File is not opened'); 
end 
  

fprintf(fid, '%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\r\n', tc_'); 
fclose(fid);

