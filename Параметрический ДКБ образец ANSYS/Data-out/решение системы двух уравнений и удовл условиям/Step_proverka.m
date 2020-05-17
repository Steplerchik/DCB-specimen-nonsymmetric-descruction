clc
clear all
close all

t=10;
Rvnesh=t/2
Rvnutr=t/4

[r, t] = meshgrid(Rvnutr:((Rvnesh-Rvnutr)/9):Rvnesh, -180:2:180); 
X=r.*cos(t./180*pi);
Y=r.*sin(t./180*pi);

a=size(r,1)
b=size(r,2)
C=0.208;
K2_K1=-0.0426;
mu=0.33;
K1=1;

sx=r;
sy=r;
for i=1:a
        for j=1:b
            sx(i,j)=0;
            sy(i,j)=0;
        end
end

for i=1:a
        for j=1:b
            f1x=cos(t(i,j)/180*pi/2)*(1-sin(t(i,j)/180*pi/2)*sin(3/2*t(i,j)/180*pi));
            f2x=-sin(t(i,j)/180*pi/2)*(2+cos(t(i,j)/180*pi/2)*cos(3/2*t(i,j)/180*pi));
            sx(i,j)=K1*(f1x/sqrt(2*pi*r(i,j))+C+K2_K1*f2x/sqrt(2*pi*r(i,j)));
            f1y=cos(t(i,j)/180*pi/2)*(1+sin(t(i,j)/180*pi/2)*sin(3/2*t(i,j)/180*pi));
            f2y=sin(t(i,j)/180*pi/2)*cos(t(i,j)/180*pi/2)*cos(3/2*t(i,j)/180*pi);
            sy(i,j)=K1*(f1y/sqrt(2*pi*r(i,j))+K2_K1*f2y/sqrt(2*pi*r(i,j)));
        end
end
   

figure
contour(X,Y,sx)
grid on
title('Изолинии поля напряжений sx');
xlabel('x');
ylabel('y');

figure
contour(X,Y,sy)
grid on
title('Изолинии поля напряжений sy');
xlabel('x');
ylabel('y');