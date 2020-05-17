clc
clear all
close all
b=150;
%Достали информацию по С

fid = fopen('C_0.4_matl.txt', 'r'); 
C(:,1)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_0.6_matl.txt', 'r'); 
C(:,2)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_0.8_matl.txt', 'r'); 
C(:,3)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.0_matl.txt', 'r'); 
C(:,4)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.2_matl.txt', 'r'); 
C(:,5)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.4_matl.txt', 'r'); 
C(:,6)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.6_matl.txt', 'r'); 
C(:,7)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_1.8_matl.txt', 'r'); 
C(:,8)= fscanf(fid, '%f'); 
fclose(fid);

fid = fopen('C_2.0_matl.txt', 'r'); 
C(:,9)= fscanf(fid, '%f'); 
fclose(fid);

% Интерполируем

a=(0.05:0.05:0.95);
for i=1:9
p4=plot(a,C(:,i),'.');
hold on
end

ai=(0:0.001:1);

for i=1:9
    yi(i,:)=spline(a,C(:,i),ai);
    plot(ai,yi(i,:))
    dyi(i,:)=diff(yi(i,:));
end
title('C(a/W)');
xlabel('a/W');
ylabel('C, [1/sqrt(mm)]');

for j=1:9
p5=plot(a(6),C(6,j),'s','MarkerEdgeColor','b','Markerfacecolor',[.1*j .1*j .1*j]);
end

for i=1:(length(ai)-1)
ai_1(i)=ai(i);
end

for j=1:9
    k=0;
for i=1:length(ai_1)
    if (k==0)&&(dyi(j,i)<0)
        k=1;
        yv(j)=yi(j,i);
        av(j)=ai_1(i);    
    end
    if (k==1)&&(dyi(j,i)>0)
        k=0;
        yn(j)=yi(j,i);
        an(j)=ai_1(i);
    end
end
end


p1=plot(av,yv,'o','Markeredgecolor','k','Markerfacecolor','r');
p2=plot(an,yn,'o','Markeredgecolor','k','Markerfacecolor','g');



for j=1:9
    k=0;
for i=1:length(ai_1)
    if (k==0)&&(yi(j,i)>0)
        k=1;
        yco(j)=yi(j,i);
        aco(j)=ai_1(i);    
    end
end
end

p3=plot(aco,yco,'o','Markeredgecolor','k','Markerfacecolor','b');
grid on

legend([p4 p3 p1 p2],'Полученные данные C из ANSYS','C=0','1-я точка экстремума С','2-я точка экстремума С')

figure

for i=1:length(yv)
    hvn(i)=0.1*(i+1);
end

for i=1:length(yco)
    hco(i)=0.1*(i+1);
end

for i=1:length(ai)
b1(i)=1.25341*ai(i)^3-3.35759*ai(i)^2+2.10440*ai(i)+0.00152;
b2(i)=0.2240*ai(i)^3-2.56509*ai(i)^2+2.20573*ai(i)+0.16180;
end

hold on
plot(ai,b1,'k','linewidth',2)
plot(ai,b2,'k','linewidth',2)

title('Поле стабильности роста трещины. (Легенду см. при увеличении a/W)');
xlabel('a/W');
ylabel('h/W');

grid on

for i=1:length(hvn)
    aa1=aco(i):0.025:av(i);
    aa2=av(i):0.025:an(i);
    aa3=an(i):0.025:1;
    plot(aa1,hvn(i),'r.')
    plot(aa2,hvn(i),'g.')
    plot(aa3,hvn(i),'r.')
end
 p2=plot(aa1(1),hvn(length(hvn)),'r.');
 p3=plot(aa2(1),hvn(length(hvn)),'g.');

for i=(length(hvn)+1):length(hco)
    aa4=aco(i):0.025:1;
    plot(aa4,hco(i),'r.')
end

for i=1:length(hco)
    aa5=0:0.025:aco(i);
    plot(aa5,hco(i),'b.')
end
p1=plot(aa5(1),hco(length(hco)),'b.');

plot(aco,hco,'b','linewidth',1)
plot(av,hvn,'r','linewidth',1)
plot(an,hvn,'g','linewidth',1)
plot(aco,hco,'o','Markeredgecolor','k','Markerfacecolor','b')
plot(av,hvn,'o','Markeredgecolor','k','Markerfacecolor','r')
plot(an,hvn,'o','Markeredgecolor','k','Markerfacecolor','g')

for j=1:9
    t=0;
    for i=1:length(yi(j,:))
if yi(j,i)>0.001
     
     razai=1-ai(i);
     r_const=0.0037062591178230937623/(yi(j,i)^2)/b;
     if (t==0)&&(r_const<razai)
      ar_(j)=ai(i);
      t=1;
     end
     
end
    end
end

p6=plot(ar_,hco,'o','Markeredgecolor','y');

legend([p1 p2 p3 p6],'Предельное значение rc равно бесконечности','Предельное значение rc уменьшается','Предельное значение rc увеличивается','Предельное значение rc/W равно 1-(a/W)')

figure
hold on

for j=1:9
    k(j)=0;
    for i=1:length(yi(j,:))
if yi(j,i)>0.025
     k(j)=k(j)+1;
     ar(j,k(j))=ai(i);
     r(j,k(j))=0.0037062591178230937623/(yi(j,i)^2)/b;
end
    end
end

title('Предельное значение rc/W (a/W)');
xlabel('a/W');
ylabel('rc/W');


for j=1:9
    arr=[];
    rr=[];
    arro=[];
    rro=[];
    for i=1:k(j)
        arr(i)=ar(j,i);
        rr(i)=r(j,i);
    end
   plot(arr,rr,'b')
   
    arro=(arr(1):0.1:1);
    rro=spline(arr,rr,arro);
    pp(j)=plot(arro(2),rro(2),'s','MarkerEdgeColor','b','Markerfacecolor',[.1*j .1*j .1*j]); 
    end

legend(pp,'h/W=0.2','h/W=0.3','h/W=0.4','h/W=0.5','h/W=0.6','h/W=0.7','h/W=0.8','h/W=0.9','h/W=1.0',0)

for i=1:length(yv)
    rv(i)=0.0037062591178230937623/(yv(i)^2)/b;
    rn(i)=0.0037062591178230937623/(yn(i)^2)/b;
end

plot(av,rv,'o','Markeredgecolor','k','Markerfacecolor','r')
plot(an,rn,'o','Markeredgecolor','k','Markerfacecolor','g')



% Мы посмотрели, как характерно влияет С на rc, а теперь глянем случаи,
% когда предельное значение rc будет меньше оставшейся длины материала на
% продолжении трещины.



