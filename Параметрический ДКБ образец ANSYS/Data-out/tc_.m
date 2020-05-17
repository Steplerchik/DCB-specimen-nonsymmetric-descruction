clc
clear all
close all

fid = fopen('tc_.txt', 'r'); 
if fid == -1 
    error('File is not opened'); 
end 

ss= fscanf(fid, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f', [19 6]); 
fclose(fid);
ss=ss';

ai=0.05:0.05:0.95;
for i=1:length(ai)
b1(i)=1.25341*ai(i)^3-3.35759*ai(i)^2+2.10440*ai(i)+0.00152;
b2(i)=0.2240*ai(i)^3-2.56509*ai(i)^2+2.20573*ai(i)+0.16180;
%b3(i)=sqrt((3*ai(i)^3-6*ai(i)^2+3*ai(i))/(ai(i)+2));
end

hold on
plot(ai,b1,'k','linewidth',2)
plot(ai,b2,'k','linewidth',2)
%plot(ai,b3,'k--','linewidth',2)

for j=1:6
    yi(j)=0.1*(j+1);
end

for j=1:6
for i=1:19
    if (ss(j,i)==0.0000)
    plot(ai(i),yi(j),'ko','linewidth',1)
    else
    plot(ai(i),yi(j),'k.','linewidth',2)    
    end
end
end
    
xlabel('a/W');
ylabel('h/W');



