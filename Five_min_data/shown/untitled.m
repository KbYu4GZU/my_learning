clear;clc;
close all;

load("..\data_origin.mat")
data=data(:,3001:end);
data=data';
M=5120;

figure;
imagesc(abs(data));

xticks([1024,2048,3072,4096,5120]);
xticklabels([60,120,180,240,300]);
xlabel("Time (sec)")
yticks([1000,2200,3400,4600,5800,7000])
yticklabels([500,1100,1700,2300,2900,3500]);
ylabel("Range (km)")
title("Complete Data")
set(gca,"YDir","normal","FontName","Times New Roman","FontSize",18);

figure;
imagesc(abs(data(400:700,:)));

xticks([1024,2048,3072,4096,5120]);
xticklabels([60,120,180,240,300]);
xlabel("Time (sec)")
yticks([0,60,120,180,240,300]);
yticklabels([200,230,260,290,320,350]);
ylabel("Range (km)")
title("Enlarge Section")
set(gca,"YDir","normal","FontName","Times New Roman","FontSize",18);

single_x=zeros(1,M);
for each_m=1:M
    single_x(each_m)=find(abs(data(:,each_m))>3,1);
end
figure;
plot(single_x)
axis([0,5120,400,600])
xticks([1024,2048,3072,4096,5120]);
xticklabels([60,120,180,240,300]);
xlabel("Time (sec)")
yticks([440,480,520,560,600]);
yticklabels([220,240,260,280,300]);
ylabel("Range (km)")
title("Start List")
set(gca,"FontName","Times New Roman","FontSize",18)