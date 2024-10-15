clear;clc;
close all;
addpath("..\..\相位校正\")

fc=430e6;       % 载频
fs=0.5e6;       % 

data=load("data_origin.mat").data;

data=data(:,3001:end);
[M,N]=size(data);

x=1:M;
% y=load("对齐方法\参数对齐\二阶_\二阶\y.mat").y;
refe_scale=[0:N/2-1,-N/2:-1]/N*fs;


fit_y=load("对齐方法\非参数对齐\加窗相邻相关对齐\fit_y.mat").fit_y;

para_y=load("对齐方法\参数对齐\二阶\find_y.mat").find_y;
rela_y=load("对齐方法\非参数对齐\加窗相邻相关对齐\find_y.mat").find_y;
entr_y=load("对齐方法\非参数对齐\单信号最小熵对齐\find_y.mat").find_y;
pso_y=load("直接成像\二阶\find_y.mat").find_y';
ephe_y=load("y_m.mat").y_m-fit_y;
ephe_y=ephe_y-min(ephe_y);


delta_r1=load("对齐方法\非参数对齐\加窗相邻相关对齐\delta_r.mat").delta_r;
delta_r2=load("对齐方法\非参数对齐\单信号最小熵对齐\delta_r.mat").delta_r;

% para_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*(fit_y+para_y)/3e8),[],2);
% rela_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*(fit_y+rela_y)/3e8),[],2);
% entr_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*(fit_y+entr_y)/3e8),[],2);
% pso_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*(fit_y+pso_y)/3e8),[],2);
% ephe_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*(fit_y+ephe_y)/3e8),[],2);


% score_his=load("对齐方法\参数对齐\二阶\score_his.mat").score_his;
% pso_data=mean(abs(rela_data),1);
% tmp_entropy2=-sum((pso_data/sum(pso_data)).*log2(pso_data/sum(pso_data)));
% 
% pso_data=mean(abs(entr_data),1);
% tmp_entropy3=-sum((pso_data/sum(pso_data)).*log2(pso_data/sum(pso_data)));
% 
% score1=ones([1,70])*tmp_entropy2;
% score2=ones([1,70])*tmp_entropy3;
% 
% figure;plot(1:70,score_his,1:70,score1,1:70,score2);
% xlabel("Times")
% ylabel("Entropy")
% set(gca,fontsize=28,fontname="Times New Roman")

% y1=y-fit_y;

figure;
colororder(["black","blue"])
yyaxis left;
plot(x,para_y/1000,x,rela_y/1000,x,entr_y/1000)
grid on
ylabel("Range (km)")

yyaxis right;
plot(x,pso_y/1000,x,ephe_y/1000)
grid on
ylabel("Range (km)")

xlabel("Time (sec)")
legend("F1","F2","F3","F4","F5")
xlim([0,5121])
xticks([1024,2048,3072,4096,5120])
xticklabels([60,120,180,240,300])
set(gca,fontsize=28,fontname="Times New Roman")


figure;imagesc(abs(para_data(:,1:end-1000)'));
xlabel("Time (sec)")
ylabel("Range (km)")
xticks([1024,2048,3072,4096,5120])
xticklabels([60,120,180,240,300])
ylim([1290,1390])
yticks([1310,1330,1350,1370,1390])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);

figure;imagesc(abs(rela_data(:,1:end-1000)'));
xlabel("Time (sec)")
ylabel("Range (km)")
xticks([1024,2048,3072,4096,5120])
xticklabels([60,120,180,240,300])
ylim([1290,1390])
yticks([1310,1330,1350,1370,1390])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);

figure;imagesc(abs(entr_data(:,1:end-1000)'));
xlabel("Time (sec)")
ylabel("Range (km)")
xticks([1024,2048,3072,4096,5120])
xticklabels([60,120,180,240,300])
ylim([1290,1390])
yticks([1310,1330,1350,1370,1390])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);


figure;imagesc(abs(fftshift(fft(para_data(:,1:end-1000)',[],2),2)));
xlabel("Range (km)")
ylabel("Range (km)")
xticks([1280,2560,3840,5120])
xticklabels([-1228.8,0,1228.8,2457.6])
yticks([1000,2000,3000,4000,5000,6000])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);

figure;imagesc(abs(fftshift(fft(rela_data(:,1:end-1000)',[],2),2)));
xlabel("Range (km)")
ylabel("Range (km)")
xticks([1280,2560,3840,5120])
xticklabels([-1228.8,0,1228.8,2457.6])
yticks([1000,2000,3000,4000,5000,6000])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);

figure;imagesc(abs(fftshift(fft(entr_data(:,1:end-1000)',[],2),2)));
xlabel("Range (km)")
ylabel("Range (km)")
xticks([1280,2560,3840,5120])
xticklabels([-1228.8,0,1228.8,2457.6])
yticks([1000,2000,3000,4000,5000,6000])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);




[para_data,~,~]=Normal_pga(para_data,15);
[rela_data,~,~]=Normal_pga(rela_data,15);
[entr_data,~,~]=Normal_pga(entr_data,15);

figure;imagesc(abs(para_data(:,1:end-1000)'));
xlabel("Range (km)")
ylabel("Range (km)")
xticks([1280,2560,3840,5120])
xticklabels([-1228.8,0,1228.8,2457.6])
yticks([1000,2000,3000,4000,5000,6000])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);clim([0,(max(abs(para_data(:)))/10)])

figure;imagesc(abs(rela_data(:,1:end-1000)'));
xlabel("Range (km)")
ylabel("Range (km)")
xticks([1280,2560,3840,5120])
xticklabels([-1228.8,0,1228.8,2457.6])
yticks([1000,2000,3000,4000,5000,6000])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);clim([0,(max(abs(rela_data(:)))/10)])

figure;imagesc(abs(entr_data(:,1:end-1000)'));
xlabel("Range (km)")
ylabel("Range (km)")
xticks([1280,2560,3840,5120])
xticklabels([-1228.8,0,1228.8,2457.6])
yticks([1000,2000,3000,4000,5000,6000])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);clim([0,(max(abs(entr_data(:)))/10)])

figure;imagesc(abs(para_data(2930:3360,1380:1810)'));
xlabel("Range (km)")
ylabel("Range (km)")
% xticks([1280,2560,3840,5120])
% xticklabels([-1228.8,0,1228.8,2457.6])
% yticks([1000,2000,3000,4000,5000,6000])
% yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);colormap("gray");clim([0,(max(abs(para_data(:)))/10)])

figure;imagesc(abs(rela_data(2370:2800,1380:1810)'));
xlabel("Range (km)")
ylabel("Range (km)")
% xticks([1280,2560,3840,5120])
% xticklabels([-1228.8,0,1228.8,2457.6])
% yticks([1000,2000,3000,4000,5000,6000])
% yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);colormap("gray");clim([0,(max(abs(rela_data(:)))/10)])

figure;imagesc(abs(entr_data(2880:3310,1380:1810)'));
xlabel("Range (km)")
ylabel("Range (km)")
% xticks([1280,2560,3840,5120])
% xticklabels([-1228.8,0,1228.8,2457.6])
% yticks([1000,2000,3000,4000,5000,6000])
% yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);colormap("gray");clim([0,(max(abs(entr_data(:)))/10)])