clear;clc;
close all;

fc=430e6;       % 载频
fs=0.5e6;  % 采样率 20e6

load("y.mat")
y=y-y(1);

data=load('..\..\data_origin.mat').data;
data=data(:,3001:end-1000);
[M,N]=size(data);
refe_scale=[0:N/2-1,-N/2:-1]/N*fs;

single_x=zeros(1,M);
for each_m=1:M
    single_x(each_m)=find(abs(data(each_m,:))>3,1);
end

x=1:M;
fit_p=polyfit(x,single_x,2);
fit_p(end)=0;
fit_y=polyval(fit_p,x)*(3e8/fs);

Y=fit_y'-y;

data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*y/3e8),[],2);

figure;plot(fit_y/1000);
% title("Y");
xlabel("Time (sec)")
ylabel("Range (km)")
xticks([1024,2048,3072,4096,5120])
xticklabels([60,120,180,240,300])
xlim([0,5121])
ylim([-90,0])
set(gca,fontname="Times New Roman",fontsize=28);

figure;plot(Y/1000)
% title("y");
xlabel("Time (sec)")
ylabel("Range (km)")
xticks([1024,2048,3072,4096,5120])
xticklabels([60,120,180,240,300])
xlim([0,5121])
ylim([0,6])
set(gca,fontname="Times New Roman",fontsize=28);

figure;imagesc(abs(data'));
% title("Range-Time")
xlabel("Time (sec)")
ylabel("Range (km)")
xticks([1024,2048,3072,4096,5120])
xticklabels([60,120,180,240,300])
yticks([1000,2000,3000,4000,5000,6000])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);

figure;imagesc(abs(fftshift(fft(data',[],2),2)));
% title("Range-Doppler")
xlabel("Range (km)")
ylabel("Range (km)")
xticks([1280,2560,3840,5120])
xticklabels([-1228.8,0,1228.8,2457.6])
yticks([1000,2000,3000,4000,5000,6000])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);colormap("gray");


figure;imagesc(abs(data'));
% title("Range-Time")
xlabel("Time (sec)")
ylabel("Range (km)")
ylim([500,700])
xticks([1024,2048,3072,4096,5120])
xticklabels([60,120,180,240,300])
yticks([550,600,650,700])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);


figure;imagesc(abs(fftshift(fft(data',[],2),2)));
% title("Range-Doppler")
xlabel("Range (km)")
ylabel("Range (km)")
xlim([0,1000])
ylim([2500,5500])
xticks([200,400,600,800,1000])
xticklabels([-2265.6,-2073.6,-1881.6,-1689.6,-1497.6])
yticks([3000,3500,4000,4500,5000,5500])
yticklabels(yticks*600/1000)
set(gca,"YDir","normal",fontname="Times New Roman",fontsize=28);colormap("gray");










