clear;clc;
close all;

fc=430e6;       % 载频
fs=0.5e6;       % 

addpath("..\..\相位校正\")

% % 参数化1
% find_y=load("对齐方法\参数对齐\二阶\find_y.mat").find_y;
% fit_y=load("对齐方法\参数对齐\二阶\fit_y.mat").fit_y;
% y=fit_y+find_y;


% % 最小熵
% find_y=load("对齐方法\参数对齐\二阶\find_y.mat").find_y;
% fit_y=load("对齐方法\参数对齐\二阶\fit_y.mat").fit_y;
% y=fit_y+find_y;

% % 相邻相关
% find_y=load("对齐方法\非参数对齐\加窗相邻相关对齐\find_y.mat").find_y;
% fit_y=load("对齐方法\非参数对齐\加窗相邻相关对齐\fit_y.mat").fit_y;
% y=fit_y+find_y;

% 相邻相关
find_y=load("对齐方法\非参数对齐\单信号最小熵对齐\find_y.mat").find_y;
fit_y=load("对齐方法\非参数对齐\单信号最小熵对齐\fit_y.mat").fit_y;
y=fit_y+find_y;

% 直接成像



load('data_origin.mat')
data=data(:,3001:end);
[M,N]=size(data);

figure;
subplot(121)
imagesc(abs(data'));title("Lunar dara")
xticks((1:5)*1000)
xticklabels((1:5)*60)
xlabel("Time (sec)")
yticklabels(yticks/1000*300)
ylabel("Range")
set(gca,"YDir","normal","FontName","Times New Roman","FontSize",24);

% figure;
subplot(122)
imagesc(abs(data(:,351:750)'));title("Lunar dara")
xticks((1:5)*1000)
xticklabels((1:5)*60)
xlabel("Time (sec)")
yticklabels((yticks+350)/1000*300)
ylabel("Range")
set(gca,"YDir","normal","FontName","Times New Roman","FontSize",24);


% refe_scale=[0:N/2-1,-N/2:-1]/N*fs;
% 
% data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*y/3e8),[],2);
% data=data(:,501:end-500);
% % figure;plot(fit_y);title('粗对齐曲线')
% % figure;plot(y-fit_y);title("补偿曲线")
% figure;imagesc(abs(data'));set(gca,"YDir","normal");title("距离时间域")
% figure;imagesc(abs(fftshift(fft(data',[],2),2)));set(gca,"YDir","normal");title("距离多普勒域")
% [data,entro_his,~]=Normal_pga(data,10);
% 
% 
% data_entropy=-sum(sum(abs(data)./sum(abs(data(:))).*log2(abs(data)./sum(abs(data(:))))))
% 
data_mean=mean(abs(data(:)));
tmp_data=(abs(data)-data_mean).^2;

data_contrast=(mean(tmp_data(:)).^(1/2))/data_mean
% 
% figure;imagesc(abs(data'));set(gca,"YDir","normal");title("相位校正结果")
