clear;clc;close all;
addpath("..\..\相位校正\")

fc=430e6;       % 载频
fs=0.5e6;       % 

data=load("data_origin.mat").data;

data=data(:,3001:end);
[M,N]=size(data);

x=1:5120;
refe_scale=[0:N/2-1,-N/2:-1]/N*fs;

f_size=18;



% % 比较图像
% fit_y=load("对齐方法\非参数对齐\加窗相邻相关对齐\fit_y.mat").fit_y;
% y1=load("对齐方法\非参数对齐\加窗相邻相关对齐\find_y.mat").find_y;
% y2=load('对齐方法\非参数对齐\单信号最小熵对齐\find_y.mat').find_y;
% y3=load("对齐方法\参数对齐\二阶\final_y.mat").final_y;
% data1=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*(fit_y+y1)/3e8),[],2);
% data2=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*(fit_y+y2)/3e8),[],2);
% data3=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*y3/3e8),[],2);
% data1=data1(:,501:end-500);
% data2=data2(:,501:end-500);
% data3=data3(:,501:end-500);
% [data1,~,~]=Normal_pga(data1,15);
% [data2,~,~]=Normal_pga(data2,15);
% [data3,~,~]=Normal_pga(data3,15);
% data1=abs(flip(data1,1));
% data2=abs(flip(data2,1));
% data3=abs(flip(data3,1));
% figure
% subplot(321)
% imagesc(data3(2760:2920,950:1310)')
% xticks((0:4)*52)
% xticklabels((xticks+2760-2560)*(8.33/2560))
% yticks((0:4)*90)
% yticklabels((yticks+950)/1000*300)
% xlabel("Frequency (Hz)")
% ylabel("Range (km)")
% clim([0,max(max(data3(2760:2920,950:1310)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(323)
% imagesc(data1(3319:3479,950:1310)')
% xticks((0:4)*52)
% xticklabels((xticks+3319-2560)*(8.33/2560))
% yticks((0:4)*90)
% yticklabels((yticks+950)/1000*300)
% xlabel("Frequency (Hz)")
% ylabel("Range (km)")
% clim([0,max(max(data1(3319:3479,950:1310)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(325)
% imagesc(data2(2809:2969,950:1310)')
% xticks((0:4)*52)
% xticklabels((xticks+2809-2560)*(8.33/2560))
% yticks((0:4)*90)
% yticklabels((yticks+950)/1000*300)
% xlabel("Frequency (Hz)")
% ylabel("Range (km)")
% clim([0,max(max(data2(2809:2969,950:1310)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% 
% subplot(322)
% imagesc(data3(1950:2050,950:1150)')
% xticks((0:4)*25)
% xticklabels((xticks+1950-2560)*(8.33/2560))
% yticks((0:4)*50)
% yticklabels((yticks+950)/1000*300)
% xlabel("Frequency (Hz)")
% ylabel("Range (km)")
% clim([0,max(max(data3(1950:2050,950:1150)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(324)
% imagesc(data1(2509:2609,950:1150)')
% xticks((0:4)*25)
% xticklabels((xticks+2509-2560)*(8.33/2560))
% yticks((0:4)*50)
% yticklabels((yticks+950)/1000*300)
% xlabel("Frequency (Hz)")
% ylabel("Range (km)")
% clim([0,max(max(data1(2509:2609,950:1150)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(326)
% imagesc(data2(1999:2099,950:1150)')
% xticks((0:4)*25)
% xticklabels((xticks+1999-2560)*(8.33/2560))
% yticks((0:4)*50)
% yticklabels((yticks+950)/1000*300)
% xlabel("Frequency (Hz)")
% ylabel("Range (km)")
% clim([0,max(max(data2(1999:2099,950:1150)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")
% 
% % subplot(337)
% % imagesc(data3(2150:2350,4650:5050)')
% % xticks()
% % xticklabels()
% % yticks()
% % yticklabels()
% % xlabel("Frequency (Hz)")
% % ylabel("Range (km)")
% % clim([0,max(max(data3(1949:2149,4650:5050)))/2])
% % set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% % 
% % 
% % 
% % subplot(338)
% % imagesc(data1(2709:2909,4650:5050)')
% % xticks()
% % xticklabels()
% % yticks()
% % yticklabels()
% % xlabel("Frequency (Hz)")
% % ylabel("Range (km)")
% % clim([0,max(max(data1(2709:2909,4650:5050)))/2])
% % set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% % 
% % 
% % 
% % subplot(339)
% % imagesc(data2(2199:2399,4650:5050)')
% % xticks()
% % xticklabels()
% % yticks()
% % yticklabels()
% % xlabel("Frequency (Hz)")
% % ylabel("Range (km)")
% % clim([0,max(max(data2(2199:2399,4650:5050)))/2])
% % set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% % colormap("gray")





% % 后面两种方法成像结果
% fit_y=load("对齐方法\非参数对齐\加窗相邻相关对齐\fit_y.mat").fit_y;
% y1=load("对齐方法\非参数对齐\加窗相邻相关对齐\find_y.mat").find_y;
% r1=load("对齐方法\非参数对齐\加窗相邻相关对齐\delta_r.mat").delta_r;
% y2=load('对齐方法\非参数对齐\单信号最小熵对齐\find_y.mat').find_y;
% r2=load("对齐方法\非参数对齐\单信号最小熵对齐\delta_r.mat").delta_r;
% 
% y1=y1+fit_y;
% y2=y2+fit_y;
% r1=r1'+fit_y;
% r2=r2'+fit_y;
% 
% figure;
% subplot(121)
% plot(x,y1,x,r1)
% xlim([0,5120])
% 
% subplot(122)
% plot(x,y2,x,r2)
% xlim([0,5120])
% 
% figure;
% subplot(131)
% plot(x,fit_y/1000,'Color',"black","LineWidth",2)
% xlim([0,5120]);
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Time (sec)","(a)"})
% ylabel(" Range (km)")
% grid on;
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(132)
% plot(x,y1,x,r1)
% xlim([0,5120]);
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Time (sec)","(b)"})
% ylabel(" Range (m)")
% grid on;
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(133)
% plot(x,y2,x,r2)
% xlim([0,5120]);
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Time (sec)","(c)"})
% ylabel(" Range (m)")
% grid on;
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% figure;
% data1=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*(fit_y+y1)/3e8),[],2);
% data2=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*(fit_y+y2)/3e8),[],2);
% data1=data1(:,501:end-500);
% data2=data2(:,501:end-500);
% % subplot(231)
% % imagesc(abs(fftshift(fft(data1',[],2),2)))
% % xticks((1:4)*(5120/4))
% % xticklabels(((1:4)-2)*(8.33/2))
% % yticks((1:4)*1500)
% % yticklabels((1:4)*1500/1000*300)
% % xlabel({"Frequency (Hz)","(c)"})
% % ylabel(" Range (km)")
% % set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% % 
% % subplot(234)
% % imagesc(abs(fftshift(fft(data2',[],2),2)));xticks((1:4)*(5120/4))
% % xticklabels(((1:4)-2)*(8.33/2))
% % yticks((1:4)*1500)
% % yticklabels((1:4)*1500/1000*300)
% % xlabel({"Frequency (Hz)","(c)"})
% % ylabel(" Range (km)")
% % set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% [data11,~,~]=Normal_pga(data1,15);
% [data22,~,~]=Normal_pga(data2,15);
% data11=abs(flip(data11,1));
% data22=abs(flip(data22,1));
% 
% subplot(221)
% imagesc(data11');
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(c)"})
% ylabel(" Range (km)")
% clim([min(data11(:)),max(data11(:))/10])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% 
% subplot(222)
% imagesc(data22');
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(c)"})
% ylabel(" Range (km)")
% clim([min(data22(:)),max(data22(:))/10])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(223)
% imagesc(data11');
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(c)"})
% ylabel(" Range (km)")
% clim([min(data11(:)),max(data11(:))/10])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(224)
% imagesc(data22');
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(c)"})
% ylabel(" Range (km)")
% clim([min(data22(:)),max(data22(:))/10])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")

% % 图8，直接成像
% f_size=18;
% y=load("直接成像\二阶\final_y.mat").final_y;
% tmp_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*y/3e8),[],2);
% tmp_data=tmp_data(:,501:end-500);
% figure;
% subplot(221)
% plot(y/1000,Color="black",LineWidth=2);
% xlim([0,5120]);
% xticks((0:5)*1000)
% xticklabels((0:5)*60)
% xlabel({"Time (sec)","(a"})
% ylabel(" Range (km)")
% grid on;
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(222)
% imagesc(abs(tmp_data(:,1:200)'));
% xticks((0:5)*1000)
% xticklabels((0:5)*60)
% yticks((0:4)*50)
% yticklabels((0:4)*50*300/1000)
% xlabel({"Time (sec)","(b)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% tmp_data=fftshift(fft(tmp_data',[],2),2);
% subplot(223)
% imagesc(abs(tmp_data));
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels(yticks/1000*300)
% xlabel({"Frequency (Hz)","(c)"})
% ylabel(" Range (km)")
% clim([min(abs(tmp_data(:))),max(abs(tmp_data(:)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% 
% subplot(224)
% imagesc(abs(tmp_data));
% xlim([4320,5120]);
% ylim([2450,4050]);
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% yticks((1:4)*400+2450)
% yticklabels(((1:4)*400+2450)/1000*300)
% xlabel({"Frequency (Hz)","(d)"})
% ylabel(" Range (km)")
% clim([min(abs(tmp_data(:))),max(abs(tmp_data(:)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);

% % 图7 方法成像
% f_size=18;
% y=load("对齐方法\参数对齐\二阶_\二阶\y.mat").y;
% tmp_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*y/3e8),[],2);
% tmp_data=tmp_data(:,501:end-500);
% 
% figure
% subplot(231)
% plot(y/1000,Color="black",LineWidth=2);
% xlim([0,5120]);
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Times (sec)","(a)"})
% ylabel(" Range (km)")
% grid on;
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(234)
% imagesc(abs(tmp_data(:,1:200)'));
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% yticks((1:4)*50)
% yticklabels((1:4)*50*300/1000)
% xlabel({"Times (sec)","(b)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(232)
% imagesc(abs(fftshift(fft(tmp_data',[],2),2)));
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(c)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")
% clim([0,max(max(abs(fftshift(fft(tmp_data',[],2),2))))/5])
% 
% subplot(235)
% imagesc(abs(fftshift(fft(tmp_data',[],2),2)));
% xlim([1020,4100])
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(d)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")
% clim([0,max(max(abs(fftshift(fft(tmp_data',[],2),2))))/5])
% 
% [fft_data,~,~]=Normal_pga(tmp_data,15);
% fft_data=flip(fft_data,1);
% 
% subplot(233)
% imagesc(abs(fft_data'));
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(e)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")
% clim([0,max(abs(fft_data(:)))/10])
% 
% subplot(236)
% imagesc(abs(fft_data'));
% xlim([2300,3100])
% ylim([2450,4050])
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*400+2450)
% yticklabels(((1:4)*400+2450)/1000*300)
% xlabel({"Frequency (Hz)","(f)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")
% clim([0,max(abs(fft_data(:)))/10])

% % 图7 方法成像 2
% f_size=18;
% y=load("对齐方法\参数对齐\二阶_\二阶\y.mat").y;
% tmp_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*y/3e8),[],2);
% tmp_data=tmp_data(:,501:end-500);
% 
% figure
% subplot(321)
% plot(y/1000,Color="black",LineWidth=2);
% xlim([0,5120]);
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Times (sec)","(a)"})
% ylabel(" Range (km)")
% grid on;
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(322)
% imagesc(abs(tmp_data(:,1:200)'));
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% yticks((1:4)*50)
% yticklabels((1:4)*50*300/1000)
% xlabel({"Times (sec)","(b)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(323)
% imagesc(abs(fftshift(fft(tmp_data',[],2),2)));
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(c)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")
% clim([0,max(max(abs(fftshift(fft(tmp_data',[],2),2))))/5])
% 
% subplot(324)
% imagesc(abs(fftshift(fft(tmp_data',[],2),2)));
% xlim([1020,4100])
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(d)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")
% clim([0,max(max(abs(fftshift(fft(tmp_data',[],2),2))))/5])
% 
% [fft_data,~,~]=Normal_pga(tmp_data,15);
% fft_data=flip(fft_data,1);
% 
% subplot(325)
% imagesc(abs(fft_data'));
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*1500)
% yticklabels((1:4)*1500/1000*300)
% xlabel({"Frequency (Hz)","(e)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")
% clim([0,max(abs(fft_data(:)))/10])
% 
% subplot(326)
% imagesc(abs(fft_data'));
% xlim([2300,3100])
% ylim([2450,4050])
% xticks((1:4)*(5120/4))
% xticklabels(((1:4)-2)*(8.33/2))
% yticks((1:4)*400+2450)
% yticklabels(((1:4)*400+2450)/1000*300)
% xlabel({"Frequency (Hz)","(f)"})
% ylabel(" Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% colormap("gray")
% clim([0,max(abs(fft_data(:)))/10])


% % 图6 星历的曲线差
% % y=load("对齐方法\参数对齐\二阶_\二阶\y.mat").y;
% f_size=18;
% y=load("y_m.mat").y_m;
% y=y-y(1);
% 
% y1=polyval(polyfit(x,y,1),x)-y';
% y2=polyval(polyfit(x,y,2),x)-y';
% y3=polyval(polyfit(x,y,3),x)-y';
% y4=polyval(polyfit(x,y,4),x)-y';
% y5=polyval(polyfit(x,y,5),x)-y';
% 
% figure;
% subplot(231);
% plot(y/1000,Color="black",LineWidth=2);
% grid on
% xlim([0,5120]);
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Times (sec)","(a)"})
% ylabel("Range (km)")
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(232);
% plot(y1,Color="black",LineWidth=2)
% xlim([0,5120])
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Times (sec)","(b)"})
% ylabel("Range (m)")
% grid on
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(233);
% plot(y2,Color="black",LineWidth=2)
% xlim([0,5120])
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Times (sec)","(c)"})
% ylabel("Range (m)")
% grid on
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(234);
% plot(y3,Color="black",LineWidth=0.1)
% xlim([0,5120])
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Times (sec)","(d)"})
% ylabel("Range (m)")
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(235);
% plot(y4,Color="black",LineWidth=0.1)
% xlim([0,5120])
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Times (sec)","(e)"})
% ylabel("Range (m)")
% set(gca,fontname="Times New Roman",fontsize=f_size)
% 
% subplot(236);
% plot(y5,Color="black",LineWidth=0.1)
% xlim([0,5120])
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Times (sec)","(f)"})
% ylabel("Range (m)")
% set(gca,fontname="Times New Roman",fontsize=f_size)


% % 图5 星历成像结果
% y=load("y_m.mat").y_m;
% y=y-y(1);
% ephe_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*y/3e8),[],2);
% ephe_data=ephe_data(:,501:end-500);
% 
% f_size=18;
% 
% figure;
% subplot(221)
% plot(y/1000,'LineWidth',2,"Color","black")
% xlim([0,5120])
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Time (sec)","(a)"})
% ylabel("Range (km)")
% grid on;
% set(gca,fontname="Times New Roman",fontsize=f_size);
% 
% subplot(222)
% imagesc(abs(ephe_data(:,1:200)'))
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% yticklabels(yticks/1000*300)
% xlabel({"Time (sec)","(b)"})
% ylabel("Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% ephe_data=fftshift(fft(ephe_data',[],2),2);
% subplot(223)
% imagesc(abs(ephe_data))
% xticks([0,1280,2560,3840,5120])
% xticklabels([-8.33,-4.165,0,4.165,8.33])
% yticklabels(yticks/1000*300)
% xlabel({"Frequency (Hz)","(c)"})
% ylabel("Range (km)")
% colormap("gray");clim([min(abs(ephe_data(:))),max(abs(ephe_data(:)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(224)
% imagesc(abs(ephe_data))
% xlim([3300,3900])
% ylim([2450,4050])
% xticklabels(((1:4)*200+3300-2560)*(8.33/2560))
% yticklabels(((1:4)*400+2450)/1000*300)
% xlabel({"Frequency (Hz)","(d)"})
% ylabel("Range (km)")
% colormap("gray");clim([min(abs(ephe_data(:))),max(abs(ephe_data(:)))/5])
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);

% % 图4 数据展示
% f_size=22;
% figure;
% subplot(131)
% imagesc(abs(data'))
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% yticks((1:4)*1500)
% yticklabels(yticks/1000*300)
% xlabel({"Time (sec)","(a)"})
% ylabel("Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% subplot(132)
% imagesc(abs(data(:,401:600)'))
% xticks((1:5)*1000)
% xticklabels(((1:5))*60)
% yticks(((1:4))*50)
% yticklabels((yticks+400)/1000*300)
% xlabel({"Time (sec)","(b)"})
% ylabel("Range (km)")
% set(gca,"YDir","normal",fontname="Times New Roman",fontsize=f_size);
% 
% single_x=zeros(1,M);
% for each_m=1:M
%     single_x(each_m)=find(abs(data(each_m,:))>3,1);
% end
% 
% subplot(133)
% plot(single_x/1000*300)
% xlim([0,5120])
% ylim([250,350])
% xticks((1:5)*1000)
% xticklabels((1:5)*60)
% xlabel({"Time (sec)","(c)"})
% ylabel("Range (km)")
% set(gca,fontname="Times New Roman",fontsize=f_size)