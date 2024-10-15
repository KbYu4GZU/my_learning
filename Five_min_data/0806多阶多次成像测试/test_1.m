clear;clc;
close all;
% delete(gcp('nocreate'));
% parpool(feature('numCores')-1);
addpath("..\..\..\相位校正\")


fc=430e6;       % 载频
fs=0.5e6;  % 采样率 20e6

data=load('..\data_origin.mat').data;
data=data(:,3001:end);
[M,N]=size(data);
start_flag=zeros(1,M);
for each_m=1:M
    start_flag(each_m)=find(abs(data(each_m,:))>3,1);
end
refe_scale=[0:N/2-1,-N/2:-1]/N*fs;
x=1:M;

Max_test=10;
Max_power=4;
RD_entro_list=zeros(Max_power,Max_test);
RT_entro_list=zeros(Max_power,Max_test);


for each_power=2:5
    fit_p=polyfit(x,start_flag,each_power);
    fit_p(end)=0;
    p_order=floor(log10(abs(fit_p)));
    fit_y=polyval(fit_p,x)*(3e8/fs);
    fit_y=fit_y';
    power_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*fit_y/3e8),[],2);
    tmp_entropy=-sum((abs(power_data)/sum(abs(power_data))).*log2(abs(power_data)/sum(abs(power_data))));
    [~,entro_his,~]=Normal_pga(power_data,10);
    RT_entro_list(each_power-1,1)=tmp_entropy;
    RD_entro_list(each_power-1,1)=min(entro_his);



    for each_test=1:Max_test

        wolf_dim=each_power;
        wolf_agents=8;
        search_times=80;
        ub=5*power(10,p_order(1:end-1));
        lb=-5*power(10,p_order(1:end-1));
        wolfs=(rand(wolf_agents,wolf_dim)-0.5).*(ub-lb);wolfs(1,:)=0;
        a_wolf=0; b_wolf=0; d_wolf=0;a_score=inf; b_score=inf; d_score=inf;

        for each_iter=1:search_times
            a=2*(1-each_iter/search_times);
            for each_wolf=1:wolf_agents
                tmp_p=[wolfs(each_wolf,:) 0];
                y=polyval(tmp_p,x);
                y=y'*(3e8/fs);
                
                tmp_data=mean(abs(ifft(fft(power_data,[],2).*exp(2j*pi*(fc+refe_scale).*y/3e8),[],2)),1); % 距离时域fft计算以后再ifft返回
                tmp_entropy=-sum((tmp_data/sum(tmp_data)).*log2(tmp_data/sum(tmp_data)));

                % 将最优的三个狼作为最优狼，后面其他狼的移动围绕他们三个
                if tmp_entropy<a_score
                    d_score=b_score;b_score=a_score;a_score=tmp_entropy;
                    d_wolf=b_wolf;b_wolf=a_wolf;a_wolf=wolfs(each_wolf,:);
                elseif tmp_entropy>a_score && tmp_entropy<b_score
                    d_score=b_score;b_score=tmp_entropy;
                    d_wolf=b_wolf;b_wolf=wolfs(each_wolf,:);
                elseif tmp_entropy>b_score && tmp_entropy<d_score
                    d_score=tmp_entropy;d_wolf=wolfs(each_wolf,:);
                end
            end

            X1=a_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*a_wolf-wolfs);
            X2=b_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*b_wolf-wolfs);
            X3=d_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*d_wolf-wolfs);
            wolfs=(X1+X2+X3);

            flag4lb=wolfs<lb;
            flag4ub=wolfs>ub;
            while sum(sum(flag4lb+flag4ub))~=0
                wolfs=wolfs.*~(flag4ub+flag4lb)+flag4lb.*(2*lb-wolfs)+flag4ub.*(2*ub-wolfs);
                flag4lb=wolfs<lb;
                flag4ub=wolfs>ub;
            end
        end     % end for each_iter
        RT_entro_list(each_power-1,each_test+1)=a_score;
        find_y=polyval([a_wolf,0],x)*(3e8/fs);
        power_data=ifft(fft(power_data,[],2).*exp(2j*pi*(fc+refe_scale).*find_y'/3e8),[],2);
        tmp_data=power_data(:,501:end-500);
        [~,entro_his,~]=Normal_pga(tmp_data,10);
        RD_entro_list(each_power-1,each_test+1)=min(entro_his);
    end
end


x=(1:11)-1;

figure;
plot(x,RT_entro_list(1,:),x,RT_entro_list(2,:),x,RT_entro_list(3,:),x,RT_entro_list(4,:));
legend("2th-order","3th-order","4th-order","5th-order")
xlabel("Testing times")
ylabel("Entropy")
set(gca,fontname="Times New Roman",fontsize=18)

figure;
plot(x,RD_entro_list(1,:),x,RD_entro_list(2,:),x,RD_entro_list(3,:),x,RD_entro_list(4,:));
xticks()
legend("2th-order","3th-order","4th-order","5th-order")
xlabel("Testing times")
ylabel("Entropy")
set(gca,fontname="Times New Roman",fontsize=18)

save('RD_entro_list.mat','RD_entro_list');
save('RT_entro_list.mat','RT_entro_list');
