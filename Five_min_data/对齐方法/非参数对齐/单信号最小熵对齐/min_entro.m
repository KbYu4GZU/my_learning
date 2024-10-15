close all; 
clear;clc;

fc=430e6;
fs=0.5e6;

data=load("..\..\..\data_origin.mat").data;
data=data(:,3001:end);
[M,N]=size(data);

start_flag=zeros(1,M);
for each_m=1:M
    start_flag(each_m)=find(abs(data(each_m,:))>3,1);
end

refe_scale=[0:N/2-1,-N/2:-1]/N*fs;

x=1:M;
fit_p=polyfit(x,start_flag,2);
fit_p(end)=0;
fit_y=polyval(fit_p,x)*(3e8/fs);
fit_y=fit_y';
data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*fit_y/3e8),[],2);
save('fit_p.mat','fit_p');save('fit_y.mat','fit_y');

delta_r=zeros([M,1]);

wolf_dim=1;
wolf_agents=5;
search_times=80;

wolf_score=zeros([M,search_times]);
wolf_posit=zeros([M,search_times]);

lb=-20;ub=20;
com_data=zeros([1,N]);

tic
for each_m=2:M
    
    com_data=com_data+abs(data(each_m-1,:));

    % 每一次大循环需要更新一个狼群参数，包括狼群、三只头狼以及他们对应的分数
    wolfs=(rand([wolf_agents,wolf_dim])-0.5)*(ub-lb);
    a_wolf=0;b_wolf=0;d_wolf=0;a_score=inf; b_score=inf; t_score=inf;

    for each_iter=1:search_times
        a=2*(1-each_m/search_times);

        for each_wolf=1:wolf_agents
            wolf_data=ifft(fft(data(each_m,:)).*exp(2j*pi*(fc + refe_scale)*(wolfs(each_wolf)+delta_r(each_m-1))/3e8));
            tmp_data=(com_data+abs(wolf_data))/each_m;
            tmp_entropy=0-sum(tmp_data./sum(tmp_data).*log2(tmp_data./sum(tmp_data)));

            % 将最优的三个狼作为最优狼，后面其他狼的移动围绕他们三个
            if tmp_entropy<a_score
                t_score=b_score;b_score=a_score;a_score=tmp_entropy;
                d_wolf=b_wolf;b_wolf=a_wolf;a_wolf=wolfs(each_wolf);
            elseif tmp_entropy>a_score && tmp_entropy<b_score
                t_score=b_score;b_score=tmp_entropy;
                d_wolf=b_wolf;b_wolf=wolfs(each_wolf);
            elseif tmp_entropy>b_score && tmp_entropy<t_score
                t_score=tmp_entropy;d_wolf=wolfs(each_wolf);
            end % end if
        end % end for each_wolf=1:wolf_agents
        
        % 优化寻优和跟新规则
        X1=a_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*a_wolf-wolfs);
        X2=b_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*b_wolf-wolfs);
        X3=d_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*d_wolf-wolfs);
        % 在搜索阶段添加针对当前最优值的莱维飞行，在包围阶段添加权重

        wolfs=(X1+X2+X3)/3;

        flag4lb=wolfs<lb;
        flag4ub=wolfs>ub;
        while sum(sum(flag4lb+flag4ub))~=0
            wolfs=wolfs.*~(flag4ub+flag4lb)+flag4lb.*(2*lb-wolfs)+flag4ub.*(2*ub-wolfs);
            flag4lb=wolfs<lb;
            flag4ub=wolfs>ub;
        end

    end % end for each_times

    delta_r(each_m)=delta_r(each_m-1)+a_wolf;
    data(each_m,:)=ifft(fft(data(each_m,:)).*exp(2j*pi*(fc + refe_scale)*delta_r(each_m)/3e8));

    disp("i:"+string(each_m)+"。a_score："+string(a_score)+"。delta_r："+string(delta_r(each_m)))
end % end for each_m=2:M
toc

clear wolf_dim;clear wolf_agents;clear search_times;clear wolfs;
clear a_wolf;clear a_score;clear b_wolf;clear b_score;clear d_wolf;clear t_score;
clear each_m;clear each_iter;clear each_wolf;clear X1;clear X2;clear X3;
clear wolf_data;clear sum_data;clear tmp_entropy;clear com_data;

x=1:M;
delta_r=delta_r';
find_p=polyfit(x,delta_r,2);
find_p(end)=0;
find_y=polyval(find_p,x);     % 计算平滑曲线
find_y=find_y';
next_r=find_y-delta_r';        % 由于寻优时已经进行过一次计算，此次计算是计算差值，将该差值带入重新计算一遍得到最后结果


data=ifft(fft(data,[],2).*exp(2j*pi*(fc + refe_scale).*next_r/3e8),[],2);

figure;imagesc(abs(data));set(gca,"YDir","normal");title("包络对齐结果")

final_y=fit_y+find_y;
figure;plot(x,delta_r,x,find_y);
save('delta_r.mat','delta_r');
save('final_y.mat','final_y');
save("find_y.mat",'find_y');
save("find_p.mat","find_p")
