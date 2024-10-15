close all; 
clear;clc;
%% 参数确立和数据读取

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
p_order=floor(log10(abs(fit_p)));
fit_y=polyval(fit_p,x)*(3e8/fs);
fit_y=fit_y';
data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*fit_y'/3e8),[],2);
save('fit_p.mat','fit_p');save('fit_y.mat','fit_y');


com_data=zeros([1,N]);
tmp_data=zeros([1,N]);

wolf_dim=1;
wolf_agents=5;
search_times=80;
lb=-20;ub=20;

% 记录距离和结果
delta_r=zeros([1,M]);
wolf_score=zeros([1,M]);

tic
for each_m=2:M
    % 刷新算法初始值
    wolfs=(rand([wolf_agents,wolf_dim])-0.5)*(ub-lb);
    a_wolf=zeros(1,wolf_dim);b_wolf=zeros(1,wolf_dim);d_wolf=zeros(1,wolf_dim);
    a_score=0; b_score=0; d_score=0;
    
    if each_m<11
        com_data=mean(abs(data(1:each_m-1,:)),1);
    else
        com_data=mean(abs(data(each_m-10:each_m-1,:)),1);
    end

    for each_iter=1:search_times
        % a是线性收敛控制变量，用于狼群位置变化的地方
        a=2*(1-each_iter/search_times);

        % 这个for循环计算每个灰狼的适应度，即每个灰狼得到的相关性
        for each_wolf=1:wolf_agents

            % 代价函数部分，下面三行代码计算每只狼的适应度
            tmp_data=abs(ifft(fft(data(each_m,:)).*exp(2j*pi*(fc + refe_scale)*(wolfs(each_wolf)+delta_r(each_m-1))/3e8)));
            tmp_score=corrcoef(tmp_data,com_data);  % 相关性计算
            tmp_score=tmp_score(1,2);

            % 将最优的三个狼作为最优狼，后面其他狼的移动围绕他们三个
            if tmp_score>a_score
                d_score=b_score;b_score=a_score;a_score=tmp_score;
                d_wolf=b_wolf;b_wolf=a_wolf;a_wolf=wolfs(each_wolf);
            elseif tmp_score<a_score && tmp_score>b_score
                d_score=b_score;b_score=tmp_score;
                d_wolf=b_wolf;b_wolf=wolfs(each_wolf);
            elseif tmp_score<b_score && tmp_score>d_score
                d_score=tmp_score;d_wolf=wolfs(each_wolf);
            end % end if

        end % end for each_wolf=1:wolf_agents
        
        % 优化寻优和跟新规则
        X1=a_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*a_wolf-wolfs);
        X2=b_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*b_wolf-wolfs);
        X3=d_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*d_wolf-wolfs);

        wolfs=(X1+X2+X3)/3;

        flag4lb=wolfs<lb;
        flag4ub=wolfs>ub;
        while sum(sum(flag4lb+flag4ub))~=0
            wolfs=wolfs.*~(flag4ub+flag4lb)+flag4lb.*(2*lb-wolfs)+flag4ub.*(2*ub-wolfs);
            flag4lb=wolfs<lb;
            flag4ub=wolfs>ub;
        end

    end % end for ii=1:serch_times
    
    wolf_score(each_m)=a_score;
    delta_r(each_m)=a_wolf+delta_r(each_m-1);

    data(each_m,:)=ifft(fft(data(each_m,:)).*exp(2j*pi*(fc + refe_scale)*delta_r(each_m)/3e8));
    
    disp("i:"+string(each_m)+"。a_score："+string(a_score)+"。delta_r："+string(delta_r(each_m)))
end
toc

figure;plot(delta_r);title("各个信号需要移动的距离情况");

%% 由于月球与两个雷达的运动是平稳的，所以上面的delta_r需要进行平滑来得到更准确的距离信息，并重新将这个距离信息用于信号计算

% 使用拟合函数来平滑距离曲线
x=1:M;
delta_r=delta_r';
find_p=polyfit(x,delta_r,2);
find_p(end)=0;
find_y=polyval(find_p,x);     % 计算平滑曲线
find_y=find_y';
next_r=find_y-delta_r;        % 由于寻优时已经进行过一次计算，此次计算是计算差值，将该差值带入重新计算一遍得到最后结果
data=ifft(fft(data,[],2).*exp(2j*pi*(fc + refe_scale).*next_r'/3e8),[],2);

figure;imagesc(abs(data));set(gca,"YDir","normal");title("包络对齐结果")

final_y=fit_y+find_y;
save('final_y.mat','final_y');
save("find_y.mat",'find_y');
save("find_p.mat","find_p")
save("wolf_score.mat",'wolf_score')
save('delta_r','delta_r')
save("next_r",'next_r')


