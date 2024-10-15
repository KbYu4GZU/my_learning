
close all; 
clear;clc;
%% 参数确立和数据读取

fc=430e6;       % 载频
f_sample=0.5e6;  % 采样率 20e6

data=load("..\data_origin.mat").data;
[M,N]=size(data);

refe_single_scale=zeros([1,N]);     % 信号包络对齐时变换的刻度点
refe_single_scale(1:ceil(N/2))=0:(f_sample/N):(f_sample/2-f_sample/N);
refe_single_scale(ceil(N/2)+1:N)=(-f_sample/2):(f_sample/N):(0-f_sample/N);


%% 
% 记录每次距离的情况，这个在每次迭代时候会使用前一个来相加得到正确范围的距离。
delta_r=zeros([M,1]);
entropy_com=zeros([M,2]);   % 记录每次信号"对齐"前后的熵，作为对比

% 狼群三个主要参数
wolf_dim=1;
wolf_agents=4;
search_times=80;

% 记录每一个信号的大搜索里面最优灰狼得到的熵值和位置
wolf_score=zeros([M,search_times]);
wolf_posit=zeros([M,search_times]);

% 搜索范围，上下界，信号向一边偏，并且可以估计相邻信号位移量，
% 因此可以比较精准定义上下界，
lb=-50;ub=50;
tmp_data=zeros([1,N]);
tic
for each_m=2:M
    % 前面信号的和，用以加当前信号计算熵
    tmp_data=tmp_data+abs(data(each_m-1,:));

    % 每一次大循环需要更新一个狼群参数，包括狼群、三只头狼以及他们对应的分数
    wolfs=(rand([wolf_agents,wolf_dim])-0.5)*(ub-lb);
    a_wolf=0;b_wolf=0;d_wolf=0;a_score=inf; b_score=inf; t_score=inf;

    % 计算没找到最优值时候的熵值并记录
    entropy_com(each_m,1)=0-sum((tmp_data+abs(data(each_m,:)))./sum(tmp_data+abs(data(each_m,:))).*log2((tmp_data+abs(data(each_m,:)))./sum(tmp_data+abs(data(each_m,:)))));

    % 灰狼算法搜索循环，这个循环内计算完当前信号最优值
    for each_iter=1:search_times
        
        a=2*(1-each_m/search_times);    % 搜索控制参数  

        for each_wolf=1:wolf_agents
            wolf_data=ifft(fft(data(each_m,:)).*exp(2j*pi*(fc + refe_single_scale)*(wolfs(each_wolf)+delta_r(each_m-1))/3e8));
            tmp_entropy=0-sum((abs(wolf_data)+tmp_data)/sum(abs(wolf_data)+tmp_data).*log2((abs(wolf_data)+tmp_data)/sum(abs(wolf_data)+tmp_data)));

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
        end % end for iii=1:wolf_agents
        
        % 优化寻优和跟新规则
        X1=a_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*a_wolf-wolfs);
        X2=b_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*b_wolf-wolfs);
        X3=d_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*d_wolf-wolfs);
        % 在搜索阶段添加针对当前最优值的莱维飞行，在包围阶段添加权重
        if a>1
            flag4levi=rand(wolf_agents,1)<0.4;   
            levi=flag4levi.*(0.05*(wolfs-a_wolf).*(rand(wolf_agents,1)./(rand(wolf_agents,1).^(rand(wolf_agents,1)*2))));
            wolfs=(X1+X2+X3)/3+levi;
        else
    %         wolfs=(X1/a_score+X2/b_wolf+X3/t_score)/(1/a_score+1/b_score+1/t_score); % 权重计算按照三个分数在总分数中的情况自动调整
            wolfs=(X1+X2+X3)/3;
        end
    
        flag4lb=wolfs<lb;
        flag4ub=wolfs>ub;
        while sum(sum(flag4lb+flag4ub))~=0
            wolfs=wolfs.*~(flag4ub+flag4lb)+flag4lb.*(2*lb-wolfs)+flag4ub.*(2*ub-wolfs);
            flag4lb=wolfs<lb;
            flag4ub=wolfs>ub;
        end
        % 记录当前最优狼的参数。
        wolf_score(each_m,each_iter)=a_score;wolf_posit(each_m,each_iter)=a_wolf;

    end % end for ii=1:search_times

    entropy_com(each_m,2)=a_score;
    delta_r(each_m)=delta_r(each_m-1)+a_wolf;
    data(each_m,:)=ifft(fft(data(each_m,:)).*exp(2j*pi*(fc + refe_single_scale)*delta_r(each_m)/3e8));

    disp("i:"+string(each_m)+"。a_score："+string(a_score)+"。delta_r："+string(delta_r(each_m)))
end % end for i=2:M
toc

clear wolf_dim;clear wolf_agents;clear search_times;clear wolfs;
clear a_wolf;clear a_score;clear b_wolf;clear b_score;clear d_wolf;clear t_score;
clear each_m;clear each_iter;clear each_wolf;clear X1;clear X2;clear X3;
clear wolf_data;clear sum_data;clear tmp_entropy;clear tmp_data;

%% 由于月球与两个雷达的运动是平稳的，所以上面的delta_r需要进行平滑来得到更准确的距离信息，并重新将这个距离信息用于信号计算

% 使用拟合函数来平滑距离曲线
x=1:M;
delta_r=delta_r/600;
p=polyfit(x,delta_r',2);
y=polyval(p,x);     % 计算平滑曲线
next_r=y-delta_r';        % 由于寻优时已经进行过一次计算，此次计算是计算差值，将该差值带入重新计算一遍得到最后结果
next_r=next_r-min(y);
delta_r=delta_r';
y=y-min(y);
y=y';

save("y.mat",'y');
save("p.mat","p")
save("entropy_com.mat",'entropy_com')
save('delta_r','delta_r')
save("next_r",'next_r')





