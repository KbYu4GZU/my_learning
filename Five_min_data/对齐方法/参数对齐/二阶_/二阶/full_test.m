close all; 
clear;clc;
%% 参数确立和数据读取

fc=430e6;       % 载频
f_sample=0.5e6;  % 采样率 20e6

data=load("..\..\data_origin.mat").data;
[M,N]=size(data);
% figure;imagesc(abs(data));set(gca,"YDir","normal");colorbar;
% title("原始数据")
% data=gpuArray(data);

single_x=zeros(1,M);
for each_m=1:M
    single_x(each_m)=find(abs(data(each_m,:))>3,1);
end
x=1:M;
single_x(2660)=single_x(2659);
p=polyfit(x,single_x,2);
p(end)=0;
log_range=floor(log10(abs(p(1:end-1))));

refe_scale=[0:N/2-1,-N/2:-1]/N*fs;

tmp_data=zeros([M,N]);
tmp_p=p;

%% 灰狼算法的主要部分

wolf_dim=2;
wolf_agents=8;
search_times=100;

lb=-5*(10.^(log_range));
ub=5*(10.^(log_range));
wolfs=(rand(wolf_agents,wolf_dim)-0.5).*(ub-lb);       % 由于实际上拟合得到的曲线就已经比较好了，所以将第一个灰狼置零用于初始曲线。
a_wolf=0; b_wolf=0; d_wolf=0;a_score=inf; b_score=inf; d_score=inf;

% 记录算法每次得到的最优熵和最优熵对应的位置
score_his=zeros([1,search_times]);
posi_his=zeros([wolf_dim,search_times]);

% 打算使用差分变异里面的东西，这里面定义三个差分个体得到三个解的情况。
flag_4_wolf_a=zeros([1,search_times]);
flag_4_wolf_b=zeros([1,search_times]);
flag_4_wolf_d=zeros([1,search_times]);
% 记录每一次的结果并排序
wolf_score_list=zeros(1,wolf_agents);

tic
for each_iter=1:search_times
    a=2*(1-each_iter/search_times);    % 搜索控制参数
    % 对每一个值，计算一次曲线和对应数据的熵，取最优熵
    % wolf_score=zeros(wolf_agents,1);

    for each_wolf=1:wolf_agents
        tmp_p=p+[wolfs(each_wolf,:) 0];
        y=polyval(tmp_p,x);
        y=y*(3e8/f_sample);
        y=y';

        tmp_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_single_scale).*y/3e8),[],2); % 距离时域fft计算以后再ifft返回
        tmp_entropy=-sum(sum(abs(tmp_data))/sum(abs(tmp_data(:))).*log2(sum(abs(tmp_data))/sum(abs(tmp_data(:)))));
        wolf_score_list(each_wolf)=tmp_entropy;
        % 将最优的三个狼作为最优狼，后面其他狼的移动围绕他们三个
        if tmp_entropy<a_score
            d_score=b_score;b_score=a_score;a_score=tmp_entropy;
            d_wolf=b_wolf;b_wolf=a_wolf;a_wolf=wolfs(each_wolf,:);
            if each_wolf==wolf_agents
                flag_4_wolf_a(each_iter)=1;
            elseif each_wolf==wolf_agents-1
                flag_4_wolf_b(each_iter)=1;
            elseif each_wolf==wolf_agents-2
                flag_4_wolf_d(each_iter)=1;
            end

        elseif tmp_entropy>a_score && tmp_entropy<b_score
            d_score=b_score;b_score=tmp_entropy;
            d_wolf=b_wolf;b_wolf=wolfs(each_wolf,:);
            if each_wolf==wolf_agents
                flag_4_wolf_a(each_iter)=2;
            elseif each_wolf==wolf_agents-1
                flag_4_wolf_b(each_iter)=2;
            elseif each_wolf==wolf_agents-2
                flag_4_wolf_d(each_iter)=2;
            end
        elseif tmp_entropy>b_score && tmp_entropy<d_score
            d_score=tmp_entropy;d_wolf=wolfs(each_wolf,:);
            if each_wolf==wolf_agents
                flag_4_wolf_a(each_iter)=3;
            elseif each_wolf==wolf_agents-1
                flag_4_wolf_b(each_iter)=3;
            elseif each_wolf==wolf_agents-2
                flag_4_wolf_d(each_iter)=3;
            end
        end % end if

    end     % end for each_wolf
    
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
        wolfs=(X1+X2+X3)/3;
    end

    [~,idx]=sort(wolf_score_list);
    % 差分算法部分
    wolfs(idx(end),:)=(b_wolf-d_wolf)/2+a_wolf;
    wolfs(idx(end-1),:)=(a_wolf-d_wolf)/2+b_wolf;
    wolfs(idx(end-2),:)=(a_wolf-b_wolf)/2+d_wolf;

    flag4lb=wolfs<lb;
    flag4ub=wolfs>ub;
    while sum(sum(flag4lb+flag4ub))~=0
        wolfs=wolfs.*~(flag4ub+flag4lb)+flag4lb.*(2*lb-wolfs)+flag4ub.*(2*ub-wolfs);
        flag4lb=wolfs<lb;
        flag4ub=wolfs>ub;
    end
    score_his(each_iter)=a_score;
    posi_his(:,each_iter)=a_wolf;
    disp("iter："+string(each_iter)+"。a："+string(a)+"。score："+string(a_score)+"。wolfs："+string(a_wolf(:)));

end     % end for each_iter
toc
p=p+[a_wolf 0];
y=polyval(p,x);
y=y*(3e8/f_sample);
y=y';

save("p.mat",'p')
save("y.mat","y")
save("score_his.mat",'score_his')
save("posi_his.mat",'posi_his')
