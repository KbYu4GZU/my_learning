close all; 
clear;clc;
%% 参数确立和数据读取

fc=430e6;       % 载频
fs=0.5e6;  % 采样率 20e6

data=load("..\..\data_origin.mat").data;
data=data(:,3001:end);
[M,N]=size(data);

single_x=zeros(1,M);
for each_m=1:M
    single_x(each_m)=find(abs(data(each_m,:))>3,1);
end

x=1:M;
fit_p=polyfit(x,single_x,2);
fit_p(end)=0;
fit_y=polyval(fit_p,x)*(3e8/fs);
save("fit_p.mat",'fit_p');save('fit_y.mat','fit_y')

refe_scale=[0:N/2-1,-N/2:-1]/N*fs;
data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*fit_y'/3e8),[],2);

log_range=floor(log10(abs(fit_p(1:end-1))));
refe_scale=[0:N/2-1,-N/2:-1]/N*fs;
tmp_data=zeros([M,N]);


wolf_dim=2;
wolf_agents=10;
search_times=70;

lb=-5*(10.^(log_range));
ub=5*(10.^(log_range));
wolfs=(rand(wolf_agents,wolf_dim)-0.5).*(ub-lb);       % 由于实际上拟合得到的曲线就已经比较好了，所以将第一个灰狼置零用于初始曲线。
wolfs(1,:)=0;
a_wolf=0; b_wolf=0; d_wolf=0;a_score=inf; b_score=inf; d_score=inf;

% 记录算法每次得到的最优熵和最优熵对应的位置
score_his=zeros([1,search_times]);
posi_his=zeros([wolf_dim,search_times]);

tic
for each_iter=1:search_times
    a=2*(1-each_iter/search_times);    % 搜索控制参数

    for each_wolf=1:wolf_agents
        find_y=polyval([wolfs(each_wolf,:) 0],x)*(3e8/fs);
        tmp_data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*find_y'/3e8),[],2);
        tmp_data=abs(fft(tmp_data,[],1));        
        tmp_entropy=-sum((abs(tmp_data(:))./sum(abs(tmp_data(:)))).*log2((abs(tmp_data(:))./sum(abs(tmp_data(:))))));
        % 将最优的三个狼作为最优狼，后面其他狼的移动围绕他们三个
            if tmp_entropy<a_score
                d_score=b_score;b_score=a_score;a_score=tmp_entropy;
                d_wolf=b_wolf;b_wolf=a_wolf;a_wolf=wolfs(each_wolf,:);
            elseif tmp_entropy>a_score && tmp_entropy<b_score
                d_score=b_score;b_score=tmp_entropy;
                d_wolf=b_wolf;b_wolf=wolfs(each_wolf,:);
            elseif tmp_entropy>b_score && tmp_entropy<d_score
                d_score=tmp_entropy;d_wolf=wolfs(each_wolf,:);
            end % end if score
    end
    
    
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
    score_his(each_iter)=a_score;
    posi_his(:,each_iter)=a_wolf;
    disp("iter："+string(each_iter)+"。a："+string(a)+"。score："+string(a_score)+"。wolfs："+string(a_wolf(:)));

end     % end for each_iter
toc

find_p=[a_wolf 0];
find_y=polyval(find_p,x)*(3e8/fs);

data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*find_y'/3e8),[],2);
figure;imagesc(abs(data));set(gca,"YDir","normal");title("RT域")
figure;imagesc(abs(fft(data,[],1)));set(gca,"YDir","normal");title("RD域")

final_y=fit_y'+find_y';
final_p=find_p+fit_p;


save("find_p.mat",'find_p')
save("find_y.mat","find_y")
save("final_p.mat",'final_p')
save("final_y.mat","final_y")
save("score_his.mat",'score_his')