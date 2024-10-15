clear;clc;
close all;
% delete(gcp('nocreate'));
% parpool(feature('numCores')-1);


fc=430e6;       % 载频
fs=0.5e6;  % 采样率 20e6

data=load('..\..\..\data_origin.mat').data;
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
data=ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*fit_y/3e8),[],2);
save('fit_p.mat','fit_p');save('fit_t.mat','fit_y');


tic
wolf_dim=2;
wolf_agents=8;
search_times=70;

lb=zeros(1,2);
ub=zeros(1,2);
lb(1)=-5*power(10,p_order(1));
ub(1)=5*power(10,p_order(1));
lb(2)=-5*power(10,(p_order(2)-1));
ub(2)=5*power(10,(p_order(2)-1));

wolfs=(rand(wolf_agents,wolf_dim)-0.5).*(ub-lb);wolfs(1,:)=0;
a_wolf=0; b_wolf=0; d_wolf=0;a_score=inf; b_score=inf; d_score=inf;

score_his=zeros([1,search_times]);
posi_his=zeros([wolf_dim,search_times]);

flag_4_wolf_a=zeros([1,search_times]);
flag_4_wolf_b=zeros([1,search_times]);
flag_4_wolf_d=zeros([1,search_times]);

wolf_score_list=zeros(1,wolf_agents);

tic
for each_iter=1:search_times
    a=2*(1-each_iter/search_times);
    
    for each_wolf=1:wolf_agents
        tmp_p=[wolfs(each_wolf,:) 0];
        y=polyval(tmp_p,x);
        y=y'*(3e8/fs);

        % tic
        % for each_time=0:(read_time/2-1)
        %     tmp_data=tmp_data+sum(abs(ifft(fft(data((each_time*125+1):((each_time+1)*125),:),[],2).*exp(2j*pi*(fc+refe_scale).*y((each_time*125+1):((each_time+1)*125))/3e8),[],2)),1);
        %     disp(each_time)
        % end
        % toc
        
        tmp_data=mean(abs(ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*y/3e8),[],2)),1); % 距离时域fft计算以后再ifft返回
        % tmp_data=tmp_data/M;
        
        tmp_entropy=-sum((tmp_data/sum(tmp_data)).*log2(tmp_data/sum(tmp_data)));

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
    
    [~,idx]=sort(wolf_score_list);
    X1=a_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*a_wolf-wolfs);
    X2=b_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*b_wolf-wolfs);
    X3=d_wolf- (rand(wolf_agents,wolf_dim).*2.*a-a).*abs(rand(wolf_agents,wolf_dim).*2.*d_wolf-wolfs);

    % flag4levi=rand(wolf_agents,1)<0.9;   
    % levi=flag4levi.*(0.05*(wolfs-a_wolf).*(rand(wolf_agents,1)./(rand(wolf_agents,1).^(rand(wolf_agents,1)*2))));
    % wolfs=(X1+X2+X3)/3+levi;
    wolfs=(X1+X2+X3)/3;

    % wolfs(idx(end),:)=(b_wolf-d_wolf)/2+a_wolf;
    % wolfs(idx(end-1),:)=(a_wolf-d_wolf)/2+b_wolf;
    % wolfs(idx(end-2),:)=(a_wolf-b_wolf)/2+d_wolf;

    flag4lb=wolfs<lb;
    flag4ub=wolfs>ub;
    while sum(sum(flag4lb+flag4ub))~=0
        wolfs=wolfs.*~(flag4ub+flag4lb)+flag4lb.*(2*lb-wolfs)+flag4ub.*(2*ub-wolfs);
        flag4lb=wolfs<lb;
        flag4ub=wolfs>ub;
    end
    score_his(each_iter)=a_score;
    posi_his(:,each_iter)=a_wolf;
    disp("iter："+string(each_iter)+"。a："+string(a)+"。score："+string(a_score)+"。wolfs："+string(a_wolf));

end     % end for each_iter
toc



find_p=[a_wolf 0];
find_y=polyval(find_p,x)*3e8/fs;
find_y=find_y';

tmp_data=(ifft(fft(data,[],2).*exp(2j*pi*(fc+refe_scale).*find_y/3e8),[],2));
figure;imagesc(abs(fftshift(fft(tmp_data',[],2),2)));set(gca,"YDir","normal");colorbar;

final_y=fit_y+find_y;

save('find_p.mat','find_p')
save('find_y.mat','find_y');
save('final_y.mat','final_y');
save('score_his.mat','score_his');
save('flag_4_wolf_a.mat','flag_4_wolf_a');
save('flag_4_wolf_b.mat','flag_4_wolf_b');
save('flag_4_wolf_d.mat','flag_4_wolf_d');

