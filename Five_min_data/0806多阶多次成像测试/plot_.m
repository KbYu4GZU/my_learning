load('RD_entro_list.mat');
load("RT_entro_list.mat");
x=1:10;
figure;
% subplot(1,2,1)

plot(x,RT_entro_list(2,:),x,RT_entro_list(3,:),x,RT_entro_list(4,:),x,RT_entro_list(5,:));
legend(["2th-order","3th-order","4th-order","5th-order"])
xlabel("")
ylabel("Entropy")
set(gca,fontname="Times New Roman",fontsize=18)

figure;
% subplot(1,2,2)
plot(x,RD_entro_list(2,:),x,RD_entro_list(3,:),x,RD_entro_list(4,:),x,RD_entro_list(5,:));
legend(["2th-order","3th-order","4th-order","5th-order"])
xlabel("")
ylabel("Entropy")

