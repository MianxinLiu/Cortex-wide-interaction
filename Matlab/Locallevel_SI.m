clear
load('./Data/label_hie_all.mat')

% load SSp-bfd regional signal
mouse={'2017 Apr 30\M1533M','2017 Apr 30\M1697F','2017 Aug 01\M2242M','2017 Jan 25\M1532M','2017 Jan 25\M1533M','2017 Jan 31\M1533M','2017 Jul 31\M2242M' };

SSP=zeros(301,700);
for id=1:7
    datapath=['D:\HKBU\mouse\brain stimuli project\results\' mouse{id} '\right\'];
    load([datapath 'voltage_trial_hemis_SSP.mat']);
    SSP(:,(1:50)+50*(id-1))=squeeze(Voltage_all(4,:,:));
end

for id=1:7
    datapath=['D:\HKBU\mouse\brain stimuli project\results\' mouse{id} '\left\'];
    load([datapath 'voltage_trial_hemis_SSP.mat']);
    SSP(:,(1:50)+50*(id-1)+350)=squeeze(Voltage_all(3,:,:));
end
save('./Data/SSP.mat','SSP')

load('./Data/SSP.mat')

% compute activation level
init=mean(SSP(75:100,:));
resp=mean(SSP(105:110,:));
SSP_act=resp-init;

% plot activation level of each type
for i=1:3
    boxplot(SSP_act(label==i),'Positions',i);
    hold on;
end
xlim([0.5,3.5])
ylim([-0.002,0.007])
xticks([1,2,3])
xticklabels({'Type 1', 'Type 2', 'Type 3'})
ylabel('SSp-bfd activation level')
set(gca,'FontSize',15)

plot([1.2,1.8],[5.2*10^-3,5.2*10^-3],'k-','LineWidth',2)
plot([2.2,2.8],[5.2*10^-3,5.2*10^-3],'k-','LineWidth',2)
plot([1.2,2.8],[6*10^-3,6*10^-3],'k-','LineWidth',2)

% statistics
p1=ranksum(SSP_act(label==1),SSP_act(label==2));
p2=ranksum(SSP_act(label==2),SSP_act(label==3));
p3=ranksum(SSP_act(label==1),SSP_act(label==3));

p=[p1,p2,p3];
[~,~,p]=fdr(p);

text(1.3,5.5*10^-3,['p=' num2str(p(1),'%.3f')])
text(2.3,5.5*10^-3,['p=' num2str(p(2),'%.3f')])
text(1.8,6.3*10^-3,['p=' num2str(p(3),'%.3f')])

% show computational illustrations
plot(SSP(75:120,1:50))
hold on;
fill([1,1,25,25],[-3e-3,3e-3,3e-3,-3e-3],'b')
alpha(.1)
fill([30,30,35,35],[-3e-3,3e-3,3e-3,-3e-3],'r')
alpha(.1)