clear;
load('./Data/mask1D.mat')
load('./Data/Voltage_datamatrix.mat');
load('./Data/zscore_datamatrix.mat');
load('./Data/label_hie_all.mat')

% TTV in each conditions and each mice
TTV=zeros(700,150);
for id=1:14
    for t=51:200
        temp=squeeze(Vol(:,t,(1:50)+50*(id-1)));
        dis=pdist(temp');
        dis=squareform(dis);
        TTV((1:50)+50*(id-1),t-50)=mean(dis);
    end
end

subplot(2,2,1)
for id=1:7
    plot1=plot(mean(TTV((1:50)+50*(id-1),:)),'b');
    plot1.Color(4) = 0.3;
    hold on
end

for id=8:14
    plot1=plot(mean(TTV((1:50)+50*(id-1),:)),'r');
    plot1.Color(4) = 0.3;
    hold on
end

TTV=zeros(700,150);
for side=1:2
    for t=51:200
        temp=squeeze(Vol(:,t,(1:350)+350*(side-1)));
        dis=pdist(temp');
        dis=squareform(dis);
        TTV((1:350)+350*(side-1),t-50)=mean(dis);
    end
end

plot(mean(TTV(1:350,:)),'b','LineWidth',3);
plot(mean(TTV(351:end,:)),'r','LineWidth',3);
title('Dist_{amp}');
xticklabels([-50,0,50,100]/150*1000)
xlabel('Time from onset (ms)')
ylabel('Exp-wise mean STP-TTV')
set(gca,'FontSize',15)

TTV2=zeros(700,150);
for id=1:14
    for t=51:200
        temp=squeeze(Vol_z(:,t,(1:50)+50*(id-1)));
        dis=pdist(temp');
        dis=squareform(dis);
        TTV2((1:50)+50*(id-1),t-50)=mean(dis);
    end
end

subplot(2,2,3)
for id=1:7
    plot1=plot(mean(TTV2((1:50)+50*(id-1),:)),'b');
    plot1.Color(4) = 0.3;
    hold on
end


for id=8:14
    plot1=plot(mean(TTV2((1:50)+50*(id-1),:)),'r');
    plot1.Color(4) = 0.3;
    hold on
end

TTV2=zeros(700,150);
for side=1:2
    for t=51:200
        temp=squeeze(Vol_z(:,t,(1:350)+350*(side-1)));
        dis=pdist(temp');
        dis=squareform(dis);
        TTV2((1:350)+350*(side-1),t-50)=mean(dis);
    end
end
plot(mean(TTV2(1:350,:)),'b','LineWidth',3);
plot(mean(TTV2(351:end,:)),'r','LineWidth',3);
title('Dist_{z-score}');
xticklabels([-50,0,50,100]/150*1000)
xlabel('Time from onset (ms)')
ylabel('Exp-wise mean STP-TTV')
set(gca,'FontSize',15)

subplot(2,2,2)
p=zeros(1,50);
for t=51:100
    p(t-50)=signrank(TTV(:,50),TTV(:,t),'tail','left');
end
[~,~,p]=fdr(p);
boxplot(TTV(:,40:100), 'Symbol','b');
hold on
for t=1:50
    if p(t)<0.01
       plot(11+t,0.13,'k*') 
    end
end
ylim([0,0.15]);
xticks([1,11,36,61])
xticklabels([-10,0,25,50]/150*1000)
title('Dist_{amp}');
xlabel('Time from onset (ms)')
ylabel('STP-TTV')
set(gca,'FontSize',15)

subplot(2,2,4)
p=zeros(1,50);
for t=51:100
    p(t-50)=signrank(TTV2(:,50),TTV2(:,t),'tail','right');
end
[~,~,p]=fdr(p);
boxplot(TTV2(:,40:100), 'Symbol','b');
hold on
for t=1:50 
    if p(t)<0.01
       plot(11+t,57,'k*') 
    end
end
xticks([1,11,36,61])
xticklabels([-10,0,25,50]/150*1000)
title('Dist_{z-score}');
xlabel('Time from onset (ms)')
ylabel('STP-TTV')
set(gca,'FontSize',15)


% response type specific
TTV_type=zeros(6,150);
TTV2=zeros(700,150);
for c=1:3
    for t=51:200
        temp=squeeze(Vol_z(:,t,label(1:350)==c));
        dis=pdist(temp');
        dis=squareform(dis);
        TTV2(label(1:350)==c,t-50)=mean(dis);
        
        temp=squeeze(Vol_z(:,t,350+find(label(351:end)==c)));
        dis=pdist(temp');
        dis=squareform(dis);
        TTV2(350+find(label(351:end)==c),t-50)=mean(dis);
    end
end

for c=1:3
    TTV_type(c,:)=mean(TTV2(label(1:350)==c,:));
end
for c=1:3
    TTV_type(3+c,:)=mean(TTV2(350+find(label(351:end)==c),:));
end

subplot(2,1,1)
for c=1:3
    plot(mean(TTV2(label(1:350)==c,:)),'LineWidth',2);
    hold on;
end
xticklabels([-50,0,50,100]/150*1000)
xlabel('Time from onset (ms)')
ylabel('Type-wise mean TTV')
set(gca,'FontSize',15)

subplot(2,1,2)
for c=1:3
    plot(mean(TTV2(350+find(label(351:end)==c),:)),'LineWidth',2);
    hold on;
end

xticklabels([-50,0,50,100]/150*1000)
xlabel('Time from onset (ms)')
ylabel('Type-wise mean TTV')
set(gca,'FontSize',15)

% correlating state space volume change with TTV change
load('./Data/volume.mat')
interval=55:80;

for i=1:6
    [r,p]=corr(Volume(i,interval)',TTV_type(i,interval)')
end