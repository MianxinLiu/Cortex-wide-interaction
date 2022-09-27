clear;
load('./Data/mask1D.mat')
load('./Data/zscore_datamatrix.mat');
load('./Data/label_hie_all.mat');

% PCA on zscore sequences
datamatrix=Vol_z(:,51:200,:);
datamatrix=reshape(datamatrix,[1236,105000]);
[V,projection,~,~,explained,~]=pca(datamatrix','NumComponents',9);

% visualize the PCs
AA=zeros(66,57,9);
[xx,yy]=find(mask_resize);

for i=1:1236
    AA(xx(i),yy(i),:)=V(i,:);
end

signalGrid=mask_resize;
signalGrid(signalGrid==0)=nan;
figure;
for i=1:9
    subplot(3,3,i)
    h=imagesc(squeeze(AA(:,:,i)),[-0.10,0.10]);
    set(h,'alphadata',~isnan(signalGrid))
    colormap('parula')
    title(['PC' num2str(i) ' (' num2str(explained(i),4) '% variances)'])
    axis off
    set(gca,'FontSize',13)
end

% compute projection scores of the whole sequence
score=zeros(301,700,3);

for trial=1:700
    for dim=1:3
        score(:,trial,dim)=squeeze(Vol_z(:,:,trial)')*V(:,dim);
    end
end

% State space atraction analysis
% plot scatter cloud changes
t=[85,90,95,100,105,110,115,120,125];

for tt=1:9
    subplot(3,3,tt)
    for c=1:3
        pos=find(label==c);
        scatter3(score(t(tt),pos,1),score(t(tt),pos,2),score(t(tt),pos,3),10,'filled')
        alpha(.4)
        hold on
    end

    pos=find(label==1);
    plot3(mean(score(t(tt),pos,1),2),mean(score(t(tt),pos,2),2),mean(score(t(tt),pos,3),2),'ks','MarkerSize',6, 'LineWidth',2)
    hold on
    pos=find(label==2);
    plot3(mean(score(t(tt),pos,1),2),mean(score(t(tt),pos,2),2),mean(score(t(tt),pos,3),2),'kd','MarkerSize',6, 'LineWidth',2)
    pos=find(label==3);
    plot3(mean(score(t(tt),pos,1),2),mean(score(t(tt),pos,2),2),mean(score(t(tt),pos,3),2),'kp','MarkerSize',6, 'LineWidth',2)
    
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    set(gca,'FontSize',15)
    title(['t=' num2str((t(tt)-100)/150*1000,'%.2f') 'ms'],'Fontsize',14)
    xlim([-40,40])
    ylim([-35,35])
    zlim([-30,30])
    grid off
    view([33.2098, 28.5669])
end


% trajectory (whole)
for c=1:3
pos=find(label==c);
tr=squeeze(mean(score(75:150,pos,1:3),2));
quiver3(tr(1:(end-1),1),tr(1:(end-1),2),tr(1:(end-1),3),diff(tr(:,1)),diff(tr(:,2)),diff(tr(:,3)),'LineWidth',1);
hold on
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
set(gca,'FontSize',15)
pos=find(label==1);
scatter3(mean(score(100,pos,1),2),mean(score(100,pos,2),2),mean(score(100,pos,3),2),200,'b+', 'LineWidth',2)
pos=find(label==2);
scatter3(mean(score(100,pos,1),2),mean(score(100,pos,2),2),mean(score(100,pos,3),2),200,'r+', 'LineWidth',2)
pos=find(label==3);
scatter3(mean(score(100,pos,1),2),mean(score(100,pos,2),2),mean(score(100,pos,3),2),200,'k+', 'LineWidth',2)

view([174.2177, 58.6569])

% trajectory (pre-stimulus)
for c=1:3
pos=find(label==c);
tr=squeeze(mean(score(75:100,pos,1:3),2));
quiver3(tr(1:(end-1),1),tr(1:(end-1),2),tr(1:(end-1),3),diff(tr(:,1)),diff(tr(:,2)),diff(tr(:,3)),'LineWidth',1);
hold on
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
set(gca,'FontSize',15)
pos=find(label==1);
scatter3(mean(score(100,pos,1),2),mean(score(100,pos,2),2),mean(score(100,pos,3),2),200,'b+', 'LineWidth',2)
pos=find(label==2);
scatter3(mean(score(100,pos,1),2),mean(score(100,pos,2),2),mean(score(100,pos,3),2),200,'r+', 'LineWidth',2)
pos=find(label==3);
scatter3(mean(score(100,pos,1),2),mean(score(100,pos,2),2),mean(score(100,pos,3),2),200,'k+', 'LineWidth',2)
view([174.2177, 58.6569])

% compute volume
Volume=zeros(6,150);

for t=1:150
    for c=1:3
        Volume(c,t)=std(score(50+t,label(1:350)==c,1))*std(score(50+t,label(1:350)==c,2))*std(score(50+t,label(1:350)==c,3));
    end
end
for t=1:150
    for c=1:3
        Volume(3+c,t)=std(score(50+t,350+find(label(351:end)==c),1))*std(score(50+t,350+find(label(351:end)==c),2))*std(score(50+t,350+find(label(351:end)==c),3));
    end
end
save('./Data/volume.mat','Volume')

subplot(2,1,1)
plot(Volume(1:3,:)'/100,'LineWidth',2)
xticklabels([-50,0,50,100]/150*1000)
xlabel('Time from onset (ms)')
ylabel('Volume in PC space (/100)')
set(gca,'FontSize',15)
subplot(2,1,2)
plot(Volume(4:end,:)'/100,'LineWidth',2)
xticklabels([-50,0,50,100]/150*1000)
xlabel('Time from onset (ms)')
ylabel('Volume in PC space (/100)')
set(gca,'FontSize',15)

Volume=zeros(2,150);

for t=1:150
    Volume(1,t)=std(score(50+t,1:350,1))*std(score(50+t,1:350,2))*std(score(50+t,1:350,3));
    Volume(2,t)=std(score(50+t,351:end,1))*std(score(50+t,351:end,2))*std(score(50+t,351:end,3));
end

figure;
plot(Volume'/100,'LineWidth',2)

xticklabels([-50,0,50,100]/150*1000)
xlabel('Time from onset (ms)')
ylabel('Volume in PC space (/100)')
set(gca,'FontSize',15)

%draw difference at three PCs
for dim=1:3
subplot(3,1,dim)
p=zeros(151,3);
dis12=zeros(1,151);dis13=zeros(1,151);dis23=zeros(1,151);
for t=50:200
    cen1=mean(score(t,label==1,dim));
    cen2=mean(score(t,label==2,dim));
    cen3=mean(score(t,label==3,dim));
    dis12(t-49)=cen1-cen2;
    dis13(t-49)=cen1-cen3;
    dis23(t-49)=cen2-cen3;
    [~,p(t-49,1)]=ttest2(score(t,label==1,dim),score(t,label==2,dim));
    [~,p(t-49,2)]=ttest2(score(t,label==1,dim),score(t,label==3,dim));
    [~,p(t-49,3)]=ttest2(score(t,label==2,dim),score(t,label==3,dim));
end
plot(((50:200)-100)/150*1000,dis12,'LineWidth',2)
hold on;
plot(((50:200)-100)/150*1000,dis13,'LineWidth',2)
plot(((50:200)-100)/150*1000,dis23,'LineWidth',2)

ylabel(['Distance PC' num2str(dim)])

if dim==1
    legend('Type 1 vs 2','Type 1 vs 3','Type 2 vs 3','autoupdate','off','box','off')
end
if dim==3
   xlabel('Time from onset (ms)') 
end
hold on

p=p*150;
% for c=1:3
%     [~,p(:,c),~]=fdr(p(:,c));
% end
top=max([dis12,dis13,dis23]);
bot=min([dis12,dis13,dis23]);
for t=50:200
    for c=1:3
        if p(t-49,c)<0.05
            plot((t-100)/150*1000,top+1+2.2*c,'*','Color',color(c,:));
        end
    end
end
set(gca,'FontSize',12)
xlim([-333.3,666.7])
% ylim([bot-1,top+1+7])
end

% Departure analysis
% single trial

tr=3;
quiver3(score(51:100,tr,1),score(51:100,tr,2),score(51:100,tr,3),diff(score(51:101,tr,1)),diff(score(51:101,tr,2)),diff(score(51:101,tr,3)),'k','LineWidth',1);
hold on;
quiver3(score(101:150,tr,1),score(101:150,tr,2),score(101:150,tr,3),diff(score(101:151,tr,1)),diff(score(101:151,tr,2)),diff(score(101:151,tr,3)),'r','LineWidth',1);
axis equal
set(gca,'view',[ -48.7696,36.7318])
set(gca,'FontSize',15)

xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

% compute statisics - distance to preongoing mean
meanprepos=squeeze(mean(score(51:100,:,1:3),1));
transdis=zeros(700,150);

for t=51:200
    for tr=1:700
        transdis(tr,t-50)=norm(squeeze(score(t,tr,1:3))-meanprepos(tr,:));
    end
end

% draw for distance
color=[    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250];
for c=1:3
plot(((51:200)-100)/150*1000,mean(transdis(label==c,:),1),'LineWidth',2,'Color',color(c,:))
hold on
end
legend('Type 1','Type 2', 'Type 3','AutoUpdate','off','box','off');
for c=1:3
for t=51:150
    temp = transdis(label==c,1:50);
    temp = reshape(temp,1,[]);
    p = ranksum(transdis(label==c,t),temp','Tail','right');
    if p<(0.05/length(51:150))
       plot((t-50)/150*1000,50+c,'*','Color',color(c,:))
    end
end
end
set(gca,'FontSize',15);
xlabel('Time from onset (ms)')
ylabel('Mean distance to pre-ongoing centroid')

% compute statisics - speed
vel=diff(score(51:200,:,1:3),1,1);
velamp=zeros(149,700);
for t=1:149
    for tr=1:700
      velamp(t,tr)=norm(squeeze(vel(t,tr,:))); 
    end
end

% draw for speed
for c=1:3
plot(((52:200)-100)/150*1000,mean(velamp(:,label==c),2),'LineWidth',2,'Color',color(c,:))
hold on
end

for c=1:3
for t=51:149
    temp = velamp(1:49,label==c);
    temp = reshape(temp,1,[]);
    p = ranksum(velamp(t,label==c),temp,'Tail','right');
    if p<(0.05/length(51:149))
       plot((t-49)/150*1000,7+c*0.1,'*','Color',color(c,:))
    end
end
end
set(gca,'FontSize',15);
xlabel('Time from onset (ms)')
ylabel('Mean speed')

% compute statisics - angle changes
velangdiff=zeros(148,700);
for t=1:148
    for tr=1:700
      velangdiff(t,tr)=atan2(norm(cross(squeeze(vel(t,tr,:)),squeeze(vel(t+1,tr,:)))),dot(squeeze(vel(t,tr,:)),squeeze(vel(t+1,tr,:))));
    end
end

% draw for angle
for c=1:3
plot(((53:200)-100)/150*1000,mean(velangdiff(:,label==c),2),'LineWidth',2,'Color',color(c,:))
hold on
end

for c=1:3
for t=50:148
    temp = velangdiff(1:48,label==c);
    temp = reshape(temp,1,[]);
    p = ranksum(velangdiff(t,label==c),temp,'Tail','left');
    if p<(0.05/length(50:148))
       plot((t-48)/150*1000,1.1+c*0.02,'*','Color',color(c,:))
    end
end
end
set(gca,'FontSize',15);
xlabel('Time from onset (ms)')
ylabel('Mean angle changes')

% plot distribution

for dim=1:3
    subplot(1,3,dim)
    rest=squeeze(score(1:100,:,dim));
    rest=reshape(rest,1,[]);
    e1=squeeze(score(116,label==1,dim));
    e1=reshape(e1,1,[]);
    e2=squeeze(score(116,label==2,dim));
    e2=reshape(e2,1,[]);
    e3=squeeze(score(116,label==3,dim));
    e3=reshape(e3,1,[]);

    boxplot(rest,'Positions',1, 'Colors',[0.5,0.5,0.5], 'BoxStyle','filled', 'Symbol','+')
    hold on
    boxplot(e1,'Positions',2, 'Colors',[0,0.4470,0.7410], 'BoxStyle','filled',  'Symbol','+')
    boxplot(e2,'Positions',3, 'Colors',[0.85,0.325,0.098], 'BoxStyle','filled',  'Symbol','+')
    boxplot(e3,'Positions',4, 'Colors',[0.929,0.694,0.125], 'BoxStyle','filled',  'Symbol','+')
    xlim([0.5,4.5])

    ylabel(['PC' num2str(dim)])
    set(gca, 'xtick', [1,2,3,4]);
    set(gca, 'xticklabel', {'Spon', 'Type 1', 'Type 2', 'Type 3'});

    p1=ranksum(rest,e1)
    p2=ranksum(rest,e2)
    p3=ranksum(rest,e3)

    sigstar({[1,2],[1,3],[1,4]},[p1,p2,p3]);
end

% compute statisics - overlapping rate

propotion=zeros(200,3);
predata=score(1:100,:,:);
predata=reshape(predata,[100*size(predata,2),3]);

pre_centroid=mean(predata,1);

dist=pdist([pre_centroid; predata]);
dist=dist(1:size(predata,1));

% check radius distribution
histogram(dist,'Normalization','probability')
xlabel('Distance to ongoing centroid')
ylabel('Frequency')
radius=prctile(dist,80);
hold on;
plot([radius,radius],[0,0.03],'r-','Linewidth',2)


for c=1:3
    
    postdata=score(1:200,label==c,:);
    for t=1:200
        thecount=0;
        count=0;
        for tr=1:length(find(label==c))
            dist = pdist([squeeze(postdata(t,tr,:))';pre_centroid]);
            if dist<=radius
                count=count+1;
            end
            thecount=thecount+1;
        end
        propotion(t,c)=count/length(find(label==c));
    end
end

% draw for overlapping rate
for c=1:3
    plot([-50:100]/150*1000, propotion(50:end,c))
    hold on
end
xlabel('Time from onset (ms)')
ylabel('Proportion of trials within spontaneous space')

hold on

p=zeros(100,3);
for t=100:200
    for c=1:3
        p(t-99,c)=sum(propotion(1:100,c)<propotion(t,c))/100;
    end
end

for t=101:200
    for c=1:3
        if p(t-99,c)<(0.05/100)
            plot((t-100)/150*1000,0.9+0.05*c,'*','Color',color(c,:));
        end
    end
end
