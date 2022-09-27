clear;

% Process to obtain ME
% can skip with directly loading the MoEn_ts.mat
path='/media/user/4TB/facial/download/';
cd(path)

mouse={'2017 Apr 30/M1533M/M1533M/','2017 Apr 30/M1697F/M1697F/', '2017 Aug 01/M2242M/', ...
    '2017 Jan 25/M1532M/M1532M/', '2017 Jan 25/M1533M/M1533M/', '2017 Jan 31/M1533M/', '2017 Jul 31/M2242M/M2242M/'};

data=zeros(74,90,301,700);

for i=1:7

if i==1||i==2
    onset=675;
    low=onset-100; high=onset+200;
else
    onset=750;
    low=onset-100; high=onset+200;
end

files=dir([path mouse{i} 'Exp001_*']);
for tr=1:50
    imageFullFileName=[path mouse{i} files(tr).name];
    info = imfinfo(imageFullFileName);
    numberOfPages = length(info);
    count=1;
    for k = low : high
        % Read the kth image in this multipage tiff file.
        thisPage = imread(imageFullFileName, k);
        thisPage = imresize(double(thisPage),1/2);
        if i==4
            thisPage=thisPage(21:end,7:end);
        end
        data(:,:,count,50*(i-1)+tr)=double(thisPage);
        count=count+1;
    end	
end

files=dir([path mouse{i} 'Exp002_*']);
for tr=1:50
    imageFullFileName=[path mouse{i} files(tr).name];
    info = imfinfo(imageFullFileName);
    numberOfPages = length(info);
    count=1;
    for k = low : high
        % Read the kth image in this multipage tiff file.
        thisPage = imread(imageFullFileName, k);
        thisPage = imresize(double(thisPage),1/2);
        if i==4
            thisPage=thisPage(21:end,7:end);
        end
        data(:,:,count,350+50*(i-1)+tr)=double(thisPage);
        count=count+1;
    end	
end
end

%% compute Motion Energy time series
MoEn=zeros(300,700);
for i=1:7
    temp=abs(diff(squeeze(data(:,:,:,(1:50)+(i-1)*50)),1,3));
    dim=size(temp);
    temp=reshape(temp,[dim(1)*dim(2),dim(3),dim(4)]);
    MoEn(:,(1:50)+(i-1)*50)=squeeze(mean(temp,1));
    
    temp=abs(diff(squeeze(data(:,:,:,(1:50)+(i-1)*50+350)),1,3));
    dim=size(temp);
    temp=reshape(temp,[dim(1)*dim(2),dim(3),dim(4)]);
    MoEn(:,(1:50)+(i-1)*50+350)=squeeze(mean(temp,1));
end

MoEn_z=zscore(MoEn,[],1);
%% ME in each mouse
subplot(1,2,1)
for i=1:7
    ts=squeeze(mean(MoEn_z(:,50*(i-1)+(1:50)), 2));
    peak=max(ts);
    plot((-99:200)/150*1000,ts);
    hold on
    ylabel('ME')
    xlabel('Time from onset (ms)')
end

save('./Data/MoEn_ts.mat','MoEn','MoEn_z');

subplot(1,2,2)
for i=1:7
%     subplot(3,3,i)
    ts=squeeze(mean(MoEn_z(:,350+50*(i-1)+(1:50)), 2));
    peak=max(ts);
    plot((-99:200),ts);
    hold on
    ylabel('ME')
    xlabel('frame from onset')
end

% Statistical comparison
load('./Data/label_hie_all.mat');

%exp001
% pos1=find(label(1:350)==1);
% pos2=find(label(1:350)==2);
% pos3=find(label(1:350)==3);

% exp002
% pos1=350+find(label(351:end)==1);
% pos2=350+find(label(351:end)==2);
% pos3=350+find(label(351:end)==3);

% all
pos1=find(label==1);
pos2=find(label==2);
pos3=find(label==3);

ts=squeeze(mean(MoEn_z(:,pos1), 2));
peak=max(ts);
plot((-50:100)/150*1000,ts(50:200), 'LineWidth', 2);
hold on
ts=squeeze(mean(MoEn_z(:,pos2), 2));
peak=max([ts',peak]);
plot((-50:100)/150*1000,ts(50:200), 'LineWidth', 2);
ts=squeeze(mean(MoEn_z(:,pos3), 2));
peak=max([ts',peak]);
plot((-50:100)/150*1000,ts(50:200), 'LineWidth', 2);
ylabel('zscored ME')
xlabel('Time from onset (ms)')
set(gca,'FontSize',15)
xlim([-50,100]/150*1000);
ylim([-1,4])

p=zeros(50,3);
for t=101:150
        p(t-100,1)=ranksum(squeeze(MoEn_z(t,pos1)),squeeze(MoEn_z(t,pos2)));
        p(t-100,2)=ranksum(squeeze(MoEn_z(t,pos1)),squeeze(MoEn_z(t,pos3)));
        p(t-100,3)=ranksum(squeeze(MoEn_z(t,pos2)),squeeze(MoEn_z(t,pos3)));
end

for i=1:3
    [~,p(:,i),~]=fdr(p(:,i));
end

for t=101:150
    if p(t-100,1)<0.05
       plot((t-100)/150*1000,3.6,'k*'); 
    end
    if p(t-100,2)<0.05
       plot((t-100)/150*1000,3.7,'*','Color',[0.5,0.5,0.5]);
    end
    if p(t-100,3)<0.05
       plot((t-100)/150*1000,3.8,'*','Color',[0.75,0.75,0.75]); 
    end
end

% Show type average
% Need to load MoEn 4D sequence (MoEn_1.mat or MoEn_2.mat)
Load(./Data/MoEn_2.mat)
figure;
k=300;
for c=1:3
    subplot(3,1,c)
%     pos=find(label(1:50)==c);
    pos2=k+find(label((351:400)+k)==c);

    aveimg=squeeze(mean(MoEn2(:,:,:,pos2), 4));
    
    for i=1:10
        subplot(3,10,i+10*(c-1))
        h=imagesc(squeeze(aveimg(:,:,100+3*(i-1))), [0,6]);
        colormap('jet')
        axis off
        if c==1
            title(['t=' num2str((3*(i-1))/150*1000) 'ms'],'Fontsize',12)
        end
    end
end