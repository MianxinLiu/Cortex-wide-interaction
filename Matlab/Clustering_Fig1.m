clear;
load('./Data/Voltage_datamatrix.mat');

% compute zscore sequences from voltage sequences
% or just load('./Data/zscore_datamatrix.mat');
Vol_z=zeros(size(Vol));
for t=1:301
    for tr=1:700
        Vol_z(:,t,tr)=zscore(squeeze(Vol(:,t,tr)));
    end
end

% Hausdorff distance 
% EXP 1&2 seperated
disM_L=zeros(350,350);
for i=1:350
    for j=1:350
        if i~=j
            mindis=zeros(1,11);
            for t=110:120
                set=[squeeze(Vol_z(:,t,i)),squeeze(Vol_z(:,110:120,j))];
                temp=pdist(set');
                mindis(t-109)=min(temp(1:11));
            end
            disM_L(i,j)=max(mindis);
        end
    end
end

for i=1:350
    for j=(i+1):350
        themax=max(disM_L(i,j),disM_L(j,i));
        disM_L(i,j)=themax;
        disM_L(j,i)=themax;
    end
end

disM_R=zeros(350,350);
for i=1:350
    for j=1:350
        if i~=j
            mindis=zeros(1,11);
            for t=110:120
                set=[squeeze(Vol_z(:,t,i+350)),squeeze(Vol_z(:,110:120,j+350))];
                temp=pdist(set');
                mindis(t-109)=min(temp(1:11));
            end
            disM_R(i,j)=max(mindis);
        end
    end
end

for i=1:350
    for j=(i+1):350
        themax=max(disM_R(i,j),disM_R(j,i));
        disM_R(i,j)=themax;
        disM_R(j,i)=themax;
    end
end

save('./Data/HD_sides.mat','disM_L','disM_R')

% EXP 1&2 combined
load('./Data/zscore_flip.mat');

disM=zeros(700,700);
for i=1:700
    for j=1:700
        if i~=j
            mindis=zeros(1,11);
            for t=110:120
                set=[squeeze(Vol_z(:,t,i)),squeeze(Vol_z(:,110:120,j))];
                temp=pdist(set');
                mindis(t-109)=min(temp(1:11));
            end
            disM(i,j)=max(mindis);
        end
    end
end

for i=1:700
    for j=(i+1):700
        themax=max(disM(i,j),disM(j,i));
        disM(i,j)=themax;
        disM(j,i)=themax;
    end
end

save('./Data/HD_all.mat','disM')

% Principal Coordinate analysis and hierarchical clustering
[V, eigvals] = cmdscale(disM);
[eigvals eigvals/max(abs(eigvals))]

Z = linkage(disvec,'ward');
figure;
dendrogram(Z)
label = cluster(Z,'maxclust',3);

clunum=zeros(1,3);
for i=1:3
    scatter3(V(label==i,1),V(label==i,2),V(label==i,3),'filled') 
    clunum(i)=sum(label==i);
    hold on
 end
 
save('./Data/label_hie_all.mat','label')

% measure distance preservation
[V_L, eigvals] = cmdscale(disM_L);
dis1=squareform(disM_L);
dis2=pdist(V_L(:,:));
pos=(dis1==0);
dis1(pos)=[];
dis2(pos)=[];
[r,p]=corr(dis1',dis2')
scatter(dis1,dis2)

% k selection criteria

% clustering tree (recommended cutoff)
[V, eigvals] = cmdscale(disM_L);
Z = linkage(V,'ward');
figure;
dendrogram(Z)

[H,T,outperm]=dendrogram(Z, 'ColorThreshold','default');

% select k elbow
WSS=zeros(1,15);
for k=1:15
    [idx,C,sumD,D]= kmeans(V_L(:,1:3),k);
    WSS(k)=sum(sumD);
end
plot(WSS)
xlabel('k')
ylabel('WSS')

% select k silhouette
s=zeros(1,15);
for k=1:15
    [idx,C,sumD,D]= kmeans(V_L(:,1:3),k);
    s(k)=mean(silhouette(V_L(:,1:3),idx));
end
plot(1:15, s)
xlabel('k')
ylabel('silhouette')

% select gap statisics
evaluation = evalclusters(V_R(:,1:3),"kmeans","gap","KList",1:15)
subplot(3,2,4)
plot(evaluation)
