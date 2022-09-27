clear
load('./Data/vol_flip.mat');
load('./Data/zscore_flip.mat');
load('./Data/label_hie_all.mat');

% poststimulus
datamatrix=Vol;
average_voltage=zeros(1236,301,3);
for c=1:3
    average_voltage(:,:,c)=mean(datamatrix(:,:,find(label==c)),3);
end

AA=zeros(66,57,301,3);
[xx,yy]=find(mask_resize);

for i=1:1236
    AA(xx(i),yy(i),:,:)=average_voltage(i,:,:);
end

signalGrid=mask_resize;
signalGrid(mask_resize==0)=nan;
figure;
for c=1:3
    rsData=squeeze(AA(:,:,:,c));
    for i=1:11
        subplot(3,11,i+11*(c-1))
        h=imagesc(squeeze(rsData(:,:,100+3*(i-1))),[-0.003,0.003]);
        set(h,'alphadata',~isnan(signalGrid))
        colormap(darkb2r(-0.003,0.003))
        axis off
        if c==1
            title(['t=' num2str((3*(i-1))/150*1000) 'ms'],'Fontsize',12)
        end
    end
end

datamatrix=Vol_z;
average_voltage=zeros(1236,301,3);
for c=1:3
    average_voltage(:,:,c)=mean(datamatrix(:,:,find(label==c)),3);
end

AA=zeros(66,57,301,3);
[xx,yy]=find(mask_resize);

for i=1:1236
    AA(xx(i),yy(i),:,:)=average_voltage(i,:,:);
end

signalGrid=mask_resize;
signalGrid(mask_resize==0)=nan;
figure;
for c=1:3
    rsData=squeeze(AA(:,:,:,c));
    for i=1:11
        subplot(3,11,i+11*(c-1))
        h=imagesc(squeeze(rsData(:,:,100+3*(i-1))), [-3,3]);
        set(h,'alphadata',~isnan(signalGrid))
        colormap(darkb2r(-3,3))
        axis off
        if c==1
            title(['t=' num2str((3*(i-1))/150*1000) 'ms'],'Fontsize',12)
        end
    end
end

% prestimulus
load('./Data/Voltage_datamatrix');
load('./Data/zscore_datamatrix');

datamatrix=Vol;
average_voltage=zeros(1236,301,3);
for c=1:3
    average_voltage(:,:,c)=mean(datamatrix(:,:,find(label==c)),3);
end

AA=zeros(66,57,301,3);
[xx,yy]=find(mask_resize);

for i=1:1236
    AA(xx(i),yy(i),:,:)=average_voltage(i,:,:);
end

signalGrid=mask_resize;
signalGrid(mask_resize==0)=nan;
figure;
for c=1:3
    rsData=squeeze(AA(:,:,:,c));
    for i=1:10
        subplot(3,10,i+10*(c-1))
        h=imagesc(squeeze(rsData(:,:,70+3*(i-1))),[-0.001,0.001]);
        set(h,'alphadata',~isnan(signalGrid))
        colormap(darkb2r(-0.001,0.001))
        axis off
        if c==1
            title(['t=' num2str((3*(i-1)-30)/150*1000) 'ms'],'Fontsize',12)
        end
    end
end


datamatrix=Vol_z;
average_voltage=zeros(1236,301,3);
for c=1:3
    average_voltage(:,:,c)=mean(datamatrix(:,:,find(label==c)),3);
end

AA=zeros(66,57,301,3);
[xx,yy]=find(mask_resize);

for i=1:1236
    AA(xx(i),yy(i),:,:)=average_voltage(i,:,:);
end

signalGrid=mask_resize;
signalGrid(mask_resize==0)=nan;
figure;
for c=1:3
    rsData=squeeze(AA(:,:,:,c));
    for i=1:10
        subplot(3,10,i+10*(c-1))
        h=imagesc(squeeze(rsData(:,:,70+3*(i-1))), [-1,1]);
        set(h,'alphadata',~isnan(signalGrid))
        colormap(darkb2r(-1,1))
        axis off
        if c==1
            title(['t=' num2str((3*(i-1)-30)/150*1000) 'ms'],'Fontsize',12)
        end
    end
end

