clear;
load('D:\HKBU\mouse\brain stimuli project\results\mask1D.mat')
load('D:\HKBU\mouse\brain stimuli\results\Voltage_datamatrix.mat');
load('D:\HKBU\mouse paper\To Neuron\Cell reports revision\codes\matlab\label_hie_all.mat')

% remove mean
datamatrix=Vol;

average_voltage=zeros(1236,301,3);
for c=1:3
    average_voltage(:,:,c)=mean(datamatrix(:,:,(label(1:350)==c)),3);
%     for t=1:301
%         average_voltage(:,t,c)=average_voltage(:,t,c)/norm(squeeze(average_voltage(:,t,c)));
%     end
end

corr_rec_rmres=zeros(700,50);
corr_rec_before=zeros(700,50);

for id=1:7
    for trial=1:50
        onset=squeeze(datamatrix(:,101,trial+50*(id-1)));
        for t=1:50
        rm_res=squeeze(datamatrix(:,101+t,trial+50*(id-1)))-average_voltage(:,101+t, label(trial+50*(id-1)));
        back=squeeze(datamatrix(:,101-t,trial+50*(id-1)));        
        corr_rec_rmres(trial+50*(id-1),t)=corr2(rm_res,onset);        
        corr_rec_before(trial+50*(id-1),t)=corr2(back,onset);
        end
    end
end

average_voltage=zeros(1236,301,3);
for c=1:3
    average_voltage(:,:,c)=mean(datamatrix(:,:,(label(351:end)==c)),3);
%     for t=1:301
%         average_voltage(:,t,c)=average_voltage(:,t,c)/norm(squeeze(average_voltage(:,t,c)));
%     end
end

for id=8:14
    for trial=1:50
        onset=squeeze(datamatrix(:,101,trial+50*(id-1)));
        for t=1:50
        rm_res=squeeze(datamatrix(:,101+t,trial+50*(id-1)))-average_voltage(:,101+t, label(trial+50*(id-1)));
        back=squeeze(datamatrix(:,101-t,trial+50*(id-1)));        
        corr_rec_rmres(trial+50*(id-1),t)=corr2(rm_res,onset);        
        corr_rec_before(trial+50*(id-1),t)=corr2(back,onset);
        end
    end
end

figure;
for c=1:3
    subplot(3,2,2*c-1)
    plot((1:50)/150*1000, median(corr_rec_rmres(label==c,:),1), 'LineWidth',2)
    hold on;
    plot((1:50)/150*1000, median(corr_rec_before(label==c,:),1), 'LineWidth',2)
    xlim([0,51/150*1000])
    set(gca,'FontSize',15)
    xlabel('Time from onset (ms)')
    ylabel('Median correlation')
end

for c=1:3
    subplot(3,2,2*c)
    pos=(label==c);
    for t=1:50
        p(t) = signrank(corr_rec_rmres(pos,t),corr_rec_before(pos,t));
    end
    boxplot(corr_rec_rmres(pos,1:50),'Boxstyle','filled','Positions',(1:2:100)/150*1000,'Colors','b','Symbol','b')
    set(gca,'XTick',[])
    set(gca,'XTickLabel',{})
    hold on;
    boxplot(corr_rec_before(pos,1:50),'Boxstyle','filled','Positions',(2:2:100)/150*1000,'Colors','r','Symbol','r')
    set(gca,'XTickLabel',{})
    %[~,~,p]=fdr(p);
    for t=1:50
        if p(t)<0.05
           scatter((t*2-0.5)/150*1000,1,'k*'); 
        end
    end
    set(gca,'XTick',([1,(10:10:100)]+0.5)/150*1000)
    set(gca,'XTickLabel',round(([1,5:5:50])/150*1000))
    xlabel('Time from onset (ms)')
    ylabel('Median correlation')
    set(gca,'FontSize',15)
end

