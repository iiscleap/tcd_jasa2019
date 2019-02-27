


close all; clear all; clc;
data_path = './data/PLDA/output/';
fname = 'offline_results.csv';

spkrs = [374,2843,5456,7447,7505];

data = readtable([data_path fname],'Delimiter',',');

spkr_wise_indx = cell(5,5);
for i = 1:length(data.FNAME)
    fname = data.FNAME{i};
    indx = strfind(fname,'_');
    spkr_1 = str2double(fname(indx(1)+1:indx(2)-1));
    spkr_2 = str2double(fname(indx(6)+1:indx(7)-1));
    tchange = str2double(fname(indx(10)+1:indx(11)-1));
    
    for j = 1:length(spkrs)
        if spkr_1 == spkrs(j)
            for k = 1:length(spkrs)
                if spkr_2 == spkrs(k)
                    spkr_wise_indx{j,k} = [spkr_wise_indx{j,k} i];
                end
            end
        end
    end
end

spkr_wise_rt = cell(5,5);
spkr_wise_hit = zeros(5,5);
spkr_wise_miss = zeros(5,5);
spkr_wise_fa = zeros(5,5);
tot_stimuli_hit = zeros(5,5);
tot_stimuli_fa = zeros(5,5);

RT = zeros(200,3);
cnt = 0;
indx_hit = [];
indx_miss = [];
indx_fa = [];
hit_window = 750; % in msec
for i = 1:size(spkr_wise_indx,1)
    for j = 1:size(spkr_wise_indx,2)
        for k = 1:length(spkr_wise_indx{i,j})
            dindx = spkr_wise_indx{i,j}(k);
            fname = data.FNAME{dindx};
            indx = strfind(fname,'_');
            spkr_1 = str2double(fname(indx(1)+1:indx(2)-1));
            spkr_2 = str2double(fname(indx(6)+1:indx(7)-1));
            tchange = str2double(fname(indx(10)+1:indx(11)-1));
%             [spkr_1 spkr_2]
            if spkr_1~=spkr_2
                indx_hit = [indx_hit dindx];
                rchange = data.TCHANGE(dindx)*1000;
                RT(indx_hit(end),1) = tchange;
                RT(indx_hit(end),2) = rchange; 
                RT(indx_hit(end),3) = rchange - tchange; 
                spkr_wise_rt{i,j} = [spkr_wise_rt{i,j} RT(indx_hit(end),2)]; 
                
                if (RT(indx_hit(end),3)>(-hit_window)) && (RT(indx_hit(end),3)<2000) % hit detection window
                    spkr_wise_hit(i,j) = spkr_wise_hit(i,j)+1;
                    tot_stimuli_hit(i,j) = tot_stimuli_hit(i,j)+1;
                elseif (RT(indx_hit(end),3)>2000) || (data.NSPKRS(dindx) == 1) || (RT(indx_hit(end),2) == 0) % miss
                    if (RT(indx_hit(end),2) == 0)
                        RT(indx_hit(end),2) = 12500;
                    end
                    indx_miss = [indx_miss dindx];
                    spkr_wise_miss(i,j) = spkr_wise_miss(i,j)+1;
                    tot_stimuli_hit(i,j) = tot_stimuli_hit(i,j)+1;
                else % fa
                    indx_fa = [indx_fa dindx];
                    spkr_wise_fa(i,j) = spkr_wise_fa(i,j)+1;                    
                    tot_stimuli_fa(i,j) = tot_stimuli_fa(i,j)+1;
                end
            else 
                if data.NSPKRS(dindx)>1 %fa
                    spkr_wise_fa(i,j) = spkr_wise_fa(i,j)+1;
                end
                tot_stimuli_fa(i,j) = tot_stimuli_fa(i,j)+1;
            end
        end
    end
end

% ----- make plots
close all;

FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 3;

clr = {[0.0 0.5 0.0],[0.0 0.5 1],'r','c'};
% plot(RT(:,1),RT(:,2),'o','color',clr{1},'MarkerSize',4,'MarkerFaceColor',clr{1},'MarkerEdgeColor','k');

% ----- plot change instant versus response instant 
if 1
RT(:,1) = RT(:,1)/1000;     
RT(:,2) = RT(:,2)/1000;     
figure; 
plot(RT(indx_hit,1),RT(indx_hit,2),'o','color',clr{1},'MarkerSize',1.5,'MarkerFaceColor',clr{1}); hold on;
h(1) = plot(RT(indx_hit(1),1),RT(indx_hit(1),2),'o','color',clr{1},'MarkerSize',1.5,'MarkerFaceColor',clr{1});
hold on;
plot(RT(indx_miss,1),RT(indx_miss,2),'o','color',clr{2},'MarkerSize',1.5,'MarkerFaceColor',clr{2});
% h(2) = plot(RT(indx_miss(1),1),RT(indx_miss(1),2),'o','color',clr{3},'MarkerSize',1.5,'MarkerFaceColor',clr{2});
hold on;
plot(RT(indx_fa,1),RT(indx_fa,2),'o','color',clr{3},'MarkerSize',1.5,'MarkerFaceColor',clr{3});
h(2) = plot(RT(indx_fa(1),1),RT(indx_fa(1),2),'o','color',clr{3},'MarkerSize',1.5,'MarkerFaceColor',clr{3});
hold on;
plot([5 10],[5 10],'-','color','k',LW,LWval);    
% hold on;
% plot([5 10],[5+.225 10+.225],'-.','color',[0.7 0.7 0.7],LW,LWval);    
hold on;
plot([5 10],[5-hit_window/1000 10-hit_window/1000],'-.','color',[0.7 0.7 0.7],LW,LWval);    
hold on;
plot([5 10],[5+2 10+2],'-.','color',[0.7 0.7 0.7],LW,LWval);    
end
% hlegend = legend(h,{'HIT','FA'});
% rect = [0.85 0.25 .001 .05]; %[left bottom width height]
% set(hlegend,'Position',rect);

xlim([4500 10000+300]/1000);
ylim([0 12800]/1000)
yticks([0:2000:12000]/1000);
yticklabels({'0','2','4','6','8','10','12'});

grid on;
set(gca,FS,FSval,'box','on');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );

if 1
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3 3]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/offlinePLDA_tchange_rchange.fig']);
print(['./figures/offlinePLDA_tchange_rchange.eps'],'-depsc','-r300');  
print(['./figures/offlinePLDA_tchange_rchange.png'],'-dpng','-r300');    
end

return;
% ----- make plots for subjectwise MISS, FA, HIT 

mu_hit = sum(spkr_wise_hit(:))/sum(tot_stimuli_hit(:));
mu_miss = sum(spkr_wise_miss(:))/sum(tot_stimuli_hit(:));
mu_fa = sum(spkr_wise_fa(:))/(sum(tot_stimuli_hit(:))+sum(tot_stimuli_fa(:)));

% ----- store mean hit/ miss/ fa/ mat file
if 0
store_path = './data/results/';
fileID = fopen([store_path 'offlinePLDA_mean_hit_miss_fa.csv'],'w');
fprintf(fileID,'MU_HIT,MU_MISS,MU_FA\n');
fprintf(fileID,'%.4f,%.4f,%.4f\n',mu_hit,mu_miss,mu_fa);
fclose(fileID);
end

figure
h(1) = plot([1 1],[0 mu_hit*100],'o-.','color','k',LW,LWval,'MarkerSize',MSval,'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor',[0 0.5 0]);
hold on;
h(2) = plot([1 1],[0 mu_miss*100],'o-.','color','k',LW,LWval,'MarkerSize',MSval,'MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[0 0.5 1]);
h(3) = plot([1 1],[0 mu_fa*100],'o-.','color','k',LW,LWval,'MarkerSize',MSval,'MarkerFaceColor','r','MarkerEdgeColor','r');
grid on;
ylabel('PERCENTAGE');
xlabel('MACHINE INDEX');
% xticks([1:17])
% xticklabels({'1','2','3','4','5','6','7','8','9',...
%              '10','11','12','13','14','15','16','17'});
yticks([0:20:100])
xlim([0 18]);
ylim([0 105]);
grid on;
set(gca,FS,FSval,'box','on');
hlegend = legend(h,{'HIT','MISS','FALSE ALRAM'});
rect = [0.6 0.5 .25 .25]; %[left bottom width height]
set(hlegend,'Position',rect);

if 0
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 8 4]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/offlinePLDA_hit_miss_fa.fig']);
    print(['./figures/offlinePLDA_hit_miss_fa.eps'],'-depsc','-r300');    
end           

% ----- plot speaker wise FA
% get the speaker wise false alarm matrix
mu_fa = zeros(5,1);
for i = 1:5
    mu_fa(i,1) = sum(spkr_wise_fa(i,:))/(sum(tot_stimuli_hit(i,:))+sum(tot_stimuli_fa(i,:)));
end
MU_FA = mu_fa;

% ----- store speaker-wise fa mat file
if 0
store_path = './data/results/';
fileID = fopen([store_path 'offlinePLDA_talkerwise_fa.csv'],'w');
fprintf(fileID,'Talker_ID, MU_FA\n');
for i = 1:length(spkrs)
    fprintf(fileID,'%d,%.4f\n',spkrs(i), mu_fa(i));
end
fclose(fileID);
end

if 0
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 3;
figure;
stem(1:5,MU_FA*100,'o','color','k','MarkerSize',MSval,'MarkerFaceColor','r','MarkerEdgeColor','r');
xlim([0 6]);
xticks(1:1:5)
xticklabels({'T1','T2','T3','T4','T5'});
grid on;
set(gca,FS,FSval,'box','on');
xlabel('TALKER INDEX');
ylabel('FALSE ALARM %')

if 0
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 5 5]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/offlinePLDA_talker_wise_fa.fig']);
    print(['./figures/offlinePLDA_talker_wise_fa.eps'],'-depsc','-r300');    
end           

end



