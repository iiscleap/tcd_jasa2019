
addpath('../plotting/')
addpath('../plotting/cbrewer/')
data_path = './data/results/';



% ----- plot hit/mis/fa
data_1 = readtable([data_path 'human_mean_hit_miss_fa.csv'],'Delimiter',',');
data_2 = readtable([data_path 'ibmWatson_mean_hit_miss_fa.csv'],'Delimiter',',');
data_3 = readtable([data_path 'offlinePLDA_mean_hit_miss_fa.csv'],'Delimiter',',');
data_4 = readtable([data_path 'onlinePLDA_mean_hit_miss_fa.csv'],'Delimiter',',');
data_5 = readtable([data_path 'textLSTM_mean_hit_miss_fa.csv'],'Delimiter',',');

% ----- make plot
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 3;

close all;
figure;
Y = [data_1.MU_HIT data_1.MU_MISS data_1.MU_FA;data_2.MU_HIT data_2.MU_MISS data_2.MU_FA;...
    data_3.MU_HIT data_3.MU_MISS data_3.MU_FA; data_4.MU_HIT data_4.MU_MISS data_4.MU_FA;...
    data_5.MU_HIT data_5.MU_MISS data_5.MU_FA];
h = bar(Y*100);
ylim([0 105])
cmap = cbrewer('seq','Blues',100);
colormap(cmap);
grid on;
l = cell(1,3);
l{1} = 'HIT';
l{2} = 'MISS';
l{3} = 'FA';
hlegend = legend(h,l);
rect = [0.58 0.75 .025 .15]; %[left bottom width height]
set(hlegend,'Position',rect);
xticks([1 2 3 4 5]);
xticklabels({'HUMAN','WATSON','OFF-PLDA','ON-PLDA','TEXT'});
xlabel('APPROACH')
ylabel('PERCENTAGE')
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
    set(gcf, 'PaperPosition', [0 0 8 4]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/final_approach_hit_ma_fa_comparison.fig']);
    print(['./figures/final_approach_hit_ma_fa_comparison.eps'],'-depsc','-r300');    
end           


% ----- plot talkerwise fa
data_1 = readtable([data_path 'human_talkerwise_fa.csv'],'Delimiter',',');
data_2 = readtable([data_path 'ibmWatson_talkerwise_fa.csv'],'Delimiter',',');
data_3 = readtable([data_path 'offlinePLDA_talkerwise_fa.csv'],'Delimiter',',');
data_4 = readtable([data_path 'onlinePLDA_talkerwise_fa.csv'],'Delimiter',',');
data_5 = readtable([data_path 'textLSTM_talkerwise_fa.csv'],'Delimiter',',');


% ----- make plot
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 3;

figure;
Y = [data_1.MU_FA'; data_2.MU_FA';data_3.MU_FA'; data_4.MU_FA';data_5.MU_FA' ];
h = bar(Y*100);
cmap = cbrewer('seq','Blues',100);
colormap(cmap);
grid on;
l = cell(1,3);
l{1} = 'T1';
l{2} = 'T2';
l{3} = 'T3';
l{4} = 'T4';
l{5} = 'T5';

hlegend = legend(h,l);
rect = [0.6 0.7 .025 .15]; %[left bottom width height]
set(hlegend,'Position',rect);
xticks([1 2 3 4 5]);
xticklabels({'HUMAN','WATSON','OFF-PLDA','ON-PLDA','TEXT'});
xlabel('APPROACH')
ylabel('PERCENTAGE')
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
    set(gcf, 'PaperPosition', [0 0 8 4]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/final_approach_talkerwise_fa_comparison.fig']);
    print(['./figures/final_approach_talkerwise_fa_comparison.eps'],'-depsc','-r300');    
end           



