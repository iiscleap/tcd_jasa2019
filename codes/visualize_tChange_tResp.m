
clear all;
%close all;

%data_path = '/home/neeks/mnt/work/neeks/spkrChange/dBase/libriSpeech/test-listening/M_5_spkrs_8_pairs_v0/';
% feats_path = '/home/neeks/mnt_1/work/neeks/tools/lre17_baseline_system/base/exp/ivectors/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/';
feats_path{1} = './data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part/mfcc/';
feats_path{2} = './data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part/fbank/';
feats_path{3} = './data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part/corti/';
feats_path{4} = './data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part/f0/';
feats_path{5} = './data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part/meanAmp/';
feats_path{6} = './data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part/engy/';


flist = dir([feats_path{1} '*.h5']);
nfiles = length(flist);

data_path = './data/rt_feats/';
% data_std = readtable([data_path 'subject_wise_mean_rt.csv'],'Delimiter',',');
data_std = readtable([data_path 'subject_wise_noise_mean_rt.csv'],'Delimiter',',');

ptcpnt_id = {'CMU-1','cmu2','CMU_3','CMU-4','cmu_5','CMU_6','CMV-7','CMU_8','cmu_9','CMU_10','cmu_11','cmu12',...
    'CMU_13','CMU-14','cmu_15','CMU_16','CMU_17','iisc_01','iisc_02','iisc_03','iisc_04','iisc_05','iisc_06',...
    'iisc_07','iisc_08','iisc_09','iisc_10','iisc_11','iisc_12'};

sel_subjs =  [13:17 18:25 26:29];

indx_subj = cell(length(ptcpnt_id),1);
for i = 1:length(ptcpnt_id)
    indx_subj{i} = [];
    for j = 1:nfiles
        if length(strfind(flist(j).name,ptcpnt_id{i}))
            indx_subj{i} =[indx_subj{i} j]; 
        end
    end
end

rt = cell(length(ptcpnt_id),1);
feats_dist = cell(length(ptcpnt_id),1);
flag_same_spk = cell(length(ptcpnt_id),1);
flag_resp = cell(length(ptcpnt_id),1);
label = cell(length(ptcpnt_id),1);
b = cell(length(ptcpnt_id),1);
mu_MAE = cell(length(ptcpnt_id),1);
std_MAE = cell(length(ptcpnt_id),1);

store_path = './data/rt_feats/';

RT = [];
indx_hit = [];
indx_miss = [];
indx_fa = [];
indx_junk_fa = [];

cnt = 1;
loop_cnt = 0;
for loop = sel_subjs
    loop_cnt = loop_cnt+1; 
    rt{loop} = zeros(1,length(indx_subj{loop}));
    flag_same_spk{loop} = zeros(1,length(indx_subj{loop}));
    flag_resp{loop} = zeros(1,length(indx_subj{loop}));
    label{loop} = zeros(1,length(indx_subj{loop}));
    b{loop} = zeros(2,1);
    feats_dist{loop} = zeros(length(indx_subj{loop}),6);
    for i = 1:length(indx_subj{loop})
        [fpath,fname,fext] = fileparts([feats_path{1} flist(indx_subj{loop}(i)).name]);
        indx = strfind(fname,'_');
        str_1 = 'tChange_';
        str_2 = '_ms_rstamp';

        indx_1 = strfind(fname,str_1)+length(str_1);
        indx_2 = strfind(fname,str_2)-1;
        tChange = str2double(fname(indx_1:indx_2));
        
        str_1 = 'rstamp_';
        str_2 = '_ms_flag';

        indx_1 = strfind(fname,str_1)+length(str_1);
        indx_2 = strfind(fname,str_2)-1;
        rChange = str2double(fname(indx_1:indx_2));
        temp = rChange - tChange;
        
        str_1 = '_lid';
        indx_1 = strfind(fname,str_1)-1;
        flag_resp{loop}(i) = str2double(fname(indx_1));
        
        rt{loop}(i) = temp;
        
        if 1
        spk_1 = str2double(fname(indx(1)+1:indx(2)-1));
        spk_2 = str2double(fname(indx(6)+1:indx(7)-1));
        if spk_1 == spk_2
            flag_same_spk{loop}(i) = 1;
        else
            flag_same_spk{loop}(i) = 0;
        end
        if (temp>225 && temp < 2000 && flag_resp{loop}(i) && ~flag_same_spk{loop}(i))
            label{loop}(i) = 3; % hit
        elseif ((temp>2000 || rChange ==0) && ~flag_resp{loop}(i) && ~flag_same_spk{loop}(i))
            label{loop}(i) = 2; % miss
        elseif ((temp<226 && ~flag_resp{loop}(i) && ~flag_same_spk{loop}(i)))
            label{loop}(i) = 1; % false alarm
        elseif ((temp<225 && temp>0 && ~flag_resp{loop}(i) && ~flag_same_spk{loop}(i)))
            label{loop}(i) = 0.5; % junk false alarm
        else
            label{loop}(i) = 0; % junk
        end
                
        if (((label{loop}(i) == 3) || (label{loop}(i) == 2) || (label{loop}(i) == 1) || (label{loop}(i) == 0.5)))% || (label{loop}(i) == 0)))% && rChange>0)
            RT(cnt,:) = [tChange rChange 0 0]; 
            FNAME{cnt} = fname;
            if (label{loop}(i) == 3)
                indx_hit = [indx_hit cnt];
                RT(cnt,3) = (RT(cnt,2)-RT(cnt,1));% - (data_std.meanRT_ms(loop_cnt));
                RT(cnt,4) = 1*(log10(RT(cnt,2)-RT(cnt,1)) - log10(data_std.meanRT_ms(loop_cnt)));
%                 RT(cnt,4) = 1./(RT(cnt,2)-RT(cnt,1) - data_std.meanRT_ms(loop_cnt));

%                 RT(cnt,3) = (RT(cnt,2)-RT(cnt,1)) - (data_std.medRT_ms(loop_cnt));
%                 RT(cnt,4) = log10(RT(cnt,2)-RT(cnt,1)) - log10(data_std.medRT_ms(loop_cnt));
            elseif (label{loop}(i) == 2)
                indx_miss = [indx_miss cnt];
                if (rChange == 0)
                    RT(cnt,2) = 12500;
                end
            elseif (label{loop}(i) == 1 && rChange>0)
                indx_fa = [indx_fa cnt];
            else
                indx_junk_fa = [indx_junk_fa cnt];
            end
            cnt = cnt+1;
        end
        end
    end
end
[length(indx_hit) length(indx_miss) length(indx_fa)]

if 1
close all;
% ----- make plots
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 5; 
LWval = 0.5;
MSval = 1;

clr = {[0.5 0.5 0.5],[0.0 0.5 0.0],[0.0 0.5 1],'r','c'};
% plot(RT(:,1),RT(:,2),'o','color',clr{1},'MarkerSize',4,'MarkerFaceColor',clr{1},'MarkerEdgeColor','k');

% ----- plot change instant versus response instant 
if 1
RT(:,1) = RT(:,1)/1000;     
RT(:,2) = RT(:,2)/1000;     
figure;
% plot(RT(:,1),RT(:,2),'o','color',clr{1},'MarkerSize',1.5,'MarkerFaceColor',clr{1});
% hold on;
plot(RT(indx_hit,1),RT(indx_hit,2),'o','color',clr{2},'MarkerSize',1.5,'MarkerFaceColor',clr{2}); hold on;
h(1) = plot(RT(indx_hit(1),1),RT(indx_hit(1),2),'o','color',clr{2},'MarkerSize',1.5,'MarkerFaceColor',clr{2});
hold on;
plot(RT(indx_miss,1),RT(indx_miss,2),'+','color',clr{3},'MarkerSize',3.5,'MarkerFaceColor',clr{3});
h(2) = plot(RT(indx_miss(1),1),RT(indx_miss(1),2),'+','color',clr{3},'MarkerSize',3.5,'MarkerFaceColor',clr{3});
hold on;
plot(RT(indx_fa,1),RT(indx_fa,2),'x','color',clr{4},'MarkerSize',3.5,'MarkerFaceColor',clr{4});
h(3) = plot(RT(indx_fa(1),1),RT(indx_fa(1),2),'x','color',clr{4},'MarkerSize',3.5,'MarkerFaceColor',clr{4});
hold on;
% plot(RT(indx_junk_fa,1),RT(indx_junk_fa,2),'o','color',clr{5},'MarkerSize',1.5,'MarkerFaceColor',clr{5});
% h(4) = plot(RT(indx_junk_fa(1),1),RT(indx_junk_fa(1),2),'o','color',clr{5},'MarkerSize',1.5,'MarkerFaceColor',clr{5});
% hold on;
plot([5 10],[5 10],'-','color',[0.7 0.7 0.7],LW,LWval);    
hold on;
plot([5 10],[5+.225 10+.225],'-.','color',[0.7 0.7 0.7],LW,LWval);    
hold on;
plot([5 10],[5+2 10+2],'-.','color',[0.7 0.7 0.7],LW,LWval);    

% uncomment below for legend
hlegend = legend(h,{'HIT','MISS','FA','DISCARD'});
rect = [0.65 0.25 .001 .05]; %[left bottom width height]
set(hlegend,'Position',rect);

xlim([4500 10000+300]/1000);
ylim([0 12800]/1000)
ylabel('RESPONSE INSTANT, tr [in s]');
xlabel('GROUND TRUTH CHANGE INSTANT, tc [in s]');
yticks([0:2000:12000 12500]/1000);
yticklabels({'0','2','4','6','8','10','12','NP'});

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
set(gcf, 'PaperPosition', [0 0 5 5]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/tchange_rchange.fig']);
print(['./figures/tchange_rchange.eps'],'-depsc','-r300');    
end

end
return;

if 0
% ----- plot change instant versus hit 
rt_hit_unique = unique(RT(indx_hit,1));
rt_hit_count = zeros(length(rt_hit_unique),1);
for i = 1:length(indx_hit)
    for j = 1:length(rt_hit_unique)
        if (rt_hit_unique(j)==RT(indx_hit(i)))
           rt_hit_count(j) = rt_hit_count(j)+1;  
        end
    end
end

rt_hit_tot_count = zeros(length(rt_hit_unique),1);
temp = [indx_hit indx_miss];
for i = 1:length(temp)
    for j = 1:length(rt_hit_unique)
        if (rt_hit_unique(j)==RT(temp(i),1))
           rt_hit_tot_count(j) = rt_hit_tot_count(j)+1;  
        end
    end
end
figure;
plot(rt_hit_unique,100*rt_hit_count./rt_hit_tot_count,'o-.','color','k',LW,LWval,'MarkerSize',MSval,'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor',[0 0.5 0]);
ylim([70 102]);
xlim([4500 10000+500]);
yticks([70:5:100]);
xlabel('CHANGE INSTANT [in msec]');
ylabel('HIT PERCENTAGE ');
grid on;
set(gca,FS,FSval,'box','on');

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );

if do_fig_save
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 3]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/hit_rchange.fig']);
print(['./figures/hit_tchange.eps'],'-depsc','-r300');    
end

end

% ----- plot change instant versus miss 
if 0
rt_miss_unique = unique(RT(indx_miss,1));
rt_miss_count = zeros(length(rt_miss_unique),1);
for i = 1:length(indx_miss)
    for j = 1:length(rt_miss_unique)
        if (rt_miss_unique(j)==RT(indx_miss(i)))
           rt_miss_count(j) = rt_miss_count(j)+1;  
        end
    end
end

rt_miss_tot_count = zeros(length(rt_miss_unique),1);
temp = [indx_hit indx_miss];
for i = 1:length(temp)
    for j = 1:length(rt_miss_unique)
        if (rt_miss_unique(j)==RT(temp(i),1))
           rt_miss_tot_count(j) = rt_miss_tot_count(j)+1;  
        end
    end
end
end


% ----- plot change instant versus reaction time 
if 0 

temp = unique(RT(indx_hit,1));

for i = 1:length(temp)
    cnt = 1;
    rt_same_indx{i} = [];
    for j = 1:length(indx_hit)
        if (RT(indx_hit(j),1)-temp(i))==0
            rt_same_indx{i}(cnt) = j;
            cnt = cnt+1;
        end
    end
end

figure;
for i = 1:length(temp)
%     plot(RT(indx_hit(rt_same_indx{i}),1),(RT(indx_hit(rt_same_indx{i}),2)-RT(indx_hit(rt_same_indx{i}),1)),'o','color',clr{2},'MarkerSize',2,'MarkerFaceColor',clr{2});
    tc(i) = RT(indx_hit(rt_same_indx{i}(1)),1);  
    mu_temp(i) = median(RT(indx_hit(rt_same_indx{i}),3));
    std_temp(i) = std(RT(indx_hit(rt_same_indx{i}),3));
end

plot(tc,mu_temp,'-o','color',clr{3},'MarkerSize',8,'MarkerFaceColor',clr{3});
hold on;
plot(tc,mu_temp+std_temp,'-','color',clr{2},'MarkerSize',8,'MarkerFaceColor',clr{2});
plot(tc,mu_temp-std_temp,'-','color',clr{2},'MarkerSize',8,'MarkerFaceColor',clr{2});

X = [ones(1,length(tc))' tc'];
Y = mu_temp';
b = pinv(X)*Y;
Y1 = X*b;
plot(tc,Y1,'-o','color',clr{1},'MarkerSize',8,'MarkerFaceColor',clr{1});
ylabel('REACTION TIME [in msec]');
xlabel('CHANGE INSTANT [in msec]');
grid on;
set(gca,FS,FSval,'box','on');

if 0
figure;
plot(tc,std_temp,'-o','color',clr{2},'MarkerSize',8,'MarkerFaceColor',clr{2});
ylabel('STD DEV. IN REACTION TIME [in msec]');
xlabel('CHANGE INSTANT [in msec]');
grid on;
set(gca,FS,FSval,'box','on');
end

end

if 1
% get chi square score
tmp_data = 10*randn(1000,1)+3;
pd_1 = fitdist(RT(indx_hit,3),'Normal');
pd_2 = fitdist(RT(indx_hit,4),'Normal');
pd_3 = fitdist(tmp_data,'Normal');

[tmp_1,tmp_1_1] = chi2gof(RT(indx_hit,3),'CDF',pd_1);
[tmp_2,tmp_2_1] = chi2gof(RT(indx_hit,4),'CDF',pd_2);
[tmp_3,tmp_3_1] = chi2gof(tmp_data,'CDF',pd_3);

[tmp_1 tmp_1_1; tmp_2 tmp_2_1;tmp_3 tmp_3_1]
% ----- plot RT distribution

figure;
h1 = histogram(RT(indx_hit,3),25);
ylabel('COUNT');
xlabel('REACTION TIME [in msec]');
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
set(gcf, 'PaperPosition', [0 0 4 3]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/hist_rt.fig']);
print(['./figures/hist_rt.eps'],'-depsc','-r300');    
end


figure;
h2 = histogram(RT(indx_hit,4),25);
ylabel('COUNT');
xlabel('LOG REACTION TIME [in msec]');
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
set(gcf, 'PaperPosition', [0 0 4 3]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/hist_log_rt.fig']);
print(['./figures/hist_log_rt.eps'],'-depsc','-r300');    
end

% plot RT normplot
figure;figure; h = normplot([RT(indx_hit,3)]);
h(1).MarkerSize = MSval;
xlabel('REACTION TIME [in msec]');
title('')
grid on;
set(gca,FS,FSval-2,'box','on');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      );

if 1
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4 5]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/normplot_rt.fig']);
print(['./figures/normplot_rt.eps'],'-depsc','-r300');    
end


figure; h = normplot([RT(indx_hit,4)]);
h(1).MarkerSize = MSval;
xlabel('LOG REACTION TIME [in msec]');
title('')
grid on;
set(gca,FS,FSval-2,'box','on');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      );

if 1
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4 5]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/normplot_log_rt.fig']);
print(['./figures/normplot_log_rt.eps'],'-depsc','-r300');    
end

end

end
