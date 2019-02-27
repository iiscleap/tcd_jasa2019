

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


flist = dir([feats_path{1} '*.htk']);
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
RT = cell(length(ptcpnt_id),1);
indx_hit = cell(length(ptcpnt_id),1);
indx_fa = cell(length(ptcpnt_id),1);
indx_miss = cell(length(ptcpnt_id),1);
indx_junk_fa = cell(length(ptcpnt_id),1);
flag_same_spk = cell(length(ptcpnt_id),1);
flag_resp = cell(length(ptcpnt_id),1);
label = cell(length(ptcpnt_id),1);

store_path = './data/rt_feats/';

cnt = 1;
loop_cnt = 0;
CONST = 10000;
for loop = sel_subjs
    loop_cnt = loop_cnt+1; 
    rt{loop} = zeros(1,length(indx_subj{loop}));
    RT{loop} = CONST*ones(length(indx_subj{loop}),4);
    flag_same_spk{loop} = zeros(1,length(indx_subj{loop}));
    flag_resp{loop} = zeros(1,length(indx_subj{loop}));
    label{loop} = zeros(1,length(indx_subj{loop}));
    indx_hit{loop} = CONST*ones(1,length(indx_subj{loop}));
    indx_miss{loop} = CONST*ones(1,length(indx_subj{loop}));
    indx_fa{loop} = CONST*ones(1,length(indx_subj{loop}));
    indx_junk_fa{loop} = CONST*ones(1,length(indx_subj{loop}));
    
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
        elseif (temp>2000 && ~flag_resp{loop}(i) && ~flag_same_spk{loop}(i))
            label{loop}(i) = 2; % miss
        elseif ((temp<0 && ~flag_resp{loop}(i) && ~flag_same_spk{loop}(i)))
            label{loop}(i) = 1; % false alarm
        elseif ((temp<225 && temp>0 && ~flag_resp{loop}(i) && ~flag_same_spk{loop}(i)))
            label{loop}(i) = 0.5; % junk false alarm
        else
            label{loop}(i) = 0; % junk
        end
                
        if (((label{loop}(i) == 3) || (label{loop}(i) == 2) || (label{loop}(i) == 1) || (label{loop}(i) == 0.5)) && rChange>0 && tChange<7000)
            RT{loop}(i,:) = [tChange rChange 0 0]; 
            FNAME{i} = fname;
            if (label{loop}(i) == 3)
                indx_hit{loop}(i) = 1;
                RT{loop}(i,3) = (RT{loop}(i,2)-RT{loop}(i,1));% - (data_std.meanRT_ms(loop_cnt));
                RT{loop}(i,4) = log10(RT{loop}(i,2)-RT{loop}(i,1));% - log10(data_std.meanRT_ms(loop_cnt));
%                 RT(cnt,3) = (RT(cnt,2)-RT(cnt,1)) - (data_std.medRT_ms(loop_cnt));
%                 RT(cnt,4) = log10(RT(cnt,2)-RT(cnt,1)) - log10(data_std.medRT_ms(loop_cnt));
            elseif (label{loop}(i) == 2)
                indx_miss{loop}(i) = 1;
            elseif (label{loop}(i) == 1)
                indx_fa{loop}(i) = 1;
                ['FA ' FNAME{i}]
            else
                indx_junk_fa{loop}(i) = 1;
                ['JFA ' FNAME{i}]
            end
            cnt = cnt+1;
        end
        end
    end
end

do_fig_save = 0;

if 1
close all;
% ----- make plots
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 3;

clr = rand(length(ptcpnt_id),3);
% clr = repmat([0.5 0.5 0.5],length(ptcpnt_id),1);


% ----- plot change instant versus subjectwise reaction time 
if 1
X2 = log10(linspace(5000,7000,2000)');
Y2 = zeros(size(X2));
X2 = [ones(length(X2),1) X2];
for loop = sel_subjs
    indx = [(find(indx_fa{loop}<CONST)+1) (find(indx_junk_fa{loop}<CONST)+1)];
    indx_hit{loop}(indx) = CONST;
    indx = find(indx_hit{loop}<CONST)';
    tmp_1 = log10(RT{loop}(indx,1));
    tmp_2 =log10(RT{loop}(indx,3));
    [tmp_1,indx] = sort(tmp_1);
    tmp_2 = tmp_2(indx);
    
%     tmp_2 = movmean(tmp_2,5);
    plot(tmp_1,tmp_2,'o','color',clr(loop,:),'MarkerSize',MSval,'MarkerFaceColor',clr(loop,:));
    hold on;
    
    indx = find(tmp_1<10000);
    tmp_1 = tmp_1(indx);
    tmp_2 = tmp_2(indx);
    Y = tmp_2;
    X = [ones(length(tmp_1),1) tmp_1];
    b = pinv(X)*Y;
    Y1 = X*b;
    Y2 = Y2+X2*b;
    b(2)
    plot(tmp_1,Y1,'-o','color',clr(loop,:),'MarkerSize',MSval,'MarkerFaceColor',clr(loop,:));
end
hold on;
plot(X2(:,2),Y2/length(sel_subjs),'-o','color','k','MarkerSize',MSval,'MarkerFaceColor','k');
end

% ----- plot change instant versus pooled reaction time 
if 0
tmp_1 = [];
tmp_2 = [];

for loop = sel_subjs
    indx = find(indx_hit{loop}<CONST)';
    tmp_1 = [tmp_1; RT{loop}(indx,1)];
    tmp_2 = [tmp_2; log10(RT{loop}(indx,3))];
end

figure;
[tmp_1,indx] = sort(tmp_1);
tmp_2 = tmp_2(indx);
plot(tmp_1,tmp_2,'-o','color',clr(loop,:),'MarkerSize',MSval,'MarkerFaceColor',clr(loop,:));
hold on;
    
indx = find(tmp_1<10000);
tmp_1 = tmp_1(indx);
tmp_2 = tmp_2(indx);
Y = tmp_2;
X = [ones(length(tmp_1),1) tmp_1];
b = pinv(X)*Y;
Y1 = X*b;
b(2)
plot(tmp_1,Y1,'-o','color',clr(loop,:),'MarkerSize',MSval,'MarkerFaceColor',clr(loop,:));
end

% ----- plot change instant versus false alarm count
if 0
    figure;
for loop = sel_subjs
    indx = find(indx_fa{loop}<CONST)';
    tmp_1 = (RT{loop}(indx,1));
    tmp_2 = indx_fa{loop}(indx);
    [tmp_1,indx] = sort(tmp_1);
    tmp_2 = tmp_2(indx);
    tmp_2 = cumsum(tmp_2);
    
%     tmp_2 = movmean(tmp_2,5);
    plot(tmp_1,tmp_2,'-o','color',clr(loop,:),'MarkerSize',MSval,'MarkerFaceColor',clr(loop,:));
    hold on;
end
end


return;
if 1
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

if do_fig_save
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

if do_fig_save
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4 3]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/hist_log_rt.fig']);
print(['./figures/hist_log_rt.eps'],'-depsc','-r300');    
end
end

end