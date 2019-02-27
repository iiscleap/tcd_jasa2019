

% analyze reaction time data

clearvars;
%close all;
%clc;
% path = './data/libriSpeech/spkrChange/listTest/M_5_spkrs_8_pairs_v0/listTest_collectedData/19Nov2017/';
% path = './data/libriSpeech/spkrChange/listTest/M_5_spkrs_8_pairs/listTest_collectedData/18Nov2017/';
% path = './data/libriSpeech/spkrChange/listTest/M_5_spkrs_8_pairs_v0/listTest_collectedData/from_29Jan_2018/';
addpath('../plotting/')
path = './data/libriSpeech/spkrChange/listTest/dataSheet/';

data = readtable([path 'accum_Nov_2017_19Jan_2018_17Apr_2018_mod.csv']);


ptcpnt_id = {'CMU-1','cmu2','CMU_3','CMU-4','cmu_5','CMU_6','CMV-7','CMU_8','cmu_9','CMU_10','cmu_11','cmu12',...
    'CMU_13','CMU-14','cmu_15','CMU_16','CMU_17','iisc_01','iisc_02','iisc_03','iisc_04','iisc_05','iisc_06',...
    'iisc_07','iisc_08','iisc_09','iisc_10','iisc_11','iisc_12'};

rt = cell(length(ptcpnt_id),1);
indx_Delta = cell(length(ptcpnt_id),1);

% for indx_pids = [1:3]%:length(ptcpnt_id)
for indx_pids = [13:29]%:length(ptcpnt_id)
    cnt = 1;
    flag = 0;
    for i = 1:length(data.ParticipantPublicID)
    if strcmp(data.ParticipantPublicID{i},ptcpnt_id{indx_pids})
        indx_ptcpnt_id(cnt) = i;
        cnt = cnt+1;
    end
    end
    
    display(['Nos: ' num2str(length(indx_ptcpnt_id))])

    cnt_break = 1;
    for i = 1:length(indx_ptcpnt_id)
    if length(strfind(data.display{indx_ptcpnt_id(i)},'rt_noise'))
       indx_display_1 = indx_ptcpnt_id(i);
    end  
    if length(strfind(data.display{indx_ptcpnt_id(i)},'rt_tone'))
       indx_display_2 = indx_ptcpnt_id(i);
    end
    if length(strfind(data.display{indx_ptcpnt_id(i)},'rt_speech_train'))
       if ~flag
       indx_display_3 = indx_ptcpnt_id(i);
       flag = 1;
       end
    end
    if length(strfind(data.display{indx_ptcpnt_id(i)},'break'))
       indx_break(cnt_break) = indx_ptcpnt_id(i);
       cnt_break = cnt_break +1;
    end
    
    if length(strfind(data.display{indx_ptcpnt_id(i)},'rt_speech_test'))
       indx_display_4 = indx_ptcpnt_id(i);
    end
    end
    
    indx_display_5 = indx_ptcpnt_id(end);
%     return;
    
    % ----- obtain rt_noise and rt_tone results
    if 0
    cnt = 1;
    cnt_1 = 1;
    for i = indx_display_1+1:indx_display_2-1
       if length(strfind(data.ZoneType{i},'response_keyboard_single'))
           r_tstamp_pretest{1}(cnt) =  str2double(data.ReactionTime(i));
           pos = strfind(data.sound{i},'_')+1;
           a_tstamp{1}(cnt) = str2double(data.sound{i}(pos(end):(end-4)));
           rt_temp = r_tstamp_pretest{1}(cnt) - a_tstamp{1}(cnt);
           if ((rt_temp > 225) & (rt_temp < 2000))
               rt{indx_pids}{1}(cnt_1) =rt_temp;
               cnt_1 = cnt_1+1;
           end
           cnt = cnt+1;
       end
    end

    cnt = 1;
    Delta = [10 25 50 100];
    nos = 5*3;
    indx_Delta{indx_pids} = zeros(length(Delta),nos);
    cnt_Delta = ones(length(Delta),1);
    for i = indx_display_2+1:indx_display_3-1
        if length(strfind(data.ZoneType{i},'response_keyboard_single'))
            r_tstamp_pretest{2}(cnt) = str2double(data.ReactionTime(i));
            pos = strfind(data.sound{i},'_')+1;
            a_tstamp{2}(cnt) = str2double(data.sound{i}(pos(end):(end-4)));
            rt_temp = r_tstamp_pretest{2}(cnt) - a_tstamp{2}(cnt);
            if ((rt_temp > 225) & (rt_temp < 2000))
                rt{indx_pids}{2}(cnt) = rt_temp;
                for j = 1:length(Delta)
                    if length(strfind(data.sound{i},['Delta_' num2str(Delta(j)) '_']))
                    indx_Delta{indx_pids}(j,cnt_Delta(j)) = cnt;
                    cnt_Delta(j) = cnt_Delta(j)+1;
                    end
                end
                cnt = cnt+1;
            end
       end
    end
    end
    
    if 1
    % ----- obtain rt_speech_test 
    indx_test = setdiff(indx_display_4+1:indx_display_5,indx_break);
    
    r_tstamp = str2double(data.ReactionTime(indx_test));

    file_id = data.sound(indx_test);
    zone_id = data.ZoneName(indx_test);
    corr_id = str2double(data.Correct(indx_test));
    gnd_corr_id = str2double(data.answer(indx_test));

    talker_str = zeros(length(file_id),1);
    for i = 1:length(file_id)
    npos = strfind(file_id{i},'_');
    talker_str(i) = str2double(file_id{i}((npos(1)+1):(npos(2)-1)));
    end
    talker_str = unique(talker_str);
    
    talker_str = sort(talker_str);
    
    for j = 1:length(talker_str)
        for k = 1:length(talker_str)
            cnt = 1;
            for l = 1:length(file_id)
                temp = file_id{l};
                temp(1:2) = '--';
                if length(strfind(temp, ['--_' num2str(talker_str(j))]))
                    if length(strfind(temp, ['SWT_LS_' num2str(talker_str(k))]))
%                         display(temp)
                        if ~(length(strfind(zone_id{l},'fixation')))
                            t1_t2_indices{j,k}(cnt) = l;
                            t1_t2_corr{j,k}(cnt) = corr_id(l); % 1/0
                            t1_t2_gnd_corr{j,k}(cnt) = gnd_corr_id(l);
                            
                            t1_t2_tstamp{j,k}(cnt) = r_tstamp(l);
                            if length(strfind(zone_id{l},'Zone2'))
                                t1_t2_zone{j,k}(cnt) = 0; % no key pressed
                            else
                                t1_t2_zone{j,k}(cnt) = 1; % key pressed
                            end
                            
                            npos = strfind(temp,'_');
                            t1_t2_gtstamp{j,k}(cnt) = ...
                            str2double(temp((npos(end-1)+1):(npos(end)-1)));
                            cnt = cnt +1;
%                             t1_t2_zone{j,k}
%                             t1_t2_corr{j,k}
                        end
                    end                   
                end
            end
        end
    end
   
    
    % measure reaction time
    pool_rt{indx_pids} = zeros(length(talker_str),length(talker_str)*8);
    
    CH_mu_rt{indx_pids} = cell(length(talker_str),length(talker_str));
    
    CH_hit{indx_pids} = zeros(length(talker_str),length(talker_str));
    CH_miss{indx_pids} = zeros(length(talker_str),length(talker_str));
    CH_fa{indx_pids} = zeros(length(talker_str),length(talker_str));
    CH_cr{indx_pids} = zeros(length(talker_str),length(talker_str));
    tot_stimuli_hit{indx_pids} = zeros(length(talker_str),length(talker_str));
    tot_stimuli_fa{indx_pids} = zeros(length(talker_str),length(talker_str));
    for j = 1:length(talker_str)
        for k = 1:length(talker_str)
            fa_miss_cnt = 0;
            hit = 0; miss = 0; fa = 0; cr =0;
            fa_cnt = 0;
            false_miss_cnt = 0;
            cnt = 1;
            for l = 1:length(t1_t2_zone{j,k})
               rt_speech = (t1_t2_tstamp{j,k}(l) - t1_t2_gtstamp{j,k}(l));
               
               if (t1_t2_gnd_corr{j,k}(l)==1) % there is a change in ground truth
                   if (t1_t2_zone{j,k}(l)==0) % no key press
                       miss = miss+1;
                   else
                       if rt_speech<225 && rt_speech>0
                           false_miss_cnt = false_miss_cnt + 1;
                           fa_cnt = fa_cnt + 1;
                           fa = fa+1;
                       elseif rt_speech<0
                           fa_cnt = fa_cnt + 1;
                           fa = fa+1;
                       elseif rt_speech>225
                           if rt_speech > 2000
                               miss = miss+1;
                           else
                               cr = cr+1;
                               hit = hit+1;
                               CH_mu_rt{indx_pids}{j,k}(cnt) = rt_speech;
                               cnt = cnt+1;
                           end
                       end
                   end
               else % there is no change in ground truth
                   fa_cnt = fa_cnt+1;
                   if (t1_t2_zone{j,k}(l)==0) % no key press
                       cr = cr+1;
                   else
                       fa = fa+1;
                   end
               end
%                 [hit miss fa cr]
            end
%             tot_stimuli_hit{indx_pids}(j,k) = 8-false_miss_cnt-fa_cnt;
            tot_stimuli_hit{indx_pids}(j,k) = 8-fa_cnt;
            tot_stimuli_fa{indx_pids}(j,k) = fa_cnt;
            if j~=k
                CH_hit{indx_pids}(j,k) = hit;
                CH_miss{indx_pids}(j,k) = miss;
                CH_fa{indx_pids}(j,k) = fa;
                CH_cr{indx_pids}(j,k) = cr;                
            else
                CH_fa{indx_pids}(j,k) = fa;
                CH_cr{indx_pids}(j,k) = cr;
            end
        end
%                        return;
    end
    end
end

return;
% ----- make plots
close all;

FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 4;

% ----- make plots for subjectwise MISS, FA, HIT 
pids = [13:17 18:29];
loop = 1;
mu_sub_miss = [];
mu_sub_hit = [];
mu_sub_fa = [];
for indx_pids = pids
    % get all non-diagonal entries
    mu_sub_hit(loop) = sum(CH_hit{indx_pids}(:))/sum(tot_stimuli_hit{indx_pids}(:));
    mu_sub_miss(loop) = sum(CH_miss{indx_pids}(:))/sum(tot_stimuli_hit{indx_pids}(:));
    mu_sub_fa(loop) = sum(CH_fa{indx_pids}(:))/sum(tot_stimuli_hit{indx_pids}(:)+tot_stimuli_fa{indx_pids}(:));
    loop = loop+1;
end

% ---- store subjectwise hit/miss/fa CSV file
if 1
store_path = './data/results/';
fileID = fopen([store_path 'human_subjectwise_hit_miss_fa.csv'],'w');
fprintf(fileID,'SID,HIT,MISS,FA\n');
for i = 1:length(pids)
fprintf(fileID,'%s,%.4f,%.4f,%.4f\n',ptcpnt_id{pids(i)},mu_sub_hit(i),mu_sub_miss(i),mu_sub_fa(i));
end
fclose(fileID);
end

if 1
store_path = './data/results/';
fileID = fopen([store_path 'human_mean_hit_miss_fa.csv'],'w');
fprintf(fileID,'MU_HIT,MU_MISS,MU_FA\n');
fprintf(fileID,'%.4f,%.4f,%.4f\n',mean(mu_sub_hit),mean(mu_sub_miss),mean(mu_sub_fa));
fclose(fileID);
end

% ----- plot mean hit/ miss/ fa/ CSV file

figure;
h(1) = plot(1:17,mu_sub_hit*100,'o-.','color','k',LW,LWval,'MarkerSize',MSval,'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor',[0 0.5 0]);
hold on;
h(2) = plot(1:17,mu_sub_miss*100,'o-.','color','k',LW,LWval,'MarkerSize',MSval,'MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[0 0.5 1]);
h(3) = plot(1:17,mu_sub_fa*100,'o-.','color','k',LW,LWval,'MarkerSize',MSval,'MarkerFaceColor','r','MarkerEdgeColor','r');
grid on;
ylabel('PERCENTAGE');
xlabel('SUBJECT INDEX');
xticks([1:17])
xticklabels({'1','2','3','4','5','6','7','8','9',...
             '10','11','12','13','14','15','16','17'});
yticks([0:5:20 85:5:100])
xlim([0 18]);
ylim([0 103]);
grid on;
% set(gca,FS,FSval,'box','on');
set(gca,FS,12,'box','on');
hlegend = legend(h,{'HIT','MISS','FALSE ALRAM'});
rect = [0.6 0.5 .25 .25]; %[left bottom width height]
set(hlegend,'Position',rect);
breakyaxis([20 80]);
legend('show')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );

if 1 % but i saved it offline to work best
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 9 6]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/subjectwise_hit_miss_fa.fig']);
    print(['./figures/subjectwise_hit_miss_fa.eps'],'-depsc','-r300');    
end   
return;
% ----- get the speaker pair-wise average RT matrix
indx_pids = 13:29;
mu_rt = cell(length(indx_pids),1);
for loop = indx_pids
    for i = 1:5
        for j = 1:5
            if i~=j
               mu_rt{loop}(i,j) = sum(CH_mu_rt{loop}{i,j})/tot_stimuli_hit{loop}(i,j);
            end
        end
    end
end

MU_RT = zeros(5,5);
for i = indx_pids
   MU_RT = MU_RT + mu_rt{i}; 
end
MU_RT = MU_RT/length(indx_pids);

% get the speaker wise false alarm matrix
indx_pids = 13:29;
mu_fa = cell(length(indx_pids),1);
for loop = indx_pids
    for i = 1:5
        mu_fa{loop}(i,1) = sum(CH_fa{loop}(i,:))/sum(tot_stimuli_hit{loop}(i,:)+tot_stimuli_fa{loop}(i,:));
    end
end

MU_FA = zeros(5,1);
for i = indx_pids
   MU_FA = MU_FA + mu_fa{i}; 
end
MU_FA = MU_FA/length(indx_pids);
% ----- store speaker-wise fa mat file
if 1
store_path = './data/results/';
fileID = fopen([store_path 'human_talkerwise_fa.csv'],'w');
fprintf(fileID,'Talker_ID, MU_FA\n');
for i = 1:length(talker_str)
    fprintf(fileID,'%d,%.4f\n',talker_str(i), MU_FA(i));
end
fclose(fileID);
end



% get the speaker wise miss matrix
indx_pids = 13:29;
mu_miss = cell(length(indx_pids),1);
for loop = indx_pids
    for i = 1:5
        for j = 1:5
            if i~=j
                mu_miss{loop}(i,j) = CH_miss{loop}(i,j)/tot_stimuli_hit{loop}(i,j);
            end
        end
    end
end

MU_MISS = zeros(5,5);
for i = indx_pids
   MU_MISS = MU_MISS + mu_miss{i}; 
end
MU_MISS = MU_MISS/length(indx_pids);


% ----- plot speaker wise FA
if 1
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 3;
figure;
stem(1:5,MU_FA*100,'o','color',[0.6 0.6 0.6],'MarkerSize',MSval,'MarkerFaceColor','r','MarkerEdgeColor','r');
xlim([0 6]);
xticks(1:1:5)
xticklabels({'T1','T2','T3','T4','T5'});
grid on;
set(gca,FS,FSval,'box','on');
xlabel('TALKER INDEX');
ylabel('FALSE ALARM %')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );

if 1
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 5 5]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/talker_wise_fa.fig']);
    print(['./figures/talker_wise_fa.eps'],'-depsc','-r300');    
end           

end
