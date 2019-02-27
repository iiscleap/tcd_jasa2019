


clear all;
%close all;

feats_path{1} = './data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part/mfcc/';

flist = dir([feats_path{1} '*.htk']);
nfiles = length(flist);

data_path = './data/rt_feats/';
% data_std = readtable([data_path 'subject_wise_mean_rt.csv'],'Delimiter',',');
data_std = readtable([data_path 'subject_wise_noise_mean_rt.csv'],'Delimiter',',');
data_words = readtable([data_path 'file_wise_tChange_nWords.csv'],'Delimiter',',');

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
    RT{loop} = CONST*ones(length(indx_subj{loop}),7);
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
        
        % gets number of words before tChange
        cwords = 0;
        str_1 = 'rstamp_';
        indx_1 = strfind(fname,str_1);
        for j = 1:length(data_words.CWORDS)
            if length(strfind(data_words.FNAME{j},fname(1:indx_1-2)))
                cwords = data_words.CWORDS(j);
                wdur = data_words.CWDUR(j);
                wrate = data_words.CWRATE(j);
            end
        end
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
        
        RT{loop}(i,:) = [tChange rChange 0 0 cwords wdur wrate];
        if (((label{loop}(i) == 3) || (label{loop}(i) == 2) || (label{loop}(i) == 1) || (label{loop}(i) == 0.5)) && rChange>0 && tChange<7800)
            RT{loop}(i,:) = [tChange rChange 0 0 cwords wdur wrate]; 
            FNAME{i} = fname;
            if (label{loop}(i) == 3)
                indx_hit{loop}(i) = 1;
                RT{loop}(i,3) = (RT{loop}(i,2)-RT{loop}(i,1));% - (data_std.meanRT_ms(loop_cnt));
                RT{loop}(i,4) = log10(RT{loop}(i,2)-RT{loop}(i,1)) - log10(data_std.meanRT_ms(loop_cnt));
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

%return;

do_fig_save = 0;

% ----- plot change instant versus subjectwise reaction time 
if 1    
    
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 1.5;


str{1} = 'dur';
str{5} = 'words';
str{6} = 'wdur';
str{7} = 'word_rate';

x_label{1} = 'SIGNAL DURATION BEFORE CHANGE [in msec]';
x_label{5} = 'WORDS BEFORE CHANGE';
x_label{6} = 'SPEECH DURATION BEFORE CHANGE [in msec]';
x_label{7} = 'WORDS-PER-SEC BEFORE CHANGE';
y_label{1} = 'log_{10} RT [in msec]';

for l = [1 5 6 7]

    C = zeros(size(indx_hit{loop}));
    Z = zeros(size(indx_hit{loop}));
    Z_std = zeros(size(indx_hit{loop}));
    TEMP = zeros(length(ptcpnt_id),length(indx_hit{loop}));
    for loop = sel_subjs

        X1 = CONST*ones(size(indx_hit{loop}));
        indx = find(indx_hit{loop}<CONST)';
        X1(indx) = indx_hit{loop}(indx);
        indx = [(find(indx_fa{loop}<CONST)+1) (find(indx_junk_fa{loop}<CONST)+1)];
        X1(indx) = CONST;
        indx = find(X1<CONST)';
    %     C(indx) = C(indx)+1;
    %     Z(indx) = Z(indx) + (RT{loop}(indx,3))'; 
        TEMP(loop,indx) = (RT{loop}(indx,4))';
    end
    
    % pool all, prune unique, and find mean std
    x = (RT{loop}(:,l));
    uniq_x = unique(x);
    tmp = cell(length(uniq_x),1);
    z95 = 1.96;
    for i = 1:length(uniq_x)
        tmp{i} = [];
        tmp_5 = [];
        tmp_5_std = [];
        for j = 1:length(x)
            if x(j) == uniq_x(i)
                tmp{i} = [tmp{i}; TEMP(:,j)];
            end
        end
    end
    
    % prune the indices which did not get any hit
    cnt = 1;
    mu = [];
    mu_std = [];
    uniq_y = [];
    for i = 1:length(tmp)
       indx = find(tmp{i}>0);
       if length(indx)~=0 
           uniq_y(cnt) = uniq_x(i);
           tmp_1 = tmp{i}(indx);
           mu(cnt) = mean(tmp_1);
           mu_std(cnt) = z95*std(tmp_1)/sqrt(length(tmp_1));
           cnt = cnt+1;
       end
    end
    [val indx] = sort(uniq_y);
    mu = mu(indx);
    mu_std = mu_std(indx);
    uniq_x = uniq_x(indx);

    tmp_3 = uniq_x;
    tmp_4 = mu';
    tmp_4_std = mu_std';
    
    figure;
    plot(tmp_3,tmp_4,'-o','color',[0.6 0.6 0.6],LW,LWval,'MarkerSize',MSval,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6]);
    hold on;
%     plot(tmp_3,tmp_4-tmp_4_std,'--');
%     plot(tmp_3,tmp_4+tmp_4_std,'--');

    % fit linear
    X = [ones(size(tmp_3(:))) tmp_3(:)];
    Y = tmp_4;
    b = pinv(X)*Y;
    Y1 = X*b;
    mdl = fitlm(X(:,2),Y);
    p_value = mdl.Coefficients.pValue(2:end)
    r_square = mdl.Rsquared.Ordinary
    plot(tmp_3,Y1,'-o','color','k',LW,LWval,'MarkerSize',MSval,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on;
    text(min(X(:,2))+0.25*(max(X(:,2))-min(X(:,2))),min(Y)+0.05*(max(Y)-min(Y)),['R-SQUARE = ' num2str(r_square)],'FontSize',6);
    text(min(X(:,2))+0.25*(max(X(:,2))-min(X(:,2))),min(Y)+0.05*(max(Y)-min(Y))+.02,['P-VALUE = ' num2str(p_value)],'FontSize',6);hold on;

    X = [ones(size(tmp_3(:))) tmp_3(:)];
    Y = tmp_4_std;
    b = pinv(X)*Y;
    Y2 = X*b;
%     plot(tmp_3,Y1+Y2,'--','color','k'); hold on;
%     plot(tmp_3,Y1-Y2,'--','color','k'); hold on;

    ylabel(y_label{1});
    xlabel(x_label{l});
    axis tight
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
        set(gcf, 'PaperPosition', [0 0 7 5]); %x_width=10cm y_width=sel_nharmcm        
        saveas(gcf,['./figures/rt_fit_' str{l} '.fig']);
        print(['./figures/rt_fit_' str{l} '.eps'],'-depsc','-r300');    
    end           
end
end
