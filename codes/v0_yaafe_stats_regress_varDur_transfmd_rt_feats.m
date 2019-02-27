
clearvars;
close all;

feats_type = {'mfcc', 'mfcc_d1','mfcc_d2','engy','envp','lsf','loudness','mel','psharp','pspread', ...
         'sflat','sflux','sroll','sshape','sslope','tshape','zcr','f0'};
     
feats_path = './data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part/';

flist = dir([feats_path feats_type{1} '/' '*.h5']);
nfiles = length(flist);

ptcpnt_id = {'CMU-1','cmu2','CMU_3','CMU-4','cmu_5','CMU_6','CMV-7','CMU_8','cmu_9','CMU_10','cmu_11','cmu12',...
    'CMU_13','CMU-14','cmu_15','CMU_16','CMU_17','iisc_01','iisc_02','iisc_03','iisc_04','iisc_05','iisc_06',...
    'iisc_07','iisc_08','iisc_09','iisc_10','iisc_11','iisc_12'};

sel_subjs =  [13:17 18:25 26:29];

data_path = './data/rt_feats/';
% data_std = readtable([data_path 'subject_wise_mean_rt.csv'],'Delimiter',',');
data_std = readtable([data_path 'subject_wise_noise_mean_rt.csv'],'Delimiter',',');

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
feats_dist = cell(length(ptcpnt_id),4);
flag_same_spk = cell(length(ptcpnt_id),1);
flag_resp = cell(length(ptcpnt_id),1);
label = cell(length(ptcpnt_id),1);
fname_id = cell(length(ptcpnt_id),1);

feats_type_merge = {'F0','lsf','mel','mfcc','mfcc_{d1}','mfcc_{d2}',...
    'temp','percp','spect','ivec','corti'};
feats = cell(length(feats_type_merge),1);

cnt = 1;

for loop = sel_subjs
    rt{loop} = zeros(1,length(indx_subj{loop}));
    flag_same_spk{loop} = zeros(1,length(indx_subj{loop}));
    flag_resp{loop} = zeros(1,length(indx_subj{loop}));
    label{loop} = zeros(1,length(indx_subj{loop}));
    feats_dist{loop} = zeros(length(indx_subj{loop}),length(feats_path));
    fname_id{loop} = cell(length(indx_subj{loop}),1);
    for i = 1:length(indx_subj{loop})
        [fpath,fname,fext] = fileparts([feats_path feats_type{1} '/' flist(indx_subj{loop}(i)).name]);
        fname_id{loop}{i} = fname;
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
        elseif ((~flag_resp{loop}(i) && flag_same_spk{loop}(i)) || ...
                (temp<225 && ~flag_resp{loop}(i) && ~flag_same_spk{loop}(i)))
            label{loop}(i) = 1; % false alarm
        else
            label{loop}(i) = 0; % junk
        end
        
        feats_set = {{'f0'},{'lsf'},{'mel'},{'mfcc'},{'mfcc_d1'},{'mfcc_d2'},...
            {'engy'},{'loudness','psharp','pspread'},...
            {'sflat','sflux','sroll','sshape','sslope'},...
            {'ivector'},{'corti'}};
        
        % ----- read feats
        for j = 1:length(feats_set)
            feats{j} = [];
            if j == 10 % ivector ivec_ac_segwise
                temp = hdf5read([feats_path feats_set{j}{1} '/' fname],'ivec_ac_segwise');
                feats{j} = [feats{j};temp];
                feats{j} = feats{j}.';
            elseif j == 11
                temp = hdf5read([feats_path feats_set{j}{1} '/' fname],'corticogram');
                feats{j} = [feats{j};temp];
            else
                for k = 1:length(feats_set{j})
                    if j == 7 % take rate of change in energy as a feature
                        if k == 1
                        temp = hdf5read([feats_path feats_set{j}{k} '/' fname],[feats_set{j}{k} '_framewise']);
                        feats{j} = [feats{j};[diff(temp) 0]];
                        else
                        temp = hdf5read([feats_path feats_set{j}{k} '/' fname],[feats_set{j}{k} '_framewise']);
                        feats{j} = [feats{j};temp];
                        end                        
                    else
                    temp = hdf5read([feats_path feats_set{j}{k} '/' fname],[feats_set{j}{k} '_framewise']);
                    feats{j} = [feats{j};temp];
                    end
                end
                if k>1
                    feats{j} = (feats{j}-mean(feats{j},2))./sqrt(var(feats{j},0,2));
                end
                if j == 1
                    feats{j} = feats{j}.';
                end
%                 if j == 3 % uncomment to take first derivative of MEL (shows little improvement)
%                     feats{j} = diff(feats{j},1,2);
%                 end                
            end
        end
        % ----- compute distance fraction wise
        frac = [25 50 75 100];
        for j = 1:length(feats)
            % get voiced indices
            indx_unvoiced = [];%find(feats{1}<10);
            flag_all_frames = 1;
            if label{loop}(i) == 3
                for l = 1:4
                    if j == 10 % read ivector
                        str_2 = '_rstamp';
                        indx_2 = strfind(fname,str_2)-1;
                        temp = hdf5read([feats_path feats_set{j}{1} '/' fname(1:indx_2)],'ivec_bc_segwise');
                        temp = temp.';
%                         feats_dist{loop,l}(i,j) = norm(temp(:,l)/sqrt(norm(temp(:,l)))-feats{j}(:,2)/sqrt(norm(feats{j}(:,2))));
                        val = pdist([temp(:,l).'; feats{j}(:,2).'],'cosine');
                        feats_dist{loop,l}(i,j) = val;
                    elseif j == 11
                        str_2 = '_ms_rstamp';
                        indx_2 = strfind(fname,str_2)-1;

                        feats_3_bc = hdf5read([feats_path feats_set{j}{1} '/' fname(1:indx_2) '_ms_percent_bc_' num2str(frac(l))],'corticogram');

                        % make a matrix for before change data
                        temp_1 = zeros(size(feats_3_bc,1)*size(feats_3_bc,2),size(feats_3_bc,3));
                        for k = 1:size(feats_3_bc,3)
                            temp_0 = feats_3_bc(:,:,k);
                            temp_1(:,k) = temp_0(:);
                        end
                        
                        % make a matrix for after change data
                        temp_2 = zeros(size(feats{j},1)*size(feats{j},2),size(feats{j},3));
                        for k = 1:size(feats{j},3)
                            temp_0 = feats{j}(:,:,k);
                            temp_2(:,k) = temp_0(:);
                        end
                        mu = temp_1;
                        centroid = mean(mu,2);
                        
                        y = temp_2;
                        mu = mean(y,2);
                        feats_dist{loop,l}(i,j) = norm(centroid-mu);
                    else
                        clen = fix(tChange/10)-10;
                        slen = fix((l/4)*clen);
                        rlen = fix(rChange/10);

                        % zero out the unvoiced regions
                        feats{j}(:,indx_unvoiced) = 0;

                        mu = feats{j}(:,clen-slen+1:clen);
                        if ~flag_all_frames
                            tmp = diag(mu'*mu);
                            thr = .001*mean(tmp);
                            tmp_indx = find(diag(mu'*mu)>thr);
                            if length(tmp_indx)==0
                                centroid = zeros(size(mu,1),1);
                            else
                                centroid = mean(mu(:,tmp_indx),2);
                            end
                        else
                            centroid = mean(mu,2);
                        end
                        y = feats{j}(:,clen+1:rlen);
                        if ~flag_all_frames
                            tmp = diag(y'*y);
                            thr = .001*mean(tmp);
                            tmp_indx = find(diag(y'*y)>thr);
                            if length(tmp_indx)==0
                                mu = zeros(size(y,1),1);
                            else
                                mu = mean(y(:,tmp_indx),2);
                            end
                        else
                            mu = mean(y,2);
                        end

                        feats_dist{loop,l}(i,j) = norm(centroid-mu);
                    end
                    if isnan(feats_dist{loop,l}(i,j))
                        disp('Check for NaN');
                        return;
                    end
               end
            else
               feats_dist{loop,l}(i,j) = 0;
            end
        end
    end
    for j = 1:length(feats)
        for k = 1:4
            feats_std = std(feats_dist{loop,k}((feats_dist{loop,k}(:,j)>0),j));
            feats_dist{loop,k}(:,j) = feats_dist{loop,k}(:,j)/feats_std;
        end
        
    end
end

return;



addpath('../plotting/cbrewer/')
close all;
if 1
% ----- make plots

feats_type = {'F0','LSF','MEL','MFCC','MFCC-D','MFCC-DD',...
    'ENGY-D','PLOUD','SPECT'};%,'I-VEC','CORTI'};
steps = 1000;
nfeats = length(feats);
r_square = zeros(nfeats+1,17,4);
p_value =  zeros(nfeats,17,4);
p_all = cell(4,1);
data_y_ip = cell(4,1);
data_y_op = cell(4,1);
fname_id_all = {};
data_indx_all = [];

for i = 1:4
    p_all{i} = zeros(17,nfeats+1);
    data_y_ip{i} = [];
    data_y_op{i} = [];
end
b = zeros(length(sel_subjs),nfeats,2);


for indx_feats = 1:nfeats+1
    cnt = 1;
    for indx_subj = sel_subjs
        % select indices with hit
        tmp_indx = find(label{indx_subj}==3);
        % pool filenames
        for i = 1:length(tmp_indx)
            if indx_feats == nfeats+1 
                data_indx_all = [data_indx_all tmp_indx(i)];
                fname_id_all{end+1} = fname_id{indx_subj}{tmp_indx};
            end
        end
        % loop on the modeling RT for different pre-change duration
        for indx_dur = 1:4
            if indx_feats<(nfeats+1)
                n = 2;
                 X = [ones(1,length(tmp_indx))' (feats_dist{indx_subj,indx_dur}(tmp_indx,indx_feats))];
%                  X = [ones(1,length(tmp_indx))' (feats_dist{indx_subj,indx_dur}(tmp_indx,indx_feats))];
            else
                n = nfeats+2;
                 X = [ones(1,length(tmp_indx))' (feats_dist{indx_subj,indx_dur}(tmp_indx,1:nfeats))];
%                 X = [ones(1,length(tmp_indx))' (feats_dist{indx_subj,indx_dur}(tmp_indx,1:nfeats))];
            end
        %         Y = sqrt(rt{indx_subj}(tmp_indx))';% - sqrt(data_std.meanRT_ms(cnt));
            Y = log10(rt{indx_subj}(tmp_indx))';% - log10(data_std.meanRT_ms(cnt));
%            Y = (rt{indx_subj}(tmp_indx))';% - log10(data_std.meanRT_ms(cnt));

            if 0 % regression using cvx version

                cvx_begin
                    variable w(n)
                    minimize( norm(X*w-Y,2)+0.0*norm(w,1))
                cvx_end
                Y1 = X*w;
                b(cnt,indx_feats,:) = w;
            elseif 0 % regression using LS
                if indx_feats<(nfeats+1)
                b(cnt,indx_feats,:) = pinv(X)*Y;
                Y1 = X*squeeze(b(cnt,indx_feats,:));
                else
                b_all = pinv(X)*Y;
                Y1 = X*b_all;
                end
                sst = sum((Y -mean(Y)).^2); % total variance
                ssr = sum((Y1-Y).^2); % residual error
                r_square(indx_feats,cnt,indx_dur) =1- ssr/sst;                
            else % regression using fitlm inbuilt function. it also uses LS.
                if indx_feats<(nfeats+1)
                    b(cnt,indx_feats,:) = pinv(X)*Y;
                    mdl = fitlm(X(:,2),Y);
                    Y1 = X*squeeze(b(cnt,indx_feats,:));
                    p_value(indx_feats,cnt,indx_dur) =  mdl.Coefficients.pValue(2);
                else
%                  2      3     4      5  6         7    8      9
%                 mfcc, fbank, corti, f0, meanAmp, engy, sfm, spectm    
                    indx_range = [2:length(feats_type)+1];
                    mdl = fitlm(X(:,indx_range),Y);
                    p_all{indx_dur}(cnt,indx_range) = mdl.Coefficients.pValue(2:end);
                    data_y_ip{indx_dur,1} = [data_y_ip{indx_dur,1};Y];
                    data_y_op{indx_dur,1} = [data_y_op{indx_dur,1};mdl.Fitted];
                 end    
                 r_square(indx_feats,cnt,indx_dur) = mdl.Rsquared.Ordinary;
            end
            
            if 0
                if (indx_dur == 4) && (indx_feats == 10)
                figure;
                plot(X(:,2),Y1,'o','color','k');
                hold on;
                plot(X(:,2),Y,'o','color','r');
                end
            end
        end
        cnt = cnt+1;
    end
end
end

FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 3;

% scatter plot of RT and Predicted RT from 'ALL' category
figure;
MSize = 20;
scatterhist(10.^(data_y_ip{4}),10.^(data_y_op{4}),'Marker','.','Color',[0 0.5 1],MS,MSval);
% colorbar;
% cmap = cbrewer('seq','Blues',100);
% colormap(flipud(cmap));
hold on;
% scatter(temp_x_0,temp_y_0,'Marker','.');
% do piecewise regression on partioned data
temp_x_0 = 10.^(data_y_ip{4});
sel_indx = find((temp_x_0<1000) & (temp_x_0>225));
temp_x_1 = temp_x_0(sel_indx);
temp_y_0 = 10.^(data_y_op{4});
temp_y_1 = temp_y_0(sel_indx);
hold on;
scatter_kde(temp_x_0,temp_y_0,'Marker','.');%,'DisplayStyle','tile','ShowEmptyBins','on','NumBins',[100 100]);

% mdl0 = fitlm(temp_x_0,temp_y_0);
% temp_y_0_pred = mdl1.Fitted;
% mdl1 = fitlm(temp_x_1,temp_y_1);
% temp_y_1_pred = mdl1.Fitted;

[temp_x_0, indx_sort] = sort(temp_x_0);
temp_y_0 = temp_y_0(indx_sort);
[temp_x_1, indx_sort] = sort(temp_x_1);
temp_y_1 = temp_y_1(indx_sort);


sst = sum((temp_y_0 - mean(temp_y_0)).^2); % total variance
ssr = sum((temp_x_0-temp_y_0).^2); % residual error
r_square_all_fit = 1 - ssr/sst;             
r_corr = corr(temp_x_0,temp_y_0)
r_corr_part = corr(temp_x_1,temp_y_1)
% plot(temp_x_0,temp_x_0,'-','color',[0.6 0.6 0.6],LW,LWval,'MarkerSize',0.25,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot([0 2000],[0 2000],'-','color',[0.6 0.6 0.6],LW,LWval,'MarkerSize',0.25,'MarkerFaceColor','k','MarkerEdgeColor','k');
% plot(temp_x_1,temp_y_1_pred,'-o','color','k',LW,LWval,'MarkerSize',0.25,'MarkerFaceColor','k','MarkerEdgeColor','k');
% plot(temp_x_2,temp_y_2_pred,'-o','color','k',LW,LWval,'MarkerSize',0.25,'MarkerFaceColor','k','MarkerEdgeColor','k');
% plot([200 2000],[200 2000],'-','color',[0.5 0.5 0.5],LW,LWval);
set(gca,FS,FSval-1,'box','on');
xlabel('RT [in msec]');
ylabel('PREDICTED RT [in msec]');
% xlim([200 2000])
% ylim([200 2000])
xlim([0 2000])
ylim([0 2000])
text(500,1800,['corr. coeff = ' num2str(r_corr,'%.2f')],'FontSize',FSval);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on',... 
  'XGrid'       , 'on',...
  'Yminorgrid'  , 'on',... 
  'Xminorgrid'  , 'on');

if 1
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 8 8]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/rt_pred_rt_scatter_hist.fig']);
    print(['./figures/rt_pred_rt_scatter_hist.eps'],'-depsc','-r300');    
    print(['./figures/rt_pred_rt_scatter_hist.png'],'-dpng','-r300');    
end

% do zoom
if 0
figure;
MSize = 20;
scatter_kde(temp_x_1,temp_y_1,'Marker','.');%,'DisplayStyle','tile','ShowEmptyBins','on','NumBins',[100 100]);
% cmap = cbrewer('seq','Blues',100);
% colormap(flipud(cmap));
hold on;
% scatter(temp_x_0,temp_y_0,'Marker','.');
% do piecewise regression on partioned data

sst = sum((temp_y_1 - mean(temp_y_1)).^2); % total variance
ssr = sum((temp_x_1-temp_y_1).^2); % residual error
r_square_all_fit_zoomed = 1 - ssr/sst               
r_corr_part = corr(temp_x_1,temp_y_1)
plot(temp_x_1,temp_x_1,'-','color','k',LW,LWval,'MarkerSize',0.25,'MarkerFaceColor','k','MarkerEdgeColor','k');
% plot(temp_x_1,temp_y_1_pred,'-o','color','k',LW,LWval,'MarkerSize',0.25,'MarkerFaceColor','k','MarkerEdgeColor','k');
% plot(temp_x_2,temp_y_2_pred,'-o','color','k',LW,LWval,'MarkerSize',0.25,'MarkerFaceColor','k','MarkerEdgeColor','k');
set(gca,FS,FSval-2,'box','on');
xlabel('RT [in msec]');
ylabel('PREDICTED RT [in msec]');
xlim([0 1000])
ylim([0 1000])

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on',... 
  'XGrid'       , 'on',...
  'Yminorgrid'  , 'on',... 
  'Xminorgrid'  , 'on');

if 1
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 12 12]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/rt_pred_rt_scatter_hist.fig']);
    print(['./figures/zoomed_rt_pred_rt_scatter_hist.eps'],'-depsc','-r300');    
    print(['./figures/zoomed_rt_pred_rt_scatter_hist.png'],'-dpng','-r300');    
end           
end
% return;
% ----- make r-square plot
figure;
clr = {[0.25 0.75 0.851],[0.5 .85 .25],[0.85 0.25 0.5],[0 0 1]};
% clr = {[0 0.1 1],[0 0.2 1],[0 0.3 1],[0 .4 1],[0 0.5 1]};
indx_sel = [1 2 3 4 5 6 7 8 9 12];
for i = 1:indx_dur
    mu_feats = 100*mean(r_square(indx_sel,:,i),2);
    h(i) = plot(1:length(mu_feats),mu_feats,'o-.','color',clr{i},LW,LWval,'MarkerSize',MSval,'MarkerFaceColor',clr{i},'MarkerEdgeColor',clr{i});
    hold on;
end
grid on;
hlegend = legend(h,{'25%','50%','75%','100%'});
rect = [0.5 0.85 .001 .05]; %[left bottom width height]
set(hlegend,'Position',rect);
set(gca,FS,FSval,'box','on');

feats_type_chosen = {};
for i = 1:length(feats_type)
    feats_type_chosen{i} = feats_type{i};
end
% {'F0','lsf','mel','mfcc','mfcc_{d1}','mfcc_{d2}',...
%     'temp','percp','spect'};
all_feats_type = feats_type_chosen;
all_feats_type{length(all_feats_type)+1} = 'ALL';

xticks(1:length(all_feats_type));
xticklabels(all_feats_type);
xtickangle(45)
xlim([0.5 length(all_feats_type)+.5])
xlabel('FEATURE');
ylabel('EXPLAINED VARIANCE (%)');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );

if 1
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 8 6]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/rt_feats_r_square.fig']);
    print(['./figures/rt_feats_r_square.eps'],'-depsc','-r300');    
end           

% ----- make slope
clr = {[0 0.5 0],[0 0.5 0],[0 0.5 0],[0 0.5 0],[0 0.5 0],[0 0.5 0]};
figure;
bmin = min(min(min(b(:,:,2))));
bmax = max(max(max(b(:,:,2))));

for indx_feats = 1:length(feats_type)
    plot((1:17)+(indx_feats-1)*20,squeeze(b(:,indx_feats,2)),'o','color',clr{1},'MarkerSize',MSval-1.5,'MarkerFaceColor',clr{1},'MarkerEdgeColor',clr{1});
    hold on;
    plot([20*(indx_feats) 20*(indx_feats)],[bmin bmax],'-.','color',[0.5 0.5 0.5],LW,LWval);
    plot(length(sel_subjs)/2+(indx_feats-1)*20,mean(squeeze(b(:,indx_feats,2)),1),'o','color','k','MarkerSize',MSval+1,'MarkerFaceColor','k','MarkerEdgeColor','k');
end
ylabel('LINEAR REGRESSION SLOPE');
xlim([0 length(feats_type)*20]);
ylim([bmin*1.05 bmax*1.05])
set(gca,FS,FSval,'box','on');
xticks([10:20:length(feats_type)*20]);
xticklabels(feats_type);
xtickangle(45);
xlabel('FEATURE');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );

if 1
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 8 6]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/rt_feats_lreg_slope.fig']);
    print(['./figures/rt_feats_lreg_slope.eps'],'-depsc','-r300');    
end           


% ----- make p-value plot
clr = repmat([0 0 1],length(sel_subjs),1);
h = [];
P_ALL_MAT = p_all{4}(:,[2:length(feats_type)+1]);
P_ALL_MAT(P_ALL_MAT>0.05) = 1;
P_ALL_MAT(P_ALL_MAT==0.05) = 1;
% P_ALL_MAT(P_ALL_MAT<0.05) = 0;

figure;imagesc(P_ALL_MAT);
cmap = cbrewer('seq','Blues',1);
colormap(cmap([1 3],:));
cbh = colorbar ; %Create Colorbar
cbh.Ticks = [0.25 0.75] ; %Create 8 ticks from zero to 1
% cbh.TickLabels = num2cell([0 1]) ;
cbh.TickLabels = {'p<0.05','p>0.05'} ;
set(gca,FS,FSval,'box','on');

% xticks(1:1:length(sel_subjs));
% % xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'});
% yticks(1:1:length(feats_type));
% yticklabels(feats_type_chosen);
% hold on;

yticks(1:1:length(sel_subjs));
yticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'});
xticks(1:1:length(feats_type));
xticklabels(feats_type_chosen);
xtickangle(45);
hold on;

% make the grid lines
M = size(P_ALL_MAT,1);
N = size(P_ALL_MAT,2);

for k = -1.5:1:M
    x = [0 N+1];
    y = [k k];
    plot(x,y,'Color','w','LineStyle','-');
    plot(x,y,'Color',[0.8 0.8 0.8],'LineStyle','-');
end

for k = -1.5:1:N
    x = [k k];
    y = [0 M+1];
    plot(x,y,'Color','w','LineStyle','-');
    plot(x,y,'Color',[0.8 0.8 0.8],'LineStyle','-');
end
hold off
ylabel('SUBJECT INDEX');
xlabel('FEATURE');

% position the color bar
x1=get(gca,'position'); % position = [x0 y0 width height]
x=get(cbh,'Position');
x(3)=0.02;
x(1) = x(1)+0.075;
x(2) = x(2)+0.2;
x(4) = x(4)-0.2;
set(cbh,'Position',x)
set(gca,'position',x1)

v = get(gca,'Position');
set(gca,'Position',[v(1)*1.4 v(2)*1.7 v(3:4)])

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      ,...
  'XGrid'       , 'off'      );

if 1
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 7 9]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/p_value_feats_mat_subjectwise.fig']);
    print(['./figures/p_value_feats_mat_subjectwise.eps'],'-depsc','-r300');    
end           
 
