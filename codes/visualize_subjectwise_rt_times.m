% make a single csv file with RT and speaker ID

clear all; close all;

data_path = './data/rt_feats/';
fname = 'subject_wise_rt.csv';

data = readtable([data_path fname],'Delimiter',',');

ptcpnt_id = {'CMU-1','cmu2','CMU_3','CMU-4','cmu_5','CMU_6','CMV-7','CMU_8','cmu_9','CMU_10','cmu_11','cmu12',...
    'CMU_13','CMU-14','cmu_15','CMU_16','CMU_17','iisc_01','iisc_02','iisc_03','iisc_04','iisc_05','iisc_06',...
    'iisc_07','iisc_08','iisc_09','iisc_10','iisc_11','iisc_12'};

sel_subjs =  [13:17 18:25 26:29];

% read the standard RT of each subject from the csv file

indx_subj = cell(length(ptcpnt_id),1);
% data_std = readtable([data_path 'subject_wise_mean_rt.csv'],'Delimiter',',');
data_std = readtable([data_path 'subject_wise_noise_mean_rt.csv'],'Delimiter',',');

rt_mat = nan(length(sel_subjs),200);
rt_mat_hit = nan(length(sel_subjs),200);
outlier = nan(length(sel_subjs),200);
rt_mat_hit_std = nan(length(sel_subjs),200);
RT = [];
group = [];

loop = 0;
h = zeros(4,length(sel_subjs));
p = zeros(4,length(sel_subjs));
for indx_pid = sel_subjs
    cnt = 1;
    cnt_hit = 1;
    loop = loop+1;
    indx_ptcpnt_id = [];
    indx_hit = [];
    for i = 1:length(data.Var1)
        if strcmp(data.Var1{i},ptcpnt_id{indx_pid})
            indx_ptcpnt_id(cnt) = i;
            cnt = cnt+1;
            if(data.Var6(i) == 3)
                indx_hit(cnt_hit) = i;
                cnt_hit = cnt_hit+1;
            end
        end
    end
    
    rt_mat(loop,1:200) = data.Var3(indx_ptcpnt_id);
    rt_mat_hit(loop,indx_hit-200*(loop-1)) = data.Var3(indx_hit);
    display(['Nos: ' num2str(length(indx_ptcpnt_id))])

    indx = find(~isnan(rt_mat_hit(loop,:)));
    
    group = [group; loop*ones(size(indx(:)))];
    RT = [RT;rt_mat(loop,indx(:)).'];
    
    if 0
    % find the chi squared fit to check gaussianity
    % to chi square testing 
    pd = fitdist(rt_mat_hit(loop,indx)','InverseGaussian');
    [h(1,loop),p(1,loop)] = chi2gof(rt_mat_hit(loop,indx),'CDF',pd);
    
    % to chi square testing 
    pd = fitdist(rt_mat_hit(loop,indx)','Lognormal');
    [h(2,loop),p(2,loop)] = chi2gof(rt_mat_hit(loop,indx),'CDF',pd);
    
    % to chi square testing 
    pd = fitdist(rt_mat_hit(loop,indx)','Gamma');
    [h(3,loop),p(3,loop)] = chi2gof(rt_mat_hit(loop,indx),'CDF',pd);

    % to chi square testing 
    pd = fitdist(rt_mat_hit(loop,indx)','Normal');
    [h(4,loop),p(4,loop)] = chi2gof(rt_mat_hit(loop,indx),'CDF',pd);
    
    % make a 2-D PDF for rt
    [H(loop,:) centers(loop,:)] = hist(rt_mat_hit(loop,indx),100);
    end
    
    % compute outliers based on mean RTs 
%     for i = 1:length(indx)
% %         [rt_mat_hit(loop,indx(i)) data_std.meanRT_ms(loop)+2*data_std.stdRT_ms(loop)] 
%         if (rt_mat_hit(loop,indx(i))> (data_std.meanRT_ms(loop)+2*data_std.stdRT_ms(loop)))
%             outlier(loop,indx(i)) = 1;
%         end
%     end
    % log transformation as it increases gaussianity 
     rt_mat_hit_std(loop,indx) = log10(rt_mat_hit(loop,indx))- log10(data_std.meanRT_ms(loop));%/(2*log10(data_std.stdRT_ms(loop)));
%     rt_mat_hit_std(loop,indx) = (rt_mat_hit(loop,indx))- (data_std.meanRT_ms(loop));%/(2*log10(data_std.stdRT_ms(loop)));
    
    % throw out outliers
%     indx_1 = find(outlier(loop,indx)==1);
%     indx(indx_1) = [];
    
%     figure;
%     normplot(rt_mat_hit_std(loop,indx));
    pd = fitdist(rt_mat_hit(loop,indx)','Normal');
    [h(1,loop),p(1,loop)] = chi2gof(rt_mat_hit(loop,indx),'CDF',pd);

    
%     figure;
%     normplot(log10((rt_mat_hit_std(loop,indx)).^(-1.5)));
    pd = fitdist(rt_mat_hit_std(loop,indx)','Normal');
    [h(2,loop),p(2,loop)] = chi2gof(rt_mat_hit_std(loop,indx),'CDF',pd);

    % make something like a box plot
    if 0
    figure(1);
    plot(loop,mean(rt_mat_hit(loop,indx)),'o','color',[0 0 0],'MarkerSize',6,'MarkerFaceColor',[0 0 0]);
    hold on;
    plot([loop loop],[min(rt_mat_hit(loop,indx)) max(rt_mat_hit(loop,indx))],'+-','color','b','MarkerSize',4,'MarkerFaceColor','r');
    hold on;

    figure(2);
    plot(loop,mean(rt_mat_hit_std(loop,indx)),'o','color',[0 0 0],'MarkerSize',6,'MarkerFaceColor',[0 0 0]);
    hold on;
    plot([loop loop],[min(rt_mat_hit_std(loop,indx)) max(rt_mat_hit_std(loop,indx))],'+-','color','b','MarkerSize',4,'MarkerFaceColor','r');
    hold on;

    figure(3);
    yyaxis left
    plot(loop,h(1,loop),'o-','color','b','MarkerSize',6,'MarkerFaceColor','b');
    hold on;
    yyaxis right
    plot(loop,p(1,loop),'o-','color','r','MarkerSize',6,'MarkerFaceColor','r');
    hold on;

    figure(4);
    yyaxis left
    plot(loop,h(2,loop),'o-','color','b','MarkerSize',6,'MarkerFaceColor','b');
    hold on;
    yyaxis right
    plot(loop,p(2,loop),'o-','color','r','MarkerSize',6,'MarkerFaceColor','r');
    hold on;
    end
    % to make equal sample set size of RT for each subject
    if 0
    indx_std{loop} = [];
    for i = indx
        if (rt_mat_hit(loop,i)<(rt_mu(loop)+2*rt_std(loop)) && rt_mat_hit(loop,i)>(rt_mu(loop)-2*rt_std(loop)))
            indx_std{loop} = [indx_std{loop} i];
        end
    end
    length(indx_std{loop});
    indx = randsample(indx_std{loop},125);
    rt_mat_std(loop,:) = rt_mat_hit(loop,indx);
    end
    
end

close all;
% ----- make plots
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 5; 
LWval = 0.25;
MSval = 2;

% box plot
boxplot(RT,group,'Notch','on','OutlierSize',2); hold on;

uniq_group = unique(group);
mu = zeros(length(uniq_group),1);
for i = 1:length(uniq_group)
    mu(i) = mean(RT(group==uniq_group(i)));
end
plot(uniq_group,mu,'o','color','k','MarkerSize',MSval-0.5,'MarkerFaceColor','k','MarkerEdgeColor','k');
% boxplot(RT,group,'Notch','on','OutlierSize',2);
xlim([0 18]);
xticks(1:1:17)
% xticklabels({'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
%              'S11','S12','S13','S14','S15','S16','S17'});
xlabel('SUBJECT INDEX')
ylabel('REACTION TIME [in msec]')
set(gca,FS,FSval,'box','on');

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6 5]); %x_width=10cm y_width=15cm        

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );

if 1
    fig_file = ['figures/subjectwise_rt_spread'];
    saveas(gcf,[fig_file '.fig']);
    print([fig_file '.eps'],'-depsc');%,'-r300');
    print([fig_file '.png'],'-dpng','-r300');
end









