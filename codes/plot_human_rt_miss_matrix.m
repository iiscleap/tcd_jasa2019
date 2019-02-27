

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

% ----- plot miss/ fa/ dprime
% plot miss
data_miss = [mu_sub_miss.']*100;
figure;
b = bar(data_miss);
b(1).FaceColor = [0 0.5 1];
grid on;
ylabel('MISS (%)');
xlabel('SUBJECT INDEX');
xticks([1:17])
xticklabels({'1','2','3','4','5','6','7','8','9',...
             '10','11','12','13','14','15','16','17'});
% yticks([0 5 10 15])
% yticklabels({'45','60','75','90','100','5','10','15','20'});
xlim([0 18]);
ylim([0 12]);
grid on;
% h = text(0,5.5,'MISS','Color',[0 0.5 1],'FontSize',FSval,'HorizontalAlignment','center');
% set(h,'Rotation',90);
set(gca,FS,FSval-2,'box','on');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );
% set(gca,'OuterPosition',[left bottom + 0.1 width height])
if 1 % move legend and again save offline
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 6 1.5]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/subjectwise_barplot_miss.fig']);
    print(['./figures/subjectwise_barplot_miss.eps'],'-depsc','-r300');    
end

% plot FA
data_fa = [mu_sub_fa.']*100;
figure;
b = bar(data_fa);
b(1).FaceColor = 'r';
grid on;
ylabel('FA (%)');
xlabel('SUBJECT INDEX');
xticks([1:17])
xticklabels({'1','2','3','4','5','6','7','8','9',...
             '10','11','12','13','14','15','16','17'});
% yticks([0 5 10 15])
% yticklabels({'45','60','75','90','100','5','10','15','20'});
xlim([0 18]);
ylim([0 15]);
grid on;
% h = text(0,5,'FA','Color','r','FontSize',FSval,'HorizontalAlignment','center');
% set(h,'Rotation',90);
set(gca,FS,FSval-2,'box','on');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );
% set(gca,'OuterPosition',[left bottom + 0.1 width height])
if 1 % move legend and again save offline
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 6 1.5]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/subjectwise_barplot_fa.fig']);
    print(['./figures/subjectwise_barplot_fa.eps'],'-depsc','-r300');    
end   

% plot dprime
d_prime = zeros(1,length(mu_sub_miss));
for i = 1:length(mu_sub_miss)
    if mu_sub_hit(i) == 1
        temp_hit = 1-0.5/sum(tot_stimuli_hit{indx_pids}(:));
    else
        temp_hit = mu_sub_hit(i);
    end
    
    if mu_sub_fa(i) == 0
        temp_fa = 0.5/sum(tot_stimuli_hit{indx_pids}(:)+tot_stimuli_fa{indx_pids}(:));
    else
        temp_fa = mu_sub_fa(i);
    end
    
    d_prime(i) = dprime_simple(temp_hit,temp_fa);
end
figure;
b = bar(d_prime);
b(1).FaceColor = 'k';
grid on;
ylabel('d-prime');
xlabel('SUBJECT INDEX');
xticks([1:17])
xticklabels({'1','2','3','4','5','6','7','8','9',...
             '10','11','12','13','14','15','16','17'});
% yticks([0 5 10 15])
% yticklabels({'45','60','75','90','100','5','10','15','20'});
xlim([0 18]);
ylim([0 6]);
grid on;
% h = text(0,3.5,'d-prime','Color',[0 0 0],'FontSize',FSval,'HorizontalAlignment','center');
% set(h,'Rotation',90);
set(gca,FS,FSval-2,'box','on');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );
% set(gca,'OuterPosition',[left bottom + 0.1 width height])
if 1 % move legend and again save offline
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 6 1.5]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/subjectwise_dprime.fig']);
    print(['./figures/subjectwise_dprime.eps'],'-depsc','-r300');    
end   


% return;



% ----- get the speaker pair-wise average RT matrix
indx_pids = 13:29;
mu_rt = cell(length(indx_pids),1);
for loop = indx_pids
    for i = 1:5
        for j = 1:5
            if i~=j
               mu_rt{loop}(i,j) = sum(CH_mu_rt{loop}{i,j})/tot_stimuli_hit{loop}(i,j);
            else
               mu_rt{loop}(i,j) = 0;
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
if 0
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


% ----- make plots
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 3;

addpath('../plotting/cbrewer/')
x = fix(MU_RT);
figure; 
imagesc(x,[500 900]);
axis image
grid on
%# grid domains
xg = 0.5:1:5;
yg = 0.5:1:5;
%# label coordinates
[xlbl, ylbl] = meshgrid(xg+0.5, yg+0.5);
%# create cell arrays of number labels
lbl = strtrim(cellstr(num2str(x(:))));
for i = 1:5
    lbl{(i-1)*5+i} = 'na';
end
text(xlbl(:), ylbl(:), lbl(:),'color','k',...
    'HorizontalAlignment','center','VerticalAlignment','middle',FS,FSval);

cmap = cbrewer('seq','Blues',10);
colormap((cmap));
xticks(1:1:5)
xticklabels({'T1','T2','T3','T4','T5'});
yticks(1:1:5)
yticklabels({'T1','T2','T3','T4','T5'});
xlabel('TALKER INDEX');
ylabel('TALKER INDEX')
set(gca,FS,FSval,'box','on');
c = colorbar;
x1=get(gca,'position'); % position = [x0 y0 width height]
x=get(c,'Position');
x(3)=0.01;
x(1) = x(1)+0.1;
x(2) = x(2)+0.05;
x(4) = x(4)-0.1;
set(c,'Position',x)
set(gca,'position',x1)

if 1
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 5 5]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/human_talker_wise_rt_matrix.fig']);
print(['./figures/human_talker_wise_rt_matrix.eps'],'-depsc','-r300');    
end

x = fix(1000*(MU_MISS))/10;
figure; 
imagesc(x,[0 20]);
axis image
grid on
%# grid domains
xg = 0.5:1:5;
yg = 0.5:1:5;
%# label coordinates
[xlbl, ylbl] = meshgrid(xg+0.5, yg+0.5);
%# create cell arrays of number labels
lbl = strtrim(cellstr(num2str(x(:))));
for i = 1:5
    lbl{(i-1)*5+i} = 'na';
end
text(xlbl(:), ylbl(:), lbl(:),'color','k',...
    'HorizontalAlignment','center','VerticalAlignment','middle',FS,FSval);

cmap = cbrewer('seq','Blues',10);
colormap((cmap));
xticks(1:1:5)
xticklabels({'T1','T2','T3','T4','T5'});
yticks(1:1:5)
yticklabels({'T1','T2','T3','T4','T5'});
xlabel('TALKER INDEX');
ylabel('TALKER INDEX')
set(gca,FS,FSval,'box','on');
c = colorbar;
x1=get(gca,'position'); % position = [x0 y0 width height]
x=get(c,'Position');
x(3)=0.01;
x(1) = x(1)+0.1;
x(2) = x(2)+0.05;
x(4) = x(4)-0.1;
set(c,'Position',x)
set(gca,'position',x1)

if 1
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 5 5]); %x_width=10cm y_width=sel_nharmcm        
saveas(gcf,['./figures/human_talker_wise_miss_matrix.fig']);
print(['./figures/human_talker_wise_miss_matrix.eps'],'-depsc','-r300');    
end

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
ylabel('FALSE ALARM (%)')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      );

if 1
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 5 4]); %x_width=10cm y_width=sel_nharmcm        
    saveas(gcf,['./figures/talker_wise_fa.fig']);
    print(['./figures/talker_wise_fa.eps'],'-depsc','-r300');    
end           

end
