
close all;

spkIDs = [374, 2843, 5456, 7447, 7505];

mu_WR = [2.76, 2.67, 3.06, 2.21, 2.2];
std_WR = [0.30, 0.29, 0.41, 0.29, 0.29];

mu_f0 = [151 102 124 101 155];
std_f0 = [25 24 35 19 41];


% ----- make plots
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 1;
MSval = 4;

figure;
for i = 1:5
plot([i,i], [mu_WR(i)-std_WR(i) mu_WR(i)+std_WR(i)],'-+','color','k',LW,LWval,MS,MSval-1,'MarkerFaceColor',[1 .6 .6]);
hold on;
plot(i,mu_WR(i),'o','color',[0 0.5 1],LW,LWval,MS,MSval-1,'MarkerFaceColor',[0 0.5 1]);
end
grid on;
xlabel('TALKER INDEX')
ylabel('WORDS-PER-SEC')
xticks(1:1:5)
xticklabels({'T1','T2','T3','T4','T5'});
xlim([0.5 5.5]);
set(gca,FS,FSval,'box','on');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3.5 5]); %x_width=10cm y_width=15cm        

if 1
    fig_file = ['figures/talker_word_rate'];
    saveas(gcf,[fig_file '.fig']);
    print([fig_file '.eps'],'-depsc');%,'-r300');
end

figure;
for i = 1:5
plot([i,i],[mu_f0(i)-std_f0(i) mu_f0(i)+std_f0(i)],'-+','color','k',LW,LWval,MS,MSval-1,'MarkerFaceColor',[1 .6 .6]);
hold on;
plot(i,mu_f0(i),'o','color',[0 0.5 1],LW,LWval,MS,MSval-1,'MarkerFaceColor',[0 0.5 1]);
end
grid on;
xlabel('TALKER INDEX')
ylabel('FREQUENCY [in Hz]')
xticks(1:1:5)
xticklabels({'T1','T2','T3','T4','T5'});
xlim([0.5 5.5]);
set(gca,FS,FSval,'box','on');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3.5 5]); %x_width=10cm y_width=15cm        

if 1
    fig_file = ['figures/talker_f0'];
    saveas(gcf,[fig_file '.fig']);
    print([fig_file '.eps'],'-depsc');%,'-r300');
end
