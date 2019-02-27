
close all;

% ----- make plots
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 6; 
LWval = 0.5;
MSval = 3;

addpath('../plotting/cbrewer/')
x = [1 .77 .76 .73 .72;...
     .77 1 .82 .77 .77;...
     .76 .82 1 .82 .82;...
     .73 .77 .82 1 .84;...
     .72 .77 .81 .84 1];
 
imagesc(x,[0.7 1]);
axis image
grid on
%# grid domains
xg = 0.5:1:5;
yg = 0.5:1:5;
%# label coordinates
[xlbl, ylbl] = meshgrid(xg+0.5, yg+0.5);
%# create cell arrays of number labels
lbl = strtrim(cellstr(num2str(x(:))));
text(xlbl(:), ylbl(:), lbl(:),'color','w',...
    'HorizontalAlignment','center','VerticalAlignment','middle',FS,FSval);

cmap = cbrewer('seq','Blues',10);
colormap(flipud(cmap));
xticks(1:1:5)
xticklabels({'D1','D2','D3','D4','D5'});
yticks(1:1:5)
yticklabels({'D1','D2','D3','D4','D5'});
xlabel('DOCUMENT INDEX');
ylabel('DOCUMENT INDEX')
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
saveas(gcf,['./figures/doc_similarity.fig']);
print(['./figures/doc_similarity.eps'],'-depsc','-r300');    
end

