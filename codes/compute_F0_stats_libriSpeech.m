




% Code to read the files in the LibriSpeech:
% a. extract the f0 variations per speaker
% Created on: August 10, 2017 by Neeks

%path = '~/Desktop/Documents/work/dBase/LibriSpeech/dev-clean/8842/304647/';
clear all; close all; clc;


direc_list = fopen('/Users/neeks/Desktop/Documents/work/dBase/LibriSpeech/dev-clean/list_direcs.txt');
dnames = textscan(direc_list,'%s');
dnos = size(dnames{1},1);

str_direc = cell(dnos,1);
for i = 1:dnos
    str_direc{i} = dnames{1}{i};
end

path = '/Users/neeks/Desktop/Documents/work/code/matlab_codes/spkrDet/data/libriSpeech/f0_estims/new/';
fnames = dir([path '*.wav']);

cnt = zeros(dnos,1);
for i = 1:dnos
    for j = 1:length(fnames)
        if (length(strfind(fnames(j).name,['_' str_direc{i} '-'])))
            cnt(i) = cnt(i)+1;
            str_spkr{i}{cnt(i)} = fnames(j).name;
        end
    end
end

% // discard 26,22,34,20
j = 0;
cnt = zeros(dnos-4,1);
for i = 1:dnos
    if ((i~=26) && (i~=22) && (i~=34) && (i~=20)) 
        j = j+1;
        direc_id{j} = str_direc{i};
        nstr_spkr{j} = str_spkr{i};
        cnt(j) = size(nstr_spkr{j},2);
    end
end

for i = 1:length(cnt)
    temp = [];
    for j = 1:cnt(i)
        [x,Fs] = audioread([path nstr_spkr{i}{j}]);
        temp = [temp; x(:,2)*1000];
    end
    pitch{i}(1,1) = mean(temp((temp>80) & (temp<400)));
    pitch{i}(1,2) = std(temp((temp>80) & (temp<400)));
    pitch{i}(1,3) = min(temp((temp>80) & (temp<400)));
    pitch{i}(1,4) = max(temp((temp>80) & (temp<400)));
end


% ----- sort according to mean F0
for i = 1:length(cnt)
    mu_pitch(i) = pitch{i}(1,1);
    std_pitch(i) = pitch{i}(1,2);
end

[val,indx] = sort(mu_pitch);
mu_pitch = mu_pitch(indx);
std_pitch = std_pitch(indx);

for i = 1:length(direc_id)
   dnames{i} = direc_id{indx(i)};
end

fileID = fopen('selected_spkrs.txt','w');

for i = 1:length(dnames)
    fprintf(fileID,'%s\n',dnames{i});
end
fclose(fileID);

% ----- sort according to std F0
for i = 1:length(cnt)
    mu_pitch(i) = pitch{i}(1,1);
    std_pitch(i) = pitch{i}(1,2);
end

[val,indx] = sort(std_pitch);
mu_pitch = mu_pitch(indx);
std_pitch = std_pitch(indx);

for i = 1:length(direc_id)
   dnames{i} = direc_id{indx(i)};
end


return;

fileID = '/Users/neeks/Desktop/Documents/work/dBase/LibriSpeech/SPEAKERS.TXT';
spkr_list = fopen(fileID);
spkr_list = textscan(spkr_list,'%s');

for i = 1:length(spkr_list{1})
    for j = 1:length(dnames)
        if ((strfind(spkr_list{1}{i},[' ' dnames{j} ' '])))
            display(spkr_list{1}{i});
        end
    end
end


% ----- get the speaker details
fileID = '/Users/neeks/Desktop/Documents/work/dBase/LibriSpeech/SPEAKERS.TXT';
% spkr_list = fopen(fileID);
% spkr_list = textscan(spkr_list,'%s');

fid = fopen(fileID);
tline = fgetl(fid);
indx = 0;
while ischar(tline)
%    disp(tline)
    indx = indx +1;
    tline = fgetl(fid);
    str{indx} = tline;
end
fclose(fid);

for i = 1:length(dnames)
    for j = 1:length(str)
        temp = [dnames{i} ' '];
        if (strncmp(str{j},temp,length(temp)))
        display([num2str(i) ' :' str{j}]);
        fpos = regexp(str{j},'\|');
        gID(i) = str{j}(fpos(1)+2);
        end
    end
end

% ----- make gender array
Mindx = find(gID=='M');
Findx = find(gID=='F');


% ----- plot with color code of Male and Female
close all;
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth'; JL = 'jumpline';
FSval = 14; 
LSval = 12;
LWval = 0.5;
MSval = 4;

figure;
plot(1:length(Mindx),mu_pitch(Mindx),'-o','color',[0 0.5 1],LW,LWval); hold on;
plot(1:length(Mindx),mu_pitch(Mindx)+std_pitch(Mindx),'-','color',[0.5 0.5 0.5],LW,LWval); hold on;
plot(1:length(Mindx),mu_pitch(Mindx)-std_pitch(Mindx),'-','color',[0.5 0.5 0.5],LW,LWval); hold on;

plot((1:length(Findx))+length(Mindx),mu_pitch(Findx),'-o','color',[0 0.5 0]); hold on;
plot((1:length(Findx))+length(Mindx),mu_pitch(Findx)+std_pitch(Findx),'-','color',[0.5 0.5 0.5],LW,LWval); hold on;
plot((1:length(Findx))+length(Mindx),mu_pitch(Findx)-std_pitch(Findx),'-','color',[0.5 0.5 0.5],LW,LWval); hold on;

text(1,275,'MALE', 'color',[0 0.5 1])
text(20,275,'FEMALE', 'color',[0 0.5 0])
grid on;
xlabel('speaker index');
ylabel('mean F0');

set(gca,FS,FSval,'box','on');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm        
if (~ isdir('./figures'))
    mkdir ('./figures');
end
if 1
    saveas(gcf,['figures/libriSpeech_dev_clean_F0_stats', '.fig']);
    print('figures/libriSpeech_dev_clean_F0_stats.eps','-depsc');
    saveas(gcf,'figures/libriSpeech_dev_clean_F0_stats.png');
end




% some extra code I had to use for input file renaming
% some filename fixing
if 0 
path = '/Users/neeks/Desktop/Documents/work/code/matlab_codes/spkrDet/data/libriSpeech/f0_estims/';

fnames = dir([path '*.wav']);

file_list = fopen([path 'all_files_to_change.txt']);
f = textscan(file_list,'%s');
nfiles1 = size(f{1},1);


file_list = fopen([path 'all_files_new_names.txt']);
g = textscan(file_list,'%s');
nfiles2 = size(g{1},1);

for i = 1:nfiles1
    [x,Fs] = audioread([path f{1}{i}]);
    audiowrite([path 'new/' g{1}{i}],x,Fs);
end

end
