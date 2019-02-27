


% collect aal files of good speakers from the LibriSpeech dataset
clear all; close all; clc
path = '/Users/neeks/Desktop/Documents/work/dBase/LibriSpeech/train-clean-100/filePool/';
filelist = dir([path '*.wav']);

direc_list = fopen('/Users/neeks/Desktop/Documents/work/dBase/LibriSpeech/train-clean-100/filePool/good100files_train_spkrs_list.txt');
direc = textscan(direc_list,'%s');
fclose(direc_list);
dnos = size(direc{1},1);

dnames = cell(dnos,1);
for i = 1:dnos
    dnames{i} = direc{1}{i};
end

cnt = zeros(dnos,1);
for i = 1:dnos
    for j = 1:length(filelist)
        if (length(strfind(filelist(j).name,['LS_' dnames{i} '-'])))
            cnt(i) = cnt(i)+1;
            fnames{i}{cnt(i)} = filelist(j).name;
        end
    end
end
display(['No. of speakers: ' num2str(dnos)]);


minD = 5;
maxD = 10;

nos_files = zeros(length(cnt),1);

for i = 1:length(cnt) % number of distinct first speaker/talker (T1)
    for j = 1:cnt(i)
        [x,Fs] = audioread([path fnames{i}{j}]);
        x = x(:,1);
        x = x./max(abs(x));
        
        dur{i}(j) = fix(length(x)/Fs);
        if dur{i}(j)>minD % duration of T1
            % remove any noise below 200 Hz
            fc_hpf= 200;
            L = length(x);
            h_ia = fir1(L,[fc_hpf 3.4e3]*2 /Fs);
            temp = freq_filtering(x,h_ia,2);
            tmaxD = min(maxD,dur{i}(j));
            
            temp = abs(hilbert(temp(fix(minD*Fs):fix(tmaxD*Fs))));

            temp = abs(hilbert(temp));
            fc_ia = 100;
            L = length(temp);
            h_ia = fir1(L,fc_ia*2 /Fs);
            temp = freq_filtering(temp,h_ia,2);
                            
            indx = find(temp<0.01*max(temp));
            if length(indx)>1
                nos_files(i) = nos_files(i)+1;
                y = x(1:fix(minD*Fs)+indx(1));
                % add a roll-off of 500 msec at start and end
                win = hamming(fix(500e-3*Fs));
                win = win-min(win);
                win = win./max(abs(win));
                y(1:length(win)/2) = win(1:length(win)/2).*y(1:length(win)/2);
                y(end-length(win)/2+1:end) = win(length(win)/2+1:end).*y(end-length(win)/2+1:end);

    %             figure;
    %             plot(x); hold on;
    %             plot(y);

                disp(['Selected ' num2str(num2str(fix(length(y)/Fs))) ' s of ' num2str(dur{i}(j)) ' s.']);
                display([num2str(i) ': ' num2str(j)]);
                store_path = './data/libriSpeech/f0_estims/train/dur_filt/'; 
                if (~ isdir([store_path]))
                    mkdir ([store_path]);
                end
                audiowrite([store_path fnames{i}{j}],y,Fs);
            end
        end
    end
end

% 
