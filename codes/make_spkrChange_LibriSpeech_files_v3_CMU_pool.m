



% making speaker change files using LibriSpeech train-set
clear all ;close all; clc;
addpath('/Users/neeks/Desktop/Documents/work/code/matlab_codes/rev_ffs/')

seed = 10;
rng(seed);
path = './data/libriSpeech/f0_estims/train/dur_filt/';
filelist = dir([path '*.wav']);


direc_list = fopen('/Users/neeks/Desktop/Documents/work/dBase/LibriSpeech/train-clean-100/filePool/morethan81_good100files_train_spkrs_list.txt');
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
        if (length(strfind(filelist(j).name,['_' dnames{i} '-'])))
            cnt(i) = cnt(i)+1;
            fnames{i}{cnt(i)} = filelist(j).name;
        end
    end
end
display(['No. of speakers: ' num2str(dnos)]);

% ----- obtain gender
fileID = '/Users/neeks/Desktop/Documents/work/dBase/LibriSpeech/SPEAKERS.TXT';
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
        %display([num2str(i) ' :' str{j}]);
        fpos = regexp(str{j},'\|');
        gID(i) = str{j}(fpos(1)+2);
        end
    end
end
fprintf('No. of males: %d\n', length(find(gID=='M')));
fprintf('No. of females: %d\n', length(find(gID=='F')));


% ----- make a file set for all male speakers
temp_cnt = 1;
% sel_IDs = {'4898','5456','7447','2843','8425','374','8226','307','8425','7505'};
% sel_IDs = {'7447','2843','374','7505','5322','4898','5456','8425'}; used on Nov. 18-19
sel_IDs = {'7447','2843','374','5456','7505','5322','4898','8425','5322'};

mindx = zeros(length(sel_IDs),1);
for i = 1:length(sel_IDs)
    for j = 1:length(dnames)
        temp = str2double(sel_IDs{i})-str2double(dnames{j});
        if ~temp
            mindx(temp_cnt) = j;
            temp_cnt = temp_cnt+1;
        end
    end
end

% ----- make a file set for all male speakers
tot_spkr = length(mindx);
N_spkr = 5;
N_train_spkr = 3;
% uncomment below to randomize
% temp_indx = randsample(tot_spkr,N_spkr+N_train_spkr);
temp_indx = 1:(N_spkr+N_train_spkr);
indx_spkr = temp_indx(1:N_spkr);

N_files_spkr_per_spkr = 8;
N_files_spkr = N_spkr*N_files_spkr_per_spkr;

spkr_files_set = cell(N_spkr,1);

for i = 1:N_spkr % number of distinct first speaker/talker (T1)
    spkr_files_set{i} = cell(N_files_spkr,1);
    tot_file = cnt(mindx(indx_spkr(i)));
    indx_files = randsample(tot_file,N_files_spkr);
    
    for j = 1:N_files_spkr
        spkr_files_set{i}{j} = fnames{mindx(indx_spkr(i))}{indx_files(j)};
    end
end
 
tstamps_array = [];
% ----- read the above file set and make test stimuli
for i = 1:N_spkr
    for j = 1:N_spkr
        for k = 1:N_files_spkr_per_spkr
            [x,Fs] = audioread([path spkr_files_set{i}{(j-1)*N_files_spkr_per_spkr+k}]);
            [y,Fs] = audioread([path spkr_files_set{j}{N_files_spkr-(j-1)*N_files_spkr_per_spkr-k+1}]);
            
            Dur_T1 = fix((length(x)/Fs/1e-3)); % in msec
            Dur_T2 = fix(4*Fs); % in samples
            % choose 4 sec of y
            temp_y = y(1:Dur_T2);
            tChange = 0;

            % remove any noise below 200 Hz
            fc_hpf= 200;
            L = length(temp_y);
            h_ia = fir1(L,[fc_hpf 3.4e3]*2 /Fs);
            temp_yf = freq_filtering(temp_y,h_ia,2);
            % compute hilbert
            temp_yh = abs(hilbert(temp_yf));
            fc_ia = 500;
            L = length(temp_y);
            h_ia = fir1(L,fc_ia*2 /Fs);
            temp_ya = freq_filtering(temp_yh,h_ia,2);
            % compute the firxt index of voice sample
            temp_ya(1:100) = 0; % zero out first 100 samples for abrupt spikes 
            indx = find(temp_ya>(.04*max(abs(temp_ya))));
            tChange = fix((indx(1)-1)/Fs/1e-3); % round to msec
            display(num2str(tChange));
            if ~tChange
                error('File for T2 begins with a spike!');
            end
            tstamps_array = [tstamps_array; tChange];
            % join the segments
            z = [0.8*x/max(abs(x)); 0.8*temp_y/max(abs(temp_y))];
            tChange = Dur_T1 + fix((indx(1)-1)/Fs/1e-3); % round to msec
            
            % ----- store the sound file
            
            store_path = ['./data/libriSpeech/spkrChange/' ...
                gID(mindx(indx_spkr(i))) '_' num2str(N_spkr) '_spkrs_' num2str(N_files_spkr_per_spkr) ...
                '_pairs_v0/' dnames{mindx(indx_spkr(i))} '-'];
            store_path = [store_path dnames{mindx(indx_spkr(j))} '/'];
            %display([gID(mindx(i)) '-' gID(mindx(k))])
            if (~ isdir([store_path]))
            mkdir ([store_path]);
            end
            
            file_name = [spkr_files_set{i}{(j-1)*N_files_spkr_per_spkr+k}(1:end-4) '-SWT-' spkr_files_set{j}{N_files_spkr-(j-1)*N_files_spkr_per_spkr-k+1}(1:end-4) '_tChange_' num2str(tChange) '_ms.wav'];
            file_name = strrep(file_name,'-','_');
            audiowrite([store_path file_name],z,Fs);                            

        end        
    end
end

old_N_spkr = N_spkr;
old_N_files_spkr_per_spkr = N_files_spkr_per_spkr;
indx_train_spkr = temp_indx(N_spkr+1:end);
indx_spkr = indx_train_spkr;
N_spkr = N_train_spkr;
N_files_spkr_per_spkr = 2;
N_files_spkr = N_spkr*N_files_spkr_per_spkr;

spkr_files_set = cell(N_spkr,1);

for i = 1:N_spkr % number of distinct first speaker/talker (T1)
    spkr_files_set{i} = cell(N_files_spkr,1);
    tot_file = cnt(mindx(indx_spkr(i)));
    indx_files = randsample(tot_file,N_files_spkr);
    
    for j = 1:N_files_spkr
        spkr_files_set{i}{j} = fnames{mindx(indx_spkr(i))}{indx_files(j)};
    end
end

% ----- read the above files and make training stimuli
for i = 1:N_spkr
    for j = 1:N_spkr
        for k = 1:N_files_spkr_per_spkr
            [x,Fs] = audioread([path spkr_files_set{i}{(j-1)*N_files_spkr_per_spkr+k}]);
            [y,Fs] = audioread([path spkr_files_set{j}{N_files_spkr-(j-1)*N_files_spkr_per_spkr-k+1}]);
            
            Dur_T1 = fix((length(x)/Fs/1e-3)); % in msec
            Dur_T2 = fix(4*Fs); % in samples
            % choose 4 sec of y
            temp_y = y(1:Dur_T2);

            % remove any noise below 200 Hz
            fc_hpf= 200;
            L = length(temp_y);
            h_ia = fir1(L,[fc_hpf 3.4e3]*2 /Fs);
            temp_yf = freq_filtering(temp_y,h_ia,2);
            % compute hilbert
            temp_yh = abs(hilbert(temp_yf));
            fc_ia = 500;
            L = length(temp_y);
            h_ia = fir1(L,fc_ia*2 /Fs);
            temp_ya = freq_filtering(temp_yh,h_ia,2);
            % compute the firxt index of voice sample
            temp_ya(1:100) = 0; % zero out first 100 samples for abrupt spikes 
            indx = find(temp_ya>(.04*max(abs(temp_ya))));
            tChange = fix((indx(1)-1)/Fs/1e-3); % round to msec
            display(num2str(tChange));
            if ~tChange
                error('File for T2 begins with a spike!');
            end
%             tstamps_array = [tstamps_array; tChange];
            % join the segments
            z = [0.8*x/max(abs(x)); 0.8*temp_y/max(abs(temp_y))];
            tChange = Dur_T1 + fix((indx(1)-1)/Fs/1e-3); % round to msec
            
            % ----- store the sound file
            store_path = ['./data/libriSpeech/spkrChange/' ...
                gID(mindx(indx_spkr(i))) '_' num2str(old_N_spkr) '_spkrs_' num2str(old_N_files_spkr_per_spkr) ...
                '_pairs_v0/training/' dnames{mindx(indx_spkr(i))} '-'];
            store_path = [store_path dnames{mindx(indx_spkr(j))} '/'];
            %display([gID(mindx(i)) '-' gID(mindx(k))])
            if (~ isdir([store_path]))
            mkdir ([store_path]);
            end
            
            file_name = [spkr_files_set{i}{(j-1)*N_files_spkr_per_spkr+k}(1:end-4) '-SWT-' spkr_files_set{j}{N_files_spkr-(j-1)*N_files_spkr_per_spkr-k+1}(1:end-4) '_tChange_' num2str(tChange) '_ms.wav'];
            file_name = strrep(file_name,'-','_');
            audiowrite([store_path file_name],z,Fs);                            

        end        
    end
end





    
    











