

function add_csv_training(sound_path,store_path,fileID)

files = dir(sound_path);
dFlags = [files.isdir];
dnames = files(dFlags);

j = 1;
for i = 3:length(dnames)
    str = dnames(i).name;
    fpos = strfind(str,'-');
    if (length(strfind(str(1:fpos-1),str(fpos+1:end))))
        spkID{j} = str(1:fpos-1);
        j = j+1;
    end
end

cnt = 0;
for i = 1:2%length(spkID)
    for j = 1:2%length(spkID)
        fnames = dir([sound_path spkID{i} '-' spkID{j} '/*.wav']);
        if i ==j
            flag = 0;
        else
            flag = 1;
        end
        if ~cnt
            fnames1 = dir([sound_path spkID{1} '-' spkID{1} '/*.wav']);
            fnames2 = dir([sound_path spkID{2} '-' spkID{2} '/*.wav']);
            fnames3 = dir([sound_path spkID{3} '-' spkID{3} '/*.wav']);
        
            fnames4 = dir([sound_path spkID{3} '-' spkID{2} '/*.wav']);
            fnames5 = dir([sound_path spkID{1} '-' spkID{2} '/*.wav']);
            fnames6 = dir([sound_path spkID{2} '-' spkID{3} '/*.wav']);
            
            copyfile(strcat(sound_path, spkID{1}, '-' ,spkID{1}, '/' ,fnames1(1).name),store_path);
            copyfile(strcat(sound_path, spkID{2}, '-' ,spkID{2}, '/' ,fnames2(1).name),store_path);
            copyfile(strcat(sound_path, spkID{3}, '-' ,spkID{3}, '/' ,fnames3(1).name),store_path);

            copyfile(strcat(sound_path, spkID{3}, '-', spkID{2}, '/', fnames4(1).name),store_path);
            copyfile(strcat(sound_path, spkID{1}, '-', spkID{2}, '/', fnames5(1).name),store_path);
            copyfile(strcat(sound_path, spkID{2}, '-', spkID{3}, '/', fnames6(1).name),store_path);

            % ----- append to the csv file
            % see below for field reference
            %fprintf(fileID,'randomise_blocks,randomise_trials,display,sound,hc_ans,sound_nc_1,sound_nc_3,sound_nc_2,sound_c_1,sound_c_2,sound_c_3,progress,answer\n');
            
            fprintf(fileID,',,rt_speech_train,,,%s,%s,%s,%s,%s,%s,\n',fnames1(1).name,fnames2(1).name,fnames3(1).name,fnames4(1).name,...
            fnames5(1).name,fnames6(1).name);
        end
        for k = 1:1%length(fnames)
            copyfile(strcat(sound_path, spkID{i}, '-', spkID{j}, '/', fnames(k).name),store_path);
            % ----- append to the csv file
            fprintf(fileID,',,trials,%s,,,,,,,,,%d\n',fnames(k).name,flag);
            cnt = cnt+1;
        end
    end
end


