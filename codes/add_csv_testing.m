

function add_csv_testing(sound_path,store_path,fileID)

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


progress = 1;
cnt = 0;
for i = 1:length(spkID) %:2
    for j = 1:length(spkID) %:2
        fnames = dir([sound_path spkID{i} '-' spkID{j} '/*.wav']);
        if i ==j
            flag = 0;
        else
            flag = 1;
        end
        if ~cnt
            fprintf(fileID,',,rt_speech_test,,\n');
        end
        
        for k = 1:length(fnames)
            if ((~rem(cnt,25)) && (cnt>0))
                disp(rem(cnt,25))
                disp(cnt)
                fprintf(fileID,',,break,,\n');
            end
            
            copyfile(strcat(sound_path, spkID{i}, '-', spkID{j}, '/', fnames(k).name),store_path)
            fprintf(fileID,',1,trials,%s,,,,,,,,%d,%d\n',fnames(k).name,progress,flag);
            cnt = cnt+1;
        end
    end
end


