

function add_csv_pretest(sound_path,store_path,fileID)

files = dir([sound_path '*.wav']);
nfiles = length(files);

ids{1} = 500;
ids{2} = [0 10 25 50 100];

cnt = 1;
for i = 1:nfiles
    for j = 1:length(ids{1})
        if (length(strfind(files(i).name,['noise_tone_F0_' num2str(ids{1}(j))])))
        pre_test(cnt).name = files(i).name;
        flag(cnt) = 1;
        cnt = cnt+1;
        end
    end
end

cnt_1 = cnt-1; 

for i = 1:nfiles
    for j = 1:length(ids{2})
        if (length(strfind(files(i).name,['tone_tone_F0_' num2str(ids{2}(j))])))
        pre_test(cnt).name = files(i).name;
        if (length(strfind(files(i).name,'Delta_0_')))
        flag(cnt) = 0;
        else
        flag(cnt) = 1;
        end    
        cnt = cnt+1;
        end
    end
end
cnt = cnt-1;

% ----- append to the csv file
% see below for field reference
%fprintf(fileID,'randomise_blocks,randomise_trials,display,sound,hc_ans,sound_nc_1,sound_nc_3,sound_nc_2,sound_c_1,sound_c_2,sound_c_3,progress,answer\n');

indx = randsample(1:cnt_1,cnt_1);
nrep = 3;
indx = repmat(indx,1,nrep);
indx
fprintf(fileID,',,rt_noise,%s,\n',pre_test(indx(3)).name);
copyfile(strcat(sound_path, pre_test(indx(3)).name),store_path);
for i = 1:nrep*cnt_1
    fprintf(fileID,',,trials,%s,,,,,,,,,%d\n',pre_test(indx(i)).name,flag(indx(i)));
    copyfile(strcat(sound_path, pre_test(indx(i)).name),store_path);
end


indx = randsample((cnt_1+1):cnt,cnt-(cnt_1+1)+1);
nrep = 3;
indx = repmat(indx,1,nrep);
indx
fprintf(fileID,',,rt_tone,%s,\n',pre_test(indx(3)).name);
copyfile(strcat(sound_path, pre_test(indx(3)).name),store_path);

for i = 1: length(indx)
    fprintf(fileID,',,trials,%s,,,,,,,,,%d\n',pre_test(indx(i)).name,flag(indx(i)));
    copyfile(strcat(sound_path,pre_test(indx(i)).name),store_path);
end
end

