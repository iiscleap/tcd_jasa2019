
function add_csv_headphone_check(sound_path,store_path,fileID)

file(1).name = 'antiphase_HC_IOS.mp3';
file(2).name = 'antiphase_HC_ISO.mp3';
file(3).name = 'antiphase_HC_OIS.mp3';

for i = 1:length(file)
copyfile(strcat(sound_path, file(i).name),store_path);
end

fprintf(fileID,',,headphone_check,%s,%d\n',file(1).name,3);
fprintf(fileID,',,headphone_check,%s,%d\n',file(2).name,2);
fprintf(fileID,',,headphone_check,%s,%d\n',file(3).name,3);

end
