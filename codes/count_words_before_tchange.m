

annot_path = './data/libriSpeech/spkrChange/listTest/M_5_spkrs_8_pairs_v0/annotations/';

files = dir([annot_path '*.csv']);

TCHANGE = zeros(length(files),1); % change instant
CWORDS = zeros(length(files),1); % words till change
TWORDS = zeros(length(files),1); % total number of words in file
CWDUR = zeros(length(files),1); % word duration till change
CWRATE = zeros(length(files),1); % word-rate till change instant

for i = 1:length(files)
    data = readtable([annot_path files(i).name],'Delimiter',',');
    
    [fpath,fname,fext] = fileparts([annot_path files(i).name]);
    str_1 = 'tChange_';
    str_2 = '_ms';

    indx_1 = strfind(fname,str_1)+length(str_1);
    indx_2 = strfind(fname,str_2)-1;
    tChange = str2double(fname(indx_1:indx_2));
    TCHANGE(i) = tChange;
    
    data.Var2 = fix(1000*data.Var2);
    data.Var3 = fix(1000*data.Var3);
    
    indx = find(data.Var2<(tChange+1));
    CWORDS(i) = length(indx);
    CWRATE(i) = length(indx)/TCHANGE(i)*1000;
    TWORDS(i) = length(data.Var3);
    
    dur = 0;
    for j = 1:length(indx)
        dur = dur+data.Var3(j)-data.Var2(j);
    end
    CWDUR(i) = dur;
end


% store the noise to tone change detection RT as baseline for each subject
if 1
store_path = './data/rt_feats/';
fileID = fopen([store_path 'file_wise_tChange_nWords.csv'],'w');
fprintf(fileID,'FNAME,TCHANGE,CWORDS,TWORDS,CWDUR,CWRATE\n');
for i = 1:length(files)
    [fpath,fname,fext] = fileparts([annot_path files(i).name]);
    fprintf(fileID,'%s,%d,%d,%d,%d,%.2f',fname,TCHANGE(i),CWORDS(i),TWORDS(i),CWDUR(i),CWRATE(i));
    fprintf(fileID,'\n');
end
fclose(fileID);
end
[min(CWORDS) max(CWORDS) mean(CWORDS)]








