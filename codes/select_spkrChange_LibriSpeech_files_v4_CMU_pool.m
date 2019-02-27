
% select files randomly with following constraints:
% a. choose each speaker out of 20.
% b. choose 10 files from its own combination
% c. choose 10 files from different speaker combination

clear all;

% ----- make the training file entries into the CVS file
%folderID = 'telugu_M_4_spkrs_8_pairs_v0';

folderID = 'M_5_spkrs_8_pairs_v0';
store_path = ['./data/libriSpeech/spkrChange/listTest/' folderID '/'];

if (~ isdir([store_path]))
mkdir ([store_path]);
end
% ----- make the cvs file
fileID = fopen(['./data/libriSpeech/spkrChange/listTest/' folderID '/gorilla_list.csv'],'w');
% ----- add the header
fprintf(fileID,'randomise_blocks,randomise_trials,display,sound,hc_ans,sound_nc_1,sound_nc_3,sound_nc_2,sound_c_1,sound_c_2,sound_c_3,progress,answer\n');
if 1
% ----- add instructions
fprintf(fileID,',,instructions,,,,,,,,,,\n');
% ----- add volume calibration
% sound_path = '/Users/neeks/Desktop/Documents/work/code/matlab_codes/spkrDet/data/libriSpeech/spkrChange/listTest/headphoneCheck/mp3/';
sound_path = ['./data/libriSpeech/spkrChange/' folderID '/training/'];
add_csv_volume_calibration(sound_path,store_path,fileID);
% ----- add headphone check
sound_path = '/Users/neeks/Desktop/Documents/work/code/matlab_codes/spkrDet/data/libriSpeech/spkrChange/listTest/headphoneCheck/mp3/';
add_csv_headphone_check(sound_path,store_path,fileID);
% ----- add pre-test
sound_path = '/Users/neeks/Desktop/Documents/work/code/matlab_codes/syntDet/data/pre_test_stimuli/';
add_csv_pretest(sound_path,store_path,fileID);
% ----- training
sound_path = ['./data/libriSpeech/spkrChange/' folderID '/training/'];
add_csv_training(sound_path,store_path,fileID);
end
% ----- testing
sound_path = ['./data/libriSpeech/spkrChange/' folderID '/'];
add_csv_testing(sound_path,store_path,fileID);

fprintf(fileID,',,end,,,,,\n');
fclose(fileID);

