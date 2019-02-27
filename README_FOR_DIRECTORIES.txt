
Following are the details about the paths and files used in the project.
Edit: 30 Aug 2018 (Neeraj)
Last Edit: 24 Feb 2019 (Neeraj)

###################################################
LISTENING TEST DATA:
###################################################
Listening Test WAV files: 
./data/audioFiles/listTest/M_5_spkrs_8_pairs_v0/
Listening Test WAV files with subjectwise RT time-stamp: 
./data/libriSpeech/spkrChange/listTest_resp/M_5_spkrs_8_pairs_v0/added_part/

Listening Test WAV file annotations:
./data/libriSpeech/spkrChange/listTest/M_5_spkrs_8_pairs_v0/annotations/

Listening Test WAV file transcipts with change label $$:
./data/libriSpeech/spkrChange/listTest/M_5_spkrs_8_pairs_v0/books/change_trans.txt

Listening Test WAV files F0 trajectories:
./data/fEATS/libriSpeech/test-listening/M_5_spkrs_8_pairs_v0/

Listening Test Colleted Data CSV file:
./data/libriSpeech/spkrChange/listTest/dataSheet/accum_Nov_2017_19Jan_2018_17Apr_2018_mod.csv

###################################################
LISTENING TEST WAVE FILES with SUBJECT-WISE RT time stamps:
###################################################
Create using CODES: a. feats_extract_yaafe.py
Output will be stored in:
./data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part

###################################################
CODES:
###################################################
a. Feature Extraction:
feats_extract_yaafe.py
(has dependency on YAAFE toolkit: http://yaafe.sourceforge.net/)

b. Talkers:
word speaking rate: get_word_rate.py

c. Using IBMWatson
To call IBM Watson API from terminal:
curl -X POST -u T -u {username}:{password} --header "Content-Type: audio/flac" --data-binary @examp_intvw_plant_brain.flac  "https://stream.watsonplatform.net/speech-to-text/api/v1/recognize?model=en-US_NarrowbandModel&speaker_labels=true"

To parse watson output to get reaction time:
watson_annot_parse.py


d. Linear regression on RT and features:
v0_yaafe_stats_regress_varDur_transfmd_rt_feats.m


e. Plotting:-
Visualizing response instant versus change instant:
1. Humans (pool all): visualize_tChange_tResp.m
2. Humans (subjectwise) RT times: visualize_subjectwise_rt_times.m
3a. Humans (subjectwise): visualize_tChange_tResp_subjectwise.m
3b. Humans (subjectwise): visualize_wChange_tResp_subjectwise.m
4. Talker word rate and f0: plot_word_rate_f0_feats.m
4. Human machine comparison: get_human_machine_comparison.m
5. Document similarity matrix: document_similarity_matrix.m
6. Talker RT matrix, miss and FA: plot_human_rt_miss_matrix.m
7. Counts words before change: count_words_before_tchange.m
8. Utitlity for making scatter plot with histogram: scatter_kde.m
9. Make DET curve: plot_DET.m

f. Human Listening test data analysis:-
1. Speech test: get_listener_wise_speechtest.m
2. Noise and tone test: get_listner_wise_pretest.m


g. Machine test data analysis:-
The following uses: {ibmWatson_listest_outputs.csv, offline_results.csv, progressive_results.csv, LSTM_text_output.csv}
1. IBM Watson: get_ibmWatson_speechtest.m 
2. offlinePLDA: get_offlinePLDA_speechtest.m
3. onlinePLDA: get_onlinePLDA_speechtest.m
4. LSTM_Text: get_LSTMtext_speechtest.m


h. Stimuli concatenation:
1. compute_F0_stats_libriSpeech.m
2. get_dur_filtered_train_libriSpeech.m
3. select_spkrChange_LibriSpeech_files_v4_CMU_pool.m
4. make_spkrChange_LibriSpeech_files_v3_CMU_pool.m


i. Make spreadsheet for Gorilla:
1. make_spreadsheet_gorilla.m
2. add_csv_headphone_check.m
3. add_csv_training.m
4. add_csv_pretest.m
5. add_csv_testing.m
6. add_csv_volume_calibration.m

j. Pretest stimuli:
1. make_tone_tone.m
2. make_noise_tone.m
3. make_sil_tone.m

k. Confi file for HTK feature extraction
1. fbank_config.cfg
2. mfcc_config.cfg

l. Some SOX scripts:
audio_duration.sh
wav_to_mp3.sh

m. Details on code usgae for machine system to be updated here.