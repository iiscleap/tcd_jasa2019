
% simple signals for change detection
% a. SIL- TONE
clear all; clc; close all;

store_path = './data/calibration_signals/noise_tone/';

store_path = './data/pre_test_stimuli/';

seed = 1;
rng(seed);

if (~ isdir(store_path))
mkdir (store_path);
end

f0 = [500];

dur = 6;
tdur = 2;
Fs = 48e3;

tchange_intv = [3 5];
nos = 5;

for i = 1:length(f0)
   for j = 1:nos
       sil_nsamp = randsample(fix(tchange_intv(1)*Fs):fix(tchange_intv(2)*Fs),1);
%        x = zeros(sil_nsamp,1);
       x = 0.1*(rand(sil_nsamp,1)-0.5);
       y = sin(2*pi*f0(i)*(0:1/Fs:tdur-1/Fs)).';
       z = [x;y];
       tChange = fix(sil_nsamp/Fs/1e-3);
       audiowrite([store_path 'noise_tone_F0_' num2str(f0(i)) '_Hz_' num2str(j) '_tChange_' num2str(tChange) '.wav'],0.7*z,Fs);
   end
end    