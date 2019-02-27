
% simple signals for change detection
% a. SIL- TONE
clear all; clc; close all;

store_path = './data/calibration_signals/tone_tone/';

seed = 1;
rng(seed);

if (~ isdir(store_path))
mkdir (store_path);
end

f0 = 500;
Delta = [0 10 25 50 100];

dur = 6;
tdur = 2;
Fs = 48e3;

tchange_intv = [3 5];
nos = 3;

t = 0:1/Fs:7-1/Fs;
alpha = 25;
traj = @(t,t0) 1./(1+exp(-alpha*(t-t0)));

for i = 1:length(Delta)
   for j = 1:nos
       t0 = (2*rand(1))+3;
       ftraj = f0+Delta(i)*traj(t,t0);
       phase = cumsum(ftraj*1/Fs);
       z = sin(2*pi*phase);
%        figure; plot(t,ftraj);
       audiowrite([store_path 'tone_tone_F0_' num2str(f0) '_Hz_Delta_' num2str(Delta(i)) '_' num2str(j) '_tChange_' num '.wav'],0.7*z,Fs);
   end
end    
    
    
    






