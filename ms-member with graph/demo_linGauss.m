% This is a demo script for a single run of MS-CPHD, MS-TCPHD and MS-MeMBer
% filters in UKF and SMC implementations.
%
% Contact: augustin.saucan@mail.mcgill.ca

close all
clear all
clc
global flag
flag=0; %默认为0,没消除杂波集Kp的影响。
mcnt=50;
%restoredefaultpath       % clear previous path
addpath(genpath(pwd));   % append all folders to current path

nb_sensors = 3;         % specify number of sensors
pD = 0.85 * ones(1,nb_sensors);   % sensor prob. of detection

%% Generate model parameters (kinematic and observation) 
model = gen_model_linGauss(pD);
%% Generate ground truth tracks
% read_truth;
truth=gen_truth(model);
%% Generate filtering parameters
filter_params = gen_filters(model);

for k=1:mcnt

%% Generate Measurements
meas{k}  =  gen_meas_linGauss(model,truth);
% Plot True tracks
% handles = plot_tracks(model,truth);
%% Run filters
disp(k);
% % 
%%% MS multi-Bernoulli Kalman filter with graph model and S-D assignment
est_ms_member2{k} = filter_ms_member_kf2(filter_params.ms_member_kf, model, meas{k}, truth);                         % run filter
fprintf('\n MS-MeMBer2 Time-average OSPA %d',mean(est_ms_member2{k}.ospa));
fprintf('\n MS-MeMBer2 Run Time  %d sec \n', est_ms_member2{k}.all_time);
% handles = plot_results(model, truth, est_ms_member2{k}, 0); 


% 
% %%% MS multi-Bernoulli Kalman filter with graph model
est_ms_member1{k} = filter_ms_member_kf1(filter_params.ms_member_kf, model, meas{k}, truth);                         % run filter
fprintf('\n MS-MeMBer1 Time-average OSPA %d',mean(est_ms_member1{k}.ospa));
fprintf('\n MS-MeMBer1 Run Time  %d sec \n', est_ms_member1{k}.all_time);
% handles = plot_results(model, truth, est_ms_member1{k}, 0); 
% 
%%% origion MS multi-Bernoulli Kalman filter
est_ms_member{k} = filter_ms_member_kf(filter_params.ms_member_kf, model, meas{k}, truth);                         % run filter
fprintf('\n MS-MeMBer0 Time-average OSPA %d',mean(est_ms_member{k}.ospa));
fprintf('\n MS-MeMBer0 Run Time  %d sec \n', est_ms_member{k}.all_time);
% handles = plot_results(model, truth, est_ms_member{k}, 0); 

end

plot_figure;
%% plot ospa with differemt sensor_n
% plot_ave_ospa;

