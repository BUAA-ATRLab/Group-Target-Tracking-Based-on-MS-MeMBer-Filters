function filter_params = gen_filters(model)


%% MS-MeMBer Gaussian mixture (Kalman filter) filter parameters
filter_params.ms_member_kf.title              = 'MS MeMBer KF';
filter_params.ms_member_kf.comp_max           = 100;         % maximum number of components
filter_params.ms_member_kf.nb_comp_tg         = 4;           % nb of components per target, used for pruning
filter_params.ms_member_kf.comp_threshold     = 1e-3;        % pruning threshold for components/hypotheses
filter_params.ms_member_kf.W_max              = 4;           % maximum number of subsets
filter_params.ms_member_kf.P_max              = 4;           % maximum number of partitions
filter_params.ms_member_kf.SD_max              = 4;
filter_params.ms_member_kf.window_len = 3;                   %MDA滑窗长度3,6
filter_params.ms_member_kf.run_flag           = 'silence';   %'disp' or 'silence' for on the fly output
filter_params.ms_member_kf.epsilonx=15;  
filter_params.ms_member_kf.epsilonv=10;

filter_params.ms_member_kf.P_G= 0.9999999;                           %gate size in percentage
filter_params.ms_member_kf.gamma= chi2inv(filter_params.ms_member_kf.P_G,model.obs.z_dim);

filter_params.ms_member_kf.sGrasp = struct( ...
        'alpha', 0.05, ... % relax factor, 0 < alpha < 1         %基于最大权重独立集问题的全局假设生成时的参数设置
        'ntup', 10, ... % maximum number of tuples 
        'niter', 2, ... % number of randomized iteration for each tuple
        'useSearch', true, ... % ture or false
        'sortType', 'weight'); % 'weight-degree' or 'weight'



end