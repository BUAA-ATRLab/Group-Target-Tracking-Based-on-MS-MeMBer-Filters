function model = gen_model_linGauss( pD )



% basic parameters
model.kin.x_dim = 4;   %dimension of state vector
model.obs.z_dim = 2;   %dimension of observation vector
model.kin.v_dim = 4;   %dimension of process noise
model.obs.w_dim = 2;   %dimension of observation noise


%% dynamical model parameters (CT model)
% state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
model.kin.T = 1;                         %sampling period
model.kin.F = [eye(2,2) model.kin.T*eye(2,2); zeros(2,2) eye(2,2)]; % state transition matrix
model.kin.F = [eye(2,2) model.kin.T*eye(2,2); zeros(2,2) eye(2,2)]; % group state transition matrix
model.kin.sigma_v = 3;
model.kin.Q = model.kin.sigma_v^2 * [(model.kin.T^3/3)*eye(2,2) (model.kin.T^2/2)*eye(2,2); (model.kin.T^2/2)*eye(2,2) (model.kin.T)*eye(2,2)];   
model.kin.sqrt_Q = sqrtm(model.kin.Q);

model.kin.prop_fn = @(x,v) gen_newstate_fn(model,x,v);                     % define the target state transition kernel, i.e., f_{k+1}(x_{k+1} | x_k ) 

% survival/death parameters
model.kin.pS = .99;                                                        % define p_S as a constant 
model.kin.compute_ps = @(x) compute_pS(model, x);                          % define p_S as a function of the target state x

% kinematical model used for simmulating the true tracks
model.kin.sim.T = 1;                         %sampling period
model.kin.sim.F = [eye(2,2) model.kin.T*eye(2,2); zeros(2,2) eye(2,2)]; % state transition matrix
model.kin.sim.sigma_v = 0.25;
model.kin.sim.Q = model.kin.sim.sigma_v^2 * [(model.kin.T^3/3)*eye(2,2) (model.kin.T^2/2)*eye(2,2); (model.kin.T^2/2)*eye(2,2) (model.kin.T)*eye(2,2)];   
model.kin.sim.sqrt_Q = sqrtm(model.kin.sim.Q);


%% birth parameters 
model.birth.J_gamma = 11;                                        % number of births
model.birth.r_gamma = 0.1*ones(1,model.birth.J_gamma);                        % weights of birth components
model.birth.P_gamma = diag([100 100 100 100]);                      % birth covariances
model.birth.mix_w_gamma = ones(1,model.birth.J_gamma);                            % mixture weights for birth pdf (here we only take one Gaussian per component)
model.birth.sqrt_P_gamma = sqrtm(model.birth.P_gamma);                                                                % squared root of birth covariance matrix
model.birth.m_gamma(:,1)=[ -200+30;2000-20;0;0];
model.birth.m_gamma(:,2)=[ -200+30; 2000+20; 0; 0];
model.birth.m_gamma(:,3)=[ -200; 2000; 0; 0];
model.birth.m_gamma(:,4)=[ -200+30; 2000+60; 0; 0];
model.birth.m_gamma(:,5)=[ 0; 500; 0; 0];
model.birth.m_gamma(:,6)=[ -1000; 300; 0; 0];
model.birth.m_gamma(:,7)=[ -1040; 300; 0; 0];
model.birth.m_gamma(:,8)=[ -1020; 320; 0; 0];
model.birth.m_gamma(:,9)=[ 300; 2000; 0; 0];
model.birth.m_gamma(:,10)=[ 300+30; 2000-30; 0; 0];
model.birth.m_gamma(:,11)=[ 300; 2000-260; 0; 0];
%% observation model prams
model.obs.nb_sensors = size(pD, 2);                                        % number of sensors
model.obs.xrange = [-1500 1500];                                           % monitoring region dimensions in meters
model.obs.yrange = [0 3000];
model.obs.H = [eye(model.obs.z_dim) zeros(model.obs.z_dim,model.obs.z_dim)];   % observation matrix      
model.obs.sigma_r   = 6;                                                  % observation noise std in meters
model.obs.R = model.obs.sigma_r^2 * eye(model.obs.z_dim);                  % measurement noise covariance
model.obs.R1 = (1*model.obs.sigma_r)^2 * eye(model.obs.z_dim);             %¡ø≤‚¬À≥˝√≈œﬁ

% detection parameters
model.obs.pD = pD;                                                         % define p_D as a constant

model.obs.gen_sen_obs = @(s_idx, x, w)  gen_observation_fn_linGauss(model, s_idx, x, w);   % define the observation function
model.obs.lik_fn      = @(Z, Hx, R)  log( mvnpdf((Z - Hx)', zeros(1, model.obs.z_dim), R ) )';


% clutter parameters
model.obs.lambda_c = 100;                                   %poisson average rate of clutter (per scan) and identical for all sensors (this could be made sensor dependent in future releases)
model.obs.range_c  = [ -1500 1500; 0 3000 ];          %uniform clutter on XY
model.obs.pdf_c    = 1/prod(model.obs.range_c(:,2)-model.obs.range_c(:,1)); %uniform clutter density

% Plot parameters
model.plot.xrange = [-1500 1500]; 
model.plot.yrange = [0 3000];

% OSPA parameters
model.ospa_cutoff = 100;
model.ospa_order  = 1;


end


