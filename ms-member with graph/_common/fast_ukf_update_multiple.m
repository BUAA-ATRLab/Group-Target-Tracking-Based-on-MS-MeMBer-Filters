function [qz_update, m_update, P_update, eta_update, K_update, S_update] = fast_ukf_update_multiple(Z, gen_sensor_obs, m, P, model, ukf_params)
% Perform UKF updates on multiple Gaussian components. The update equation
% involves sensor coordinates specified model and the current sensor index s

nb_comp = size(m,2);
nb_meas = size(Z,2);

qz_update  = zeros(nb_meas, nb_comp);
m_update   = zeros(model.kin.x_dim, nb_comp*nb_meas);
eta_update = zeros(model.obs.z_dim, nb_comp);
S_update   = zeros(model.obs.z_dim, model.obs.z_dim, nb_comp);
K_update   = zeros(model.kin.x_dim, model.obs.z_dim, nb_comp);
P_update   = zeros(model.kin.x_dim, model.kin.x_dim, nb_comp);

% process each Gaussian component 
for i = 1:nb_comp
    % Single component update

    [eta_comp, K_comp, P_comp, S_comp] = ukf_update_single( gen_sensor_obs, m(:,i), P(:,:,i), model, ukf_params);
        
    eta_update(:,i) = eta_comp;
    K_update(:,:,i) = K_comp;
    S_update(:,:,i) = S_comp;
    P_update(:,:,i) = P_comp;
    
    
    % Evaluate the likelihood of all pairings between measurements versus all predicted componets 
    qz_update(:,i) = exp( model.obs.lik_fn(Z, repmat(eta_update(:,i),1,nb_meas), S_update(:,:,i)) )';
    
    
    % form all pairings of updated component means for all measurements versus all predicted componets 
    innov = model.obs.angle_wrap(Z - repmat(eta_comp,1,nb_meas));
    m_update(:,i + (0:nb_meas-1)*nb_comp) = repmat(m(:,i),1,nb_meas) + K_comp * innov;
    
end

end

function [eta, K, P, S] = ukf_update_single( gen_sensor_obs, m_comp, P_comp, model, ukf)
% Standard unscented 

% Generate the sigma points
% Each column in ukf_m is the concatenation of a target state vector
% and a dimMeax1 column vector of zeros
ukf_m = [m_comp; zeros(model.obs.z_dim, 1)];
ukf_P = blkdiag(P_comp, model.obs.R);

temp = sqrtm( (ukf.L + ukf.lambda) * ukf_P );
ukf_y = [ukf_m, bsxfun(@plus, ukf_m, temp), bsxfun(@minus, ukf_m, temp)];

% We split ukf_y into two matrices
% ukf_x is the x_dim x noSigma matrix that contains the state
% vector part of each sigma point
% ukf_x thus represents our sigma points
% ukf_e is the z_dim x noSigma matrix where each column
% contains the measurement noise for the corresponding sigma
% point
ukf_x = ukf_y(1:model.kin.x_dim, :);
ukf_e = ukf_y(model.kin.x_dim + 1:end,:);

% Construct one set of sigma points for each possible target
ukf_z = gen_sensor_obs(ukf_x) + ukf_e;

% Compute predicted measurement for component i
eta = ukf_z * ukf.Wm';  

% Computed covariance matrix
tempDif = bsxfun(@minus, ukf_z, eta);
temp = bsxfun(@times, ukf.Wc, tempDif);
S = temp*tempDif' + model.obs.R;

% Compute inverse of curr_S(:,:,i)
inverseS = S \ eye(model.obs.z_dim);

temp = bsxfun(@times, ukf.Wc, bsxfun(@minus, ukf_x, m_comp));
G = temp*tempDif';
K = G*inverseS;
P = P_comp - G*inverseS*G';
          
end



