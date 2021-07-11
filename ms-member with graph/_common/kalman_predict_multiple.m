function [m_predict,P_predict] = kalman_predict_multiple(model, m, P)      
% Kalman predict a collection of Gaussian mixtures specified as the means m
% and covariances P. 

p_length = size(m,2);

m_predict = zeros(size(m));
P_predict = zeros(size(P));

for idxp = 1:p_length
    [m_temp, P_temp] = kalman_predict_single(model.kin.F, model.kin.Q, m(:,idxp), P(:,:,idxp));
    m_predict(:,idxp) = m_temp;
    P_predict(:,:,idxp) = P_temp;
end
end

function [m_predict,P_predict] = kalman_predict_single(F,Q,m,P)

m_predict = F * m;
P_predict = Q + F*P*F'; 

end

