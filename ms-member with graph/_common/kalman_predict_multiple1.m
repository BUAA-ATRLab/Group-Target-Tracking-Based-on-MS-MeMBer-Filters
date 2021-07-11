function [m_predict,P_predict] = kalman_predict_multiple1(model, m, P, G_id, Gr)      
% Kalman predict a collection of Gaussian mixtures specified as the means m
% and covariances P. 

p_length = size(m,2);

m_predict = zeros(size(m));
P_predict = zeros(size(P));

for idxp = 1:p_length
    if G_id(:,idxp)==0
        [m_temp, P_temp] = kalman_predict_single(model.kin.F, model.kin.Q, m(:,idxp), P(:,:,idxp));
    else
%         ng = length(Gr(G_id(:,idxp)).index);
%         F = [eye(2,2) model.kin.T*eye(2,2); zeros(2,2) eye(2,2)]; % group state transition matrix
%         B = [zeros(2,2) model.kin.T*eye(2,2); zeros(2,2) zeros(2,2)]; 
%         u = Gr(G_id(:,idxp)).center-m(:,idxp);
%         [m_temp, P_temp] = kalman_predict_single1(F, model.kin.Q, m(:,idxp), P(:,:,idxp), B, u);
        
        ng = length(Gr(G_id(:,idxp)).index);
        F = [eye(2,2) model.kin.T/ng*eye(2,2); zeros(2,2) eye(2,2)]; % group state transition matrix
        B = [zeros(2,2) model.kin.T*eye(2,2); zeros(2,2) zeros(2,2)]; 
        u = Gr(G_id(:,idxp)).center-1/ng*m(:,idxp);
        [m_temp, P_temp] = kalman_predict_single1(F, model.kin.Q, m(:,idxp), P(:,:,idxp), B, u);
        
%         ng = length(Gr(G_id(:,idxp)).index);
%         F = eye(4,4); % group state transition matrix
%         B = [0 0 model.kin.T 0; 0 0 0 model.kin.T;0 0 0 0;0 0 0 0]; 
%         u = Gr(G_id(:,idxp)).center;
%         [m_temp, P_temp] = kalman_predict_single1(F, model.kin.Q, m(:,idxp), P(:,:,idxp), B, u);
    end
    m_predict(:,idxp) = m_temp;
    P_predict(:,:,idxp) = P_temp;
end
end

function [m_predict,P_predict] = kalman_predict_single(F,Q,m,P)

m_predict = F * m;
P_predict = Q + F*P*F'; 

end

function [m_predict,P_predict] = kalman_predict_single1(F,Q,m,P,B,u)
% ori_v=m([3 4],:);
% m([3 4],:)=center([3 4],:);
% m_predict = F * m;
% m_predict([3 4],:)=ori_v;
% P_predict = Q + F*P*F'; 
m_predict = F*m+ B*u;
P_predict = Q + F*P*F'; 

end