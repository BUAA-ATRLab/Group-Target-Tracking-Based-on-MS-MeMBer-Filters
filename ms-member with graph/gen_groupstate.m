function X=gen_groupstate(Xd, model)
    X= zeros(size(Xd));
    %-- short hand
    L= size(Xd,2);
    T= model.kin.T; 
    omega= Xd(5,:);
    %-- pre calcs
    sin_omega_T= sin(omega*T);
    cos_omega_T= cos(omega*T);
    %--  x/y pos/
    vmean=mean(Xd([2 4],:),2);
    X([1 3],:)= Xd([1 3],:)+repmat(vmean,1,L)*T;
%     X([2 4],:)= Xd([2 4],:);
    X(2,:)= cos_omega_T.*Xd(2,:)- sin_omega_T.*Xd(4,:);
    X(4,:)= sin_omega_T.*Xd(2,:)+ cos_omega_T.*Xd(4,:);
    %-- turn rate
    X(5,:)= Xd(5,:);
end