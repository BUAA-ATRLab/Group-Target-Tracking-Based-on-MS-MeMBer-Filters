function meas = gen_meas_linGauss(model, truth, varargin)


%variables
meas.K = truth.K;
meas.Z = cell(truth.K, 1);

% set up the Mersenne Twister generator for the repeatability of the simulations
if (~isempty(varargin)) && (~isempty(varargin{1}))
    stream = RandStream('mt19937ar', 'seed', varargin{1});
    RandStream.setGlobalStream(stream);
end

% generate measurements
for k = 1:truth.K   
    
    X = truth.X{k};
    Z = cell(1,model.obs.nb_sensors);
    nb_targ = size(X,2);
    
    Zt={};
    Zc={};
    
    for s = 1:model.obs.nb_sensors        
        % Generate target measurements
        if ~isempty(X)
            Wl = sqrtm(model.obs.R) * randn(model.obs.z_dim, nb_targ);
%             Wl = sqrtm(model.obs.R) * zeros(model.obs.z_dim, nb_targ);
            Zl = model.obs.gen_sen_obs( s, X, Wl);
   
            idx = find(rand(1,nb_targ) <= model.obs.pD(s));
            Z{s} = [];
            if ~isempty(idx)
                Z{s} = Zl(:,idx);
            end
        end
        
        % Generate clutter for sensor s
        N_c       = poissrnd(model.obs.lambda_c);                                                                  % number of clutter points
        C         = repmat(model.obs.range_c(:,1),[1 N_c])+ diag(model.obs.range_c*[ -1; 1 ])*rand(model.obs.z_dim, N_c);  % clutter generation
        Zc{s}= C;
        Zt{s}= Z{s};
        % append clutter measurements to target measurements
        Zord = [ Z{s} C ];
        % Radomize order of target-originated and clutter measuremnts
        nb_meas = size(Zord,2); 
        if nb_meas == 0
            Z{s} = [];
        else 
            Zrand = Zord(:, randperm(nb_meas));
            Z{s} = Zrand;
        end      
    end % for s 
    meas.Z{k} = Z;
    Ztarget{k} = Zt;
    Zclutter{k} = Zc;
    
end % for k
    

% figure;
% for k=1:meas.K
%     if ~isempty(Zclutter{k}{1})
%         plot(Zclutter{k}{1}(1,:),Zclutter{k}{1}(2,:),'.','Color',0.6*ones(1,3));hold on
%     end
% end
% for k=1:meas.K
%     if ~isempty(Ztarget{k}{1})
%         plot(Ztarget{k}{1}(1,:),Ztarget{k}{1}(2,:),'k.','MarkerSize',5);hold on
%     end
% end
% axis equal
% axis([-1500 1500 0 3000]);
% xlabel('x(m)'),ylabel('y(m)');




end
    
