function X = gen_newstate_fn(model, Xd, V)

% Linear state space equation (NCV model)
if ~isnumeric(V)
    if strcmp(V,'noise')
        V = model.kin.sqrt_Q * randn(model.kin.x_dim, size(Xd,2));
    elseif strcmp(V,'noiseless')
        V = zeros(model.kin.x_dim, size(Xd,2));
    end
end

if isempty(Xd)
    X = [];
else  % we employ a linear trasition model
    X = model.kin.F * Xd + V;                                              % modify here for user specified transition model
end