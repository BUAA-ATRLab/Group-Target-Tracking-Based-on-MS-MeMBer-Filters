function Z = gen_observation_fn_linGauss(model, s, X, W)


if isempty(X)
    Z = [];
else % modify below here for user specified measurement model
    Z = model.obs.H * X + W;
end

end