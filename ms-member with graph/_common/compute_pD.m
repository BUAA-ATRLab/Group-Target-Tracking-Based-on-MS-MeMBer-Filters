function pD = compute_pD(model,X, s_indx)

if isempty(X)
    pD = [];
else
pD = model.obs.pD(s_indx) * ones(1,size(X,2));
pD = pD(:);
end
