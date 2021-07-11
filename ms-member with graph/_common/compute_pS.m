function pS = compute_pS(model, X)

if isempty(X)
    pS = [];
else
    pS = model.kin.pS * ones(size(X,2),1);
    pS = pS(:);
end
