function [ X_tilde, W_tilde, index ] = Resampling( X, W, N_sample, choice )
%function [ X_tilde, W_tilde, index ] = Resampling( X, W, choice )
%   Performs multinomial resampling of the input distribution ( X, W )
%   The resulting distribution is given by ( X_tilde, W_tilde ), where
%   W_tilde = 1/(number of particles). The index vector serves as a LUT 
%   (Look Up Table betwen the two particle vectors X and X_tilde.)
%   The choice character string selects either the 'multinomial', or 
%   'stratified' or 'systematic' or the 'residual'  resampling schemes.   
%   All algorithms are   
%
% Implementations based on:
%   "On resampling algorithms for particle filters" by: Jeroen Hol, Thomas 
%     Schon and Frederik Gustafsson.
%   "Comparison of resampling schemes for particle filtering" by: Randal
%     Douc, Olivier Capp� and Eric Moulines
%   "Bayesian signal processing: Classical, modern and particle filtering 
%     methods" by James V. Candy

if isempty(X)
    X_tilde = [];
    W_tilde = [];
    index   = [];
else
    
    % State vector dimension and Number of Particles
    [Xdim, NbPar] = size(X);
    index = zeros(N_sample, 1);
    
    % Make sure no weights are NaN, otherwise just return the input
    % samples
    if sum(isnan(W))
        X_tilde = X;
        W_tilde = 1/NbPar * ones(1,NbPar);
        return;
    end
    
    % Non-sorted and sorted uniformily distributed random vectors
    u_tilde    = rand(1,N_sample);
    u          = zeros(1,N_sample);
    
    % Sorting the random vector
    
    if strcmpi(choice, 'multinomial')
        
        u(NbPar)   = (u_tilde(NbPar))^(1/NbPar);
        for l = NbPar-1:-1:1,
            u(l)     = (u_tilde(l))^(1/l)*u(l+1);
        end;
        
    elseif strcmpi(choice, 'stratified')
        
        u = (0:N_sample-1) + u_tilde;
        u = u / N_sample;
        
    elseif strcmpi(choice, 'systematic')
        
        u = (0:N_sample-1) + u_tilde(1);
        u = u / N_sample;
        
    else  % Residual resampling
        
        % Repetition counts
        N_det = floor(NbPar*W);
        % Remainder
        R = sum(N_det);
        % The number of particles which will be drawn stochastically
        N_rnd = NbPar - R;
        
        % Modified weights for the random resempling part
        W_rnd = (NbPar*W - N_det)/ N_rnd;
        
        % Draw the deterministic part
        i = 1;
        for j = 1:NbPar,
            for k = 1:N_det(j)
                index(i) = j;
                i = i + 1;
            end
        end
        
        % Draw the stochastic part
        
        Frep       = [0 cumsum(W_rnd)'];
        
        while i<=NbPar
            u_rnd = rand(1,1);
            j = 1;
            while ( ((u_rnd-Frep(j))*(u_rnd-Frep(j+1))>0) && (j < NbPar))
                j = j + 1;
            end
            index(i) = j;
            i = i + 1;
        end
        
        % Resampled distributions
        X_tilde = X(:,index);
        W_tilde = 1/NbPar * ones(NbPar,1);
        
        return;
    end
    
    
    % Projecting the sorted random vec tor onto the cumulative function
    % of the input distribution
    
    % Form cumulative distribution function of the input distribution
    Frep       = [0 cumsum(W)];
    
    seuil      = 1;
    par_aux    = zeros(Xdim, N_sample);
    
    for l=1:N_sample,
        n = seuil;
        while ( ((u(l)-Frep(n))*(u(l)-Frep(n+1))>0) && (n < NbPar))
            n=n+1;
        end;
        seuil = n;
        par_aux(:,l) = X(:,seuil);
        index(l) = seuil;
    end;
    
    % Resampled distributions
    X_tilde = par_aux;
    W_tilde = 1/N_sample * ones(1,N_sample);
    
end

end