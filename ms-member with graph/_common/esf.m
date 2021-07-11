function s = esf(Z)
% function s = esf(Z)
% Function that computes the elementary symmetric function using Mahler's recursive formula.
% The definition of the function can be found as eq. (11.141) at page 369
% of Mahler's book "Statistical Multisource Multitarget Information Fusion"
% Z  = input ector 
% s =  function value for the specific input

if isempty(Z)
    s= 1;
    return;
end

n_z = length(Z);
F = zeros(2,n_z);

i_n = 1;
i_nminus = 2;

for n = 1:n_z
    F(i_n,1) = F(i_nminus,1) + Z(n);
    for k = 2:n
        if k==n
            F(i_n,k) = Z(n)*F(i_nminus,k-1);
        else   
            F(i_n,k) = F(i_nminus,k) + Z(n)*F(i_nminus,k-1);
        end            
    end
    tmp = i_n;
    i_n = i_nminus;
    i_nminus = tmp;
end    

s= [1; F(i_nminus,:)']; 

end