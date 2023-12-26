function [si_trunc,sigma_trunc,phi_trunc,NumBases] = truncated_matrix(si,sigma,phi,tol)
%  TRUNCATED_MATRIX Brief summary of this function.
% 
% Detailed explanation of this function.
% Input sigma is a diagonal matrix obtained from the svd of snapshots.
    [~,c]=size(sigma);
    if c==1
        sigma=sigma(1);
    else
        sigma=diag(sigma);
    end
%     sigma=diag(sigma);
    sum_sigma=sum(sigma.^2);
    a1=sigma(1).^2;
    a= sigma(1).^2/sum_sigma;
    ii=1;
    while a<tol
        a1=a1+sigma(ii+1).^2;
        a=a1/sum_sigma;
        ii=ii+1;
    end
    sigma_trunc=sigma(1:ii);
    si_trunc=si(:,1:ii);
    phi_trunc=phi(:,1:ii);
    NumBases=ii;    % number of basis retained after truncation 
end