function [U,S,V] = svdecon(X)
% Input:
% X : m x n matrix
%
% Output:
% X = U*S*V'
%
% Description:
% Does equivalent to svd(X,'econ') but faster
%
% Vipin Vijayan (2014)

% X = bsxfun(@minus,X,mean(X,2));
[m,n] = size(X);

if  m <= n
    C = X*X';
    [U,D] = eig(C, 'vector');
    clear C;
    
    [d,ix] = sort(abs(D),'descend');
    U = U(:,ix);    
    
    if nargout > 2
        V = X'*U;
        s = sqrt(d);
        V = V./s';
        S = diag(s);
    else
        s = sqrt(d);
        S = diag(s);
    end
else
    C = X'*X;  % calculate covariance matrix - Real symmetric matrix
    [V,D] = eig(C);  % calc eig values and vectors of C
    clear C;
    
    [d,ix] = sort(abs(diag(D)),'descend');  % sort highest energy first
    V = V(:,ix);    % sort V as well
    
    U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
    s = sqrt(d); % calculate singular values of original X matrix = lambda^(1/2)
    U = U./s'; % convert orthogonal matrix to orthonormal matrix
    S = diag(s);  % makes into diag matrix form
end
