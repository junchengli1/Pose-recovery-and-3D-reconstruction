function E = estimateEmatrix(X1,X2)
% Estimate E matrix given a set of 
% pairs of matching *calibrated* points
% X1,X2: Nx2 matrices of calibrated points
%   i^th row of X1 matches i^th row of X2

% Kronecker products
% Your code goes here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [X1(:,1).*X2(:,1), X1(:,1).*X2(:,2), X1(:,1), ...
         X1(:,2).*X2(:,1), X1(:,2).*X2(:,2), X1(:,2), ...
         X2(:,1),X2(:,2), ones(size(X1,1),1)];


[U,S,V] = svd(A);
% V=V/V(end);
E =reshape(V(:,end),3,3);
[U,S,V] = svd(E);
E=U*diag([1 1 0])*V';
% Project E on the space of essential matrices



% End of your code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
