function [E, bestInliers] = estimateEmatrixRANSAC(X1,X2)
% Estimate E matrix given a set of
% pairs of matching *calibrated* points
% X1,X2: Nx2 matrices of calibrated points
%   i^th row of X1 matches i^th row of X2
%
% E: robustly estimated E matrix
% bestInliers: indices of the rows of X1 (and X2) that where in the
% largest consensus set

nIterations = 3000;
sampleSize = 8;

%fractionInliers = 0.6;
%nInliers = floor((size(X1,1) - sampleSize) * fractionInliers);
%bestError = Inf;
eps = 10^(-5);
bestNInliers = 0;




for i=1:nIterations
    indices = randperm(size(X1,1));
    sampleInd = indices(1:sampleSize);
    testInd =  indices(sampleSize+1:length(indices));
    % Your code goes here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x1_sample = X1(sampleInd,1:2);
    x2_sample = X2(sampleInd,1:2);
    E_sample = estimateEmatrix(x1_sample, x2_sample);
    X1_R = [X1(testInd, :), ones(size(testInd,2),1)];
    X2_R = [X2(testInd, :), ones(size(testInd,2),1)];
    e_hat =[0, -1, 0;
        1, 0, 0;
        0, 0, 0];  
    e2 = (E_sample*X1_R')'; % N*3
    TOP12 = sum((X2_R .* e2),2).^2; % sum(N*3) 
    BOT12 = (sum((e_hat*(e2)').^2, 1));
    d1 = TOP12'./BOT12;  
    e1 = (X2_R*E_sample); %N*3
    TOP21 = sum((X1_R .* e1),2).^2; %sum(N*3)
    BOT21 = (sum((e_hat*(e1)').^2, 1));
    d2 = TOP21'./BOT21; 
    residuals = d1 + d2;
    L = residuals < eps;
    curInliers = [sampleInd, testInd(L)];    % don't forget to include the sampleInd
    
    % End of your code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    curNInliers = length(curInliers);
    
    if curNInliers > bestNInliers
        bestNInliers = curNInliers;
        bestInliers = curInliers;
        E = E_sample;
    end
end

disp(['Best number of inliers: ' bestNInliers '/' size(X1,1)]);
