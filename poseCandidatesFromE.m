function [transfoCandidates] = poseCandidatesFromE(E)
% Return the 4 possible transformations for an input matrix E
% transfoCandidates(i).T is the 3x1 translation
% transfoCandidates(i).R is the 3x3 rotation

transfoCandidates = repmat(struct('T',[],'R',[]),[4 1]);
% Fill in the twisted pair for E and the twisted pair for -E
% The order does not matter.

% Your code goes here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[U,S,V] = svd(E);
W = [0 -1 0; 1 0 0; 0 0 1];
 t1 = U(:,3);
 t2 = -U(:,3);
 t3 = U(:,3);
 t4 = -U(:,3);
R1 = U*W*V';
R2 = U*W*V';
R3 = U*W'*V';
R4 = U*W'*V';
if det(R1) < 0
   t1=-t1; R1 = -R1;
end
if det(R2) < 0
   t2=-t2; R2 = -R2;
end
if det(R3) < 0
   t3=-t3; R3 = -R3;
end
if det(R4) < 0
   t4=-t4; R4 = -R4;
end

transfoCandidates(1).T=t1;
transfoCandidates(1).R=R1;
transfoCandidates(2).T=t2;
transfoCandidates(2).R=R2;
transfoCandidates(3).T=t3;
transfoCandidates(3).R=R3;
transfoCandidates(4).T=t4;
transfoCandidates(4).R=R4;
% [U,S,V] = svd(E);
% T_e = U(:,end);
% R_z1 = [0 -1 0;1 0 0;0 0 1];
% R_z2 = [0 1 0;-1 0 0;0 0 1];
% %%%T%%%
% transfoCandidates(1).T = T_e;
% transfoCandidates(2).T = -T_e;
% transfoCandidates(3).T = T_e;
% transfoCandidates(4).T = -T_e;
% 
% %%%R%%%
% transfoCandidates(1).R = U*R_z1'*V';
% transfoCandidates(2).R = U*R_z2'*V';
% transfoCandidates(3).R = U*R_z2'*V';
% transfoCandidates(4).R = U*R_z1'*V';
end
% 
% 
% % End of your code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

