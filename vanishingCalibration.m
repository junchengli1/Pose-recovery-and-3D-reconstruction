function [K] = vanishingCalibration(im)
% Calibrate camera from three vanishing points
% Assumption: no pixel distortion, just: 
% - focal length f
% - camera center projection: 3x1 homogeneous vector c = [u0;v0;1]

% Click on points and intersect lines to locate vanishing points
nPoints = 4; % number of points to estimate a line
nLines = 2;  % numer of lines to estimate a vanishing point

v = zeros(3,3); % three vanishing points vi = v(:,i)

L = {};
for v_num = 1:3               % loop through vanishing pts to estimate
    L{v_num} = zeros(3,nLines);
    for line_num = 1:nLines   % loop through lines intersecting at
                              % a given vanishing point
        disp(['Line ' num2str(line_num) ' going through v' num2str(v_num)]);
        L{v_num}(:,line_num) = fitLine(im,nPoints);
    end
    v(:,v_num) = findIntersection(L{v_num});
end


% Focal length estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want to solve a linear system for s = 1/f^2

% Your code goes here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% End of your code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Projection center estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute [u0;v0;1] as the  orthocenter of v1,v2,v3
% Use orthogonality equations to define a least-square problem and
% solve it

% Your code goes here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l3=cross(v(1,:),v(2,:))
l2=cross(v(1,:),v(3,:))
l1=cross(v(2,:),v(3,:))
l3=cross(v(1,:),v(2,:))/sqrt(l3(1)^2+l3(2)^2)
l2=cross(v(1,:),v(3,:))/sqrt(l2(1)^2+l2(2)^2)
l1=cross(v(2,:),v(3,:))/sqrt(l1(1)^2+l1(2)^2)
v3p=[0,-1,v(2,3);1,0,-v(1,3);-v(2,3),v(1,3),0];
v2p=[0,-1,v(2,2);1,0,-v(1,2);-v(2,2),v(1,2),0];
v1p=[0,-1,v(2,1);1,0,-v(1,1);-v(2,1),v(1,1),0];

X1=l1*v1p
X2=l2*v2p
X3=l3*v3p
X=[X1',X2',X3']
[U,S,V] = svd(X);
u = reshape(V(:,end),3,1)';
  u=[u(1)/u(3);u(2)/u(3);1]

l3=cross(v(1,:),v(2,:))
l2=cross(v(1,:),v(3,:))
l1=cross(v(2,:),v(3,:))
l3=cross(v(1,:),v(2,:))/sqrt(l3(1)^2+l3(2)^2)
l2=cross(v(1,:),v(3,:))/sqrt(l2(1)^2+l2(2)^2)
l1=cross(v(2,:),v(3,:))/sqrt(l1(1)^2+l1(2)^2)



u0=316;
v0=220;
c =[u0;v0;1]  ; % Vector [u0;v0;1]
uc=v(1,3);
vc=v(2,3);
ua=v(1,1);
va=v(2,1);
ub=v(1,2);
vb=v(2,2);
a=inv([uc-ua,vc-va;uc-ub,vc-vb])*[ub*(uc-ua)+vb*(vc-va);ua*(uc-ub)+va*(vc-vb)];
u0=a(1);
v0=a(2);
c =[u0;v0;1]  ; % Vector [u0;v0;1]
a1=(v(1,1)*v(1,2)+v(2,1)*v(2,2)-u0*(v(1,1)+v(1,2))-v0*(v(2,1)+v(2,2))+(u0^2+v0^2));
a2=(v(1,3)*v(1,2)+v(2,3)*v(2,2)-u0*(v(1,3)+v(1,2))-v0*(v(2,3)+v(2,2))+(u0^2+v0^2));
a3=(v(1,1)*v(1,3)+v(2,1)*v(2,3)-u0*(v(1,1)+v(1,3))-v0*(v(2,1)+v(2,3))+(u0^2+v0^2));

a1=[v(1,1)*v(1,2)+v(2,1)*v(2,2),(v(1,1)+v(1,2)),(v(2,1)+v(2,2)),1];
a2=[v(1,3)*v(1,2)+v(2,3)*v(2,2),(v(1,3)+v(1,2)),(v(2,3)+v(2,2)),1];
a3=[v(1,1)*v(1,3)+v(2,1)*v(2,3),(v(1,1)+v(1,3)),(v(2,1)+v(2,3)),1];
a=[a1;a2;a3]
[U,S,V] = svd(a);
B=V(:,end)
% B=[B(1),0,B(2);0,B(1),B(3);B(2),B(3),B(4)]
% L = chol(B,'lower')'
% inv(L)
Py=-B(3)/B(1)
px=-B(2)/B(1)
f=sqrt(B(4)/B(1)-(py^2+px^2))
s=a\[-1;-1;-1]
f = 1 / sqrt(s)%%%
;
% End of your code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = [f 0 c(1); 
     0 f c(2); 
     0 0 1];