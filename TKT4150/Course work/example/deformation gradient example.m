clear all
close all
% See deformation gradient example on: 
% http://www.continuummechanics.org/deformationgradient.html

% initial geometry
Xc = [-2,-2; ...
     -2,2; ...
      2,2; ...
      2,-2];

% deformed geometry  
xc = [4,2; ...
     2,6; ...
     8,6; ...
     6,2];

% deformation map
x = @(X) [X(1) + 0.25*X(1)*X(2) + 5, X(2) + 4]';
  
% deformation gradient
F = @(X)[1+0.25*X(2), 0.25*X(1); 0, 1];

disp('=================== Xa = [0,2], Xb = [2,2] =====================')
Xa = [0,2]';
Xb = [2,2]';
xa = x(Xa);
xb = x(Xb);
dx = F(Xa)*(Xb-Xa);
U = sqrt(F(Xa)'*F(Xa));
R = F(Xa)*inv(U);

figure(1)
ax = gca; 
co = get(ax, 'ColorOrder');
hold all
plot(ax,Xc(:,1),Xc(:,2),'b.','MarkerSize',15)
plot(ax,xc(:,1),xc(:,2),'r.','MarkerSize',15)
quiver(ax,Xa(1),Xa(2),Xb(1)-Xa(1),Xb(2)-Xa(2),0,'LineWidth',1,'color',co(1,:),'AutoScale','on','MaxHeadSize',0.5)
quiver(ax,xa(1),xa(2),xb(1)-xa(1),xb(2)-xa(2),0,'LineWidth',1,'color',co(2,:),'AutoScale','on','MaxHeadSize',0.5)
xlim(ax,[-4,12])
ylim(ax,[-4,10])

disp('dx = F(Xa)*dX = ')
disp(dx)
disp('F(Xa) = ')
disp(F(Xa))
disp('U(Xa)')
disp(U)
disp('R(Xa) = ')
disp(R)

disp('=================== Xa = [2,0], Xb = [2,2] =====================')
Xa = [2,0]';
Xb = [2,2]';
dX = Xb - Xa;
xa = x(Xa);
xb = x(Xb);
dx = F(Xa)*(Xb-Xa);
U = sqrt(F(Xa)'*F(Xa));
R = F(Xa)*inv(U);

quiver(ax,Xa(1),Xa(2),Xb(1)-Xa(1),Xb(2)-Xa(2),0,'LineWidth',1,'color',co(1,:),'AutoScale','on','MaxHeadSize',0.5)
quiver(ax,xa(1),xa(2),xb(1)-xa(1),xb(2)-xa(2),0,'LineWidth',1,'color',co(2,:),'AutoScale','on','MaxHeadSize',0.5)

disp('dx = F(Xa)*dX = ')
disp(dx)
disp('F(Xa) = ')
disp(F(Xa))
disp('U(Xa)')
disp(U)
disp('R(Xa) = ')
disp(R)