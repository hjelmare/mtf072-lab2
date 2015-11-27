%Main program in CFD assignment K1
clc;
clear variables;

%Declaration of scalar variables
maxDiff = 1e-3;
initialT = 5;
kFactor = 1;
gamma = 1/500;
L = 1;
H = 1;
T1 = 10;
T2 = 20;
T3 = 5;
T4 = 10;
c1 = 20;
c2 = 0.2;
BC = [2 2 0 2];

% Non-uniform grid
edgesX = dlmread('data/xc.dat')';
edgesY = dlmread('data/yc.dat')';


tic
%Initializing mesh and temperature
[T, y, x] = initializeMesh(edgesY,edgesX,T1,T2,T3,T4);
deltaX = diff(edgesX);
deltaX = [1 deltaX 1];
deltaY = diff(edgesY);
deltaY = [1 deltaY 1];

%Pre-calculating coefficients
aCoeff = CalcCoefficients(T,x,y,deltaX,deltaY,gamma,BC,kFactor);

%Gauss-Seidel loop
epsilon = inf;
while (epsilon > maxDiff)
   
    T = GaussSeidel(T,aCoeff);
    
    epsilon = CalcEpsilon(T,aCoeff);
    
end

%Calculate gradient
[dX,dY] = CalcGradient(T,x,y);

T = T(2:end-1,2:end-1);
[xMesh,yMesh] = meshgrid(x(2:end-1),y(2:end-1));

%Plotting result
figure(1);
contourf(xMesh,yMesh,T,20);
hold on
quiver(x(2:2:end-1),y(2:2:end-1),-dX(1:2:end,1:2:end),-dY(1:2:end,1:2:end),'r','AutoScaleFactor',5);
axis equal
axis([x(2) x(end-1) y(2) y(end-1)]);
colorbar;

% Plot grid points
%plot(xMesh,yMesh,'r.')


% Boundary conditions (green for heat flux (Dirichlet), red for Neumann)
plot([x(2) x(end-1)],[y(2) y(2)],'g','LineWidth',3)
plot([x(end-1) x(end-1)],[y(2) y(end-1)],'g','LineWidth',3)
plot([x(2) x(end-1)],[y(end-1) y(end-1)],'r','LineWidth',3)
plot([x(2) x(2)],[y(2) y(end-1)],'g','LineWidth',3)
hold off

time = toc;
disp([num2str(length(x)) 'x' num2str(length(y)) ' pts in ' num2str(time) ' s' ])
 
saveas(gcf,['vector' num2str(length(x)) 'x' num2str(length(y)) '.png'],'png')

