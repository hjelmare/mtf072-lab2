%Main program in CFD assignment K2
clc;
clear variables;

%Declaration of scalar variables
%grid = 'coarse';
grid = 'fine';
maxDiff = 1e-3;
kFactor = 1;
rho = 1;
gamma = 1/500;
L = 1;
H = 1;
T1 = 10;
T2 = 0;
T3 = 0;
T4 = 10;
Ta = 20;
ha = 1.97;
BC = [0 0 0 2];
%BC = [2 0 0 2];

% Loading grid and velocity data
edgesX = dlmread(['data/grid2/' grid '_grid/xc.dat'])';
edgesY = dlmread(['data/grid2/' grid '_grid/yc.dat'])';
u = dlmread(['data/grid2/' grid '_grid/u.dat'])';
v = dlmread(['data/grid2/' grid '_grid/v.dat'])';

%Reshapeing velocity vectors
Nx = length(edgesX) + 1;
Ny = length(edgesY) + 1;
u = reshape(u,Nx,Ny);
v = reshape(v,Nx,Ny);

tic
%Initializing mesh and temperature
[T, y, x] = initializeMesh(edgesY,edgesX,T1,T2,T3,T4);
T(y>ha,1) = Ta;
deltaX = diff(edgesX);
deltaX = [1 deltaX 1];
deltaY = diff(edgesY);
deltaY = [1 deltaY 1];

%Pre-calculating coefficients
aCoeff = CalcCoefficients(T,x,y,u,v,rho,deltaX,deltaY,gamma,BC,kFactor);

%Gauss-Seidel/TMDA loop
epsilon = inf;
while (epsilon > maxDiff)
   
    %T = TDMA2(T,aCoeff);
    T = GaussSeidel(T,aCoeff);
    epsilon = CalcEpsilon(T,aCoeff,y);
    
end

%Calculate gradient
[dX,dY] = CalcGradient(T,x,y);

%Deleting frame/border-values
T = T(2:end-1,2:end-1);
u = u(2:end-1,2:end-1);
v = v(2:end-1,2:end-1);
x = x(2:end-1);
y = y(2:end-1);
[xMesh,yMesh] = meshgrid(x,y);

%Plotting result
figure(1);
contourf(xMesh,yMesh,T,20);
hold on
%quiver(x(1:2:end),y(1:2:end),-dX(1:2:end,1:2:end),-dY(1:2:end,1:2:end),'r','AutoScaleFactor',5);
quiver(x,y,u,v,5)
axis equal
axis([x(1) x(end) y(1) y(end)]);
colorbar;

% Plot grid points
%plot(xMesh,yMesh,'r.')


% Boundary conditions (green for heat flux (Dirichlet), red for Neumann)
color = 'rbg';
plot([x(1) x(end)],[y(1) y(1)],color(BC(1)+1),'LineWidth',3)
plot([x(end) x(end)],[y(1) y(end)],color(BC(2)+1),'LineWidth',3)
plot([x(1) x(end)],[y(end) y(end)],color(BC(3)+1),'LineWidth',3)
plot([x(1) x(1)],[y(1) y(end)],color(BC(4)+1),'LineWidth',3)
hold off

time = toc;
disp([num2str(length(x)) 'x' num2str(length(y)) ' pts in ' num2str(time) ' s' ])
 
saveas(gcf,['vector_gs_bc1' num2str(length(x)) 'x' num2str(length(y)) '.png'],'png')


% Plot temp along wall
figure(2);
plot(y,T(:,end)')
xlabel('y','FontSize',12)
ylabel('T','FontSize',12)

saveas(gcf,['temp_gs' num2str(length(x)) 'x' num2str(length(y)) '.png'],'png')