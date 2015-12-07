%Main program in CFD assignment K2
clc;
clear variables;

%Declaration of scalar variables
%grid = 'coarse';
grid = 'fine';
maxDiff = 1e-3;
maxDiff2 = 1e-3;
kFactor = 1;
rho = 1;
c_p = 500;
gamma = 1/c_p;
L = 1;
H = 1;
T1 = 10;
T2 = 0;
T3 = 0;
T4 = 10;
Ta = 20;
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

%Finding inlet/outlet
inletIndex = find(u(:,1)>0.1);
outletIndex = find(u(:,1)<0);


%Initializing mesh and temperature
[T, y, x] = initializeMesh(edgesY,edgesX,T1,T2,T3,T4);
T(inletIndex,1) = Ta;
deltaX = diff(edgesX);
deltaX = [1 deltaX 1];
deltaY = diff(edgesY);
deltaY = [1 deltaY 1];

%Pre-calculating coefficients
aCoeff = CalcCoefficients3(T,x,y,u,v,rho,deltaX,deltaY,gamma,BC,kFactor,...
    inletIndex,outletIndex);

%Gauss-Seidel/TMDA loop
epsilon = inf;
eps_save = zeros(10000,1);
iteration = 0;

while (epsilon > maxDiff)
   
    iteration = iteration + 1;
    %T = GaussSeidel(T,aCoeff);
    T = TDMA2(T,aCoeff);
    epsilon = CalcEpsilon(T,aCoeff,y);
    eps_save(iteration) = epsilon; 
    
end
eps_save(eps_save==0) = [];
T_GS = T;

%Initializing mesh and temperature
[T, y, x] = initializeMesh(edgesY,edgesX,T1,T2,T3,T4);
T(inletIndex,1) = Ta;
deltaX = diff(edgesX);
deltaX = [1 deltaX 1];
deltaY = diff(edgesY);
deltaY = [1 deltaY 1];

epsilon = inf;
eps_save = zeros(10000,1);
iteration = 0;

while (epsilon > maxDiff2)
   
    iteration = iteration + 1;
    T = GaussSeidel(T,aCoeff);
    %T = TDMA2(T,aCoeff);
    epsilon = CalcEpsilon(T,aCoeff,y);
    eps_save(iteration) = epsilon; 
    
end
eps_save(eps_save==0) = [];
T_TDMA = T;


%Deleting frame/border-values
T_TDMA = T_TDMA(2:end-1,2:end-1);
T_GS = T_GS(2:end-1,2:end-1);
u = u(2:end-1,2:end-1);
v = v(2:end-1,2:end-1);
x = x(2:end-1);
y = y(2:end-1);
[xMesh,yMesh] = meshgrid(x,y);

T_diff = T_TDMA - T_GS;

%Plotting result 
figure(1);
contourf(xMesh,yMesh,T_diff,20);
hold on
quiver(x,y,u,v,5,'b')
axis equal
axis([x(1) x(end) y(1) y(end)]);
colorbar;




