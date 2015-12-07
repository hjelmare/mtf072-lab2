%Main program in CFD assignment K2
clc;
clear variables;

%Declaration of scalar variables
%grid = 'coarse';
grid = 'fine';
maxDiff = 1e-3;
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
tic
while (epsilon > maxDiff)
   
    iteration = iteration + 1;
    %T = TDMA2(T,aCoeff);
    T = GaussSeidel(T,aCoeff);
    epsilon = CalcEpsilon(T,aCoeff,y);
    eps_save(iteration) = epsilon; 
    
end
eps_save(eps_save==0) = [];

time = toc;
disp([num2str(length(x)) 'x' num2str(length(y)) ' pts in ' num2str(time) ' s' ])

%Calculate gradient
[dX,dY] = CalcGradient(T,x,y);

u_orig = u;
v_orig = v;
T_orig = T;

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
quiver(x,y,u,v,5,'b')
axis equal
axis([x(1) x(end) y(1) y(end)]);
colorbar;

% Boundary conditions (green for heat flux (Dirichlet), red for Neumann)
color = 'rbg';
plot([x(1) x(end)],[y(1) y(1)],color(BC(1)+1),'LineWidth',3)
plot([x(end) x(end)],[y(1) y(end)],color(BC(2)+1),'LineWidth',3)
plot([x(1) x(end)],[y(end) y(end)],color(BC(3)+1),'LineWidth',3)
plot([x(1) x(1)],[y(1) y(end)],color(BC(4)+1),'LineWidth',3)
hold off
 
%saveas(gcf,['vector_gs' num2str(length(x)) 'x' num2str(length(y)) '.png'],'png')

figure;
plot(1:iteration,eps_save);

% Plot temp along wall
% figure(2);
% plot(y,T(:,end)')
% xlabel('y','FontSize',12)
% ylabel('T','FontSize',12)

%saveas(gcf,['temp_gs' num2str(length(x)) 'x' num2str(length(y)) '.png'],'png')


% Conservation

absDiffusion = kFactor * sum(abs(deltaY(2:end-1) .* -dX(:,2)'));
diffusion = kFactor * sum(deltaY(2:end-1) .* -dX(:,2)');

absConvection = c_p*rho * sum(abs(deltaY(1:end) .* u_orig(:,1)' .* T_orig(:,2)'));
convection = c_p*rho * sum(deltaY(1:end) .* u_orig(:,1)' .* T_orig(:,2)');

error = (diffusion + convection) / (absDiffusion + absConvection)



