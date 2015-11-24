%-------------------------------------------------------------------------
%                 CFD of Turbulent Flow - MTF072
%                 Task K2 - Grid #1 - Case #3
%                 Group #43 - Behrad Gharraee
%-------------------------------------------------------------------------
clc
clear all
close all
tic

T_i = 0; %Initial temp.
T_in = 20 ; %inlet temp.
s_bar = 0; % source term
kfactor=10; %for increasing or decreasing k value
error = 0.001 ; % number of iterations
k = 1 ;
c_p = 500 ;
gamma = k*kfactor/c_p ;
rho = 1;
h_in=10; %Inlet height (in number of nodes)
h_out=10; %Outlet height (in number of nodes)
bc = 0 ; %B.C. flag: Neumann=0, Dirichlet=1
solver = 0; % Solver: G-S=0 , TDMA=1

%% Begin skeleton
% read xc
load xc.dat
nim1=length(xc);
% nim1 = ni-1 = number of grid lines. Number of cell nodes = ni
ni=nim1+1;
%
% read yc
load yc.dat
njm1=length(yc);
nj=njm1+1;
%
% read u
load u.dat
u2d=reshape(u,ni,nj);
% read v
load v.dat
v2d=reshape(v,ni,nj);
%
% xc and yc are the coordinates of the grid
%
%       o---------------o  yc(j)
%       |               |
%       |               |
%       |       P       |
%       |               |
%       |               |
%       o---------------o  yc(j-1)
%
%     xc(i-1)          xc(i)          
%
% The cell-node (P) has the coordinates xp(i), yp(j)
%
% compute the x-coordinates of the cell centres
for i=2:nim1
   xp(i)=0.5*(xc(i)+xc(i-1));
end
xp(1)=xc(1);
xp(ni)=xc(nim1);
%
% take the transpose of x
xp=xp';

% compute the y-coordinates of the cell centres
for j=2:njm1
   yp(j)=0.5*(yc(j)+yc(j-1));
end
yp(1)=yc(1);
yp(nj)=yc(njm1);
%
% take the transpose of y
yp=yp';
% End of Skeleton

%% Face U values
U_facex=zeros(nj-2,ni-1); %Velocity in x-dir at faces
for j=1:nj-2
    for i=1:ni-3
        U_facex(j,i+1)=interp1(xp,u2d(j+1,:),xc(i+1));
    end
end
U_facex(:,1)=u2d(2:nj-1,1);
U_facex(:,ni-1)=u2d(2:nj-1,ni);
U_w = U_facex(:,1:ni-2);
U_e = U_facex(:,2:ni-1);
 
V_facey=zeros(nj-1,ni-2); %Velocity in y-dir at faces
for j=1:ni-2
    for i=1:nj-3
        V_facey(i+1,j)=interp1(yp,v2d(:,j+1),yc(i+1));
    end
end
V_facey(1,:)=v2d(1,2:ni-1);
V_facey(nj-1,:)=v2d(nj,2:ni-1);
V_s = V_facey(1:nj-2,:);
V_n = V_facey(2:nj-1,:);

%% Distances between neighboring nodes
for i=2:ni-1
    deltax_WP(i-1)=xp(i)-xp(i-1);
end
for i=2:ni-1
    deltax_PE(i-1)=xp(i+1)-xp(i);
end
for i=2:nj-1
    deltay_SP(i-1)=yp(i)-yp(i-1);
end
for i=2:nj-1
    deltay_PN(i-1)=yp(i+1)-yp(i);
end

%% Values for "A"
for i=1:nj-2
    A_w(i)=yc(i+1)-yc(i);
end
A_e=A_w;
for i=1:ni-2
    A_s(i)=xc(i+1)-xc(i);
end
A_n=A_s;

%% Diffusion Terms
for i=1:ni-2
    for j=1:nj-2
        D_w(i,j)=gamma*A_w(i)./deltax_WP(j);
    end
end
    
for i=1:ni-2
    for j=1:nj-2
        D_e(i,j)=gamma*A_e(i)./deltax_PE(j);
    end
end

for j=1:nj-2
    for i=1:ni-2
        D_s(i,j)=gamma*A_s(j)./deltay_SP(i);
    end
end

for j=1:nj-2
    for i=1:ni-2
        D_n(i,j)=gamma*A_n(j)./deltay_PN(i);
    end
end

%% Convection Terms
for i=1:ni-2
    for j=1:nj-2
        F_w(i,j)=rho*U_w(i,j).*A_w(i);
    end
end

for i=1:ni-2
    for j=1:nj-2
        F_e(i,j)=rho*U_e(i,j).*A_e(i);
    end
end

for j=1:nj-2
    for i=1:ni-2
        F_s(i,j)=rho*V_s(i,j).*A_s(j);
    end
end

for j=1:nj-2
    for i=1:ni-2
        F_n(i,j)=rho*V_n(i,j).*A_n(j);
    end
end
    
%% Coefficient Matrix E,W
for i = 1:(length(xc)-1)
    xn(i) = (xc(i+1)- xc(i))/2; % half cell length x-dir
end
for j = 1:(length(yc)-1)
    yn(j) = (yc(j+1)-yc(j))/2; % half cell length x-dir
end

fex= xn./deltax_PE;
fwx= xn./deltax_WP;
fny = yn./deltay_PN;
fsy = yn./deltay_SP;

for z= 1:nj-2
     for s=1:ni-2
aW_1= D_w(z,s)+fwx(s)*F_w(z,s);
aW_2 = D_w(z,s)+F_w(z,s) ;
aw(1) = aW_1 ;
aw(2) = aW_2 ;
aw(3) = D_w(z,s) ;
 aW(z,s) = max (aw);
     end
end

for z= 1:nj-2
     for s=1:ni-2
aE_1=  D_e(z,s) -fex(s)*F_e(z,s);
aE_2 = D_e(z,s) - F_e(z,s);
ae(1) = aE_1 ;
ae(2) = aE_2 ;
ae(3) = D_e(z,s) ;
 aE(z,s) = max (ae);
     end
end

%aW(1:nj-h_in-2,1) = 0; % West wall Neumann B.C.
% aE(1:h_out,nj-2)=0; % East Wall Dirichlet
if (bc == 0) 
%aE(1:nj-2-h_out,ni-2) = 0 ; % East wall Neumann B.C.
end

%% Coefficient Matrix N,S

for z= 1:nj-2
     for s=1:ni-2
aS_1= D_s(z,s)+fsy(s)*F_s(z,s);
aS_2 = D_s(z,s)+F_s(z,s);
as(1) = aS_1 ;
as(2) = aS_2 ;
as(3) = D_s(z,s) ;
 aS(z,s) = max (as);
     end
end

for z= 1:nj-2
     for s=1:ni-2
aN_1=  D_n(z,s)-fny(s)*F_n(z,s);
aN_2 = D_n(z,s) - F_n(z,s);
an(1) = aN_1 ;
an(2) = aN_2 ;
an(3) = D_n(z,s) ;
 aN(z,s) = max (an);
     end
end

aS(1,:) = 0 ; % South wall Neumann B.C.
aN(nj-2,:) = 0 ; % North wall Neumann B.C.

%% Node Coefficient Matrix
aP = aE+aW+aN+aS;

%% T-Matrix
T=zeros(nj,ni);
if (bc==1)
    for i=h_out+2:nj
    T(i,ni)=10; %East wall dirichlet B.C.
    end
end
for i=nj-h_in:nj
    T(i,1)=T_in; %Inlet temperature
end

T(2:nj-1,2:ni-1) = T_i; % Remaining nodes initial temp.

%% Source Matrix
B=zeros(nj-2,ni-2);
for s= 1:nj-2
    for z= 1:ni-2
 B(s,z)=  s_bar*4*xn(z)*yn(s);
    end
end

%% solver
h_A = yc(nj-1)-yc(nj-h_in-2); %Inlet length
iter=0;
cnv=0.1;
%Gauss-Seidel------------------------------------------------------------
if (solver==0) 
while (cnv>error)
    for s= 1:nj-2
        for z=1:ni-2
            T(s+1,z+1)=(aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+...
                aN(s,z)*T(s+2,z+1)+aS(s,z)*T(s,z+1)+B(s,z))/(aP(s,z));                        
        end
    end

T(nj,:)=T(nj-1,:); %B.C North wall
T(1,:)=T(2,:); %B.C South Wall

if (bc == 0)
    T(:,ni)=T(:,ni-1); %B.C. East wall Neumann
end

T(2:h_out+1,ni)=T(2:h_out+1,ni-1); %Outlet B.C. Neumann
T(1:nj-h_in-1,1)=T(1:nj-h_in-1,2); %B.C. West Wall Neumann (Except Inlet)

R=0;

for s= 1:nj-2
        for z=1:ni-2
         R=(abs((aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+...
                aN(s,z)*T(s+2,z+1)+aS(s,z)*T(s,z+1)+B(s,z))-(aP(s,z)*T(s+1,z+1))))+R;
        end
end

iter=iter+1
F = rho*U_w(nj-2,1)*h_A*(T_in-T(2,ni-1));
conv_hist(iter)=R/F;
cnv=conv_hist(iter)

 end
end
%TDMA---------------------------------------------------------------------
if (solver==1)
while (cnv>error)
    T(nj,:)=T(nj-1,:); %B.C North wall
    T(1,:)=T(2,:); %B.C South Wall

    if (bc == 0)
        T(:,ni)=T(:,ni-1); %B.C. East wall Neumann
    end

    T(2:h_out+1,ni)=T(2:h_out+1,ni-1); %Outlet B.C. Neumann
    T(1:nj-h_in-1,1)=T(1:nj-h_in-1,2); %B.C. West Wall Neumann (Except Inlet)

    d_j = zeros(ni-2,nj-2);
    P = zeros(ni-2,nj-2);
    Q = zeros(ni-2,nj-2);
for s= 1:nj-2
      for z=1:ni-2
        
a_j = aP;
b_j = aN;
c_j = aS;
d_j(s,z) = aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+(B(s,z));%(aN(s,z)*T(s,z))+(aS(s,z)*T(s,z+1))+(B(s,z));
P(s,1) = b_j(s,1)/a_j(s,1);
Q(s,1) = (d_j(s,1)+ (c_j(s,1)*T(s+1,z)))/a_j(s,1);
       end
end

for p = 1:ni-2
    for q = 2:ni-2
    P(p,q) = b_j(p,q)/(a_j(p,q)-(c_j(p,q)*P(p,q-1)));
    Q(p,q) = (d_j(p,q)+(c_j(p,q)*Q(p,q-1)))/(a_j(p,q)-(c_j(p,q)*P(p,q-1)));
    end
end

for n = 1:ni-2
   for m = nj-2:-1:1
    T(n+1,m+1) = (P(n,m)*T(n+1,m+2)) + (Q(n,m));
   end
end

R=0;
for s= 1:nj-2
        for z=1:ni-2
         R=(abs((aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+...
                aN(s,z)*T(s+2,z+1)+aS(s,z)*T(s,z+1)+B(s,z))-(aP(s,z)*T(s+1,z+1))))+R;
        end
end

iter=iter+1
F = rho*U_w(nj-2,1)*h_A*(T_in-T(2,ni-1));
conv_hist(iter)=R/F;
cnv=conv_hist(iter)

end
end

%% Face Temperatures
Tf_x = zeros(nj-2,ni-1) ; % Faces in x-dir
for s= 1:nj-2
        for z=1:ni-3
            Tf_x(s,z+1)= (T(s+1,z+2)-T(s+1,z+1))*(xc(z+1)-xp(z+1))/(xp(z+2)-xp(z+1))+T(s+1,z+1) ;
        end
end
Tf_x(:,1)= T(2:nj-1,1);
Tf_x(:,ni-1)= T(2:nj-1,ni);

Tf_y = zeros(nj-1,ni-2) ; % Faces in y-dir
for z= 1:ni-2
        for s=1:nj-3
            Tf_y(s+1,z)= (T(s+2,z+1)-T(s+1,z+1))*(yc(s+1)-yp(s+1))/(yp(s+2)-yp(s+1))+T(s+1,z+1) ;
        end
end
Tf_y(1,:)= T(1,2:ni-1);
Tf_y(nj-1,:)= T(nj,2:ni-1);

%% Flux Matrix
flux_x = zeros(nj-2,ni-2);
flux_y = zeros(nj-2,ni-2);
for s= 1:nj-2
for z=1:ni-2
         flux_x(s,z) = -gamma*(Tf_x(s,z+1)-Tf_x(s,z))./(xc(z+1)-xc(z));
         flux_y(z,s) = -gamma*(Tf_y(z+1,s)-Tf_y(z,s))./(yc(z+1)-yc(z));
end
end

%% Global Conservation
glb_cnv = sum( sum(flux_x(:,1))+sum(flux_x(:,ni-2))+sum(flux_x(1,:))+sum(flux_x(nj-2,:))+...
    sum(flux_y(:,1))+sum(flux_y(:,ni-2))+sum(flux_y(1,:))+sum(flux_y(nj-2,:)));

%% Plots
[xnc_m,ync_m] = meshgrid(xp(2:(ni-1)),yp(2:(nj-1)));
figure(1) % Temperature contour plot
hold on
T_node = T(2:(nj-1),2:(ni-1));
contourf(xnc_m,ync_m,T_node)
vec= 5;
quiver(xp,yp,u2d,v2d,vec)

figure(2) % Horizontal temperature flux
hold on
plot(xp(2:ni-1),flux_x(5,:),'r')
plot(xp(2:ni-1),flux_x(49,:),'g')
plot(xp(2:ni-1),flux_x(30,:),'b')
plot(xp(2:ni-1),flux_x(20,:),'m')
legend({'Row 5','Row 49','Row 30','Row 20'},'Location', 'NorthWest');
title('Temperature') ;
xlabel({'x direction [m] '});
ylabel({'Temperature flux [C]'});

toc