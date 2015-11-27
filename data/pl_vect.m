%close all
clear all
% read xc
load ('grid2/coarse_grid/xc.dat')
load ('grid2/coarse_grid/yc.dat')
nim1=length(xc);
% nim1 = ni-1 = number of grid lines. Number of cell nodes = ni
ni=nim1+1;
%
% read yc

njm1=length(yc);
nj=njm1+1;
%
% read u
load ('grid2/coarse_grid/u.dat')
load ('grid2/coarse_grid/v.dat')
u2d=reshape(u,ni,nj);

% read v

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
%
% plot the velocity field
% length of vectors = vec
vec= 5;
quiver(xp,yp,u2d,v2d,vec)
axis('equal')
%print vectxy.ps -deps
