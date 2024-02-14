%Figure 2(b)
clear
addpath '../../useCase';
addpath '../src';
t0=0;
te=3/2;
nNodes=2^8+1;
nNodeTs=2^8+1;
MovingSurface=@FixedSurface1;
vel=@velTaylorGreenVortex2;
pfun=@pfunTriFunc;
order=[6,6,6,10];
tic
[~,pts, SDR]= fluxDR3D(MovingSurface, vel, pfun,t0, te,[nNodes,nNodeTs],order);
save(['../Figures/Data/TaylorGreenVortex2ptsSmall.mat'],'pts');
save(['../Figures/Data/TaylorGreenVortex2SDRSmall.mat'],'SDR');
toc

