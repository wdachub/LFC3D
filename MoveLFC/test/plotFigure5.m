%Figure 5
clear
addpath '../../useCase';
addpath '../src';
t0=0;
te=1;
nNodes=2^8+1;
nNodeTs=2^8+1;
MovingSurface=@FixedSurface2;
vel=@velIntersection3D;
pfun=@pfunIntersection3D;
order=[6,6,6,10];
tic
[flux,pts, SDR]= fluxDR3D(MovingSurface, vel, pfun, t0, te,[nNodes,nNodeTs],order);
save(['../Figures/Data/Intersection3DwithFixedpts.mat'],'pts');
save(['../Figures/Data/Intersection3DwithFixedSDR.mat'],'SDR');
flux
toc