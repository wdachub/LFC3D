% Figure 3
clear
addpath '../../useCase';
addpath '../src';
t0=0;
te=[0.3,0.5,0.8,1];
nNodes=256;
nNodeTs=ceil(nNodes*te);
MovingSurface=@EPSurface1;
vel=@velTaylorGreenVortex2;
pfun=@pfunTriFunc;
order=[6,6,6,10];
for i = 1:size(te,2)
tic
[flux,pts, SDR]= fluxDR3D(MovingSurface, vel, pfun,t0, te(i),[nNodes,nNodeTs(i)],order);
save(['../Figures/Data/MTaylorGreenVortex2pts' int2str(i) '.mat'],'pts')
save(['../Figures/Data/MTaylorGreenVortex2SDR' int2str(i) '.mat'],'SDR')
toc
end
