% Figure 5
clear
addpath '../../useCase';
addpath '../src';
tic
t0=0;
te=[1,2];
nNodes=256;
nNodeTs=ceil(nNodes*te);
MovingSurface=@EPSurface2;
vel=@velIntersection3D;
pfun=@pfunIntersection3D;

order=[6,6,6,10];
for i = 1:size(te,2)
tic
[~,pts, SDR]= fluxDR3D(MovingSurface, vel, pfun, t0, te(i),[nNodes,nNodeTs(i)],order);
save(['../Figures/Data/MovingIntersection3Dpts' int2str(i) '.mat'],'pts')
save(['../Figures/Data/MovingIntersection3DSDR' int2str(i) '.mat'],'SDR')
toc
end
