% Generate the up part on table 4.
clear
addpath '../../useCase';
addpath '../src';
t0=0;
te=1;
MovingSurface=@EPSurface2;
vel=@velIntersection3D;
pfun=@pfunIntersection3D;

nSegS0=32;
nSegT0=32;
nNodeSs=nSegS0*2.^[0,1,2,3,4];
nNodeTs=nSegT0*2.^[0,1,2,3,4];
orders=[2,4,6];
p1=[0;0;0];
p2=[0;1;0];
p3=[1;1;0];
exactflux = computeFlux(p1, p2, p3,t0, te, vel, pfun, MovingSurface);
for j=1:length(orders)
for i=1:length(nNodeSs)
    tic
nNodeS=nNodeSs(i);
nNodeT=nNodeTs(i);
order=[orders(j),orders(j),orders(j),10];
%     surface_SplineOrder=order(1);
%     DR_SplineOrder=order(2);
%     RKOrder = order(3);
%     cubatureOrder=order(4);
[flux,pts, SDR]= fluxDR3D(MovingSurface, vel, pfun, t0, te,[nNodeS,nNodeT],order);
toc

% plotDR(pts.u,pts.v,pts.tp,SDR)
numflux(j,i)=flux;
error(j,i)=abs(flux-exactflux)/abs(exactflux);
clear pts 

end
end

for j=1:length(orders)
result(j,1)=error(j,1);
for i=1:length(error)-1
   result(j,2*i+1)=error(j,i+1);
   result(j,2*i)=log(error(j,i)/error(j,i+1))/log(2);
end
end
result

outputDir = '../Tables/';
writeTestResultsToLatexTable(result, nSegS0, orders, 'Table4_1.txt', outputDir);
