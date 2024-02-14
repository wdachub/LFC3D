%See the contour of figure 5.
clear
close all
figs = figure('Units','centimeters','Position',[10 6 7.5 5.5]);

leftMargin=1;
BottomMargin=0.5;
width=4.5;
hight=3.9;
wDistant=0.3;
clim=[-0.55,0.15];
alphaVal = 0.8;
fontSize = 9;

%Figure 5 contour
ax1=axes('Units','centimeters','Position',[leftMargin+0.1,BottomMargin+0.4,width,hight]);
A = load('../Data/Intersection3DwithFixedpts.mat');
B = load('../Data/Intersection3DwithFixedSDR.mat');
PlotContour(A.pts.u,A.pts.v,A.pts.tp,B.SDR,alphaVal,fontSize)
xlim([-1,1])
ylim([-1,1])
exportgraphics(gcf,'./Contour.pdf','Resolution',400)