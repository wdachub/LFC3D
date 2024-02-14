%Generate figure 2.
clear
close all
figs = figure('Units','centimeters','Position',[10 5 7 4]);

leftMargin=0.8;
BottomMargin=0;
width=4.5;
hight=3.9;
wDistant=0.3;
clim=[-0.55,0.15];
alphaVal = 0.8;
fontSize = 9;

%Figure 5
ax1=axes('Units','centimeters','Position',[leftMargin+0.1,BottomMargin+0.4,width,hight]);
A = load('../Data/TaylorGreenVortex2ptsSmall.mat');
B = load('../Data/TaylorGreenVortex2SDRSmall.mat');
PaperPlotDR(A.pts.u,A.pts.v,A.pts.tp,B.SDR,alphaVal,fontSize)
view(42.8770,24.9836)
hold on
axes('Units','centimeters','Position',[leftMargin+width+0.3,BottomMargin+0.4,0.05,hight-0.4]);
c=colorbar('Location', 'eastoutside');
c.FontSize=fontSize;
caxis(clim);
% Manually flush the event queue and force MATLAB to render the colorbar
% necessary on some versions
drawnow
% Get the color data of the object that correponds to the colorbar
cdata = c.Face.Texture.CData;
% Change the 4th channel (alpha channel) to 10% of it's initial value (255)
cdata(end,:) = uint8(alphaVal * cdata(end,:));
% Ensure that the display respects the alpha channel
c.Face.Texture.ColorType = 'truecoloralpha';
% Update the color data with the new transparency information
c.Face.Texture.CData = cdata;
% Make sure that the renderer doesn't revert your changes
drawnow
c.Face.ColorBinding = 'discrete';
set(gca,'visible','off')
%PlotContour(A.pts.u,A.pts.v,A.pts.tp,B.SDR,alphaVal,fontSize);
exportgraphics(gcf,'./Figure2.pdf','Resolution',300)