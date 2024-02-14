%Generate figure 6.
clear
close all
figs = figure('Units','centimeters','Position',[10 6 12.6 5.5]);

leftMargin=1;
BottomMargin=0.5;
width=4.5;
hight=3.9;
wDistant=0.3;
clim=[-0.55,0.15];
alphaVal = 0.8;
fontSize = 9;

%Left ones
ax1=axes('Units','centimeters','Position',[leftMargin+0.1,BottomMargin+0.4,width,hight]);
A = load('../Data/MovingIntersection3Dpts1.mat');
B = load('../Data/MovingIntersection3DSDR1.mat');
PaperPlotDR(A.pts.u,A.pts.v,A.pts.tp,B.SDR,alphaVal,fontSize)
title(['$$t\in[0,1]$$'],'interpreter','latex')
set(gca,'CLim',clim)
view(-135,45)

%Right ones
ax2=axes('Units','centimeters','Position',[leftMargin+width+leftMargin+0.4,BottomMargin+0.4,width,hight]);

A = load('../Data/MovingIntersection3Dpts2.mat');
B = load('../Data/MovingIntersection3DSDR2.mat');
%triPlot(A.pts.u,A.pts.v,A.pts.tp,B.SDR,'./obj/te_equal_2/')
PaperPlotDR(A.pts.u,A.pts.v,A.pts.tp,B.SDR,alphaVal,fontSize)
title(['$$t\in[0,2]$$'],'interpreter','latex')
view(-135,45)
set(gca,'CLim',clim)
hold on
axes('Units','centimeters','Position',[leftMargin+2*width+5*wDistant,BottomMargin+0.4,0.05,hight-0.4]);
c = colorbar;
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

exportgraphics(gcf,'./Figure6.pdf','Resolution',300)