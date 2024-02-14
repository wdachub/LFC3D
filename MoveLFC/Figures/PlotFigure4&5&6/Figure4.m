%Generate figure 4 in the paper.
clear
close all
addpath '../../../useCase';

figs = figure('Units','centimeters','Position',[10 6 16.5 7.8]);

leftMargin=0.8;
BottomMargin=0.7;
width=3;
hight=3;
wDistant=0.3;
clim=[0,0.25];
alphaVal = 0.8;
fontSize = 8;

%Line 1
nNodes = 50;
MovingSurface=@EPSurface1;
pts.u=linspace(0,1,nNodes);
pts.v=linspace(0,1,nNodes);
[uup,vvp]=meshgrid(pts.u,pts.v);
uup=uup';
vvp=vvp';
t = [0.3, 0.5, 0.8, 1];
for i = 1:size(t,2)
ttp = t(i)*ones(size(uup));
sur = MovingSurface(uup,vvp,ttp);
pts.surface(:,:,1) = sur(:,1:nNodes);
pts.surface(:,:,2) = sur(:,nNodes+1:2*nNodes);
pts.surface(:,:,3) = sur(:,2*nNodes+1:3*nNodes);
ax=axes('Units','centimeters','Position',[leftMargin*i+width*(i-1),BottomMargin+3.5,width,hight]);
surf(pts.surface(:,:,1),pts.surface(:,:,2),pts.surface(:,:,3),'FaceAlpha',alphaVal)
ylabel('$y$','Interpreter','Latex','FontSize',fontSize)
xlabel('$x$','Interpreter','Latex','FontSize',fontSize)
zlabel('$z$','Interpreter','Latex','FontSize',fontSize,'rotation',0)
ax.ZLabel.HorizontalAlignment='right';
set(gca,'CLim',clim)
shading flat;  % style
light('Position',[4 4 2]);  % position of light
light('Position',[-4 -4 3]);
material shiny;
title(['$$t_e = $$' num2str(t(i))],'interpreter','latex')
xlim([-1 1])
ylim([-1 1])
zlim([0 0.25])
set(gca,'FontSize',fontSize)

%Line 2
ax=axes('Units','centimeters','Position',[leftMargin*i+width*(i-1),BottomMargin,width,hight]);
A = load(['../Data/MTaylorGreenVortex2pts' num2str(i) '.mat']);
B = load(['../Data/MTaylorGreenVortex2SDR' num2str(i) '.mat']);
PaperPlotDR(A.pts.u,A.pts.v,A.pts.tp,B.SDR,alphaVal,fontSize)
set(gca,'CLim',clim)
view(46.5739,28.2421)
xlim([-1 1])
ylim([-1 1])
zlim([0 0.25])
end

hold on
axes('Units','centimeters','Position',[15.5,BottomMargin+0.2,0.05,2*hight]);
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

exportgraphics(gcf,'./Figure4.pdf','Resolution',300)