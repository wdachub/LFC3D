% Plot the donating region.
function  PaperPlotDR(u,v,tp,SDR,alphaVal,fontSize)
[uu,vv]=meshgrid(u,v);
[~,ttp]=meshgrid(u,tp);

plotDRcomponent(uu,vv,SDR.surface,alphaVal)
hold on
plotDRcomponent(uu,vv,SDR.timeline,alphaVal)
hold on
for i=1:4
  plotDRcomponent(uu,ttp,SDR.streakline{i},alphaVal)  
end
hold off

shading flat;  % style
light('Position',[3 3 2]);  % position of light
light('Position',[-3 -3 3]);
material shiny;
ylabel('$y$','Interpreter','Latex','FontSize',fontSize)
xlabel('$x$','Interpreter','Latex','FontSize',fontSize)
zlabel('$z$','Interpreter','Latex','FontSize',fontSize,'rotation',0)
ax=gca;
ax.ZLabel.HorizontalAlignment='right';
axis tight
set(gca,'FontSize',fontSize)
grid on