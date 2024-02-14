% Plot the figure of donating region in our paper.
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
light('Position',[1 3 2]);  % light position
light('Position',[-3 -1 3]);
material shiny;
ylabel('$y$','Interpreter','Latex','FontSize',fontSize)
xlabel('$x$','Interpreter','Latex','FontSize',fontSize)
zlabel('$z$','Interpreter','Latex','FontSize',fontSize,'rotation',0)
ax=gca;
ax.ZLabel.HorizontalAlignment='right';
axis equal 
axis tight
set(gca,'FontSize',fontSize)
grid on