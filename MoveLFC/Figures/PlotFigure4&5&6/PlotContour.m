% Plot the donating region.
function  PlotContour(u,v,tp,SDR,alphaVal,fontSize)
[uu,vv]=meshgrid(u,v);
[~,ttp]=meshgrid(u,tp);

%plotDRcomponentcontour(uu,vv,SDR.surface,alphaVal)
%hold on
plotDRcomponentcontour(uu,vv,SDR.timeline,alphaVal)
hold on
for i=1:4
  plotDRcomponentcontour(uu,ttp,SDR.streakline{i},alphaVal)  
end
hold off

ylabel('$y$','Interpreter','Latex','FontSize',fontSize)
xlabel('$x$','Interpreter','Latex','FontSize',fontSize)
zlabel('$z$','Interpreter','Latex','FontSize',fontSize,'rotation',0)
ax=gca;
ax.ZLabel.HorizontalAlignment='right';
axis tight
set(gca,'FontSize',fontSize)
grid on