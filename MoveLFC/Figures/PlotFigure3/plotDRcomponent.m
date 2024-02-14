%plot one of the components of the donating region.
function  plotDRcomponent(uu,vv,szSp,alphaVal)

u=szSp.x.breaks{1};
v=szSp.x.breaks{2};
X=ppual(szSp.x,{u,v});
Y=ppual(szSp.y,{u,v});
Z=ppual(szSp.z,{u,v});
surface(X,Y,Z,'FaceAlpha',alphaVal)
return
end

