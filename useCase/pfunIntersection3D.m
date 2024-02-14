function pfun = pfunIntersection3D(x,y,z,t)
%passive function conserved by time dependent strain flow.
pfun=exp(-5*t).*(x.^2+y.^2+z.^2);
end   