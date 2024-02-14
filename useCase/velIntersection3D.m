function vel = velIntersection3D(t, xyz,T)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
% this velocity field produces winding numbers +1,0,-1,-2
a=1;
% b=3*(2+time);
b=2*pi;
vel = [a*x+b*(y+z)  -b*(x+z)+a*y a*z+b*(y-x)];

