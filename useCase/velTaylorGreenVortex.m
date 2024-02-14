function vel = velTaylorGreenVortex(t, xyz,T)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
T = 3;
vel = [2.*sin(pi.*x).^2.*sin(2.*pi.*y).*sin(2.*pi.*z).*cos(pi.*t./T) ...
       -sin(2.*pi.*x).*sin(pi.*y).^2.*sin(2.*pi.*z).*cos(pi.*t./T) ...
       -sin(2.*pi.*x).*sin(2.*pi.*y).*sin(pi.*z).^2.*cos(pi.*t./T)];