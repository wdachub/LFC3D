 function exactflux = computeFlux(p1, p2,p3, t0, te, vel, pfun, surface, T)
% Return the exact flux through line segment p0-p1 with time [t0,te]
%  for velocity field vortexShear2.
% The divergence-free condition of velocity field is required.
% The outward normal is the RHS of ray p0->p1;
%  which agrees to the convention of counter-clockwise polygon.
% Hence positive flux implies out-flow.
% integrand = 1 : the flux is \int u\cdot n
%           = 2 : the flux is \int f*u \cdot n where
%                 f = sin(2\pi x)sin(2*\pi y)

if exist('T','var')==0
    T = 3;
end


x1 = p1(1); y1 = p1(2);z1 = p1(3);
x2 = p2(1); y2 = p2(2);z2 = p2(3);
x3 = p3(1); y3 = p3(2);z3 = p3(3);

    if isequal(vel, @velTaylorGreenVortex) && isequal(pfun,@pfunTriFunc)&& isequal(surface,@FixedSurface1)
        exactflux=-sin(pi*z1)^3*(sin(pi*x2)^3-sin(pi*x1)^3)*(sin(pi*y3)^3-sin(pi*y1)^3)*(sin(pi*te/T)-sin(pi*t0/T))*4*T/(9*pi^3);
    elseif isequal(vel, @velTaylorGreenVortex) && isequal(pfun,@pfunConstant) && isequal(surface,@EPSurface1)...
            && isequal(t0,0) && isequal(te,1)
        exactflux=-0.09756854721114007620521726;
    elseif isequal(vel, @velTaylorGreenVortex) && isequal(pfun,@pfunTriFunc) && isequal(surface,@EPSurface1)...
            && isequal(t0,0) && isequal(te,1)
        exactflux=-0.00001097648149308549609143151;
        
    elseif isequal(vel, @velIntersection3D) && isequal(pfun,@pfunIntersection3D) && isequal(surface,@EPSurface2)...
            && isequal(t0,0) && isequal(te,1)
        exactflux=-0.001783095654530444584141087;
    elseif isequal(vel, @velIntersection3D) && isequal(pfun,@pfunIntersection3D) && isequal(surface,@EPSurface2)...
            && isequal(t0,0) && isequal(te,2)
        exactflux=-0.002464276943522421041523338;
    elseif isequal(vel, @velIntersection3D) && isequal(pfun,@pfunIntersection3D) && isequal(surface,@FixedSurface2)...
            && isequal(t0,0) && isequal(te,1)
        exactflux=7*(1-exp(-5))/48;
    end
    exactflux=-exactflux;%adjust the orientation
 end


function y=fluxnu(t,s,p0,p1,velCase,funt)
r=p1-p0;
n=[r(2),-r(1)];
y=zeros(size(t));
for i=1:size(t,1)
    for j=1:size(t,2)
        y(i,j)=funt(t(i,j),p0+s(i,j)*r)*n*velCase(t(i,j),p0+s(i,j)*r);
    end
end
end
