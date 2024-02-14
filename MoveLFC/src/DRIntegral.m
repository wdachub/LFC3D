function flux = DRIntegral(t0,pfun,pts,SDR,order)
% Calculate the integral $\oint{S} F(x,y,z)dy\wedge dz$ on page 13,
% (29) by Gauss quadrature formula, while $S$ is the donating region, $F$
% is anti-derivative of given function pfun relative to coordinate 
% parameter x.
%
% List of Parameters:
% =========================================================================
%   IN
%       Name  |         Description
%-------------|------------------------------------------------------------
%      t0     | The time parameter of the given scalar function pfun. In
%             | practise, pfun is the initial value of a time-related
%             | scalar function, which means that t0 is chosen to be a 
%             | constant, which is equal to the initial time.
%      -------|------------------------------------------------------------
%      pfun   | The function $f_{0}(x)$ in page 1 of our paper. It is 
%             | stored as a function handle.  
%      -------|------------------------------------------------------------
%      pts    | The given points on the donating region with parameter u
%             | and v. The points are stored in a structure as:
%             | [u,v,surface,timeline,tp,streakline]. 
%             | u,v: the value of argument u and v. In fact, pts are stored
%             | in a matrix P, P_{ij}=x(u(i),v(j)).
%             | surface: The points on the moving surface.
%             | timeline: The points on the preimage of flow map.
%             | tp: time parameters.
%             | streakline: The four boundary surface consisted of
%             | streaklines derived by the boundary of the given surface.
%      -------|------------------------------------------------------------
%      SDR    | The spline-approximated donating region stored in a
%             | structure.
%             | SDR.surface:the spline-approximated initial surface.
%             | SDR.timeline:the spline-approximated timeline.
%             | SDR.streakline:the spline-approximated streakline
%      -------|------------------------------------------------------------
%      order  | The convergence order for each steps.
%             | surface_SplineOrder=order(1);
%             | DR_SplineOrder=order(2);
%             | RKOrder = order(3);
%             | cubatureOrder=order(4);
% -------------------------------------------------------------------------
%
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
%       flux  | The lagrangian flux derived by initial value pfun and
%             | donating region pts&SDR.

DR_SplineOrder=order(2);
cubatureOrder=order(4);
cubature_type = 4;%GAUSS-LEGENDRE (TENSORIAL).
% Find a interval [a,b] such that $\min x=a, \max x=b$.
% Derive a
temp(1)=min(min(pts.surface(:,:,1)));
temp(2)=min(min(pts.timeline(:,:,1)));
for i=1:4
    temp(2+i)=min(min(pts.streakline{i}(:,:,1)));
end
xi1=min(temp);%the min of the x component of all input points
clear temp
% Derive b
temp(1)=max(max(pts.surface(:,:,1)));
temp(2)=max(max(pts.timeline(:,:,1)));
for i=1:4
    temp(2+i)=max(max(pts.streakline{i}(:,:,1)));
end
xi2=max(temp);%the max of the x component of all input points
% Get the value of auxiliary parameter xi.
xi = (xi1+xi2)/2;

clear temp
SplineDegree=DR_SplineOrder-1;
q=cubature_manager(cubatureOrder,SplineDegree,cubature_type);%generate quadrature rule.

%Quadrature on the given surface
productInt(1)=splinegauss(SDR.surface,q,xi,t0,pfun);

%Quadrature on the timeline
productInt(2)=splinegauss(SDR.timeline,q,xi,t0,pfun);

%Quadrature on the surface consists of streakline
for i=1:4
    productInt(2+i)= splinegauss(SDR.streakline{i},q,xi,t0,pfun);
end

%The flux value.
flux=productInt(1)-productInt(2)+productInt(3)+productInt(4)+productInt(5)+productInt(6);
return
