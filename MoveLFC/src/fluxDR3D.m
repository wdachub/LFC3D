function [flux,pts,SDR] = fluxDR3D(MovingSurface, vel, pfun,t0,te ,nSeg,order)
%Generate the donating region and calculate the flux via LFC algorithm.
% order: the convergence order of each steps.
% List of Parameters:
% =========================================================================
%   IN
%       Name  |         Description
%-------------|------------------------------------------------------------
%MovingSurface| The given moving surface in LFC problem. Stored in a
%             | function handle, point=MovingSurface(u,v,t).
%      -------|------------------------------------------------------------
%      vel    | The velocity field. Given by a function handle, 
%             | velocity=vel(x,y,z,t).
%      -------|------------------------------------------------------------
%      pfun   | The scalar function as the initial value for scalar 
%             | conservation law. Given by a function handle.
%      -------|------------------------------------------------------------
%      t0     | Initial time parameter.
%      -------|------------------------------------------------------------
%      te     | Final time parameter.
%      -------|------------------------------------------------------------
%      nSeg   | The number of discretized intervals.
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
%       flux  | The lagrangian flux derived by initial value pfun and
%             | donating region pts&SDR.

%%
% Discretize the space and time interval.
nNodes=nSeg(1);%the number of spacial nodes.
dt=(te-t0)/(nSeg(2)-1);
DR_SplineOrder=order(2);
RKOrder = order(3);
pts.u=linspace(0,1,nNodes);
pts.v=linspace(0,1,nNodes);
[uup,vvp]=meshgrid(pts.u,pts.v);
uup=uup';
vvp=vvp';

% Get the sample points on the moving surface.
sur = MovingSurface(uup,vvp,t0);
pts.surface(:,:,1) = sur(:,1:nNodes);
pts.surface(:,:,2) = sur(:,nNodes+1:2*nNodes);
pts.surface(:,:,3) = sur(:,2*nNodes+1:3*nNodes);

% Generate a spline-approximated surface via piecewise polynomial.
SDR.surface.x=spapi({DR_SplineOrder,DR_SplineOrder},{pts.u,pts.v},pts.surface(:,:,1));
SDR.surface.y=spapi({DR_SplineOrder,DR_SplineOrder},{pts.u,pts.v},pts.surface(:,:,2));
SDR.surface.z=spapi({DR_SplineOrder,DR_SplineOrder},{pts.u,pts.v},pts.surface(:,:,3));
    SDR.surface.x=sp2pp(SDR.surface.x);%convert to pp form 
    SDR.surface.y=sp2pp(SDR.surface.y);
    SDR.surface.z=sp2pp(SDR.surface.z);

%%
%Get the preimage of the given surface on flow map.
sur = MovingSurface(uup,vvp,te);
pts.surface(:,:,1) = sur(:,1:nNodes);
pts.surface(:,:,2) = sur(:,nNodes+1:2*nNodes);
pts.surface(:,:,3) = sur(:,2*nNodes+1:3*nNodes);

for i=1:length(pts.u)
    for j=1:length(pts.v)
        pts.timeline(i,j,:)=flowmap(squeeze(pts.surface(i,j,:))', te, t0, vel, RKOrder, -dt);%timeline points
    end
end

% Spline-approximated timeline.
SDR.timeline.x=spapi({DR_SplineOrder,DR_SplineOrder},{pts.u,pts.v},pts.timeline(:,:,1));
SDR.timeline.y=spapi({DR_SplineOrder,DR_SplineOrder},{pts.u,pts.v},pts.timeline(:,:,2));
SDR.timeline.z=spapi({DR_SplineOrder,DR_SplineOrder},{pts.u,pts.v},pts.timeline(:,:,3));
    SDR.timeline.x=sp2pp(SDR.timeline.x);%convert to pp form 
    SDR.timeline.y=sp2pp(SDR.timeline.y);
    SDR.timeline.z=sp2pp(SDR.timeline.z);

%%
% Generate the surface derived by streaklines. These streaklines start on
% the boundary of given moving surface.
%
% Notation:
% The boundary is a closed curve with anticlockwise orientation.
% (xxp==0|xxp==1|yyp==0|yyp==1)
% 1:nSeg x=0
% nSeg:nSeg:nSeg*(nSeg-1) y=1
% nSeg^2:-1:nSeg*(nSeg-1)+1 x=1
% nSeg*(nSeg-2)+1:-nSeg:nSeg y=0
% (0,0)->(0,1)->(1,1)->(1,0)->(0,0)
% In order to save computation of strealine, we store
% (0,0)->(0,1)->(1,1)->(1,0)
% streakline=zeros(4*nSeg-3,1,3);%assume that there is same grid size on u and v.

sur = MovingSurface(uup,vvp,t0);
pts.surface(:,:,1) = sur(:,1:nNodes);
pts.surface(:,:,2) = sur(:,nNodes+1:2*nNodes);
pts.surface(:,:,3) = sur(:,2*nNodes+1:3*nNodes);
pstreakline(1:nNodes,1,:)=pts.surface(1,1:nNodes,:);%clockwise of matrix
pstreakline(nNodes+1:2*nNodes-1,1,:)=pts.surface(2:nNodes,nNodes,:);
pstreakline(2*nNodes:3*nNodes-2,1,:)=pts.surface(nNodes,nNodes-1:-1:1,:);
pstreakline(3*nNodes-1:4*nNodes-4,1,:)=pts.surface(nNodes-1:-1:2,1,:);
pts.tp=linspace(t0,te,nSeg(2));
for indbbtp=2:length(pts.tp)
    sur = MovingSurface(uup,vvp,pts.tp(indbbtp));
    pts.surface(:,:,1) = sur(:,1:nNodes);
    pts.surface(:,:,2) = sur(:,nNodes+1:2*nNodes);
    pts.surface(:,:,3) = sur(:,2*nNodes+1:3*nNodes);
    pstreakline(1:nNodes,indbbtp,:)=pts.surface(1,1:nNodes,:);%clockwise of matrix
    pstreakline(nNodes+1:2*nNodes-1,indbbtp,:)=pts.surface(2:nNodes,nNodes,:);
    pstreakline(2*nNodes:3*nNodes-2,indbbtp,:)=pts.surface(nNodes,nNodes-1:-1:1,:);
    pstreakline(3*nNodes-1:4*nNodes-4,indbbtp,:)=pts.surface(nNodes-1:-1:2,1,:);
    pstreakline(:,indbbtp,:)=flowmap(squeeze(pstreakline(:,indbbtp,:)), pts.tp(indbbtp), t0, vel, RKOrder, -dt);
end
pstreakline(4*nNodes-3,:,:)=pstreakline(1,:,:);


pts.streakline{1}=pstreakline(1:nNodes,:,:);
pts.streakline{2}=pstreakline(nNodes:2*nNodes-1,:,:);
pts.streakline{3}=pstreakline(2*nNodes-1:3*nNodes-2,:,:);
pts.streakline{4}=pstreakline(3*nNodes-2:4*nNodes-3,:,:);




for i=1:4
    SDR.streakline{i}.x=spapi({DR_SplineOrder,DR_SplineOrder},{pts.u,pts.tp},pts.streakline{i}(:,:,1));
    SDR.streakline{i}.y=spapi({DR_SplineOrder,DR_SplineOrder},{pts.u,pts.tp},pts.streakline{i}(:,:,2));
    SDR.streakline{i}.z=spapi({DR_SplineOrder,DR_SplineOrder},{pts.u,pts.tp},pts.streakline{i}(:,:,3));
    SDR.streakline{i}.x=sp2pp(SDR.streakline{i}.x);%convert to pp form 
    SDR.streakline{i}.y=sp2pp(SDR.streakline{i}.y);
    SDR.streakline{i}.z=sp2pp(SDR.streakline{i}.z);
end

% Derive the flux via donating region.
flux=DRIntegral(t0,pfun,pts,SDR,order);



end

