function timeline = flowmap(p, t0, te, vel, RKOrder, dt)
% Calculate the numerical flow map of point p with velocity field vel, 
% during time interval [t0,te]
% Marked as $\phi_{t_{0}}(p,t_{e})$ in our paper. See algorithm 1 in 
% our paper for detail.
% The start position is a single point $p$.
% The flow map marks the trajectory of a point in velocity field.
% 
% List of Parameters:
% =========================================================================
%   IN
%       Name  |         Description
%-------------|------------------------------------------------------------
%      p      | Defines the initial point, stored as:
%             |        [x1 y1 z1]
%             | This vector marks a point in space $\mathbb{R}^{3}$.
%      -------|------------------------------------------------------------
%      t0     | The initial time parameter in the flow map.
%      -------|------------------------------------------------------------
%      te     | The final time parameter in the flow map.
%      -------|------------------------------------------------------------
%      vel    | The given velocity field. This parameter is stored as a
%             | function handle. Its form must be vel(t,p), t is a real
%             | number marks time parameter, p marks a point in space 
%             | $\mathbb{R}^{3}$, and vel(t,p) returns a vector in 
%             | $\mathbb{R}^{3}$, marks the velocity vector on point p at 
%             | time t.  
%      -------|------------------------------------------------------------
%      RKorder| The convergence order of the numerical ode algorithm. In
%             | our paper, we choose explicit Runge-Kutta algorithm, so
%             | this parameter marks the convergence order of RK method. 
%             | This parameter must be an integer ranges from 2 to 8.
%      -------|------------------------------------------------------------
%      dt     | The time step of the ode solver.
% -------------------------------------------------------------------------
%
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
%     timeline| Marks the position of point p at time te. Stored as:
%             |        [x1 y1 z1]
%             | It's actually the image of point p on discrete flow map.

% Partition the time interval [t0,te]
ttmp=t0:dt:te;
if ttmp(end)~=te
    ttmp=[ttmp,te];
end

% Calculate the discretized flow map by Runge-Kutta method.
ptmp=p;
for i=1:length(ttmp)-1
    ptmp=RungeKutta(ptmp, ttmp(i), ttmp(i+1), vel, RKOrder);
end
timeline=ptmp;
end