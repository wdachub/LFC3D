function productInt =splinegauss(szSp,q,xi,t0,pfun)
% Calculate the integral $\int_{S_{i}} F(x,y,z)dy\wedge dz$ on page 13,
% (29) by Gauss quadrature formula. szSp is a spline-approximated surface
% and we should calculate the numerical surface integral on it. pfun
% represents function $f$ in our paper, while $\frac{\partial F}{\partial
% x}=pfun.
%
% List of Parameters:
% =========================================================================
%   IN
%       Name  |         Description
%-------------|------------------------------------------------------------
%      szSp   | Defines the given surface $S_{i} which we need to derive 
%             | the integral formula on. It is a struct constructed with:
%             | [szSp.x,szSp.y,szSp.z],each of its components is a
%             | piecewise two-variate polynomial. Actually, these three
%             | components compose the spline-approximated surface.
%      -------|------------------------------------------------------------
%      q      | The quadrature rule of Gauss quadrature rule, consists of
%             | integral nodes and weights. This variable is determined by
%             | the convergence order of Gauss quadrature.
%      -------|------------------------------------------------------------
%      xi     | Auxiliary parameter. In our paper, we need to calculate the
%             | antiderivate of $f$ related to $x$, and we should make the 
%             | interval $[\xi,x]$ changes to $[-1,1]$ by affine map. In
%             | practise, for coordinate x varies from a to b, xi is
%             | defined to be (a+b)/2.
%      -------|------------------------------------------------------------
%      t0     | The time parameter of the given scalar function pfun. In
%             | practise, pfun is the initial value of a time-related
%             | scalar function, which means that t0 is chosen to be a 
%             | constant, which is equal to the initial time.
%      -------|------------------------------------------------------------
%      pfun   | The function $f_{0}(x)$ in page 1 of our paper. It is 
%             | stored as a function handle.  
% -------------------------------------------------------------------------
%
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
%   productInt|  The approximated value of the integral 
%             |  $\int_{S_{i}} F(x,y,z)dy\wedge dz$

% Derivatives of y(u,v) and z(u,v)
szSp.yu=fnder(szSp.y,[1,0]); 
szSp.yv=fnder(szSp.y,[0,1]);
szSp.zu=fnder(szSp.z,[1,0]);
szSp.zv=fnder(szSp.z,[0,1]);

% The quadrature rule in each square.

u=szSp.x.breaks{1};
v=szSp.x.breaks{2};


lpu=szSp.x.pieces(1);%the number of block in u direction
lpv=szSp.x.pieces(2);%the number of block in v direction
diffu=diff(u)/2;% difference in u nodes, (u(i+1)-u(i))/2
diffv=diff(v)/2;% difference in v nodes.
addu=(u(1:end-1)+u(2:end))/2; % (u(i+1)+u(i))/2
addv=(v(1:end-1)+v(2:end))/2; % (v(i+1)+v(i))/2
productInt = 0;

% Derive the surface integral block by block.
for i=1:lpu
    for j=1:lpv
        temp_nodesu=diffu(i)*q.yn+addu(i);
        temp_nodesv=diffv(j)*q.yn+addv(j);
        temp_nodesx=ppualmy(szSp.x,[temp_nodesu';temp_nodesv'],[i,j]);
        temp_nodesy=ppualmy(szSp.y,[temp_nodesu';temp_nodesv'],[i,j]);
        temp_nodesz=ppualmy(szSp.z,[temp_nodesu';temp_nodesv'],[i,j]);
        temp_diff_x=(temp_nodesx-xi)/2;
        temp_nodesxl=bsxfun(@times,repmat(temp_diff_x,1,1,length(q.xn)),shiftdim(q.xn,-2))+...
            repmat((temp_nodesx+xi)/2,1,1,length(q.xn));
        jacob=ppualmy(szSp.yu,[temp_nodesu';temp_nodesv'],[i,j]).*...
            ppualmy(szSp.zv,[temp_nodesu';temp_nodesv'],[i,j])-...
            ppualmy(szSp.yv,[temp_nodesu';temp_nodesv'],[i,j]).*...
            ppualmy(szSp.zu,[temp_nodesu';temp_nodesv'],[i,j]);
        temp_weight=(q.yw*q.yw').*(diffu(i)*diffv(j)').*jacob.*temp_diff_x;
        
        
        temp_weight=bsxfun(@times,repmat(temp_weight,1,1,length(q.xw)),shiftdim(q.xw,-2));
        temp_nodesxl=temp_nodesxl(:);
        temp_nodesy=temp_nodesy(:);
        temp_nodesy=repmat(temp_nodesy,length(q.xw),1);
        temp_nodesz=temp_nodesz(:);
        temp_nodesz=repmat(temp_nodesz,length(q.xw),1);
        temp_weight=temp_weight(:);
        
        %delete zero weights
        zerow=temp_weight==0;
        temp_nodesxl(zerow)=[];
        temp_nodesy(zerow)=[];
        temp_nodesz(zerow)=[];
        temp_weight(zerow)=[];
        %record the nodes.
        if isempty(temp_weight)
            productInt = productInt + 0;
        else
            fNodes=pfun(temp_nodesxl,temp_nodesy,temp_nodesz,t0);
            productInt = productInt + temp_weight'*fNodes;
        end
    end
end
productInt=-productInt;%adjust the orientation
% ppual function, find the antiderivative of $f$ by 1-dimension quadrature
% rule.
function v = ppualmy(pp,x,ix)

lx=length(x(1,:));%number of points
% locate the scattered data in the break sequences:
[mx,~] = size(x);

m = 2;%two dimension problem
if mx~=m, error(message('SPLINES:PPUAL:wrongx', num2str( m ))), end

localcoef=pp.coefs(1,ix(1):pp.pieces(1):ix(1)+(pp.order(1)-1)*pp.pieces(1),...
    ix(2):pp.pieces(2):ix(2)+(pp.order(2)-1)*pp.pieces(2));%choose the local coef.
localcoef=shiftdim(localcoef,1);%1*n*n to n*n

ry0=repmat(x(2,:),pp.order(1),1)-pp.breaks{2}(ix(2));
rx0=x(1,:)-pp.breaks{1}(ix(1));
temp=repmat(localcoef(:,1),1,lx);
for j=2:pp.order(2)
    temp =bsxfun(@plus,temp.*ry0,localcoef(:,j));
end

v=zeros(lx,lx);
for i=1:lx%loop of x
    rxi0=repmat(rx0(i),1,lx);
    tempi=temp;
    for j=2:pp.order(1)
        tempi(1,:) = tempi(1,:).*rxi0+tempi(j,:);
    end
    v(i,:)=tempi(1,:);
end









