function quadrature=cubature_manager(M,p_i,cubature_type)


%--------------------------------------------------------------------------1
% OBJECT:
%-----------
%
% THE TARGET OF THIS ROUTINE IS TO PROVIDE NODES AND WEIGHTS OF A
% QUADRATURE ROUTINE ON [-1,1]. THE CODES ARE DESCRIBED BY L.N. TREFETHEN
% IN HIS CLENSHAW-CURTIS PAPER AND BY WALDVOGEL (PUBLISHED BY BIT).
%
%--------------------------------------------------------------------------
% INPUTS:
%----------
%
% M : DEGREE OF PRECISION OF 2D CUBATURE RULE.
%
% p_i: SPLINE DEGREE
%
% cubature_type: IT POINTS A CUBATURE ON THE SQUARE [-1,1] x [-1,1].
%
%     [cubature_type=0]: PADUA POINTS.
%     [cubature_type=1]: FEJER 1 (TENSORIAL).
%     [cubature_type=2]: FEJER 2 (TENSORIAL).
%     [cubature_type=3]: CLENSHAW CURTIS (VIA WEIGHTS) (TENSORIAL).
%     [cubature_type=4]: GAUSS-LEGENDRE (TENSORIAL).
%     [cubature_type=5]: GAUSS-LEGENDRE-LOBATTO (TENSORIAL).
%
%----------
% OUTPUTS:
%----------
%
% nodes : M x 2 MATRIX OF NODES, IN THE INTERVAL [-1,1].
%
% weights: M x 1 COLUMN VECTOR OF WEIGHTS.
%
%--------------------------------------------------------------------------
% ADDITIONAL ROUTINES :
%----------------------
%
% 1. quadrature_rules_1D
% 2. pdcub
%
%--------------------------------------------------------------------------

switch cubature_type
    
    case 4
        
        
        
        % lower degree of a Gauss-Legendre rule with ADE equal to Nu_deg.
        Nm_deg=ceil(((M+p_i)*p_i)/2);
        Nm=Nm_deg-1;%When the input of the subroutine quadrature_rules_1D is n, the output is a quadrature rules with n+1 points.
        [quadrature.yn,quadrature.yw]=quadrature_rules_1D(Nm,cubature_type);
        
        % lower degree of a Gauss-Legendre rule with ADE equal to M.
        Nl_deg=ceil((M+1)/2);
        % In the quadrature rule, for M as input we have a rule with M+1 pts.
        Nl=Nl_deg-1;
        [quadrature.xn,quadrature.xw]=quadrature_rules_1D(Nl,cubature_type);

        
end






%--------------------------------------------------------------------------
% quadrature_rules_1D.
%--------------------------------------------------------------------------

function [nodes,weights]=quadrature_rules_1D(n,quadrature_type)

%--------------------------------------------------------------------------
% OBJECT:
%-----------
% THE TARGET OF THIS ROUTINE IS TO PROVIDE NODES AND WEIGHTS OF A
% QUADRATURE ROUTINE ON [-1,1]. THE CODES ARE DESCRIBED BY L.N. TREFETHEN
% IN HIS CLENSHAW-CURTIS PAPER AND BY WALDVOGEL (PUBLISHED BY BIT).
%
%--------------------------------------------------------------------------
% INPUTS:
%----------
%
% n : NUMBER OF NODES OF THE QUADRATURE RULE (NOT THE DEGREE!!).
%
% quadrature_type: IT POINTS A QUADRATURE RULE
%           [quadrature_type=1]: FEJER 1.
%           [quadrature_type=2]: FEJER 2.
%           [quadrature_type=3]: CLENSHAW CURTIS (VIA WEIGHTS).
%           [quadrature_type=4]: GAUSS-LEGENDRE.
%           [quadrature_type=5]: GAUSS-LEGENDRE-LOBATTO.
%           [quadrature_type=6]: COMPOSITE TRAPEZOIDAL RULE.
%
%----------
% OUTPUTS:
%----------
%
% nodes : M x 2 MATRIX OF NODES, IN THE INTERVAL [-1,1].
%
% weights: M x 1 COLUMN VECTOR OF WEIGHTS.
%
%--------------------------------------------------------------------------
% ADDITIONAL ROUTINES :
%----------------------
%
% 1. r_jacobi
% 2. gauss
% 3. lobatto_jacobi
%
% THESE ROUTINES ARE WRITTEN BY D. LAURIE AND W. GAUTSCHI, AND CAN BE FOUND
% IN W. GAUTSCHI HOMEPAGE. THEY ARE ATTACHED IN THIS FILE (SEE THE BOTTOM
% OF THE FILE). NO EXTERNAL FILE IS NEEDED.
%
%--------------------------------------------------------------------------
% EXAMPLE:
%----------
%
% [nodes,weights]=quadrature_rules_1D(5,4)
% 
% nodes =
% 
%    -0.9325
%    -0.6612
%    -0.2386
%     0.2386
%     0.6612
%     0.9325
% 
% 
% weights =
% 
%     0.1713
%     0.3608
%     0.4679
%     0.4679
%     0.3608
%     0.1713
%
%--------------------------------------------------------------------------
% RELATED PAPERS:
%-----------------
%
% [1] W. GAUTSCHI, ORTHOGONAL POLYNOMIALS AND QUADRATURE.
%
% [2] LLOYD N. TREFETHEN, IS GAUSS QUADRATURE BETTER THAN CLENSHAW?CURTIS?
%
% [3] J. WALDVOGEL, FAST CONSTRUCTION OF THE FEJER AND CLENSHAW-CURTIS
%     QUADRATURE RULES.
%
%--------------------------------------------------------------------------

switch quadrature_type
    
    case 1 % FEJER 1.
        N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';
        v0=[2*exp(i*pi*K/n)./(1-4*K.^2); zeros(l+1,1)];
        v1=v0(1:end-1)+conj(v0(end:-1:2));
        weights=ifft(v1);
        k=(1/2):(n-(1/2)); nodes=(cos(k*pi/n))';
        
    case 2 % FEJER 2.
        N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';
        v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
        v2=-v0(1:end-1)-v0(end:-1:2);
        wf2=ifft(v2); weights=[wf2;0];
        k=0:n; nodes=(cos(k*pi/n))';
        
    case 3 % CLENSHAW CURTIS.
        n=n-1;
        N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';
        g0=-ones(n,1); g0(1+l)=g0(1+l)+n; g0(1+m)=g0(1+m)+n;
        g=g0/(n^2-1+mod(n,2));
        v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
        v2=-v0(1:end-1)-v0(end:-1:2);
        wcc=ifft(v2+g); weights=[wcc;wcc(1,1)];
        k=0:n; nodes=(cos(k*pi/n))';
        
    case 4 % GAUSS LEGENDRE.
        beta=0.5./sqrt(1-(2*(1:n)).^(-2));
        T=diag(beta,1)+diag(beta,-1);
        [V,D]=eig(T);
        x=diag(D); [x,index]=sort(x); x=x';
        w=2*V(1,index).^2;
        nodes=x';
        weights=w';
        
    case 5 % GAUSS LEGENDRE LOBATTO.
        xw=lobatto_jacobi(n);
        nodes=xw(:,1);
        weights=xw(:,2);
        
    case 6 % TRAPEZOIDAL RULE.
        h=2/n;
        x=-1:h:1;
        w=h*[0.5 ones(1,length(x)-2) 0.5];
        nodes=x';
        weights=w';
        
end





%--------------------------------------------------------------------------
% ADDITIONAL ROUTINES BY W. GAUTSCHI.
%--------------------------------------------------------------------------

%---------------------------
% r_jacobi.
%---------------------------

function ab=r_jacobi(N,a,b)

nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
if N==1
    ab=[nu mu]; return
end

N=N-1;
n=1:N;
nab=2*n+a+b;
nuadd=(b^2-a^2)*ones(1,N)./(nab.*(nab+2));
A=[nu nuadd];
n=2:N;
nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
abadd=[mu; B1; B'];
ab=[A' abadd];





%---------------------------
% gauss.
%---------------------------

function xw=gauss(N,ab)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
    J(n,n-1)=sqrt(ab(n,2));
    J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);  %JV = VD
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];





%---------------------------
% lobatto_jacobi.
%---------------------------

function xw=lobatto_jacobi(N,a,b)
if nargin<2, a=0; end; if nargin<3, b=a; end
ab=r_jacobi(N+2,a,b);
ab(N+2,1)=(a-b)/(2*N+a+b+2);
ab(N+2,2)=4*(N+a+1)*(N+b+1)*(N+a+b+1)/((2*N+a+b+1)* ...
    (2*N+a+b+2)^2);
xw=gauss(N+2,ab);
