function pfun = pfunTriFunc(x,y,z,t)
pfun=sin(pi.*x).*sin(pi.*y).*sin(pi.*z);
% if size(x,1) > 1e4
%     n=1e4;  %一次计算量
%     pfun=zeros(numel(x),1);
%     nums = floor(numel(x)/n);
%     for ii=1:nums
%         pfun(n*ii-n+1:n*ii)=bsxfun(@times,bsxfun(@times,sin(pi*x(n*ii-n+1:n*ii)),sin(pi*y(n*ii-n+1:n*ii))),sin(pi*z(n*ii-n+1:n*ii)));
%         if ii == nums && n*nums ~= numel(x)
%             pfun(n*nums+1:end)=bsxfun(@times,bsxfun(@times,sin(pi*x(n*nums+1:end)),sin(pi*y(n*nums+1:end))),sin(pi*z(n*nums+1:end)));
%         end
%     end
% else
%     pfun=sin(pi.*x).*sin(pi.*y).*sin(pi.*z);
% end
end