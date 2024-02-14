function sur = EPSurface2(s1,s2,t)
v1=s1;
v2=s2;
if size(t) == [1,1]
    ttp = t*ones(size(s1));
    sur=[v1.*sin(pi*ttp/2), v2.*sin(pi*ttp/2), (v1.^2./4 + v2.^2./9)./2.*ttp.^2];
else
    sur=[v1.*sin(pi*t/2), v2.*sin(pi*t/2), (v1.^2./4 + v2.^2./9)./2.*t.^2];
    %sur=[s1.*sin(pi*t/2), s2.*sin(pi*t/2), (s1.^2./4 + s2.^2./9)./2.*t.^2];
end
end