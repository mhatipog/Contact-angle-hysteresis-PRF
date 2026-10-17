function [f,a] = fullhistory_rhs(t,y,N,nu,A,ell,xref,mu,mode)
% Both-moving then pinned-left hybrid PDE; original fout.m is UNCHANGED.
% Physical advancing pressure 1+A*sin(k*(x-xref))-mu*s*velocity.
% Artificial receding modes in probes are only an algebraic device.
if ~(mu>=0 && abs(A)<.01 && ell>0 && mode(2)==2 && any(mode(1)==[0 2]))
 error('Invalid parameter or regime');
end
if mode(1)==2
  f0=fout(t,y,N,nu,[2,2],0);
  Bl=fout(t,y,N,nu,[1,2],1)-f0;
  Br=fout(t,y,N,nu,[2,1],1)-f0;
else
  f0=fout(t,y,N,nu,[0,2],0);
  Bl=zeros(size(f0));
  Br=fout(t,y,N,nu,[0,1],1)-f0;
end
B=[Bl Br]; S=diag([-1 1]); Q=B(N+1:N+2,:);
G=A*[sin(2*pi*(y(N+1)-xref)/ell);sin(2*pi*(y(N+2)-xref)/ell)];
v=(eye(2)+mu*Q*S)\(f0(N+1:N+2)+Q*G);
kap=1+G-mu*S*v;
f=f0+B*(kap-1);
if nargout>1
 a=struct('vel',v,'kappa',kap,'Young',1+G,...
   'contact_dissipation',mu*(v'*v),'matrix_condition',cond(eye(2)+mu*Q*S));
end
end