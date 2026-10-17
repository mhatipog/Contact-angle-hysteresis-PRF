function [f,aux] = fullhistory_cox_rhs(t,y,N,nu,A,ell,xref,M,mode)
% Full local angle-cubed Cox/Voinov-type kinetic law, not linearized.
% s_i*v_i = M*(theta_app_i^3-theta_Young_i^3), with s=[-1;1].
% Pressure multiplier kappa_app=cos(theta_app)/cref.
% Two moving lines: 2x2 implicit Newton; pinned-left: same with left B=0.
% Original hydrodynamic-beam fout.m unchanged.
cref=.94;
if ~(M>0&&ell>0&&abs(A)<.01&&mode(2)==2&&any(mode(1)==[0 2]))
  error('Unsupported dynamic configuration');
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
B=[Bl Br];Q=B(N+1:N+2,:);
S=[-1;1];
g=sin(2*pi*([y(N+1);y(N+2)]-xref)/ell);
kY=1+A*g;
thetaY=acos(cref*kY);
v=f0(N+1:N+2)+Q*(kY-1);
for iter=1:20
 svel=S.*v;
 cubic=thetaY.^3+svel/M;
 if any(cubic<=0)
   error('Angle-cube closure outside physical range');
 end
 ta=cubic.^(1/3);
 kappa=cos(ta)/cref;
 dkappa=-(sin(ta)./cref).*S./(3*M*ta.^2);
 residual=v-f0(N+1:N+2)-Q*(kappa-1);
 if max(abs(residual))<1e-10
    break;
 end
 jac=eye(2)-Q*diag(dkappa);
 dv=jac\residual;
 v=v-dv;
 if iter==20
    error('Cox boundary Newton method did not converge');
 end
end
f=f0+B*(kappa-1);
diss=sum(cref*(kY-kappa).*(S.*v));
if diss< -1e-13 || any(abs(cref*kappa)>1) || any(abs(cref*kY)>1)
  error('Invalid local kinetic dissipation or contact angle');
end
if nargout>1
 aux=struct('velocity',v,'theta_app',ta,'theta_Young',thetaY,...
    'young_kappa',kY,'pressure_kappa',kappa,...
    'contact_dissipation',diss,'iterations',iter);
end
end