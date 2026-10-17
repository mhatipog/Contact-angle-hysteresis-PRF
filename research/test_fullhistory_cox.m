% Pointwise Cox-Voinov boundary validation close to equilibrium.
% An ALL-TIME nonlinear Cox simulation was attempted but the stiff early
% solver probes states beyond the physically admissible positive-angle
% cubic law. This only tests the closure at quasi-steady states.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));addpath(fileparts(mfilename('fullpath')));
nu=4;V=.2;xp0=.65;xref=.45;ell=.25;cref=.94;th0=acos(cref);
for N=[12,16,20]
 y=[V*ones(N,1);xp0-V;xp0];
 initopts=odeset('RelTol',1e-7,'AbsTol',1e-8,'Events',...
   @(t,z)EventsFcn(t,z,N,nu,[2,2],.04,acosd(1/1.04),0));
 start=ode15s(@(t,z)fout(t,z,N,nu,[2,2],.04),[0,.012],y,initopts);
 assert(~isempty(start.ie)&&start.ie(1)==1,'No initial pinning');
 xm=start.y(N+1,end);
 for A=[-.002,0,.002]
   fp=@(t,z)fullhistory_rhs(t,z,N,nu,A,ell,xref,0,[0,2]);
   sim=ode15s(fp,[start.x(end),.1],start.y(:,end),...
       odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0005));
   yeq=sim.y(:,end);
   for M=[100,300]
     [fc,aux]=fullhistory_cox_rhs(.1,yeq,N,nu,A,ell,xref,M,[0,2]);
     mu=tan(th0)/(3*M*th0^2);
     fl=fullhistory_rhs(.1,yeq,N,nu,A,ell,xref,mu,[0,2]);
     assert(aux.contact_dissipation>=-1e-13,'Negative dissipation near equilibrium');
     assert(abs(fc(N+1))<1e-12,'Left line moved while pinned');
     assert(abs(fc(N+2)-fl(N+2))<.001,...
         'Nonlinear Cox closure fails the linear tangent-limit check');
     assert(abs(fc(N+2))<.02,'Not sufficiently near stationary state');
     fprintf('COX_LOCAL N=%d M=%g A=%+.4f cox_velocity=%.8g linear_velocity=%.8g contact_diss=%.8g\n',...
         N,M,A,fc(N+2),fl(N+2),aux.contact_dissipation);
   end
 end
end
fprintf('PASS angle-cubed boundary law checked ONLY in its valid low-speed regime.\n');