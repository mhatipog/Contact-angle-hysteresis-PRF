% Full patterned startup at valid linear low-speed closure, followed by
% unlinearized Cox-Voinov cubic law ONLY after left contact line is pinned.
% Original initial meniscus speed is O(40-70), so applying quasi-steady
% Cox law there would exceed its physically admissible positive-angle range.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));addpath(fileparts(mfilename('fullpath')));
nu=4;V=.2;xp0=.65;xref=.45;ell=.25;lambda=.0293;tol=.001;
cref=.94;theta0=acos(cref);
for N=[12,16,20]
 initial=[V*ones(N,1);xp0-V;xp0];
 for M=[100,300]
  mu=tan(theta0)/(3*M*theta0^2); % exact tangent to angle-cubed law at v=0
  for A=[.002,0,-.002]
   startdyn=@(t,y)fullhistory_rhs(t,y,N,nu,A,ell,xref,mu,[2,2]);
   o1=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0001,...
     'Events',@(t,y)fullhistory_pin_event(t,y,N,nu,A,ell,xref,mu));
   early=ode15s(startdyn,[0,.012],initial,o1);
   assert(~isempty(early.ie)&&early.ie(1)==1,'Initial patterned pinning is not reachable');
   tpin=early.x(end);ypin=early.y(:,end);xm=ypin(N+1);
   [~,aux0]=fullhistory_cox_rhs(tpin,ypin,N,nu,A,ell,xref,M,[0,2]);
   assert(aux0.contact_dissipation>=0,'Negative contact-line dissipation');
   slow=@(t,y)fullhistory_cox_rhs(t,y,N,nu,A,ell,xref,M,[0,2]);
   op=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0001);
   traj=ode15s(slow,[tpin,.1],ypin,op);
   error_mass=max(abs(mean(traj.y(1:N,:),1)-V));
   assert(error_mass<4e-5,'Mass conservation error');
   eq=steady_patterned_yield(xm,V,nu,A,ell,xref);
   Gamma=left_yield_observable(traj.y(:,end),N,nu);
   assert(abs(Gamma-eq.gamma)<.00018,'Slow Cox branch equilibrium differs from exact beam');
   oevent=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.00005,...
     'Events',@(t,y)fullhistory_release_event(t,y,N,nu,A,ell,xref,lambda,tol));
   s=ode15s(slow,[tpin,.1],ypin,oevent);
   released=~isempty(s.ie);
   assert(released==(A<=0),'Pattern sign classification fails full nonlinear late-stage law');
   if released,trelease=s.x(end);else,trelease=Inf;end
   fprintf('COX_SLOW N=%d M=%g A=%+.4f tpin=%.10g xm=%.10g Gamma=%.10g eq=%.10g released=%d trelease=%g dissip=%g mass_err=%g\n',...
     N,M,A,tpin,xm,Gamma,eq.gamma,released,trelease,aux0.contact_dissipation,error_mass);
  end
 end
end
fprintf('PASS full patterned startup and later nonlinear Cox evolution.\n');