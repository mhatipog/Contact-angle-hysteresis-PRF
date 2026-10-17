% Full patterned startup + nonlinear Cox law + original hydrodynamics.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));addpath(fileparts(mfilename('fullpath')));
nu=4;V=.2;xp0=.65;xref=.45;ell=.25;lambda=.0293;tol=.001;
for N=[12,16,20]
 start=[V*ones(N,1);xp0-V;xp0];
 for M=[300,1000]
  for A=[.002,0,-.002]
   ff=@(t,y)fullhistory_cox_rhs(t,y,N,nu,A,ell,xref,M,[2,2]);
   o1=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0001,...
     'Events',@(t,y)fullhistory_cox_pin_event(t,y,N,nu,A,ell,xref,M));
   initial=ode15s(ff,[0,.012],start,o1);
   assert(~isempty(initial.ie)&&initial.ie(1)==1,'Nonlinear Cox initial pin event not reached');
   tpin=initial.x(end);yp=initial.y(:,end);xm=yp(N+1);
   assert(max(abs(mean(initial.y(1:N,:),1)-V))<4e-5,'Initial mass conservation');
   [~,detail]=fullhistory_cox_rhs(tpin,yp,N,nu,A,ell,xref,M,[2,2]);
   assert(detail.contact_dissipation>=0,'Negative Cox contact dissipation');
   xp_ref=yp(N+2);
   fp=@(t,y)fullhistory_cox_rhs(t,y,N,nu,A,ell,xref,M,[0,2]);
   o2=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0001);
   all=ode15s(fp,[tpin,.1],yp,o2);
   err=max(abs(mean(all.y(1:N,:),1)-V));
   assert(err<4e-5,'Late mass conservation');
   Gamma=left_yield_observable(all.y(:,end),N,nu);
   eq=steady_patterned_yield(xm,V,nu,A,ell,xref);
   assert(abs(Gamma-eq.gamma)<.00018,'Nonlinear dynamics wrong static limit');
   oe=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.000035,...
    'Events',@(t,y)fullhistory_release_event(t,y,N,nu,A,ell,xref,lambda,tol));
   event=ode15s(fp,[tpin,.1],yp,oe);
   released=~isempty(event.ie);
   assert(released==(A<=0),'Opposite-phase release classification not preserved under Cox dynamics');
   if released,time=event.x(end);else,time=Inf;end
   fprintf('COX_FULL_HISTORY N=%d M=%g A=%+.4f t_pin=%.10g x_minus=%.10g right_at_pin=%.10g gamma_final=%.9g gamma_exact=%.9g released=%d release_time=%g dissipation_at_pin=%g mass_error=%g\n',...
      N,M,A,tpin,xm,xp_ref,Gamma,eq.gamma,released,time,detail.contact_dissipation,err);
  end
 end
end
fprintf('PASS full local angle-cubed Cox law retains phase-selective first release.\n');