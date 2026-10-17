% Both contact lines patterned FROM THE INITIAL STATE; no prepared pin.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));addpath(fileparts(mfilename('fullpath')));
nu=4;V=.2;xp0=.65;xref=.45;ell=.25;lambda=.0293;tol=.001;CREF=.94;
for N=[12,16,20]
 init=[V*ones(N,1);xp0-V;xp0];
 a=fout(0,init,N,nu,[2,2],0);
 b=fullhistory_rhs(0,init,N,nu,0,ell,xref,0,[2,2]);
 assert(max(abs(a-b))<1e-9,'Failed exact recovery of unmodified solver');
 for mu=[0,.01]
  for A=[.002,0,-.002]
   ff=@(t,y)fullhistory_rhs(t,y,N,nu,A,ell,xref,mu,[2,2]);
   o1=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0001, ...
        'Events',@(t,y)fullhistory_pin_event(t,y,N,nu,A,ell,xref,mu));
   s=ode15s(ff,[0,.012],init,o1);
   if isempty(s.ie)||s.ie(1)~=1
      fprintf('FULL_HISTORY N=%d A=%g mu=%g FAILURE first event=%s\n',N,A,mu,mat2str(s.ie));
      error('The original first pinning scenario is not reachable');
   end
   tpin=s.x(end); yp=s.y(:,end); xm=yp(N+1);
   assert(max(abs(mean(s.y(1:N,:),1)-V))<3e-5,'Early mass conservation failed');
   [~,ex]=fullhistory_rhs(tpin,yp,N,nu,A,ell,xref,mu,[2,2]);
   assert(ex.contact_dissipation>=0 && ex.matrix_condition<1e6,'Kinetic matrix issue');
   for x=[xm,yp(N+2)]
       k=1+A*sin(2*pi*(x-xref)/ell);
       assert(CREF*k>0 && CREF*(k+lambda)<1,'Contact angle cos outside physical range');
   end
   fp=@(t,y)fullhistory_rhs(t,y,N,nu,A,ell,xref,mu,[0,2]);
   o2=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0001);
   run=ode15s(fp,[tpin,.1],yp,o2);
   assert(max(abs(mean(run.y(1:N,:),1)-V))<3e-5,'Late mass conservation failed');
   g=arrayfun(@(j)left_yield_observable(run.y(:,j),N,nu),1:length(run.x));
   gl=A*sin(2*pi*(xm-xref)/ell);
   threshold=lambda+tol*(1+lambda);
   crossing=any(g-gl>=threshold);
   eq=steady_patterned_yield(xm,V,nu,A,ell,xref);
   assert(abs(eq.gamma-g(end))<.00015,'Stationary PDE disagrees with beam BVP');
   ev=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.00003, ...
     'Events',@(t,y)fullhistory_release_event(t,y,N,nu,A,ell,xref,lambda,tol));
   piece=ode15s(fp,[tpin,.1],yp,ev);
   triggered=~isempty(piece.ie);
   assert(triggered==crossing,'Event detector and full pinned trace disagree');
   if triggered, tevent=piece.x(end);else,tevent=Inf;end
   fprintf('FULL_HISTORY N=%d mu=%.4g A=%+.4g t_pin=%.9g xm=%.9g eq_gamma=%.9g local=%.9g max_demand=%.9g released=%d release_time=%g\n', ...
     N,mu,A,tpin,xm,eq.gamma,gl,max(g-gl),triggered,tevent);
  end
 end
end
fprintf('PASS: patterned full-history first pin and first release tests.\n');