% Diagnostic: event-time convergence for negative-phase pattern.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));addpath(fileparts(mfilename('fullpath')));
V=.2;nu=4;xp0=.65;ell=.25;lambda=.0293;A=-.002;
for N=[12,16,20]
  y=[V*ones(N,1);xp0-V;xp0];
  initialopts=odeset('RelTol',1e-7,'AbsTol',1e-8,...
    'Events',@(t,z)EventsFcn(t,z,N,nu,[2,2],.04,acosd(1/1.04),0));
  s=ode15s(@(t,z)fout(t,z,N,nu,[2,2],.04),[0,.012],y,initialopts);
  assert(~isempty(s.ie)&&s.ie(1)==1);
  tpin=s.x(end);ypin=s.y(:,end);xm=ypin(N+1);
  f=@(t,z)fout(t,z,N,nu,[0,1],A*sin(2*pi*(z(N+2)-xm)/ell));
  times=[];
  for maxstep=[.0002,.00008,.00004,.00002]
    opts=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',maxstep,...
       'Events',@(t,z)left_yield_event(t,z,N,nu,lambda));
    z=ode15s(f,[tpin,.03],ypin,opts);
    assert(~isempty(z.ie),'Expected negative-phase first event');
    times(end+1)=z.x(end);
    fprintf('TIME_CONVERGENCE N=%d maxstep=%.8g t=%.12g event=%d\n',...
       N,maxstep,z.x(end),z.ie(1));
  end
  spread=max(times)-min(times);
  fprintf('TIME_CONVERGENCE N=%d variation=%.10g refined=%.10g\n',N,spread,times(end));
end
fprintf('PASS: sampled event-time sensitivity to integration step.\n');
