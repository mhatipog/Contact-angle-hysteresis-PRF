% Exact patterned BEAM + original FULL PDE + independently computed events.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));
addpath(fileparts(mfilename('fullpath')));
V=.2;nu=4;xp0=.65;ell=.25;lambda=.0293;threshold=lambda+.001*(1+lambda);
for N=[12,16,20]
  initial=[V*ones(N,1);xp0-V;xp0];
  o0=odeset('RelTol',1e-7,'AbsTol',1e-8,...
    'Events',@(t,y)EventsFcn(t,y,N,nu,[2,2],.04,acosd(1/1.04),0));
  first=ode15s(@(t,y)fout(t,y,N,nu,[2,2],.04),[0,.012],initial,o0);
  assert(~isempty(first.ie)&&first.ie(1)==1,'No initial left pin');
  tpin=first.x(end);ypin=first.y(:,end);xm=ypin(N+1);
  fprintf('TRANSITION N=%d left_pin=%.10g xm=%.10g threshold=%.10g\n',...
         N,tpin,xm,threshold);
  for A=[+.002,0,-.002]
    f=@(t,y)fout(t,y,N,nu,[0,1],A*sin(2*pi*(y(N+2)-xm)/ell));
    options=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0002);
    full=ode15s(f,[tpin,.1],ypin,options);
    gamma=left_yield_observable(full.y(:,end),N,nu);
    gammas=zeros(size(full.x));
    for k=1:length(gammas)
      gammas(k)=left_yield_observable(full.y(:,k),N,nu);
    end
    eq=steady_patterned_yield(xm,V,nu,A,ell);
    assert(abs(gamma-eq.gamma)<8e-5,'PDE stationary value disagrees with exact beam BVP');
    assert(max(abs(mean(full.y(1:N,:),1)-V))<2e-5,'Mass conservation failed');
    yield_triggered=max(gammas)>=threshold;
    expected=A<=0;
    assert(yield_triggered==expected,'Unexpected transition by first passage');
    eopts=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.00008,...
      'Events',@(t,y)left_yield_event(t,y,N,nu,lambda));
    response=ode15s(f,[tpin,.1],ypin,eopts);
    actual=~isempty(response.ie);
    assert(actual==expected,'Original discretized PDE event disagrees with stationary prediction');
    if actual
      assert(response.ie(1)==1,'Unexpected event');
      tau=response.x(end);
    else
      tau=Inf;
      assert(abs(response.x(end)-.1)<1e-7,'Undetected incomplete integration');
    end
    fprintf('TRANSITION N=%d phaseA=%+.4f eq_gamma=%.10g PDE_gamma=%.10g peak=%.10g first_event=%g expected=%d\n',...
       N,A,eq.gamma,gamma,max(gammas),tau,expected);
  end
end
fprintf('PASS: matched-extrema opposite-phase surfaces yield different first-depinning outcomes.\n');
