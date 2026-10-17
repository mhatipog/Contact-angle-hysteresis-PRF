% Exploratory long-time diagnostic using ORIGINAL (unmodified) MATLAB RHS.
% Not a classification of true infinite-time trapping and not paper N=25.
% Compare the archived Figure 11 threshold xplus0_escape~0.61117 at
% nu=4, lambda_max=0.02 with an N=16 trajectory at xplus0=0.60.
here=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'Auxillary_functions'));
N=16; nu=4; V=0.2; xp0=0.60; lambda_max=0.02; T=10;
theta_a=acosd(1/(1+lambda_max)); theta_r=0;
IC=[V*ones(N,1); xp0-V; xp0]; modes=[2,2]; tnow=0;
status='unknown'; maxerr=0; transitions=0;
fprintf('Finite-horizon original solver: N=%d nu=%g V=%g xp0=%g lambda=%g T=%g\n', ...
    N,nu,V,xp0,lambda_max,T);
for segment=1:35
  opts=odeset('RelTol',1e-7,'AbsTol',1e-8, ...
    'Events',@(t,z) EventsFcn(t,z,N,nu,modes,lambda_max,theta_a,theta_r));
  sol=ode15s(@(t,z) fout(t,z,N,nu,modes,lambda_max),[tnow,T],IC,opts);
  if isempty(sol.x) || sol.x(end)<=tnow+1e-12
    error('No progress in original solver segment %d',segment);
  end
  IC=sol.y(:,end);
  maxerr=max(maxerr,max(abs(mean(sol.y(1:N,:),1)-V)));
  tnow=sol.x(end);
  fprintf('  segment=%d t=%.9g xm=%.8g xp=%.8g mode=[%d,%d] event=', ...
    segment,tnow,IC(N+1),IC(N+2),modes(1),modes(2));
  if isempty(sol.ie)
    fprintf('none\n');
    if abs(tnow-T)<1e-7
      status='not escaped at finite horizon';
    else
      status='no terminal event';
    end
    break
  end
  fprintf('%d\n',sol.ie(1));
  if any(ismember(sol.ie,[5,6,7]))
    if any(sol.ie==5),status='escaped at right end';
    elseif any(sol.ie==6),status='reached left end';
    else,status='walls touching';end
    break
  end
  old=modes;
  modes=update_pinning_conditions(modes,sol);
  if all(old==modes)
    error('Original switching table produced no state change');
  end
  transitions=transitions+1;
end
fprintf('DIAGNOSTIC status=%s last_t=%.9g last_xm=%.9g last_xp=%.9g event_count=%d max_mass_error=%.4g\n', ...
   status,tnow,IC(N+1),IC(N+2),transitions,maxerr);
assert(maxerr<1e-4,'Poor mass conservation');
assert(~strcmp(status,'unknown'),'Unresolved diagnostic, too many switching cycles');
