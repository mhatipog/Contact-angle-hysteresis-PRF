% Parameter-free first-passage yield envelope in the ORIGINAL bendotaxis model.
% Copyright original MATLAB parts: Alex Bradley, MIT licence in upstream repo.
% We DO NOT modify the original RHS or event switching implementation.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));
addpath(fileparts(mfilename('fullpath')));
V=.2;nu=4;xu0=.65;obs=[.005,.010,.015,.020,.025,.028,.029,.0294,.0296,.0297,.030,.035,.040];
for N=[12,16,20]
  z0=[V*ones(N,1);xu0-V;xu0];
  theta_a=acosd(1/(1+.04));theta_r=0;
  opts=odeset('RelTol',1e-7,'AbsTol',1e-8,...
       'Events',@(t,z) EventsFcn(t,z,N,nu,[2,2],.04,theta_a,theta_r));
  sol0=ode15s(@(t,z)fout(t,z,N,nu,[2,2],.04),[0,.012],z0,opts);
  assert(~isempty(sol0.ie)&&sol0.ie(1)==1,'Unexpected first switching event');
  tpin=sol0.x(end);pinned=sol0.y(:,end);xm=pinned(N+1);
  opts2=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.002);
  sol=ode15s(@(t,z)fout(t,z,N,nu,[0,2],.04),[tpin,.10],pinned,opts2);
  tt=sol.x; gam=zeros(size(tt));
  for j=1:length(tt)
     v=sol.y(:,j); dz=1/N; w=v(N+2)-v(N+1);
     pxl=(7*v(1)-33*v(2)+62*v(3)-58*v(4)+27*v(5)-5*v(6))/(2*dz^4*w^5);
     uxl=(35*v(1)-35*v(2)+21*v(3)-5*v(4))/16;
     gam(j)=-pxl*uxl/(nu*w)-1;
  end
  eq=steady_pinned_equilibrium(xm,V,nu);
  [peak,ind]=max(gam);
  fprintf('GRID N=%d first_pin=%.10g xm=%.10g peak_Gamma=%.10g tpeak=%.10g Gamma_end=%.10g static=%.10g\n',...
          N,tpin,xm,peak,tt(ind),gam(end),eq.lambda_e);
  for k=1:length(obs)
     target=obs(k)+1e-3*(1+obs(k));
     iz=find(gam>=target,1,'first');
     if isempty(iz)
        fprintf('  lambda=%.3f predicted first depin: none through t=0.10\n',obs(k));
     else
        if iz==1,tdep=tpin;
        else,tdep=interp1(gam(iz-1:iz),tt(iz-1:iz),target);end
        fprintf('  lambda=%.3f predicted depin=%.10g\n',obs(k),tdep);
     end
  end
  if N==16
    target=.02+1e-3*(1+.02);
    iz=find(gam>=target,1,'first');
    assert(~isempty(iz),'Expected .02 depinning in pinned trace');
    tdep=interp1(gam(iz-1:iz),tt(iz-1:iz),target);
    assert(abs(tdep-.00647691109)<.00035,'First passage disagrees with original event-driven solver');
    assert(peak>.02 && peak<.07,'Gamma peak unreasonable');
    assert(abs(eq.lambda_e-.030746643)<2e-3,'Stationary benchmark disagrees');
  end
end
fprintf('PASS: parameter-free pinned branch predicts contact-line depinning.\n');
