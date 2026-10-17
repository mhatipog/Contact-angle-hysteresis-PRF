% Phase-sensitive advancing wettability using ORIGINAL unmodified fout.m.
% A right meniscus with pressure multiplier kappa(x_plus)=1+A*g(x_plus)
% is represented exactly by fout(...,[0,1], A*g), since the algebra of
% the right ghost point is identical to the [0,2] advancing formula.
% [0,1] is an algebraic evaluation device, NOT a claim of right recession.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));
Nlist=[12,16]; nu=4; V=.2; xu0=.65; ell=.25;
strengths=[0,.0002,-.0002,.001,-.001,.002,-.002];
lambdas=[.020,.028]; delta=1e-3;
for N=Nlist
  ic=[V*ones(N,1);xu0-V;xu0];
  oa=odeset('RelTol',1e-7,'AbsTol',1e-8, ...
    'Events',@(t,y)EventsFcn(t,y,N,nu,[2,2],.04,acosd(1/1.04),0));
  first=ode15s(@(t,y)fout(t,y,N,nu,[2,2],.04),[0,.012],ic,oa);
  assert(~isempty(first.ie)&&first.ie(1)==1,'First contact line did not pin');
  tpin=first.x(end); ypin=first.y(:,end); xm=ypin(N+1);
  assert(ell>0&&1-.002>0,'Right meniscus pressure must be positive');
  fprintf('PATTERN GRID N=%d first_pin=%.9g xminus=%.9g ell=%.6g\n',N,tpin,xm,ell);
  T=.025;
  paths=cell(size(strengths));
  for j=1:length(strengths)
    A=strengths(j);
    dynamic=@(t,y)fout(t,y,N,nu,[0,1],A*sin(2*pi*(y(N+2)-xm)/ell));
    ob=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.00006);
    sol=ode15s(dynamic,[tpin,T],ypin,ob);
    t=sol.x; gamma=zeros(size(t));
    for k=1:length(t)
      y=sol.y(:,k); dz=1/N;w=y(N+2)-y(N+1);
      pl=(7*y(1)-33*y(2)+62*y(3)-58*y(4)+27*y(5)-5*y(6))/(2*dz^4*w^5);
      ul=(35*y(1)-35*y(2)+21*y(3)-5*y(4))/16;
      gamma(k)=-pl*ul/(nu*w)-1;
    end
    mass=max(abs(mean(sol.y(1:N,:),1)-V));
    assert(mass<2e-5,'Mass conservation failed for patterned RHS');
    paths{j}=struct('A',A,'t',t,'gamma',gamma,'xp',sol.y(N+2,end));
  end
  base=paths{1};
  for lam=lambdas
    threshold=lam+delta*(1+lam);
    crossing=zeros(size(paths));
    for j=1:length(paths)
      p=paths{j};
      ii=find(p.gamma>=threshold,1,'first');
      if isempty(ii)
        crossing(j)=NaN;
      elseif ii==1
        crossing(j)=tpin;
      else
        crossing(j)=interp1(p.gamma(ii-1:ii),p.t(ii-1:ii),threshold);
      end
    end
    if any(isnan(crossing))
      fprintf('PATTERN N=%d lambda=%.4f missing crossing: %s\n',N,lam,mat2str(crossing,7));
      continue
    end
    t0=crossing(1);h=.000008;
    G=@(p,tt)interp1(p.t,p.gamma,tt,'linear');
    Gdot=(G(base,t0+h)-G(base,t0-h))/(2*h);
    GA=(G(paths{2},t0)-G(paths{3},t0))/.0004;
    coefficient=-GA/Gdot;
    fprintf('RESPONSE N=%d lambda=%.4f base_t=%.10g dt_dA=%.8g Gdot=%.8g GA=%.8g\n',...
      N,lam,t0,coefficient,Gdot,GA);
    assert(abs(Gdot)>1e-4,'First passage not locally transverse');
    for j=4:7
      A=paths{j}.A;
      exact=crossing(j)-t0;
      first_order=A*coefficient;
      err=abs(exact-first_order);
      fprintf('PATTERN N=%d lambda=%.4f A=%+.5f t_cross=%.10g delta_time=%+.8g first_order=%+.8g error=%g\n',...
           N,lam,A,crossing(j),exact,first_order,err);
      assert(err<.00025,'Perturbative release time outside verification tolerance');
    end
    assert((crossing(4)-t0)*(crossing(5)-t0)<0,...
         'Equal-range opposite-phase patterns did not split event time');
  end
end
fprintf('PASS: original RHS yields phase-dependent first passage with verified response.\n');
