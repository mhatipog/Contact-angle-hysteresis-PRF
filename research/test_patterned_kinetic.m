% Check rate-dependent advancing contact line against original PDE.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));addpath(fileparts(mfilename('fullpath')));
nu=4;V=.2;xp0=.65;ell=.25;lambda=.0293;target=lambda+.001*(1+lambda);
for N=[12,16]
  initial=[V*ones(N,1);xp0-V;xp0];
  oa=odeset('RelTol',1e-7,'AbsTol',1e-8,...
   'Events',@(t,z)EventsFcn(t,z,N,nu,[2,2],.04,acosd(1/1.04),0));
  s=ode15s(@(t,z)fout(t,z,N,nu,[2,2],.04),[0,.012],initial,oa);
  assert(~isempty(s.ie)&&s.ie(1)==1);
  tpin=s.x(end);y0=s.y(:,end);xm=y0(N+1);
  for A=[.002,-.002]
    full=fout(tpin,y0,N,nu,[0,1],A*sin(2*pi*(y0(N+2)-xm)/ell));
    alg=patterned_kinetic_rhs(tpin,y0,N,nu,xm,A,ell,0);
    assert(max(abs(full-alg))<.0001,'Not algebraically equivalent at zero kinetic drag');
  end
  for mu=[0,.002,.01]
    for A=[.002,-.002]
      f=@(t,y)patterned_kinetic_rhs(t,y,N,nu,xm,A,ell,mu);
      ob=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.00015);
      s2=ode15s(f,[tpin,.1],y0,ob);
      eq=steady_patterned_yield(xm,V,nu,A,ell);
      g=left_yield_observable(s2.y(:,end),N,nu);
      assert(abs(g-eq.gamma)<.00012,'Kinetic model equilibrium disagrees with exact static beam');
      first=NaN;
      for j=1:length(s2.x)
        yi=s2.y(:,j);gg=left_yield_observable(yi,N,nu);
        if gg>=target
          first=s2.x(j);break
        end
      end
      if A>0
        assert(isnan(first),'Positive-phase pattern should remain pinned');
      else
        assert(isfinite(first),'Negative-phase pattern should release');
      end
      fprintf('KINETIC N=%d mu=%g A=%+.4f first_t=%g final_Gamma=%.10g exact_equilibrium=%.10g\n',...
         N,mu,A,first,g,eq.gamma);
    end
  end
end
fprintf('PASS: static phase boundary survives dissipative kinetic closure.\n');
