% Real original MATLAB RHS + original hysteresis event switching, not a Python port.
% Figures 5 and paper Sec III predict a left receding event for lambda=0.02
% and continued left pinning for lambda=0.04, at nu=4 and V=0.2.
where=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(where),'Auxillary_functions'));
N=16; nu=4; V=0.2; t_end=0.03;
for lambda_max=[0.02,0.04]
    theta_a=acosd(1/(1+lambda_max)); theta_r=0;
    IC=[V*ones(N,1); 0.45; 0.65];
    modes=[2,2]; tnow=0; got_receding=false; got_first_pin=false;
    max_mass_error=0;
    fprintf('Original MATLAB dynamics: lambda_max=%.3f\n',lambda_max);
    for segment=1:18
        opts=odeset('RelTol',1e-7,'AbsTol',1e-8, ...
           'Events',@(t,z) EventsFcn(t,z,N,nu,modes,lambda_max,theta_a,theta_r));
        sol=ode15s(@(t,z) fout(t,z,N,nu,modes,lambda_max),[tnow,t_end],IC,opts);
        if isempty(sol.x)
            error('ODE solver produced no time output');
        end
        max_mass_error=max(max_mass_error,max(abs(mean(sol.y(1:N,:),1)-V)));
        new_t=sol.x(end); IC=sol.y(:,end);
        if ~isempty(sol.ie) && any(ismember(sol.ie,[5,6,7]))
            error('Unexpected terminal contact/wall event in t<=0.03 diagnostic');
        end
        if isempty(sol.ie)
            assert(abs(new_t-t_end)<1e-7,'Solver did not reach requested t_end');
            break;
        end
        event=sol.ie(1);
        previous=modes;
        modes=update_pinning_conditions(modes,sol);
        fprintf('  t=%.9g event=%d state=[%d,%d] -> [%d,%d] xm=%.9g xp=%.9g\n', ...
                new_t,event,previous(1),previous(2),modes(1),modes(2),IC(N+1),IC(N+2));
        assert(any(modes~=previous),'A switching event did not change the mode');
        assert(new_t>tnow+1e-10,'Non-progressing event loop');
        if previous(1)==2 && modes(1)==0, got_first_pin=true; end
        if previous(1)==0 && modes(1)==1, got_receding=true; end
        assert(modes(2)==2,'Right meniscus changed regime unexpectedly');
        tnow=new_t;
    end
    assert(got_first_pin,'No first left pinning event');
    if lambda_max<0.03
        assert(got_receding,'Expected left depinning for lambda=0.02');
    else
        assert(~got_receding,'Unexpected left depinning for lambda=0.04');
        assert(modes(1)==0,'Expected left to remain pinned until t=0.03');
    end
    assert(new_t>=t_end-1e-8,'Diagnostic stopped before t_end');
    assert(max_mass_error<1e-5,'Mass conservation error');
    fprintf('  PASS lambda=%.2f final_left=%d mass_error=%.4g\n',lambda_max,modes(1),max_mass_error);
end
