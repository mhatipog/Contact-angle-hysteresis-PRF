% Validate first left-contact-line pin using the ORIGINAL MATLAB RHS and EventsFcn.
% Short event diagnostic; not a full escape/trapping classification.
here=fileparts(mfilename('fullpath'));
root=fileparts(here);
addpath(fullfile(root,'Auxillary_functions'));
N=16; nu=4; lambda_max=0.02; V=0.2; xp0=0.65;
theta_a=acosd(1/(1+lambda_max)); theta_r=0;
IC=[V*ones(N,1); xp0-V; xp0]; mode=[2,2];
opts=odeset('RelTol',1e-7,'AbsTol',1e-8, ...
    'Events',@(t,z) EventsFcn(t,z,N,nu,mode,lambda_max,theta_a,theta_r));
sol=ode15s(@(t,z) fout(t,z,N,nu,mode,lambda_max), [0,0.01], IC,opts);
t_first=sol.x(end);
if ~isfield(sol,'ie') || isempty(sol.ie) || ~any(sol.ie == 1)
    error('Expected left advancing-to-pinned event 1 in original solver');
end
xm=sol.y(N+1,end); xp=sol.y(N+2,end);
masserr=max(abs(mean(sol.y(1:N,:),1)-V));
assert(abs(t_first-0.00331)<6e-4,'First pinning time outside coarse benchmark');
assert(abs(xm-0.44544)<8e-4,'First pinned position outside coarse benchmark');
assert(abs(xp-0.65577)<8e-4,'Upper meniscus position outside coarse benchmark');
assert(masserr<5e-6,'Mass conservation error');
fprintf('Original ode15s first pin: t=%.9g, xminus=%.9g, xplus=%.9g, mass_err=%g\n',t_first,xm,xp,masserr);
