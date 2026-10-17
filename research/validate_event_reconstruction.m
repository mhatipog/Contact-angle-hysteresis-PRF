% Independent original-event vs first-passage reconstruction for the same pinned state.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));
N=16;nu=4;V=.2;xu0=.65;
IC=[V*ones(N,1);xu0-V;xu0];
optsA=odeset('RelTol',1e-7,'AbsTol',1e-8,...
 'Events',@(t,y)EventsFcn(t,y,N,nu,[2,2],.04,acosd(1/1.04),0));
s0=ode15s(@(t,y)fout(t,y,N,nu,[2,2],.04),[0,.012],IC,optsA);
assert(~isempty(s0.ie) && s0.ie(1)==1,'Unexpected first mode switch');
t0=s0.x(end);y0=s0.y(:,end);
opts=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0001);
s=ode15s(@(t,y)fout(t,y,N,nu,[0,2],.04),[t0,.1],y0,opts);
t=s.x;gamma=zeros(size(t));
for j=1:length(t)
 y=s.y(:,j);dz=1/N;w=y(N+2)-y(N+1);
 pxl=(7*y(1)-33*y(2)+62*y(3)-58*y(4)+27*y(5)-5*y(6))/(2*dz^4*w^5);
 uxl=(35*y(1)-35*y(2)+21*y(3)-5*y(4))/16;
 gamma(j)=-pxl*uxl/(nu*w)-1;
end
fprintf('ORIGINAL FIRST-PASSAGE reconstruction N=16, tpin=%.10g peak=%.10g\n',t0,max(gamma));
for lam=[.010,.020,.025,.028,.029,.0294,.0296,.0297,.030]
 threshold=lam+1e-3*(1+lam);
 ix=find(gamma>=threshold,1,'first');
 if isempty(ix)
   predicted=Inf;
 elseif ix==1
   predicted=t0;
 else
   predicted=interp1(gamma(ix-1:ix),t(ix-1:ix),threshold);
 end
 theta_a=acosd(1/(1+lam));
 optsE=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.0001,...
   'Events',@(tau,y)EventsFcn(tau,y,N,nu,[0,2],lam,theta_a,0));
 se=ode15s(@(tau,y)fout(tau,y,N,nu,[0,2],lam),[t0,.1],y0,optsE);
 if isempty(se.ie)
   actual=Inf;
 else
   assert(se.ie(1)==2,'A different event interrupts pinned branch');
   actual=se.x(end);
 end
 if isfinite(predicted) && isfinite(actual)
   fprintf('CROSSCHECK lambda=%.4f predicted=%.10g original_event=%.10g difference=%.3g\n',...
        lam,predicted,actual,abs(predicted-actual));
   assert(abs(predicted-actual)<.00025,'First passage fails original event verification');
 else
   fprintf('CROSSCHECK lambda=%.4f predicted=%g original_event=%g\n',lam,predicted,actual);
   assert(isinf(predicted)&&isinf(actual),'False negative or false positive departure');
 end
end
fprintf('PASS original event solver and first-passage construction agree.\n');
