% Independent spectral gap validation using the unmodified MATLAB RHS.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));addpath(fileparts(mfilename('fullpath')));
V=.2;nu=4;xu0=.65;
for N=[12,16,20]
 ic=[V*ones(N,1);xu0-V;xu0];
 oa=odeset('RelTol',1e-7,'AbsTol',1e-8,'Events',...
   @(t,y)EventsFcn(t,y,N,nu,[2,2],.04,acosd(1/1.04),0));
 first=ode15s(@(t,y)fout(t,y,N,nu,[2,2],.04),[0,.012],ic,oa);
 assert(~isempty(first.ie) && first.ie(1)==1,'Unexpected initial pinning event');
 xmin=first.y(N+1,end);
 ob=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.001);
 sol=ode15s(@(t,y)fout(t,y,N,nu,[0,2],.04),[first.x(end),.10],first.y(:,end),ob);
 yy=sol.y(:,end);mass=mean(yy(1:N));
 q=[yy(1:N-1);yy(N+2)];
 R=@(z)pinned_branch_reduced_rhs(z,N,nu,xmin,mass);
 jac=zeros(N,N);h=3e-8;
 for j=1:N
   e=zeros(N,1);e(j)=h;
   jac(:,j)=(R(q+e)-R(q-e))/(2*h);
 end
 spectrum=sort(real(eig(jac)),'descend');
 slow=spectrum(1);second=spectrum(2);
 fprintf('SPECTRUM N=%d slow=%.10g time=%.10g second=%.10g ratio=%.8g residual=%.6g conservation=%.8g\n',...
   N,slow,-1/slow,second,abs(second/slow),max(abs(R(q))),mass);
 assert(slow< -250 && slow> -550,'Unexpected leading eigenvalue');
 assert(second < -50000,'No separated fast modes');
 assert(abs(second/slow)>100,'Insufficient spectral separation');
end
fprintf('PASS original MATLAB PDE has a grid-robust slow eigenmode.\n');
