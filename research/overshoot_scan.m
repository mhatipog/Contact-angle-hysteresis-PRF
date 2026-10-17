% Exploratory scan: peak dynamic pinned-left yield demand vs stationary beam.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Auxillary_functions'));addpath(fileparts(mfilename('fullpath')));
V=.2;N=12;
for nu=[1,2,4,6]
 for xu0=[.50,.55,.60,.65,.70,.75]
  try
   initial=[V*ones(N,1);xu0-V;xu0];
   oa=odeset('RelTol',1e-7,'AbsTol',1e-8,'Events',...
     @(t,y)EventsFcn(t,y,N,nu,[2,2],.04,acosd(1/1.04),0));
   first=ode15s(@(t,y)fout(t,y,N,nu,[2,2],.04),[0,.03],initial,oa);
   if isempty(first.ie) || first.ie(1)~=1
     fprintf('SCAN nu=%g xu0=%g SKIP (no clean left pin)\n',nu,xu0);continue
   end
   tpin=first.x(end);y0=first.y(:,end);xm=y0(N+1);
   ob=odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',.002);
   pinned=ode15s(@(t,y)fout(t,y,N,nu,[0,2],.04),[tpin,.1],y0,ob);
   g=zeros(size(pinned.x));
   for j=1:length(g)
      y=pinned.y(:,j);dz=1/N;w=y(N+2)-y(N+1);
      pl=(7*y(1)-33*y(2)+62*y(3)-58*y(4)+27*y(5)-5*y(6))/(2*dz^4*w^5);
      ul=(35*y(1)-35*y(2)+21*y(3)-5*y(4))/16;
      g(j)=-pl*ul/(nu*w)-1;
   end
   eq=steady_pinned_equilibrium(xm,V,nu);
   F=fout(pinned.x(end),pinned.y(:,end),N,nu,[0,2],.04);
   fprintf('SCAN nu=%g xu0=%g xm=%.8g peak=%.8g eq=%.8g excess=%.8g g_end=%.8g xu_final=%.8g xu_dot=%.8g\n',...
      nu,xu0,xm,max(g),eq.lambda_e,max(g)-eq.lambda_e,g(end),pinned.y(N+2,end),F(N+2));
  catch err
    fprintf('SCAN nu=%g xu0=%g SKIP solver=%s\n',nu,xu0,err.message);
  end
 end
end
fprintf('PASS: completed exploratory scan (skips documented).\n');
