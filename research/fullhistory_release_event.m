function [val,term,direction]=fullhistory_release_event(~,y,N,nu,A,ell,xref,lambda,tol)
% Local right-edge-of-pinned-interval pressure:
% kappa_R(xm)=1+lambda+A*sin(k*(xm-xref)).
gam=left_yield_observable(y,N,nu);
local=A*sin(2*pi*(y(N+1)-xref)/ell);
val=gam-lambda-local-tol*(1+lambda);
term=1;direction=1;
end