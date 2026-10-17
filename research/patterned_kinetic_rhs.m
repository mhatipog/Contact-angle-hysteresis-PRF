function f = patterned_kinetic_rhs(t,y,N,nu,xm,A,ell,mu)
% Schur-complement elimination of advancing Young-force kinetic response:
% kappa_app = 1 + A sin(k*(xplus-xm)) - mu * dxplus/dt, mu >= 0.
%
% The original unmodified fout(...,[0,1],lambda) is EXACTLY AFFINE in
% lambda on a pinned-left / right-free branch, and agrees algebraically
% with the right-advancing pressure boundary for kappa_app=1+lambda.
% Mode [0,1] is an algebraic pressure proxy only: physically advancing.
%
% Set F(y,kappa)=F0(y)+(kappa-1)B(y), where
% B=F(y,2)-F(y,1), and let v=e_+^T F.
% Then v=(v0+B_+*A*g)/(1+mu*B_+) and F=F0+B*(A*g-mu*v).
if nargin~=8 || ~(mu>=0 && abs(A)<1 && ell>0)
  error('Require mu>=0, abs(A)<1 and ell>0');
end
f0=fout(t,y,N,nu,[0,1],0);
f1=fout(t,y,N,nu,[0,1],1);
B=f1-f0;
g=sin(2*pi*(y(N+2)-xm)/ell);
den=1+mu*B(N+2);
if ~(den>0)
  error('Kinetic contact-line closure lost positive mobility');
end
velocity=(f0(N+2)+B(N+2)*A*g)/den;
f=f0+B*(A*g-mu*velocity);
end
