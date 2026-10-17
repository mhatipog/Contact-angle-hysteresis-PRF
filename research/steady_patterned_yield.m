function eq = steady_patterned_yield(xm,V,nu,A,ell)
% Exact quartic clamped/free beam equilibrium with spatially variable
% right-meniscus advancing pressure factor kappa(xplus).
if ~(xm>0 && V>0 && nu>0 && ell>0 && abs(A)<1 && xm+V<1)
  error('Invalid physical parameters');
end
init=[xm+V+.01,-nu,-.04*nu,1-.02*nu];
opts=optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
[q,res,flag]=fsolve(@(q)pattern_residual(q,xm,V,nu,A,ell),init,opts);
xp=q(1);p=q(2);c=q(3);d=q(4);w=xm-xp;
hl=p*w^4/24+c*w+d;
if (max(abs(res))>1e-8 || xp<=xm || xp>=1 || d<=0 || p>=0 || hl<=0)
  error('No admissible patterned equilibrium');
end
assert(flag~=0,'Patterned equilibrium solver failed');
eq=struct('xp',xp,'p',p,'d',d,'hl',hl,'gamma',-p*hl/nu-1,...
          'residual',max(abs(res)),'kappa',1+A*sin(2*pi*(xp-xm)/ell));
end
function res=pattern_residual(q,xm,V,nu,A,ell)
xp=q(1);p=q(2);c=q(3);d=q(4);w=xm-xp;
hl=p*w^4/24+c*w+d;hp=p*w^3/6+c;
res=[p*w^2/2-2/xm^2*(2*xm*hp-3*hl+3);
     p*w-6/xm^3*(xm*hp-2*hl+2);
     p*d+nu*(1+A*sin(2*pi*(xp-xm)/ell));
     -p*w^5/120-c*w^2/2-d*w-V];
end
