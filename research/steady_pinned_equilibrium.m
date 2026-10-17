function eq = steady_pinned_equilibrium(xl,V,nu)
% Independent steady pinned-left / advancing-right equilibrium benchmark.
% Bradley, Hewitt & Vella, Phys. Rev. Fluids 6, 114003 (2021), Eqs. (33)-(36).
% Neither the original time-dependent solver nor a published data file is used.
if nargin ~= 3 || ~(xl>0 && V>0 && nu>0 && xl+V<1)
    error('Require 0 < x_minus and volume and nu, with x_minus+volume < 1');
end
init = [xl+V+0.01, -nu, -0.04*nu, 1-0.02*nu];
opts = optimset('TolFun',1e-11,'TolX',1e-11,'Display','off');
[q, f, flag] = fsolve(@(q) equilibrium_residual(q,xl,V,nu), init, opts);
xu=q(1); p=q(2); c=q(3); d=q(4);
err=max(abs(f));
if flag <= 0 || err > 1e-8 || ~(xu>xl && xu<1 && d>0 && p<0)
    error('No verified physical equilibrium found');
end
w=xl-xu;
leftgap = p*w^4/24 + c*w + d;
if min([leftgap,d]) <= 0
    error('Channel gap nonpositive');
end
eq=struct('x_minus',xl,'x_plus',xu,'V',V,'nu',nu,...
          'pressure',p,'slope',c,'gap_plus',d,...
          'lambda_e',leftgap/d-1,'residual_max',err);
end

function F = equilibrium_residual(q,xl,V,nu)
xu=q(1);p=q(2);c=q(3);d=q(4);w=xl-xu;
hl=p*w^4/24 + c*w + d;
dhl=p*w^3/6 + c;
d2hl=p*w^2/2;
d3hl=p*w;
volume=-p*w^5/120-c*w^2/2-d*w;
F=[d2hl-2/xl^2*(2*xl*dhl-3*hl+3);
   d3hl-6/xl^3*(xl*dhl-2*hl+2);
   p*d+nu;
   volume-V];
end
