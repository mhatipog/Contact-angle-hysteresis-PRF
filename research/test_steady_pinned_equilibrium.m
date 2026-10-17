% Cross-check Figure 5 steady-pinned limit without invoking ode15s/PMAD.
where=fileparts(mfilename('fullpath'));
addpath(where);
eq=steady_pinned_equilibrium(0.45,0.2,4);
assert(abs(eq.x_plus-0.660902319126)<1e-6);
assert(abs(eq.lambda_e-0.031496059736)<1e-6);
assert(eq.residual_max<1e-8);
assert(eq.lambda_e>0.02 && eq.lambda_e<0.04);
last=-Inf;
for xl=[0.3,0.4,0.45,0.5,0.55]
    eq=steady_pinned_equilibrium(xl,0.2,4);
    assert(eq.lambda_e>last);
    assert(eq.residual_max<1e-8);
    last=eq.lambda_e;
end
fprintf('Independent stationary bendotaxis benchmark passed.\n');
