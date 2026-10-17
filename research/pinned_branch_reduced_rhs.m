function out = pinned_branch_reduced_rhs(q,N,nu,xm,mass)
% Exact conserved-volume reduction of original unmodified fout.m.
U=[q(1:N-1); N*mass-sum(q(1:N-1))];
y=[U;xm;q(N)];
f=fout(0,y,N,nu,[0,2],.04);
out=[f(1:N-1);f(N+2)];
end
