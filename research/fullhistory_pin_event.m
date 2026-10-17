function [val,term,direction]=fullhistory_pin_event(t,y,N,nu,A,ell,xref,mu)
f=fullhistory_rhs(t,y,N,nu,A,ell,xref,mu,[2,2]);
val=[f(N+1);f(N+2)];
term=[1;1];direction=[1;-1];
end