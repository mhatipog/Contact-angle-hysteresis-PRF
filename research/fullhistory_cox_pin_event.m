function [val,term,direction]=fullhistory_cox_pin_event(t,y,N,nu,A,ell,xref,M)
f=fullhistory_cox_rhs(t,y,N,nu,A,ell,xref,M,[2,2]);
val=[f(N+1);f(N+2)];term=[1;1];direction=[1;-1];
end