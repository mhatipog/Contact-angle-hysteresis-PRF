function [value,isterminal,direction] = left_yield_event(t,y,N,nu,lambda)
% Original EventsFcn pinned-left receding trigger rewritten as a single event.
tol=.001;
value=left_yield_observable(y,N,nu)-(lambda+tol*(1+lambda));
isterminal=1;
direction=+1;
end
