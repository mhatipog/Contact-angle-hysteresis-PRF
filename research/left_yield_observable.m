function gamma = left_yield_observable(y,N,nu)
dz=1/N;
w=y(N+2)-y(N+1);
p=(7*y(1)-33*y(2)+62*y(3)-58*y(4)+27*y(5)-5*y(6))/(2*dz^4*w^5);
u=(35*y(1)-35*y(2)+21*y(3)-5*y(4))/16;
gamma=-p*u/(nu*w)-1;
end
