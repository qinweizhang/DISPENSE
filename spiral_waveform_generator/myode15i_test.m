function myode15i_test
yy0 = [0; 0];
yyp0 = [0; 2];
tspan = [0:0.1: 100];

[tt, yy]  = ode15i(@my_ODE15i, tspan, yy0, yyp0);

plot(tt,yy);

 function f = my_ODE15i(t, q, dq)
        f = [dq(1)-q(2);10*dq(2) + 4*q(1)];
        
    

