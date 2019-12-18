function dydt = reducedlorenz96(t, x, Q)
lorenz96out = lorenz96(t, Q*x);
dydt = Q' * lorenz96out;
end