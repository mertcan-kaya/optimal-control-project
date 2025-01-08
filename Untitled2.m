syms s t Y

K = [0.4142 1.3522];

f = -(K(1) + 1)

F = laplace(f,t,s)

Y1 = s*Y - 1

Y2 = s*Y1 + 1

Sol = solve(Y2 + K(2)*Y - F, Y)

sol = ilaplace(Sol,s,t)

ezplot(sol,[0,10])

fplot(sol)