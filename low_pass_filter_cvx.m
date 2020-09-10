a = 1.01;
wp = 0.12*pi;
ws = 0.24*pi;
n = 30;
m = 40*n;
w = linspace(0, pi, m);
E_m = [ones(m, 1) 2*cos(w'*(1:n-1))];
x = find(w<wp);x = x(length(x));
y = find(w>ws);y = y(1);

cvx_begin
    variables r(n) R(length(w)) d
    minimize d
    subject to
        R == E_m * r
        R(y:length(R)) <= d
        1/a^2 <= R(1:x)
        R(1:x) <= a^2
        R(:) >= 0
cvx_end
d = sqrt(d);
hold on
plot(w, sqrt(R))
yline(a,"k--")
yline(1/a,"k--")
yline(d,"k--")
xlabel("w")
ylabel("|H(W)|")
legend("d = "+string(d))


        