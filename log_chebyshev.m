N = 1000;
% range of w
wa = 0.01*pi;
wb = pi;
w = logspace(log(wa), wb, N);
cvx_begin
    variables R(N) a
    minimize a
    subject to 
        inv_pos(a) - (w(:)).*R(:) <=0
        (w(:)).*R(:) <= a
        R(:)>=0
 cvx_end
 semilogx(w, log(sqrt(R)))
 xlabel("w")
 ylabel("log(H(w))")
 