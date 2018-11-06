function [] = bobEx2p3 ()

dx = power(10, -4:-1);
for i = 1:length(dx)
    x = (0:dx(i):1-dx)';
    f = sin(2*pi*x);
    exact = -(2*pi)^3 * cos(2*pi*x);
    fd = Dxxx(f, dx(i));
    
    rms(i) = sqrt(mean((fd - exact).^2));
    figure;
    plot(x,fd,x,exact)
end

rms
figure;
loglog(dx, rms)

end

function deriv = Dxxx(solution, dx)

F = zeros(length(solution) + 3, 1);
F(1:2) = solution(end-1:end);
F(end) = solution(1);
F(3:end-1) = solution;

deriv = (-F(1:end-3) + 3*F(2:end-2) - 3*F(3:end-1) + F(4:end)) / dx^3;

end
