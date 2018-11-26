function [] = bobEx4p3 ()

k = pi*(0:0.001:1);
w1 = 1 + sin(k).^2;
w2 = 1 + 2*sin(k).^2 .* cos(k).^2;

plot(2*k/pi, w1, 'linewidth', 2.0, 2*k/pi, w2, 'linewidth', 2.0)
xlabel('k / k_{max}','fontsize',15)
ylabel('w^2/f^2','fontsize',15)
set(gca,'fontsize',15)
legend('Aligned','45-degree')
xlim([0,1])

end
