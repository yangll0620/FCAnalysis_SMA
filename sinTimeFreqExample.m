f = 1;
x = 0:0.001:1/f*2;
plot(x,sin(2*pi*f*x))
hold on
plot(x,sin(2*pi*f*x-pi/2))

grid on

xticks([])
yticks([-1 0 1])

plot([1/(4*f) 1/(4*f)], [0 1], 'k--')
xlabel('time')

pos = get(gca, 'Position');
annotation(gcf,'doublearrow',[pos(1) 0.225],[0.8 0.8]);
annotation(gcf,'textbox',[0.15 0.75 0.1 0.1],'String',{'\Deltat'},'LineStyle','none');

saveas(gcf, 'sinExample', 'tif');