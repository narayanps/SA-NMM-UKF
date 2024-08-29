close all
set(0,'defaultTextInterpreter','latex'); %trying to set the default

col_y = [30 144 255]./255;
T=20000;
dt=0.001;
fs=1/dt;
t=1/fs:1/fs:T/fs;
subplot(2,1,1)
plot(t,y, 'LineWidth', 2, 'color', col_y)
xlabel('Time[s]')
ylabel('Amplitude')

subplot(2,1,2)
col_B = [255 140 0]./255;
col_b=[138 43 226]./255;
plot(t,B, 'LineWidth', 2, 'color', col_B)
hold on
plot(t,b, 'LineWidth', 2, 'color', col_b)
xlabel('Time[s]')
ylim([30 120])
ylabel('$B [mV], b [s^{-1}] $')
legend('B', 'b')

exportfig(gcf,strcat('model_output.eps'),...
    'width',24, 'height',14, ...
    'color','rgb', 'FontSize', 18)
