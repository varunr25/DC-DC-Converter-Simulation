%% Sw.m
frequency = 10000;
T_sw = 1 / frequency;
s = [0];

time = linspace(0, T_sw);
for n = 1:100
    s(n) = sw(.5, time(n));
end

plot(time,s);
axis([0,time(end), -1,2])

% Plot the triangular waveform
plot(time, s)
xlabel('Time [s]')
ylabel('x(t)')
title('sw.m')
grid on
legend('show')
disp("Complete");