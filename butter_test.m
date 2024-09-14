clear;
fc = 250;
fs = 1000;

[b,a] = butter(1,fc/(fs/2));

a
b

%% 
clear;
fc = 7;
fs = 1000;

[b,a] = butter(2,fc/(fs/2));

% freqz(b,a,[],fs)
% 
% subplot(2,1,1)
% ylim([-100 20])

w = 20;
t = linspace(0,1,100);
x = sin(2*pi*w*t);
y = filter(b,a,x);

plot(t,x);
hold on
plot(t,y)
legend('Input Data','Filtered Data')