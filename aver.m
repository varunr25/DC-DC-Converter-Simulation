% ID Number: 229,506
% ECE 31033 - Project #1
% aver.m

% The fourth file (aver.m) contains a function you create to compute the average of 
% a waveform. Specifically, the function is of the form 
%       function av = aver(x,T,dt) 
% where x is the waveform to be averaged, T is its period, and dt is the period of time 
% between samples. This function must use the last period of the input waveform to 
% calculate the average. 

function av = aver(x, T, dt)
    location = length(x);
    av = 0;
    time = 0;

    while (time <= T)
        av = av + dt * (x(location));
        time = time + dt;
        location = location - 1;
    end

    av = av / T;
end
