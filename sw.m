% ID Number: 229,506
% ECE 31033 - Project #1
% sw.m

% The first file (sw.m) contains a function (sw) that accepts the duty cycle D, and a 
% single instant of time as an input, and outputs the state (on/off) of the transistor 
% at that time instant as an output. A Fourier series-based triangle wave that you 
% create within this function should be compared with the duty cycle D to 
% determine the state of the transistor. The output of the function is a 1 if the 
% transistor is to be turned on. It is a value of 0 if it is turned off. 

function state = sw(D, t)
    T_sw = 1 / 10000;
    
    w = 2 * pi / T_sw;

    a_k = 0;
    triangle_wave = 0.5;

    N = 200; % Number of Fourier terms.    

    k = 1;
    while k <= N
        z = k * w * T_sw; % Temporary variable; to simplify code for the coefficient.

        a_k = (2 * (4 * cos(0.5 * z) - 2 * cos(z) - 2)) / (z^2);
        triangle_wave = triangle_wave + a_k * cos(k * w * t);
        k = k + 1;
    end

    if D >= triangle_wave 
        state = 1;
    else
        state = 0;
    end
end