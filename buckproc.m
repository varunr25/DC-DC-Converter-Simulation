% ID Number: 229,506
% ECE 31033 - Project #1
% buckproc.m

% The file buckproc.m first contains the circuit parameter values (i.e. L, C, fsw, time 
% step, initial conditions etc.). Only the initial value of your circuit voltages and 
% currents should be pre-established (i.e. Vload(1)=0). It then invokes buck. Finally, 
% it performs your plotting and any post-processing calculations that are done 
% using the simulated data (such as computing average values, efficiency, etc.).
%% Ideal - Given Values
V_in = 800;
V_load_avg = 400;
V_load_ripple = 10;
P_load_light = 50000;
P_load_heavy = 250000;
frequency = 10000;

ideal_boolean = 1; % = 0 if non ideal, = 1 if ideal; here, it is ON.

%% Ideal - Calculated Values
T_sw = 1 / frequency;
D = V_load_avg / V_in;  % Duty Cycle

R_load_light = (V_load_avg^2) / P_load_light;
R_load_heavy = (V_load_avg^2) / P_load_heavy;

L_crit = (R_load_light * (1 - D)) / (2 * frequency);
L = L_crit * 1.1;

C = (V_load_avg / V_load_ripple) * (T_sw^2 * (1 - D)) / (8 * L);

I_load_light = V_load_avg / R_load_light;
I_load_heavy = V_load_avg / R_load_heavy;

i_L1_light = (V_load_avg / R_load_light) - (1 - D) * T_sw * V_load_avg / (2 * L);
i_L2_light = (V_load_avg / R_load_light) + (1 - D) * T_sw * V_load_avg / (2 * L);

i_L1_heavy= (V_load_avg / R_load_heavy) - (1 - D) * T_sw * V_load_avg / (2 * L);
i_L2_heavy= (V_load_avg / R_load_heavy) + (1 - D) * T_sw * V_load_avg / (2 * L);

%% Buck Intialization - Heavy Load
% Initializing Values
k = 1;
t = 0; 
dt = 1e-7;

tend = 100 * T_sw;

% Zero Vectors (used in buck)
t_vec = [0];  
switch_state = [0];

V_L_vec = [0];
i_L_vec = [0];

V_C_vec = [0];
i_C_vec = [0];

V_load_vec = [0]; 
i_load_vec = [0];

V_switch1 = [0];
i_switch1 = [0];

V_switch2 = [0];
i_switch2 = [0];

%% Running Buck - Using R_load_heavy
R_load = R_load_heavy;
disp('Running buck for heavy load.');
buck

%% Post-processing Calculations (computing avg values, efficiency, etc)
disp("---------------------")
disp("Heavy Averages:")

V_load_avg_func_H = aver(V_load_vec, T_sw, dt);
disp("  V_load Average: " + V_load_avg_func_H);

i_load_avg_func_H = aver(i_load_vec, T_sw, dt);
disp("  i_load Average: " + i_load_avg_func_H);

V_L_func_H = aver(V_L_vec, T_sw, dt);
disp("  V_L Average: " + V_L_func_H);

i_L_func_H = aver(i_L_vec, T_sw, dt);
disp("  i_L Average: " + i_L_func_H);

V_C_func_H = aver(V_C_vec, T_sw, dt);
disp("  V_C Average: " + V_C_func_H);

i_C_func_H = aver(i_C_vec, T_sw, dt);
disp("  i_C Average: " + i_C_func_H);

V_sw1_func_H = aver(V_switch1, T_sw, dt);
disp("  V_sw1 Average: " + V_sw1_func_H);

i_sw1_func_H = aver(i_switch1, T_sw, dt);
disp("  i_sw1 Average: " + i_sw1_func_H);

V_sw2_func_H = aver(V_switch2, T_sw, dt);
disp("  V_sw2 Average: " + V_sw2_func_H);

i_sw2_func_H = aver(i_switch2, T_sw, dt);
disp("  i_sw2 Average: " + i_sw2_func_H);

P_out_H = (V_load_avg_func_H^2) / R_load;
P_in_H = V_in * i_sw1_func_H;
eff_H = P_out_H / P_in_H;

disp("Efficiency for Light Load: " + (eff_H * 100) + "%.");
disp("---------------------")

%% Plotting - Heavy Load - Transient
% Plots for the transient to steady state
figure;
sgtitle("Heavy Load: Voltage and Current Plots for All Components at Transient State");
% Plots for the Load
subplot(5,2,1);
plot(t_vec, V_load_vec);
title('Load Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,2);
plot(t_vec, i_load_vec);
title('Load Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for the Inductor
subplot(5,2,3);
plot(t_vec, V_L_vec);
title('Inductor Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,4);
plot(t_vec, i_L_vec);
title('Inductor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for the Capacitor
subplot(5,2,5);
plot(t_vec, V_C_vec);
title('Capacitor Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,6);
plot(t_vec, i_C_vec);
title('Capacitor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 1
subplot(5,2,7);
plot(t_vec, V_switch1);
title('Transistor Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,8);
plot(t_vec, i_switch1);
title('Transistor Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 2 
subplot(5,2,9);
plot(t_vec, V_switch2);
title('Diode Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,10);
plot(t_vec, i_switch2);
title('Diode Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

%% Plotting - Heavy Load - Steady State
periods_to_plot = 2; 

points_per_period = round(T_sw / dt);  % Points per period
total_periods = floor(tend / T_sw);    % Total number of periods in the simulation

start_index = max(1, (total_periods - periods_to_plot) * points_per_period + 1);
end_index = min(length(t_vec), total_periods * points_per_period);

range_to_plot = start_index:end_index;

%% Plot
figure;
sgtitle("Heavy Load: Voltage and Current Plots for All Components at Steady State");
% Plots for the Load
subplot(5,2,1);
plot(t_vec(range_to_plot), V_load_vec(range_to_plot));
title('Load Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,2);
plot(t_vec(range_to_plot), i_load_vec(range_to_plot));
title('Load Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for the Inductor
subplot(5,2,3);
plot(t_vec(range_to_plot), V_L_vec(range_to_plot));
title('Inductor Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,4);
plot(t_vec(range_to_plot), i_L_vec(range_to_plot));
title('Inductor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for the Capacitor
subplot(5,2,5);
plot(t_vec(range_to_plot), V_C_vec(range_to_plot));
title('Capacitor Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,6);
plot(t_vec(range_to_plot), i_C_vec(range_to_plot));
title('Capacitor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 1
subplot(5,2,7);
plot(t_vec(range_to_plot), V_switch1(range_to_plot));
title('Transistor Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,8);
plot(t_vec(range_to_plot), i_switch1(range_to_plot));
title('Transistor Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 2 
subplot(5,2,9);
plot(t_vec(range_to_plot), V_switch2(range_to_plot));
title('Diode Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,10);
plot(t_vec(range_to_plot), i_switch2(range_to_plot));
title('Diode Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Buck Intialization - Light Load
% Initializing Values
k = 1;
t = 0; 
dt = 1e-7;

tend = 250 * T_sw;

% Zero Vectors (used in buck)
t_vec = [0];  
switch_state = [0];

V_L_vec = [0];
i_L_vec = [0];

V_C_vec = [0];
i_C_vec = [0];

V_load_vec = [0]; 
i_load_vec = [0];

V_switch1 = [0];
i_switch1 = [0];

V_switch2 = [0];
i_switch2 = [0];

%% Running Buck - Using R_load_light
R_load = R_load_light;
disp('Running buck for light load.');
buck

%% Post-processing Calculations (computing avg values, efficiency, etc)
disp("---------------------")
disp("Light Averages:")

V_load_avg_func_L = aver(V_load_vec, T_sw, dt);
disp("  V_load Average: " + V_load_avg_func_L);

i_load_avg_func_L = aver(i_load_vec, T_sw, dt);
disp("  i_load Average: " + i_load_avg_func_L);

V_L_func_L = aver(V_L_vec, T_sw, dt);
disp("  V_L Average: " + V_L_func_L);

i_L_func_L = aver(i_L_vec, T_sw, dt);
disp("  i_L Average: " + i_L_func_L);

V_C_func_L = aver(V_C_vec, T_sw, dt);
disp("  V_C Average: " + V_C_func_L);

i_C_func_L = aver(i_C_vec, T_sw, dt);
disp("  i_C Average: " + i_C_func_L);

V_sw1_func_L = aver(V_switch1, T_sw, dt);
disp("  V_sw1 Average: " + V_sw1_func_L);

i_sw1_func_L = aver(i_switch1, T_sw, dt);
disp("  i_sw1 Average: " + i_sw1_func_L);

V_sw2_func_L = aver(V_switch2, T_sw, dt);
disp("  V_sw2 Average: " + V_sw2_func_L);

i_sw2_func_L = aver(i_switch2, T_sw, dt);
disp("  i_sw2 Average: " + i_sw2_func_L);

P_out_L = (V_load_avg_func_L^2) / R_load;
P_in_L = V_in * i_sw1_func_L;
eff_L = P_out_L / P_in_L;

disp("Efficiency for Light Load: " + (eff_L * 100) + "%.");
disp("---------------------")
%% Plotting - Light Load - Transient
% Plots for the transient to steady state
figure;
sgtitle("Light Load: Voltage and Current Plots for All Components at Transient State");
% Plots for the Load
subplot(5,2,1);
plot(t_vec, V_load_vec);
title('Load Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,2);
plot(t_vec, i_load_vec);
title('Load Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for the Inductor
subplot(5,2,3);
plot(t_vec, V_L_vec);
title('Inductor Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,4);
plot(t_vec, i_L_vec);
title('Inductor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for the Capacitor
subplot(5,2,5);
plot(t_vec, V_C_vec);
title('Capacitor Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,6);
plot(t_vec, i_C_vec);
title('Capacitor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 1
subplot(5,2,7);
plot(t_vec, V_switch1);
title('Transistor Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,8);
plot(t_vec, i_switch1);
title('Transistor Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 2 
subplot(5,2,9);
plot(t_vec, V_switch2);
title('Diode Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,10);
plot(t_vec, i_switch2);
title('Diode Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

%% Plotting - Light Load - Steady State
periods_to_plot = 2; 

points_per_period = round(T_sw / dt);  % Points per period
total_periods = floor(tend / T_sw);    % Total number of periods in the simulation

start_index = max(1, (total_periods - periods_to_plot) * points_per_period + 1);
end_index = min(length(t_vec), total_periods * points_per_period);

range_to_plot = start_index:end_index;

%% Plot
figure;
sgtitle("Light Load: Voltage and Current Plots for All Components at Steady State");
% Plots for the Load
subplot(5,2,1);
plot(t_vec(range_to_plot), V_load_vec(range_to_plot));
title('Load Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,2);
plot(t_vec(range_to_plot), i_load_vec(range_to_plot));
title('Load Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for the Inductor
subplot(5,2,3);
plot(t_vec(range_to_plot), V_L_vec(range_to_plot));
title('Inductor Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,4);
plot(t_vec(range_to_plot), i_L_vec(range_to_plot));
title('Inductor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for the Capacitor
subplot(5,2,5);
plot(t_vec(range_to_plot), V_C_vec(range_to_plot));
title('Capacitor Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,6);
plot(t_vec(range_to_plot), i_C_vec(range_to_plot));
title('Capacitor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 1
subplot(5,2,7);
plot(t_vec(range_to_plot), V_switch1(range_to_plot));
title('Transistor Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,8);
plot(t_vec(range_to_plot), i_switch1(range_to_plot));
title('Transistor Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 2 
subplot(5,2,9);
plot(t_vec(range_to_plot), V_switch2(range_to_plot));
title('Diode Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(5,2,10);
plot(t_vec(range_to_plot), i_switch2(range_to_plot));
title('Diode Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non Ideal - Given Values
V_T_on = 1;
V_D_on = 1;
R_T_on = 0.01;
R_D_on = 0.01;

ideal_boolean = 0; % = 0 if non ideal, = 1 if ideal; here, it is OFF.

%% Non Ideal - Calculated Values
R_load = R_load_heavy;
D_non_ideal = (V_D_on + V_load_avg + (R_D_on * V_load_avg / R_load));
D_non_ideal = D_non_ideal / (V_in - V_T_on + V_D_on + (R_D_on * V_load_avg / R_load) - (R_T_on * V_load_avg / R_load));

i_L1_NI = (V_load_avg / R_load) - (((1 - D_non_ideal) * T_sw * V_load_avg) / (2 * L));
i_L2_NI = (V_load_avg / R_load) + (((1 - D_non_ideal) * T_sw * V_load_avg) / (2 * L));

%% Buck Intialization - Heavy Load
% Initializing Values
k = 1;
t = 0; 
dt = 1e-7;

tend = 100 * T_sw;

% Zero Vectors (used in buck)
t_vec = [0];  
switch_state = [0];

V_L_vec = [0];
i_L_vec = [0];

V_C_vec = [0]; 
i_C_vec = [0];

V_load_vec = [0]; 
i_load_vec = [0]; 

V_switch1 = [0]; 
i_switch1 = [0]; 

V_switch2 = [0]; 
i_switch2 = [0]; 

P_switch1 = [0]; % Power loss across the transistor; new for non-ideal calculations.
P_switch2 = [0]; % Power loss across the diode; new for non-ideal calculations.

%% Running Buck - Using R_load_heavy
R_load = R_load_heavy;
disp('Running buck for heavy load and non ideal conditions.');
buck

%% Post-processing Calculations (computing avg values, efficiency, etc)
disp("---------------------")
disp("Non Ideal Averages:")

V_load_avg_func_NI = aver(V_load_vec, T_sw, dt);
disp("  V_load Average: " + V_load_avg_func_NI);

i_L_func_NI = aver(i_L_vec, T_sw, dt);
disp("  i_L Average: " + i_L_func_NI);

V_sw1_func_NI = aver(V_switch1, T_sw, dt);
disp("  V_sw1 Average: " + V_sw1_func_NI);

i_sw1_func_NI = aver(i_switch1, T_sw, dt);
disp("  i_sw1 Average: " + i_sw1_func_NI);

V_sw2_func_NI = aver(V_switch2, T_sw, dt);
disp("  V_sw2 Average: " + V_sw2_func_NI);

i_sw2_func_NI = aver(i_switch2, T_sw, dt);
disp("  i_sw2 Average: " + i_sw2_func_NI);

P_sw1_func_NI = aver(P_switch1, T_sw, dt);
disp("  P_sw1 Average: " + V_sw2_func_NI);

P_sw2_func_NI = aver(P_switch2, T_sw, dt);
disp("  P_sw2 Average: " + i_sw2_func_NI);

P_out_NI = (V_load_avg_func_NI^2) / R_load;
P_in_NI = V_in * i_sw1_func_NI;
eff_NI = P_out_NI / P_in_NI;

disp("Efficiency for Non-Ideal: " + (eff_NI * 100) + "%.");
disp("Transistor Power Loss: " + P_sw1_func_NI);
disp("Diode Power Loss: " + P_sw2_func_NI);
disp("---------------------")
%% Plotting - Non Ideal Heavy Load - Transient
% Plots for the transient to steady state
figure;
sgtitle("Non Ideal Heavy Load: Voltage and Current Plots for All Components at Transient State");
% Plots for the Load
subplot(3,2,1);
plot(t_vec, V_load_vec);
title('Load Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(3,2,2);
plot(t_vec, i_L_vec);
title('Inductor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 1
subplot(3,2,3);
plot(t_vec, V_switch1);
title('Transistor Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(3,2,4);
plot(t_vec, i_switch1);
title('Transistor Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 2 
subplot(3,2,5);
plot(t_vec, V_switch2);
title('Diode Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(3,2,6);
plot(t_vec, i_switch2);
title('Diode Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

figure;
subplot(1,2,1);
plot(t_vec, P_switch1);
title('Diode Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Power (W)');

subplot(1,2,2);
plot(t_vec, P_switch2);
title('Diode Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Power (W)');

%% Plotting - Heavy Load - Steady State
periods_to_plot = 2; 

points_per_period = round(T_sw / dt);  % Points per period
total_periods = floor(tend / T_sw);    % Total number of periods in the simulation

start_index = max(1, (total_periods - periods_to_plot) * points_per_period + 1);
end_index = min(length(t_vec), total_periods * points_per_period);

range_to_plot = start_index:end_index;

%% Plot
figure;
sgtitle("Non Ideal Heavy Load: Voltage and Current Plots for All Components at Steady State");
% Plots for the Load
subplot(3,2,1);
plot(t_vec(range_to_plot), V_load_vec(range_to_plot));
title('Load Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(3,2,2);
plot(t_vec(range_to_plot), i_L_vec(range_to_plot));
title('Inductor Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 1
subplot(3,2,3);
plot(t_vec(range_to_plot), V_switch1(range_to_plot));
title('Transistor Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(3,2,4);
plot(t_vec(range_to_plot), i_switch1(range_to_plot));
title('Transistor Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

% Plots for Switch 2 
subplot(3,2,5);
plot(t_vec(range_to_plot), V_switch2(range_to_plot));
title('Diode Switching Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');

subplot(3,2,6);
plot(t_vec(range_to_plot), i_switch2(range_to_plot));
title('Diode Switching Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');