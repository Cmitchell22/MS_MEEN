t0 = 5;
t_final = 40;


%% Time/Iteration Setup
time = zeros(1,t_final - t0);
startup_time = [t0 + 0.25, t0 + 0.5];
for t_setup = (1):t_final
    time(t_setup) = t0 + t_setup;
end
Total_Time = [startup_time, time];