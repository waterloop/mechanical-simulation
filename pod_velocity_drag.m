%% Aerodynamics Model

%Drag Force
rho = 0.0097 %density [kg/m^3]
Cd = 0.193 %drag coefficient
A = 21.927 %area [m^2]

%Rolling Force
C0 = 0.05 %rolling coefficient 1
m = 81.7; %mass [kg]
g = 9.81; %gravity [m/s^2]

sim('pod_velocity_drag_SIM',100)

figure()
plot(ans.time,ans.acceleration)
grid
title('Simulated Acceleration')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')

figure()
plot(ans.time,ans.velocity)
grid
title('Simulated Velocity')
ylabel('Velocity (m/s)')
xlabel('Time (s)')

%{
figure()
plot(ans.time,ans.displacement)
grid
title('Simulated Displacement')
xlabel('Time (s)')
ylabel('Displacement (m)')
%}
