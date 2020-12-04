clear
%Numerically solve phagocytosis model for a cylindrical particle
%David M. Richards, 23/10/2017
%Define parameters
DT = 2.5e-5; %time step [s]
DX = 0.01; %lattice step [um]
T_MAX = 60; %simulation time [s]
L = 50; %membrane length [um]
PLOT_TIMES = [3 6]; %times at which to plot p [s]
D = 1; %diffusion constant [um^2/s]
pL = 500; %bound receptor density [/um]
p0 = 50; %initial receptor density [/um]
E = 15; %ligand-receptor binding energy [-]
B = 20; %bending modulus [-]
R = 2; %target radius [um]
%Declare variables
num_steps = round(T_MAX/DT); %number of time steps [-]
num_latt_pts = round(L/DX); %number of lattice points [-]
plot_times_steps = round(PLOT_TIMES/DT); %plot times in time steps [-]
%p = zeros(num_steps,num_latt_pts); %receptor density [/um]
a = zeros(num_steps,1); %cup size [um]
p_plus = fzero(@(p_plus)p_plus/pL-log(p_plus/pL)-E+2*B/(pL*R^2)-1, ...
    [1e-10 pL]); %calculate p at edge of cup [/um]
%Initial conditions
for i = 1:num_latt_pts
    p(1,i) = p0;
end
a(1) = 0;
p_prev = p(1,:);
sol = zeros(2,num_latt_pts);
%Main loop over time
for i = 2:num_steps
    %find cup size in lattice steps
    a_latt_pts = round( a(i-1)/DX + 1 );
    %set p to pL within cup
    for j = 1:a_latt_pts-1
        p_new(j) = pL;
    end
    %impose p_plus boundary condition
    p_new(a_latt_pts) = p_plus;
    %diffuse p outside cup
    for j = a_latt_pts+1:num_latt_pts-1
        p_new(j) = p_prev(j) + ...
            D*DT/DX^2 * (p_prev(j-1)-2*p_prev(j)+p_prev(j+1)) + ...
            D*DT/(DX*j*DX) * (p_prev(j+1)-p_prev(j));
    end
    %impose no flux condition at x=L
    p_new(num_latt_pts) = p_prev(num_latt_pts) + ...
        D*DT/DX^2 * (p_prev(num_latt_pts-1)-p_prev(num_latt_pts));
    %update cup size
    p_dash_plus = (p_new(a_latt_pts+1)-p_new(a_latt_pts)) / DX;
    a_dot = D*p_dash_plus / (pL-p_plus);
    a(i) = a(i-1) + a_dot*DT;
    if i == plot_times_steps(1)
        sol(1,:) = p_new;
    elseif i == plot_times_steps(2)
        sol(2,:) = p_new;
    end
    p_prev = p_new;
end
%% Analytical Solution
g = (p0-p_plus)/(pL-p_plus);
E1 = @(x) integral(@(u) exp(-u)./u,x,inf);
alpha = fzero(@(alpha)alpha^2*exp(alpha^2)*E1(alpha^2)-g, [0 1]);
A = (p0-p_plus)/E1(alpha^2);
time = [3,6];
%%Finding a(t)
a_any = 2*alpha*sqrt(D.*time);
%%Finding p(r,t)
a_sol = zeros(2,num_latt_pts);
for i = 1:2
    for j = 1:num_latt_pts
        if j*DX < a_any(i)
            a_sol(i,j) = pL;
        else
            a_sol(i,j) = p0 - A*E1((j*DX)^2/(4*D*time(i)));
        end
    end
end
%% Plot receptor density and cup size
figure()
subplot(1,2,1);
x = 0:DX:(num_latt_pts-1)*DX;
plot(x,sol(1,:),'y-', ...
    x,sol(2,:),'g-', ...
    x,a_sol(1,:),'k--', ...
    x,a_sol(2,:),'k--','LineWidth',3);
set(gca,'FontSize',24);
legend(['t=' num2str(PLOT_TIMES(1)) 's - Numerical'], ...
    ['t=' num2str(PLOT_TIMES(2)) 's - Numerical'], ...
    ['t=' num2str(PLOT_TIMES(1)) 's - Analytical'], ...
    ['t=' num2str(PLOT_TIMES(2)) 's - Analytical']);
xlim([0 8]);
ylim([0 500]);
xlabel('Distance from centre of cup (\mu m)');
ylabel('Receptor density (\mu m^{-1})');
title('Receptor density profile during engulfment');
subplot(1,2,2);
t = 0:DT:(num_steps-1)*DT;
plot(t,a,t,2*alpha*sqrt(D*t),'LineWidth',3);
set(gca,'FontSize',24);
legend('Numerical solution','Analytic solution','Location','NorthWest');
xlabel('Time (s)');
ylabel('Cup size (\mu m)');
title('Growth of phagocytic cup');