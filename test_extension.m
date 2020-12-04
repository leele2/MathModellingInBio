clear
%Numerically solve phagocytosis model for a cylindrical particle
%David M. Richards, 23/10/2017
%Define parameters
DT = 4e-6; %time step [s]
DX = 0.01; %lattice step [um]
T_MAX = 60; %simulation time [s]
L = 50; %membrane length [um]
PLOT_TIMES = [3 6]; %times at which to plot p [s]
D0 = 1.0; %diffusion constant [um^2/s]
D1 = 10;
lambda = 50;
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
            (D0 + D1*exp(-j*DX/lambda))...
            *DT/DX^2 * (p_prev(j-1)-2*p_prev(j)+p_prev(j+1)) + ...
            (D0 + D1*exp(-j*DX/lambda))...
            *DT/(DX*j*DX) * (p_prev(j+1)-p_prev(j));
    end
    %impose no flux condition at x=L
    p_new(num_latt_pts) = p_prev(num_latt_pts) + ...
        (D0 + D1*exp(-j*DX/lambda))...
        *DT/DX^2 * (p_prev(num_latt_pts-1)-p_prev(num_latt_pts));
    %update cup size
    p_dash_plus = (p_new(a_latt_pts+1)-p_new(a_latt_pts)) / DX;
    a_dot = (D0 + D1*exp((-j+1)*DX/lambda))...
        *p_dash_plus / (pL-p_plus);
    a(i) = a(i-1) + a_dot*DT;
    if i == plot_times_steps(1)
        sol(1,:) = p_new;
    elseif i == plot_times_steps(2)
        sol(2,:) = p_new;
    end
    p_prev = p_new;
end
%% Plot receptor density and cup size
figure()
subplot(1,2,1);
x = 0:DX:(num_latt_pts-1)*DX;
plot(x,sol(1,:),'y-', x,sol(2,:),'g-','LineWidth',3)
set(gca,'FontSize',24);
legend(['t=' num2str(PLOT_TIMES(1)) 's'], ...
    ['t=' num2str(PLOT_TIMES(2)) 's'])
xlim([0 8]);
ylim([0 500]);
xlabel('Distance from centre of cup (\mu m)');
ylabel('Receptor density (\mu m^{-1})');
title('Receptor density profile during engulfment');
subplot(1,2,2);
t = 0:DT:(num_steps-1)*DT;
plot(t,a,'LineWidth',3)
set(gca,'FontSize',24);
xlim([0 60]);
ylim([0 3]);
xlabel('Time (s)');
ylabel('Cup size (\mu m)');
title('Growth of phagocytic cup');
sgtitle(['D_0 = ',num2str(D0),', D_1 = ',num2str(D1)...
    ,' and \lambda = ',num2str(lambda)],'FontSize',30)
