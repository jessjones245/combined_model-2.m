clear; clc; close all;

% Parameters from Hoore et al. (2020)
kp = 21;
kn = 0.22;
M = 6.7e7;
ds = 5.3e-9;
es = 1.9;
df = ds/30;

lambda_s = ds*M + es;
lambda_f = df*M;

% Sleep-debt and coupling parameters
a = 0.25;
tau_LA = 4.06;
Rb = 0.0012;

sleep_hours = [4 6 8];
gamma_values = [0 0.1 0.3 0.5];

% Simulation settings
years = 50;
total_days = years*365;
dt = 1/24;              % One-hour time step
steps_per_day = 24;

colors = lines(3);

figure('Position',[100 100 1600 450]);

% Repeat simulation for each gamma
for g = 1:length(gamma_values)

    gamma = gamma_values(g);

    subplot(1,4,g);
    hold on;

    % Repeat for each sleep duration
    for j = 1:length(sleep_hours)

        h = sleep_hours(j);
        awake_hours = 24-h;

        % Initial conditions
        S = Rb/lambda_s;
        F = 0;
        D = 0;

        F_results = zeros(total_days+1,1);

        % Simulate the coupled model
        for k = 1:(total_days*steps_per_day)

            hour = mod(k-1,24);

            % Awake or asleep
            if hour < awake_hours
                w = 1;
            else
                w = 0;
            end

            % A-beta production with sleep-debt coupling
            R = Rb*(1+a*(2*w-1))*(1+gamma*D);

            % Three coupled differential equations
            dS = R-kp*S*F-kn*S^2-lambda_s*S;
            dF = kp*S*F+kn*S^2-lambda_f*F;
            dD = (3*w-2-D)/tau_LA;

            % Forward Euler numerical method
            S = S+dt*dS;
            F = F+dt*dF;
            D = D+dt*dD;

            % Store fibrillar A-beta concentration daily
            if mod(k,steps_per_day) == 0
                day = k/steps_per_day;
                F_results(day+1) = F;
            end
        end

        time_years = (0:total_days)/365;

        plot(time_years,F_results, 'Color',colors(j,:), ...
            'LineWidth',2, ...
            'DisplayName',sprintf('%d h sleep',h));

        % Print final concentration
        fprintf('gamma = %.1f, sleep = %d h, F = %.6f uM\n', ...
            gamma,h,F);

    end

    title(sprintf('gamma = %.1f',gamma));
    xlabel('Time (years)');
    ylabel('Fibrillar A-beta (uM)');

    legend('Location','northwest');
    grid on;
    xlim([0 years]);
    ylim([0 0.05]);

    hold off;
end

sgtitle('Effect of Sleep Duration and Sleep-Debt Coupling');

% Save image
exportgraphics(gcf,'combined_model_gamma_05.png','Resolution',300);