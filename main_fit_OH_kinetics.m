%% ========================================================================
%  Quantitative interpretation of lateral adsorbate–intermediate
%  interactions via Butler–Volmer kinetics
%
%  Authors:  Yilin Zhao, Hengshuo Huang, Mingchuan Luo*
%  Affiliation:  School of Materials Science and Engineering,
%                Peking University, Beijing 100871, China
 
%  Corresponding:  m.luo@pku.edu.cn
%
%  This script fits the modified Butler–Volmer kinetic model 
%  to experimental ORR polarization data on Pt(111) via non-linear least squares.
%
%  !!Required data!!  origin Data-OH.xlsx (sheets: 'sita OH (E)', 'j vs E')
 
%% ========================================================================
 
 
 
 
clear; clc; close all;
 
%% Known Constants
jL = 6.10;           % Limiting diffusion current density (mA cm-2)
m = 1;               % Reaction order
E0 = 1.229;           % Equilibrium potential (V vs RHE)
b_star = -0.0519;     % Tafel slope (V/dec); 

dataE = xlsread('origin Data-OH.xlsx', 'sita OH (E)');
E_theta = dataE(:,1); % potential (V vs RHE)
theta_OH = dataE(:,2); % OH adsorption coverage (monolayer)
 
 
%% 1.3 Raw polarization data (sheet: 'j vs E')
% E_j: potential (V vs RHE); j: measured current density (mA cm-2)
dataE = xlsread('origin Data-OH.xlsx', 'j vs E');
j = dataE(:,1); % potential (V vs RHE)
E_j = dataE(:,2); % measured current density (mA cm-2)
 

% Plot interpolation results for validation
E_plot = linspace(min(E_theta), max(E_theta), 1000);
theta_plot = theta_interp(E_plot);
 
figure('Name','Spline interpolation of OH(E)','Position',[100,100,800,400]);
plot(E_theta, theta_OH, 'ro', 'MarkerSize',4, 'DisplayName','Experimental data');
hold on; grid on; box on;
plot(E_plot, theta_plot, 'b-', 'LineWidth',1.5, 'DisplayName','Cubic spline fit');
xlabel('E (V vs RHE)','FontSize',11);
ylabel('\theta_{OH} (monolayer)','FontSize',11);
legend('Location','southeast','FontSize',10);
title('OH Coverage vs Potential — Spline Interpolation','FontSize',12);
 
%% Kinetic Current j_K(E) via Koutecky–Levich 
% 3.1 Select valid data points 
% Criteria: (1) E within ?_OH data range (interpolation only); (2) retain ORR reduction current (j<0); (3) avoid near-zero denominators
 
% E_min = min(E_theta);
% E_max = max(E_theta);
 
E_min = 0.56;
E_max = 1.05;
 
% idx_valid = (E_j >= E_min) & (E_j <= E_max) & (j < 0);
idx_valid = (E_j >= E_min) & (E_j <= E_max);
 
% Extract valid data
E_valid = E_j(idx_valid);
j_valid = j(idx_valid);
 
 
%% Define the fitting model
jmodel = @(p, E) -(p(1) .* (1 - theta_interp(E)).^m ...
    .* exp( 2.303 .* (E - E0) ./ b_star ) ...
    .* exp( ( p(2) .* theta_interp(E) - log(1 + p(3) .* exp(p(2) .* theta_interp(E)))))).* jL./((p(1) .* (1 - theta_interp(E)).^m ...
    .* exp( 2.303 .* (E - E0) ./ b_star ) ...
    .* exp( ( p(2) .* theta_interp(E) - log(1 + p(3) .* exp(p(2) .* theta_interp(E)))))) + jL);
     
p0 = [10e-4, 0.05, 1e-4];  % Initial guess: [j0, eps_OH, K0CA]
lb = [1e-9, -10, 1e-20]; % Lower bounds: j0 > 0, K0CA > 0
ub = [1e-2, +10.5, 1e+5];      % Upper bounds
 
% Non-linear least-squares fitting 
options = optimoptions('lsqcurvefit', ...
    'Display','iter', ...
    'MaxIterations',1000, ...
    'MaxFunctionEvaluations',10000, ...
    'TolFun',1e-30, ...
    'TolX',1e-30);
 
[p_fit, resnorm, ~, exitflag, ~] = lsqcurvefit(jmodel, p0, E_valid, j_valid, lb, ub, options);
 
% Print fitting results
fprintf('\n Fitting Results \n');
fprintf('Convergence: exitflag = %d (1 = converged)\n', exitflag);
fprintf('Fitted j0 = %.4e mA cm-2 \n', p_fit(1));
fprintf('Fitted eps_OH = %.4f V\n', p_fit(2));
fprintf('Fitted K0CA = %.4e\n', p_fit(3));
fprintf('Residual sum of squares = %.4e\n', resnorm);
fprintf('======================================================\n');
 
% Fitted vs experimental data
j_fit = jmodel(p_fit, E_valid);
 
figure('Name','Kinetic model fit vs experiment','Position',[100,100,800,400]);
plot(E_j, j, 'bo', 'MarkerSize',4, 'DisplayName','Experimental j');
hold on; grid on; box on;
plot(E_valid, j_fit, 'r-', 'LineWidth',1.5, 'DisplayName','Model fit');
xlabel('E (V vs RHE)','FontSize',11);
ylabel('j (mA·cm^{-2})','FontSize',11);
legend('Location','northeast','FontSize',10);
title('ORR Kinetic Model — Fit vs Experiment','FontSize',12);
xlim([0.4 1])
 


