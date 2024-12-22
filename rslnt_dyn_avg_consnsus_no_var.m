

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Resilient dynamic-average consensus or Distributed Average tracking%%%%%%%
%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
n = 5;
x_ini = [-1 2 3 2.2 1]' .* rand(n,1);
z_ini =  rand(n,1);
% x_ini = zeros(n,1);
% z_ini =  zeros(n,1);
d_ini = zeros(n,1);
d_p_ini = zeros(n,1);
chi_ini =[x_ini ; z_ini ; d_ini ; d_p_ini];
tspan = [1 8];
options = odeset('RelTol',1e-6);

[t chi] = ode45(@x_position, tspan, chi_ini, options);
phi_avg = zeros(length(t),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_avg = (1/5)* (chi(:,6)+chi(:,7)+chi(:,8)+chi(:,9)+chi(:,10));
phi = (-cos(t)) + (2*sin(t+pi/2)) - 3*cos(3*t)- 4*cos(4*t) + 5*sin(5*t);
% phi = [1; -2; 2; -1 ; 3];
% phi = 0;
phi_avg = (1/5)* phi;


%%%%%%%%%%%%%%%%%%%%%%% Another plot option %%%%%%%%%%%%%%%%%
% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 17;
opts.height     = 5;
opts.fontType   = 'Times';
opts.fontSize   = 9;

% create new figure
fig = figure; clf

% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;
axis tight
% remove unnecessary white space
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

% export to png
fig.PaperPositionMode   = 'auto';
%%%%%%%%%%%%%%%%%%%%%%Plot section %%%%%%%%%%%%%%%%
% fig1 = figure;
% set(fig1,'PaperPosition',[0 0 6.5 4]);
% set(fig1,'PaperSize',[6 3.8]);
plot(t, phi_avg,'-.');
%for static consensus
% plot(phi_avg);
hold on 
plot(t, chi(:,1:5));
% plot(t, z_avg,'linewidth',1.5);
leg1 = legend('$r_{avg}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',10);
xlabel('Time (Seconds)');
ylabel('Positions of agents');
% save2pdf('dycon_attack_d_high_gains',fig,600);

function chi_dot = x_position(t, chi)
 n = 5;
 d_o = 1.1 * ones(n,1);
L = circulant([1 -1 0 0 0],1);
% L=[1 -1 0 0 0;0 1 -1 0 0;0 0 2 -1 -1;0 0 -1 2 -1;-1 0 0 -1 2]; 
% L = [2 -1 -1 0 0;
%      0 1 -1 0 0;
%      0 0 1 -1 0;
%      0 0 0 1 -1;
%      -1 0 0 0 1];
D = eye(n,n);
A = D - L;


r = [-cos(t); 2*sin(t+pi/2); - 3*cos(3*t); - 4*cos(4*t) ; 5*sin(5*t)];
v_r = [sin(t); 2*cos(t) ; 9*sin(3*t) ; 16*sin(4*t) ; 25*cos(5*t)];

% r = [1; -2; 2; -1 ; 3];
% v_r = [0; 0 ; 0 ; 0 ; 0];

% r = [-0; 0; - 0; -0; 0];            % for zero-input system
% v_r = [0; 0 ; 0 ; 0 ; 0];           % for zero-input system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F1 = -eye(n,n);
F2 = -2 * eye(n,n);
Ba = [1 3 2 1 2;
      2 -1 3 2 1;
      -2 1 1 2 3;
      2 1 -4 3 2;
      1 3 2 1 -5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% For four/five agents alpha = 90, gamma = 5, beta = 500 works well for
%%%% undirected graphs
alpha = 90;
beta = 50;


k = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : n
    x_dot(i,1) = + (alpha * (r(i) - chi(i,1))) - (alpha*beta)*(L(i,:) * chi(1:5,:)) + ( beta * ((L(i,:) * chi(6:10,:))))+ v_r(i)+ (L(i,:) * chi(11:15,:));
    z_dot(i,1) =  (alpha * (r(i) - chi(k+i,1))) -  beta * (L(i,:) * chi(1:5,:))  -  (L(i,:) * (chi(6:10,:))) + v_r(i);
end
d_dot = F1 * chi(11:15,:) + Ba * d_o;
d_p_dot = F2 * chi(16:20,:) + Ba * d_o;
chi_dot = [x_dot ; z_dot; d_dot; d_p_dot];
end