% Use integration to calculate the math expection
% If you use this, please cite 
% Statistical analysis of multichannel FxLMS algorithm for narrowband active noise control

clc;
clearvars;
close all;
% dbstop if error
% dbstop if warning

path = 'simupath';
Fs = 8e3;
I = 3;
J = 2;
K = 2;
% figure(1)
% title('SecondaryPath')
% hold on
% legend

S = cell(J,K);
S_est = cell(J,K);
if strcmp(path,'simupath')
    mu_all = 1e-4:2e-3:0.3;
    M = 12;
    tmp = zeros(M,1);
    tmp(4) = 1;
    S{1,1} = tmp;
    if K>1
        S{2,2} = tmp;
        tmp = zeros(M,1);
        tmp(5) = 0.6;
        S{1,2} = tmp;
        S{2,1} = tmp;
    end
else        % real

end
for jj = 1:J
    for  kk = 1:K
        S_est{jj,kk} = S{jj,kk};
        %                 plot(S{jj,kk},'DisplayName',['S',num2str(jj),num2str(kk)],'LineWidth',1.5)
    end
end
M_est = length(S_est{1,1});

f = [100;200;300];
f = f(1:I);
omega = 2*pi*f/Fs;
% omega=[0.1*pi,0.2*pi,0.3*pi].';

rng(0)
a_p = rand(K,I)*2-1;b_p=rand(K,I)*2-1; %K*1
norm_primary = sqrt(a_p.^2+b_p.^2);
a_p = a_p./norm_primary;
b_p = b_p./norm_primary;
% a_p=[2.0;-0.5];b_p=[-1.0;1.0]; %K*1
% a_p=[2.0,1.0,0.5];b_p=[-1.0,-0.5,0.1];

syms mu [I,1]
% mu = [0.01,0.01,0.01,0.01]/10;
% mu = mu(1:I);
alpha = cell(I,1);beta = cell(I,1);
alpha_est = cell(I,1);beta_est = cell(I,1);
C_opt = cell(I,1);
C_opt_plus = cell(I,1);
W_opt = cell(I,1);
a_opt = cell(I,1);
b_opt = cell(I,1);
for ii = 1:I
    for jj = 1:J
        for  kk = 1:K
            alpha{ii}(jj,kk) = cos((0:M-1).*omega(ii))*S{jj,kk};
            beta{ii}(jj,kk)= sin((0:M-1).*omega(ii))*S{jj,kk};
            alpha_est{ii}(jj,kk) = cos((0:M_est-1).*omega(ii))*S_est{jj,kk};
            beta_est{ii}(jj,kk) = sin((0:M_est-1).*omega(ii))*S_est{jj,kk};
        end
    end
    C_opt{ii} = [alpha{ii}.',-beta{ii}.'; beta{ii}.',alpha{ii}.'];
    if J<K
        C_opt_plus{ii} = inv( C_opt{ii}.'*C_opt{ii} )*C_opt{ii}.';
    elseif J==K
        C_opt_plus{ii} = inv( C_opt{ii});
    elseif J>K
        C_opt_plus{ii} = C_opt{ii}.'*inv( C_opt{ii}*C_opt{ii}.' );
    end
    W_opt{ii} = C_opt_plus{ii}*[a_p(:,ii);b_p(:,ii)];
    a_opt{ii} = W_opt{ii}(1:J);
    b_opt{ii} = W_opt{ii}(J+1:2*J);
end


syms x y

x_f_a = @(ii,x)  alpha{ii}.*cos(x) + beta{ii}.*sin(x);
x_f_b = @(ii,x)  -beta{ii}.*cos(x) + alpha{ii}.*sin(x);
x_f_a_est = @(ii,x) alpha_est{ii}.*cos(x) + beta_est{ii}.*sin(x);
x_f_b_est = @(ii,x) -beta_est{ii}.*cos(x) + alpha_est{ii}.*sin(x);
x_f = @(ii,x) [x_f_a(ii,x);x_f_b(ii,x)];
x_f_est = @(ii,x) [x_f_a_est(ii,x);x_f_b_est(ii,x)];

P_x_y = @(ii,x,y)  x_f(ii,x) * x_f(ii,y).';
P_xest_yest = @(ii,x,y)  x_f_est(ii,x) * x_f_est(ii,y).';
P_xest_y = @(ii,x,y)  x_f_est(ii,x) * x_f(ii,y).';
P_x_yest = @(ii,x,y)  x_f(ii,x) * x_f_est(ii,y).';

E_P_ihat_i = cell(I,1);
m_D_core = cell(I,1); % Mean sense update matrix D

%%
disp(['Mean sense D'])
% mean sense
for ii = 1:I
    int_symbol = P_xest_y(ii,x,x);
    E_once = zeros(size(int_symbol));
    for jj = 1:size(int_symbol,1)
        parfor kk = 1:size(int_symbol,2)
            Fint = int(int_symbol(jj,kk),x,[0 2*pi])/(2*pi);
            E_once(jj,kk) = Fint; % get the data
        end
    end
    m_D_core{ii} = eye(2*J) - mu(ii) *E_once;
end

%%
% mean square sense
% Mean square sense update matrix F1
disp(['Mean square sense F1'])
clear E_once
syms E_once
ms_onefre = cell(1,I);
for ii = 1:I
    int_symbol = kron((eye(2*J)-mu(ii).*P_xest_y(ii,x,x)),(eye(2*J)-mu(ii).*P_xest_y(ii,x,x)));
    for jj = 1:size(int_symbol,1)
        disp(['All is ',num2str(I*size(int_symbol,1)),', now is ',num2str((ii-1)*size(int_symbol,1)+jj )])
        for kk = 1:size(int_symbol,2)
            Fint = int(int_symbol(jj,kk),x,[0 2*pi])/(2*pi);
            E_once(jj,kk) = Fint ; % get the data
        end
    end
    ms_onefre{ii} = E_once;
end
ms_F1_core = ms_onefre;

% Mean square sense update matrix F2
disp(['Mean square sense F2'])
clear E_once
syms E_once
ms_onefre = cell(1,I);
for ii = 1:I
    int_symbol = kron(P_xest_y(ii,x,y),P_xest_y(ii,x,y));
    for jj = 1:size(int_symbol,1)
        parfor kk = 1:size(int_symbol,2)
            Fint = int(int_symbol(jj,kk),x,[0 2*pi])/(2*pi);
            Fint = int(Fint,y,[0 2*pi])/(2*pi);
            E_once(jj,kk) = Fint; % get the data
        end
    end
    ms_onefre{ii} = E_once;
end
ms_F2_core = ms_onefre;

% Mean square sense update matrix F3
disp(['Mean square sense F3'])
clear E_once
syms E_once
ms_onefre = cell(I,1);
for ii = 1:I
    int_symbol = P_xest_yest(ii,x,x);
    for jj = 1:size(int_symbol,1)
        parfor kk = 1:size(int_symbol,2)
            Fint = int(int_symbol(jj,kk),x,[0 2*pi])/(2*pi);
            E_once(jj,kk) = Fint; % get the data
        end
    end
    ms_onefre{ii} = E_once(:);
end
ms_F3_core = ms_onefre;

% Mean square sense update matrix G
disp(['Mean square sense G'])
clear E_once
syms E_once
ms_onefre = cell(1,I);
for ii = 1:I
    int_symbol = kron(x_f(ii,x).', x_f(ii,x).');
    for jj = 1:size(int_symbol,1)
        parfor kk = 1:size(int_symbol,2)
            Fint = int(int_symbol(jj,kk),x,[0 2*pi])/(2*pi);
            E_once(jj,kk) = Fint; % get the data
        end
    end
    ms_onefre{ii} = E_once;
end
ms_G_core = ms_onefre;

% save(['ms_trans_matrix_I',num2str(I),'J',num2str(J),'K',num2str(K),path,'.mat'], 'm_D_core', 'ms_F1_core', 'ms_F2_core', 'ms_F3_core','ms_G_core')

