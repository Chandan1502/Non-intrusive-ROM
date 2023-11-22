function [t,stn,nstn] = excitation_simulation(omega_min,omega_max,t_final,dt,spectrum)
t=0:dt:t_final;
% Power spectral density

S_0=0.5;
switch spectrum
    case 'KT'
        % Kanai-Tajimi
        omega_g=4; zeta_g=0.5;
        psd_my=@(omega) 2*S_0.*((1+4*(zeta_g^2).*(omega./omega_g).^2)./(((1-((omega./omega_g).^2)).^2)+(4*(zeta_g^2)).*((omega./omega_g).^2)));

    case 'CP'
        % Clough-Penzien
        omega_g=4.2; zeta_g=0.1; omega_f=2.3; zeta_f=0.1; % Parameters of CP spectrum
        psd_my=@(omega) 2*S_0*(((omega_g.^4)+(4*(zeta_g.^2)*(omega_g.^2)*(omega.^2)))./((((omega_g.^2)-(omega.^2)).^2)+(4.*(zeta_g.^2).*(omega_g.^2).*(omega.^2)))).*((omega.^4)./((((omega_f.^2)-(omega.^2)).^2)+(4.*(zeta_f.^2).*(omega_f.^2).*(omega.^2))));

    case 'WN'
        % Band limited white noise
        psd_my=@(omega) 2*S_0+(0*omega);

end
% PSD simulation

Num_PSDslice=3000;
Num_PSDpoints=Num_PSDslice+1;

omega_point=linspace(omega_min,omega_max,Num_PSDpoints);

wk1=omega_point(2:end);     % dummy variables
wk2=omega_point(1:end-1);   % dummy variables
w_k=(wk1+wk2)./2;

wk3=wk1-wk2;
A_k=sqrt(2*wk3.*psd_my(w_k));

phi_k=2*pi*rand(Num_PSDslice,1);    % Generation of uniformly distributed random numbers
cst=A_k'.*cos((w_k'+0.001).*t-phi_k);

stn=sum(cst,1); % Stationary excitation
% Genertation of nonstationary process

a1=1;
a2=0.25;
a_t=a1*t.*exp(-a2*t); % Modulating function
nstn=a_t.*stn;   % Nonstationary process

% plot(t,f_t,'-',t,f_tn,'--')
% plot(t,a_t)
end