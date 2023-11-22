clear variables
clc
% all units are in SI units
E=210*10^9;             % E= modulus of elasticity

I=1.2*10^(-6);
ro=7800;                % ro= density of beam
Ar=9*10^(-2);           % Ar= area of c/s
m=ro*Ar;

L=10;                    % L= length of beam in m

q=00*10^3;              % force in N/m

% Discretization of beam and assembly of beam elements
Nel=64;                 % Nel= Number of elements in beam
[Kglobal,Mglobal,~] = Beam_element_assembly(E,I,Ar,ro,L,Nel,q);
K=Kglobal;
M=Mglobal;

% Linear analysis
% The distributed spring is assumed to be linear
% Modification of stiffness matrix for distributed springs

% propoerties of soil
soilG=4e6; % Shear modulus of soil in Pa
soilNu=0.35; % Poisson's ratio
k_spring=2.1*soilG/(2-soilNu)*L;    % k_spring= total stiffness of sprin

k_v=k_spring/L;
%
for i=1:length(K)/2
    if i==1 || i==length(K)/2
        K(2*i-1,2*i-1)=K(2*i-1,2*i-1)+0.5*k_v*L/Nel;
    else
        K(2*i-1,2*i-1)=K(2*i-1,2*i-1)+k_v*L/Nel;
    end
end


% Eigenvalue Analysis
[modeShape,d]=eig(K,M);
funda_freq=sqrt(d(1,1)); % fundamental frequency
natural_freq=(sqrt(diag(d)));       % Numerical Result

wmax=max(max(sqrt(d)));
tcr=2/wmax;
c=zeros(length(M));
eta=0.02;
for i=1:length(M)
    c(i,i)=2*eta*sqrt(d(i,i));
end
C=M*modeShape*c*modeShape'*M;         % C=damping matrix
i_mode=1;

% Initial conditions
[Ndof,~]=size(K);
Disp_initial_HDM=zeros(Ndof,1);  % Intial displacement
% Disp_initial_HDM=modeShape(:,1);
Vel_initial_HDM=zeros(Ndof,1);   % Initial velocity

%% Generation of snapshot matrix / Training
amp_vec=2e6*modeShape(:,i_mode);
NumRealization=6000;

omega_min=5; % Hyper parameters of spectrum (should be changed)
omega_max=50;

dt=2e-4;
N=50001;
t_final=(N-1)*dt;

t=0:dt:t_final;
Nsnap=1000; % # of snapshots
Snap_interval=floor(N/Nsnap);
parfor i=1:NumRealization
    [~,~,f_t] = excitation_simulation(omega_min,omega_max,t_final,dt,'CP');   
    F=zeros(Ndof,N);
    F(:,:)=amp_vec.*f_t;

    [ ~,Disp_HDM,~,~] = Newmark(M,C,K,F,dt,N,Disp_initial_HDM,Vel_initial_HDM );
    peak_disp_BLWF(i)=max(max(abs(Disp_HDM)));
    S_big(:,:,i)=Disp_HDM(:,1:Snap_interval:Nsnap*Snap_interval);
    Force_input(:,i)=f_t(1:Snap_interval:Nsnap*Snap_interval);
end


%% Computation of POD bases
S_SVD=reshape(S_big,Ndof,[]);
[si,sigma,phi]=svd(S_SVD,"econ");
SingularValues=diag(sigma);
PODBases=si;