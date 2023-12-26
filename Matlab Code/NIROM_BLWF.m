clear all
clc
load("Results_BLWF.mat.mat")
%% Training of each component using NN
clear f_t
time_idx_start=1;
time_idx_end=501;

% Numtraining=[100 500 1000 1500 2000 2500 3000 3500 4000];
Numtraining=[4000];
time_start_training=clock;
for ii=1:length(Numtraining)
    % PCA of input force
    [u, s, v]=svd(Force_input(time_idx_start:time_idx_end,1:Numtraining(ii)));
    [~,~,~,Nf] = truncated_matrix(u,s,v,0.9999);
    basis_force=u(:,1:Nf);
    F_input=basis_force'*Force_input(time_idx_start:time_idx_end,1:Numtraining(ii));
    clear u s v
    F_input_normalize= rescale(F_input,0,1);

    minF_input=min(min(F_input));
    maxF_input=max(max(F_input));
    diff_min_max=maxF_input-minF_input;

    Xr_array_dummy=Xr_array(:,time_idx_start:time_idx_end,1:Numtraining(ii));
    Xr=reshape(Xr_array_dummy,NumBases,[]); % Converts(reshape) multidimensional array into a matrix

    for nb=1:NumBases
        each_Xr=Xr(nb,:);
        each_Xr=reshape(each_Xr,[time_idx_end-time_idx_start+1,Numtraining(ii)]);    % reshape a vector into a matrix of size(Nsnap,length(para_space))

        % PCA of output signal (training)
        [si,sigma,phi]=svd(each_Xr,"econ");

        % selection of number of singular values that contains the given accuracy
        tol=0.99999;
        [siReduced,sigmaReduced,phiReduced,Nr] = truncated_matrix(si,sigma,phi,tol);
        U_output=siReduced'*each_Xr;

        minU_output(nb)=min(min(U_output));
        maxU_output(nb)=max(max(U_output));
        diffminmaxU_output(nb)=maxU_output(nb)-minU_output(nb);
        U_output_normalize = rescale(U_output,0,1);
        % Traininig of neural network
        hiddenlayerSize=[60 45 25];
        net=fitnet(hiddenlayerSize);

        net.divideParam.trainRatio=80/100;
        net.divideParam.valRatio=10/100;
        net.divideParam.testRatio=10/100;
        net.trainFcn = 'trainscg'; %'trainscg'; %'traingdx'%,

        net.trainParam.max_fail = 10;
        net.trainParam.epochs=1500;
        net.trainParam.showWindow = 0;

        [trainednet,tr]=train(net,F_input_normalize,U_output_normalize);
        each_NN(nb).trained=trainednet; % stores NN architecture for each component
        each_NN(nb).tr=tr;
        each_NN(nb).siReduced=siReduced;

        clear siReduced each_Xr
    end
end
time_end_training=clock;
time_training=etime(time_end_training,time_start_training);
%% Testing of trained network for different realization for a given number of training points
disp("parfor started")
time_start=clock;
disp(time_start)

parfor i=1:2000
    [t,f_t,~] = excitation_simulation(omega_min,omega_max,t_final,dt,'CP');

    % HDM Solution for new force
    time_start_HDM=clock;
    F=zeros(Ndof,N);
    F(:,:)=amp_vec.*f_t;
    [ ~,Disp_HDM_test,~] = Newmark(M,C,K,F,dt,N,Disp_initial_HDM,Vel_initial_HDM );
    Disp_HDM_test_discrete=Disp_HDM_test(:,1:Snap_interval:Nsnap*Snap_interval);

    max_disp=max(max(abs(Disp_HDM_test_discrete(:,time_idx_start:time_idx_end))));
    peak_HDM(i)=max_disp;

    time_end_HDM=clock;

    time_HDM(i)=etime(time_end_HDM,time_start_HDM);


    f_t_new=f_t(1:Snap_interval:Nsnap*Snap_interval);
    f_t_new=basis_force'*f_t_new(time_idx_start:time_idx_end)';

    f_t_new_normalize= (1/diff_min_max)*(f_t_new-minF_input);
    time_start_ROM=clock;

    Xr_new=Xr_compute(each_NN,f_t_new_normalize,minU_output,diffminmaxU_output,NumBases);

    Disp_HDM_ROM=PODBases*Xr_new;
    time_end_ROM=clock;

    time_ROM(i)=etime(time_end_ROM,time_start_ROM);

    peak_ROM(i)=max(max(abs(Disp_HDM_ROM)));
    if rem(i,1001)==0
        display(i)
    end
end
time_end=clock;
time_whole=etime(time_end,time_start);

ksdensity(peak_ROM)
hold on
ksdensity(peak_HDM)
hold off
