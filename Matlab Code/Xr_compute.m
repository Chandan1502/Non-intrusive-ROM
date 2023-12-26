function Xr_new = Xr_compute(each_NN,f_t_new_normalize,minU_output,diffminmaxU_output,NumBases)
for nb=1:NumBases
    trainednet=each_NN(nb).trained;
    y_new_normalize=trainednet(f_t_new_normalize);
    y_new=minU_output(nb)+diffminmaxU_output(nb)*y_new_normalize;
    siReduced=each_NN(nb).siReduced;
    y_new_final=siReduced*y_new;
    Xr_new(nb,:)=y_new_final';
end
end