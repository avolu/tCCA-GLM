function [ REG, ADD] = perf_temp_emb_cca( X, AUX, param, flags )
% perf_temp_embed_cca uses multimodal data, here fNIRS, acceleration and auxiliary signals, to
% extract regressors for GLM-based noise reduction using temporally embedded
% CCA
% This function is an adaptation of the BLISSARD artifact rejection code by
% Alexander von Luhmann, written in 2018 at TUB
%
% INPUTS:
% X:    input nirs data with dimensions |time x channels|
% AUX:   input accel data with dimensions |time x channels|
%                    feeding in orthogonalized (PCA) data is advised.
% param.tau:  temporal embedding parameter (lag in samples)
% param.NumOfEmb: Number of temporally embedded copies
% param.ct:   correlation threshold. returns only those regressors that
%               have a canonical correlation greater than the threshold.
% flags.pcaf:  flags for performing pca of [AUX X] as preprocessing step.
%       default: [0 0]
%
% OUTPUTS:
% REG: Regressors found by temp. embedded CCA
% ADD.ccac: CCA correlation coefficients between projected fNIRS sources
%       and projected aux data
% ADD.aux_emb:   Temp. embedded accelerometer signals
% ADD.U,V:  sources in CCA space found by CCA


%% Perform PCA
if flags.pcaf(1)==true
    [coeff_afs aux_pca latent_afs] = pca(AUX);
    aux_sigs = aux_pca;
else
    aux_sigs =AUX;
end
if flags.pcaf(2)==true
    [coeff_x x_pca latent_x] = pca(X);
    nirs_sigs = x_pca;
else
    nirs_sigs = X;
end


%% Temporally embed auxiliary data
% Aux with shift left]
aux_emb = aux_sigs;
for i=1:param.NumOfEmb
    aux=circshift( aux_sigs, i*param.tau, 1);
    aux(1:2*i,:)=repmat(aux(2*i+1,:),2*i,1);
    aux_emb=[aux_emb aux];
    ADD.aux_emb=aux_emb;
end
% for i=1:param.NumOfEmb
%     aux=circshift( aux_sigs, -i*param.tau, 1);
%     aux(1:2*i,:)=repmat(aux(2*i+1,:),2*i,1);
%     aux_emb=[aux_emb aux];
%     ADD.aux_emb=aux_emb;
% end

%cut to same length of samples
s1=size(aux_emb,1);
s2=size(X,1);
if s1 ~=s2
    aux_emb=aux_emb(1:min([s1 s2]),:);
    X=X(:,1:min([s1 s2]));
end


%% Perform CCA with fNIRS signals and temporally embedded aux signals
[Au,Av,ccac,U,V] = canoncorr(X,aux_emb);
ADD.ccac =ccac;
ADD.U = U;
ADD.V = V;
ADD.Au = Au;
ADD.Av = Av;

% return auxiliary cca components that have correlation > ct
compindex=find(ccac>param.ct);
REG = V(:,compindex);

    
end

