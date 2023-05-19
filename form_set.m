% A function that creates a feature vector (standardized according to the
% VAURAS mean/std) for the ENES machine.
% Input:
% - features structure (output of a get_features_ENES() function)
% - the vector of means for standardization
% - the vector of stds for standardizadion
% Output:
% - vector of 41 standardized features in a correct order

function dataset=form_set(features,mean_stand,std_stand)



aEEG_mean_unilateral(:)=mean(features.aEEG.mean(5:6,:),1,'omitnan');
aEEG_iqr_unilateral(:)=mean(features.aEEG.iqr(5:6,:),1,'omitnan');
rEEG_mean_unilateral(:)=mean(features.rEEG.mean(5:6,:),1,'omitnan');
rEEG_iqr_unilateral(:)=mean(features.rEEG.iqr(5:6,:),1,'omitnan');
aEEG_mean_bilateral(:)=features.aEEG.mean(8,:);
aEEG_iqr_bilateral(:)=features.aEEG.iqr(8,:);
rEEG_mean_bilateral(:)=features.rEEG.mean(8,:);
rEEG_iqr_bilateral(:)=features.rEEG.iqr(8,:);

PSD_1_3_unilateral(:)=mean(features.PSD(5:6,1,:),1,'omitnan');
PSD_3_8_unilateral(:)=mean(features.PSD(5:6,2,:),1,'omitnan');
PSD_8_15_unilateral(:)=mean(features.PSD(5:6,3,:),1,'omitnan');
PSD_15_30_unilateral(:)=mean(features.PSD(5:6,4,:),1,'omitnan');
PSD_1_3_bilateral(:)=features.PSD(8,1,:);
PSD_3_8_bilateral(:)=features.PSD(8,2,:);
PSD_8_15_bilateral(:)=features.PSD(8,3,:);
PSD_15_30_bilateral(:)=features.PSD(8,4,:);

cPSD_1_3_unilateral(:)=squeeze(mean([features.cPSD(1,3,1,:) features.cPSD(2,4,1,:)],2,'omitnan'));
cPSD_3_8_unilateral(:)=squeeze(mean([features.cPSD(1,3,2,:) features.cPSD(2,4,2,:)],2,'omitnan'));
cPSD_8_15_unilateral(:)=squeeze(mean([features.cPSD(1,3,3,:) features.cPSD(2,4,3,:)],2,'omitnan'));
cPSD_15_30_unilateral(:)=squeeze(mean([features.cPSD(1,3,4,:) features.cPSD(2,4,4,:)],2,'omitnan'));
cPSD_1_3_bilateral(:)=features.cPSD(5,6,1,:);
cPSD_3_8_bilateral(:)=features.cPSD(5,6,2,:);
cPSD_8_15_bilateral(:)=features.cPSD(5,6,3,:);
cPSD_15_30_bilateral(:)=features.cPSD(5,6,4,:);

wPLI_1_3_unilateral(:)=squeeze(mean([features.wPLI(1,3,1,:) features.wPLI(2,4,1,:)],2,'omitnan'));
wPLI_3_8_unilateral(:)=squeeze(mean([features.wPLI(1,3,2,:) features.wPLI(2,4,2,:)],2,'omitnan'));
wPLI_8_13_unilateral(:)=squeeze(mean([features.wPLI(1,3,3,:) features.wPLI(2,4,3,:)],2,'omitnan'));
wPLI_13_30_unilateral(:)=squeeze(mean([features.wPLI(1,3,4,:) features.wPLI(2,4,4,:)],2,'omitnan'));
wPLI_1_3_bilateral(:)=features.wPLI(5,6,1,:);
wPLI_3_8_bilateral(:)=features.wPLI(5,6,2,:);
wPLI_8_13_bilateral(:)=features.wPLI(5,6,3,:);
wPLI_13_30_bilateral(:)=features.wPLI(5,6,4,:);
        
ASI_unilateral(:)=squeeze(mean([features.ASI(1,3,:) features.ASI(2,4,:)],2,'omitnan'));
ASI_bilateral(:)=features.ASI(5,6,:);
SC_global=features.SC;

NC_frontal_3_8(:)=mean(features.NC(1:2,1,:),1,'omitnan');
NC_frontal_8_15(:)=mean(features.NC(1:2,2,:),1,'omitnan');
NC_frontal_15_30(:)=mean(features.NC(1:2,3,:),1,'omitnan');
NC_parietal_3_8(:)=mean(features.NC(3:4,1,:),1,'omitnan');
NC_parietal_8_15(:)=mean(features.NC(3:4,2,:),1,'omitnan');
NC_parietal_15_30(:)=mean(features.NC(3:4,3,:),1,'omitnan');



    aEEG_mean_unilateral_med=median(aEEG_mean_unilateral,'omitnan');
    aEEG_iqr_unilateral_med=median(aEEG_iqr_unilateral,'omitnan');
    rEEG_mean_unilateral_med=median(rEEG_mean_unilateral,'omitnan');
    rEEG_iqr_unilateral_med=median(rEEG_iqr_unilateral,'omitnan');
    aEEG_mean_bilateral_med=median(aEEG_mean_bilateral,'omitnan');
    aEEG_iqr_bilateral_med=median(aEEG_iqr_bilateral,'omitnan');
    rEEG_mean_bilateral_med=median(rEEG_mean_bilateral,'omitnan');
    rEEG_iqr_bilateral_med=median(rEEG_iqr_bilateral,'omitnan');

    PSD_1_3_unilateral_med=median(PSD_1_3_unilateral,'omitnan');
    PSD_3_8_unilateral_med=median(PSD_3_8_unilateral,'omitnan');
    PSD_8_15_unilateral_med=median(PSD_8_15_unilateral,'omitnan');
    PSD_15_30_unilateral_med=median(PSD_15_30_unilateral,'omitnan');
    PSD_1_3_bilateral_med=median(PSD_1_3_bilateral,'omitnan');
    PSD_3_8_bilateral_med=median(PSD_3_8_bilateral,'omitnan');
    PSD_8_15_bilateral_med=median(PSD_8_15_bilateral,'omitnan');
    PSD_15_30_bilateral_med=median(PSD_15_30_bilateral,'omitnan');

    cPSD_1_3_unilateral_med=median(cPSD_1_3_unilateral,'omitnan');
    cPSD_3_8_unilateral_med=median(cPSD_3_8_unilateral,'omitnan');
    cPSD_8_15_unilateral_med=median(cPSD_8_15_unilateral,'omitnan');
    cPSD_15_30_unilateral_med=median(cPSD_15_30_unilateral,'omitnan');
    cPSD_1_3_bilateral_med=median(cPSD_1_3_bilateral,'omitnan');
    cPSD_3_8_bilateral_med=median(cPSD_3_8_bilateral,'omitnan');
    cPSD_8_15_bilateral_med=median(cPSD_8_15_bilateral,'omitnan');
    cPSD_15_30_bilateral_med=median(cPSD_15_30_bilateral,'omitnan');

    wPLI_1_3_unilateral_med=median(wPLI_1_3_unilateral,'omitnan');
    wPLI_3_8_unilateral_med=median(wPLI_3_8_unilateral,'omitnan');
    wPLI_8_13_unilateral_med=median(wPLI_8_13_unilateral,'omitnan');
    wPLI_13_30_unilateral_med=median(wPLI_13_30_unilateral,'omitnan');
    wPLI_1_3_bilateral_med=median(wPLI_1_3_bilateral,'omitnan');
    wPLI_3_8_bilateral_med=median(wPLI_3_8_bilateral,'omitnan');
    wPLI_8_13_bilateral_med=median(wPLI_8_13_bilateral,'omitnan');
    wPLI_13_30_bilateral_med=median(wPLI_13_30_bilateral,'omitnan');

    ASI_unilateral_med=median(ASI_unilateral,'omitnan');
    ASI_bilateral_med=median(ASI_bilateral,'omitnan');

    NC_frontal_3_8_med=median(NC_frontal_3_8,'omitnan');
    NC_frontal_8_15_med=median(NC_frontal_8_15,'omitnan');
    NC_frontal_15_30_med=median(NC_frontal_15_30,'omitnan');
    NC_parietal_3_8_med=median(NC_parietal_3_8,'omitnan');
    NC_parietal_8_15_med=median(NC_parietal_8_15,'omitnan');
    NC_parietal_15_30_med=median(NC_parietal_15_30,'omitnan');
    
    SC_med=median(SC_global,'omitnan');


features_vector=[aEEG_mean_unilateral_med aEEG_mean_bilateral_med aEEG_iqr_unilateral_med aEEG_iqr_bilateral_med...
    rEEG_mean_unilateral_med rEEG_mean_bilateral_med rEEG_iqr_unilateral_med rEEG_iqr_bilateral_med...
    PSD_1_3_unilateral_med PSD_1_3_bilateral_med PSD_3_8_unilateral_med  PSD_3_8_bilateral_med...
    PSD_8_15_unilateral_med PSD_8_15_bilateral_med PSD_15_30_unilateral_med PSD_15_30_bilateral_med...
    cPSD_1_3_unilateral_med cPSD_1_3_bilateral_med cPSD_3_8_unilateral_med cPSD_3_8_bilateral_med...
    cPSD_8_15_unilateral_med cPSD_8_15_bilateral_med cPSD_15_30_unilateral_med cPSD_15_30_bilateral_med...
    wPLI_1_3_unilateral_med wPLI_1_3_bilateral_med wPLI_3_8_unilateral_med  wPLI_3_8_bilateral_med ...
    wPLI_8_13_unilateral_med wPLI_8_13_bilateral_med wPLI_13_30_unilateral_med wPLI_13_30_bilateral_med...
    NC_frontal_3_8_med NC_parietal_3_8_med NC_frontal_8_15_med NC_parietal_8_15_med NC_frontal_15_30_med NC_parietal_15_30_med...
    SC_med ASI_unilateral_med ASI_bilateral_med
    ];

dataset=(features_vector-mean_stand)./std_stand;




end