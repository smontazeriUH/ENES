% The streamlined routine for ENES calculation.
% Inputs:
% aEEG, montage: F3, F4, P3, P4. The other channels are ignored. It is
% better if the sampling rate is 250 Hz or more. 
% delay between birth and beginning of the recording in minutes. So, if the
% EEG recording was started, for instance, 2 hours and 22 minutes after the birth, please
% enter 142.
%
%
% Alexey Pospelov
% University of Helsinki, Finland
% 15/05/2023


function score=ENES(name,delay)

if ~isnumeric(delay)
    error('The delay must be a nonnegative integer value and represent delay between the birth and the beginnig of the EEG recording that is used as an input')
end
if or(delay<0,round(delay)~=delay)
    error('The delay must be a nonnegative integer value and represent delay between the birth and the beginnig of the EEG recording that is used as an input')
end

sr=256; % The sampling rate at which the artefact detector works. Used as a standard sampling rate 
% EEG import
Filename=name; % full name of the file, later will be an input argument of the function
Delay=delay; % initial delay in minutes, later will be an input agrument of the function
[Ind_hdr, Ind_record] = edfread(Filename); % it is assumed that the input data is a single edf with correct montage. Only first four channels matter 

disp(['First four channels are: ' Ind_hdr.label{1} ', ' Ind_hdr.label{2} ', ' Ind_hdr.label{3} ', ' Ind_hdr.label{4} ', must be F3, F4, P3, P4']);


% test if 4-6h data and/or 22-24h data is present, report the presence

for ch=1:4
    dur(ch)=size(Ind_record(ch,:),2)/Ind_hdr.frequency(ch);% in seconds
end
if length(unique(dur))~=1
    error('Variable length of the EEG recordings in different channels, cannot proceed')
end


dur_min=unique(dur)/60; % duration of the EEG in minutes
index(1:24*60)=0;
index(Delay+1:min(Delay+dur_min,24*60))=1;

Presence_4_6=sum(index(4*60+1:6*60));
Presence_22_24=sum(index(22*60+1:24*60));

disp(['Available data: ' num2str(100*Presence_4_6/120) '% at 4-6h epoch, ' num2str(100*Presence_22_24/120) '% at 22-24h epoch'])



% resampling (goes first to handle inlikelly but possible case of different sampling rates in different channels)
for ch=1:4
    Data_temp(ch,:)=resample(Ind_record(ch,:),sr,Ind_hdr.frequency(ch));
end


% generation of differential channels 
Data(1:4,1:min(Delay*60*sr,24*60*60*sr))=NaN;

if size(Data_temp,2)<=24*60*60*sr-Delay*60*sr % the recording does not go up intil the end of the first day;
    Data=[Data Data_temp nan([4, 24*60*60*sr-Delay*60*sr-size(Data_temp,2)])];
    %disp('first')
else
    Data=[Data Data_temp(:,1:24*60*60*sr-Delay*60*sr)];
    %disp('second')
end
Data(5,:)=Data(1,:)-Data(3,:); %F3-P3 (Left)
Data(6,:)=Data(2,:)-Data(4,:); %F4-P4 (Right)
Data(7,:)=Data(1,:)-Data(2,:); %F3-F4 (Frontal)
Data(8,:)=Data(3,:)-Data(4,:); %P3-P4 (Posterior)



% filtration
Mask(1:1440*60*sr,1:size(Data,1))=0;
for chn=1:size(Data,1)
    tic
    NanMask=isnan(Data(chn,:));
    Data(chn,NanMask)=0;
    %     Data(:,chn)=eegfilt(Data(:,chn)',sr,2,0);
    %     Data(:,chn)=eegfilt(Data(:,chn)',sr,0,30);
    Data(chn,:)=my_bandpass(Data(chn,:)',[0.2 30],sr);
    Data(chn,NanMask)=NaN;
    Mask(:,chn)=NanMask';
    toc
end
for ch=1:8
    Data(ch,isnan(Data(ch,:)))=0;
end

% artefact detection
% The detection in 24h-long 8 channel recording goes in two 12h-long parts
% with a patch in the middle to remove the marginal effects. This is
% because on a regular office machine I got out of memory error. If you
% have the out of memory problem even with this design, you will have to
% increase the number of parts.
age_epoch(1,:)=240:359;
age_epoch(2,:)=1320:1439;



%patch=Data(1:size(Data,1),60*sr*sample2_duration+1-50*4*sr:60*sr*sample2_duration+50*4*sr)';




%[final_dec_patch, dec_prob_patch] = artefact_detection(patch', resnet);
    
    


threshold=0.75;
    
    
% for ch=1:8
%     for ep=1:length(dec_prob)
%         Netmask((ep-1)*sr*2+1:(ep-1)*sr*2+sr*4,ch)=dec_prob(ch,ep,1)<threshold;
%     end
% end
%     %combine the existance mask (where there is any data) and the netmask
%     %(where the network finds artifacts)
%  Mask=or(Mask,Netmask); % at this stage, "1" - artefact, "0" - clean data

% calculation of features, report if the data is insufficient (based on the artefacts detection)

EpochLenght=2*60; %in seconds
EpochOverlap=1*60; % in seconds

artifact_percentage_thresholds = ...
struct( 'ASI',      0.2, ...
        'wPLI',     0.1, ...
        'NC',       0.1, ...
        'aEEG',     0.5, ...
        'rEEG',     0.5, ...
        'PSD',      0.5, ...
        'cPSD',     0.5, ...
        'SC',       0.5);
% features = get_features_ENES(epoched_data,epoched_artefacts,artifact_percentage_thresholds, ...
% 'epoch_length_seconds', EpochLenght, 'sampling_rate', sr);



% calculation of ENES using trained SVM

mean_4_6h_coeff=[5.46845681770084	3.95017134602996	2.24431554885275	1.75524416112281	71.0588979095032	52.0442176516819	34.9656036184672	27.2998002824934	43.6674606647392	27.1048861408641	9.60928246397925	5.01145393661139	2.32246130958635	1.32168880144513	1.37121489320984	0.821774314804161	7.70806802815595	20.2284659863814	2.04186240778940	3.62043704290095	0.716908207950191	0.485542026194714	0.237807990650958	0.263274617547269	0.0501944429122198	0.0462536469071892	0.0170885658311066	0.0123685664239242	0.0120401648994853	0.00854029093916959	0.0151901674250208	0.00678566533623297	0.0932273939886459	0.0931785035831539	0.0950270284252628	0.0968471981991541	0.0973623496909681	0.0977486411279809	0.0340942024852859	5.02584949964948	4.03098428524936];
std_4_6h_coeff=[2.25537720525079	1.90651832972295	1.05323199361861	0.917722908456028	29.1652585323897	24.2305039239831	16.0264089907130	14.1037921852906	28.6791362091276	22.1543064090715	6.67690887634683	4.20634610791008	1.61752997844106	1.16970451223067	0.933726220508133	0.743440347356093	5.68070522838557	16.4746281780344	1.86926957708845	3.25571540935910	0.738717370148671	0.340894764477410	0.205371502882893	0.197350512723956	0.0656850618705244	0.0859342257692015	0.0273250660635892	0.0179924092214256	0.0192140618029446	0.0146657519462444	0.0303263495943488	0.0160648465332180	0.0163582544156614	0.0147787976157766	0.0170306521045883	0.0176310737683657	0.0178533654568003	0.0182376196575077	0.0398005769077415	1.42184089027502	1.18982340165704];

mean_22_24h_coeff=[5.40566805187259	4.13488520061947	2.07902001231696	1.72671626270776	66.4553300345256	50.6153639844883	29.7881297773875	23.4571066714771	34.1999459473346	21.1389095984343	8.86682130232586	4.54413470526470	2.28190173743729	1.46591760937336	1.19175845948406	0.795940580675981	6.66149425079273	17.9015657442184	1.93877102449189	3.65329174683221	0.791863106104130	0.497855336343879	0.231576288427550	0.270178723276787	0.0252022328483666	0.0149794920201029	0.00673098494982348	0.00378112123242685	0.00783691337106953	0.00702923035968165	0.00845074379889743	0.00380274608170098	0.100506745101342	0.0993581352445266	0.100897486682908	0.105828203412643	0.0982706513201405	0.105097994507746	0.0232906014244754	4.44008830462230	3.53129379115381];
std_22_24h_coeff=[1.85897703519250	1.77453267670948	0.767629265056801	0.838136192419123	24.6507210167533	20.4641729269902	10.5995623991299	9.39546321103515	25.0818961400791	16.3895338460142	6.30349912623249	3.65697038800793	1.40758995311461	1.36203607713251	0.783483555244832	0.667491912465856	6.19043849296858	15.8471813851839	1.68822050210795	3.37459246392845	0.693452903041026	0.361613995234480	0.186415497040232	0.184045565833715	0.0375381820231902	0.0203015999872208	0.0109191982026728	0.00633268306864877	0.0213732513753588	0.0195774976803044	0.0275469856105071	0.00902364613353722	0.0223147331430226	0.0172321826763984	0.0218697300010310	0.0218456199598436	0.0239504021863565	0.0238616543359289	0.0295696947917272	1.36487173249669	1.04387098363363];

load('SVM_trained_2epoch.mat')
load('resnet.mat');
if Presence_4_6/120>=0.5
    disp('Processing of the first period started')
    sample1=Data(1:size(Data,1),(age_epoch(1,1))*60*sr:60*sr*age_epoch(1,end))';
    [final_dec1, dec_prob1] = artefact_detection(sample1', resnet);
    for ch=1:8
        for ep=1:length(dec_prob1)
            Netmask((ep-1)*sr*2+1:(ep-1)*sr*2+sr*4,ch)=dec_prob1(ch,ep,1)<threshold;
        end
    end
    %combine the existance mask (where there is any data) and the netmask
    %(where the network finds artifacts)
    size(Mask((age_epoch(1,1))*60*sr+1:60*sr*age_epoch(1,end),:))
    size(Netmask)
    Mask1=or(Mask((age_epoch(1,1))*60*sr+1:60*sr*age_epoch(1,end),:),Netmask); % at this stage, "1" - artefact, "0" - clean data
    
    
    epoched_temp_4_6(1:EpochLenght*sr,1:8,1:118)=NaN;
    epoched_artefacts_temp_4_6(1:EpochLenght*sr,1:8,1:118)=NaN; % 0 - good data, 1 - missing data/artefact
    for ep=1:118
        epoched_temp_4_6(:,:,ep)=Data(:,(age_epoch(1,ep)-2)*EpochOverlap*sr+1:(age_epoch(1,ep)-2)*EpochOverlap*sr+EpochLenght*sr)';
        epoched_artefacts_temp_4_6(:,:,ep)=Mask1((ep-1)*EpochOverlap*sr+1:(ep-1)*EpochOverlap*sr+EpochLenght*sr,:);
    end
    epoched_data_4_6{1}=epoched_temp_4_6;
    epoched_artefacts_4_6{1}=logical(epoched_artefacts_temp_4_6);  

    features_4_6 = get_features_ENES(epoched_data_4_6,epoched_artefacts_4_6,artifact_percentage_thresholds, ...
    'epoch_length_seconds', EpochLenght, 'sampling_rate', sr);

    
    [~,~,~,posteriors_4_6]=predict(SVM_full,form_set(features_4_6,mean_4_6h_coeff,std_4_6h_coeff));
    ENES_4_6h=posteriors_4_6*[0 1 2 3]';
    score(1)=ENES_4_6h;
else
    score(1)=NaN;
    
end

if Presence_22_24/120>=0.5
    disp('Processing of the second period started')
    sample2=Data(1:size(Data,1),(age_epoch(2,1))*60*sr:60*sr*age_epoch(2,end))';
    [final_dec2, dec_prob2] = artefact_detection(sample2', resnet);
    for ch=1:8
        for ep=1:length(dec_prob2)
            Netmask((ep-1)*sr*2+1:(ep-1)*sr*2+sr*4,ch)=dec_prob2(ch,ep,1)<threshold;
        end
    end
    %combine the existance mask (where there is any data) and the netmask
    %(where the network finds artifacts)
    Mask2=or(Mask((age_epoch(2,1))*60*sr+1:60*sr*age_epoch(2,end),:),Netmask); % at this stage, "1" - artefact, "0" - clean data
    
    epoched_temp_22_24(1:EpochLenght*sr,1:8,1:118)=NaN;
    epoched_artefacts_temp_22_24(1:EpochLenght*sr,1:8,1:118)=NaN; % 0 - good data, 1 - missing data/artefact
    for ep=1:118
        epoched_temp_22_24(:,:,ep)=Data(:,(age_epoch(2,ep)-1)*EpochOverlap*sr+1:(age_epoch(2,ep)-1)*EpochOverlap*sr+EpochLenght*sr)';
        epoched_artefacts_temp_22_24(:,:,ep)=Mask2((ep-1)*EpochOverlap*sr+1:(ep-1)*EpochOverlap*sr+EpochLenght*sr,:);
    end
    epoched_data_22_24{1}=epoched_temp_22_24;
    epoched_artefacts_22_24{1}=logical(epoched_artefacts_temp_22_24);  

    features_22_24 = get_features_ENES(epoched_data_22_24,epoched_artefacts_22_24,artifact_percentage_thresholds, ...
    'epoch_length_seconds', EpochLenght, 'sampling_rate', sr);

    [~,~,~,posteriors_22_24]=predict(SVM_full,form_set(features_22_24,mean_22_24h_coeff,std_22_24h_coeff));
    ENES_22_24h=posteriors_22_24*[0 1 2 3]';
    score(2)=ENES_22_24h;
else
    score(2)=NaN;
end



clear Data Mask Mask1 Mask2 Netmask % Not sure if this is useful as the internal variables should be cleared after the function did its job anyway, but I got random "out of memory" errors when the function was run on a big set of files, so better safe than sorry. 

end 