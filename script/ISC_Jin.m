%% Code info
% Week 6 Lab
% Author: JIN KE
% Motivation and Cognition Neuroscience Lab
% MA Program in Social Sciences (MAPSS)
% University of Chicago
% Email: jinke@uchicago.edu 
% Mobile: (312)-998-5025

%% Generate toy fMRI and behavioral data
n = 20;
r = 5;
t = 50;
s = 2;
data = rand([n,r,t,s]);
behavior = rand([n,n]);

save('input.mat','data')
%% Calculate three types of ISC: temporal, dynamic, and spatial

%% Temporal ISC

% pairwise ISC
pairwise_temporal_ISC = zeros([n,n,r]);
mean_series_session = mean(data,4);

for i=1:n
    for j=1:n
        for k=1:r
            
            pairwise_temporal_ISC(i,j,k)= corr(squeeze(mean_series_session(i,k,:)),squeeze(mean_series_session(j,k,:)));
        
        end
    end
end

% leave one out ISC

loo_temporal_ISC =zeros([n,r]);
for i=1:n
    for k=1:r

        mean_series_used  = mean_series_session;
        mean_series_used(i,:,:)=[]; 
        
        loo_temporal_ISC(i,k)= corr(squeeze(mean_series_session(i,k,:)),squeeze(mean(mean_series_used(:,k,:))));
    
    end
end

%% (2) Dynamic ISC

winsize = 10;
step = 10;

loo_dynamic_ISC = zeros([n,r,(t-winsize)/step+1]);

for i=1:n
    for k=1:r
        for win = 1:(t-winsize)/step+1

        mean_series_used = mean_series_session;
        mean_series_used(i,:,:)=[]; 

        loo_dynamic_ISC(i,k,win)= corr(squeeze(mean_series_session(i,k,(win-1)*winsize+1:win*winsize)),squeeze(mean(mean_series_used(:,k,(win-1)*winsize+1:win*winsize))));
        
        end
    end
end

figure
plot(10:10:50,squeeze(loo_dynamic_ISC(1,1,:)))
xlabel('Time(s)')
ylabel('Dynamic ISC')

%% (3) Spatial ISC

loo_spatial_ISC = zeros([n,1]);

for i=1:n

        mean_series_used  = mean_series_session;
        mean_series_used(i,:,:)=[];

        loo_spatial_ISC(i)= corr(mean(mean_series_session(i,:,:),3)',mean(mean(mean_series_used(:,:,:),3))') ;

end

%% Intra-subject correlation

intrasubject_temporal_ISC = zeros([n,r]);
for i=1:n
    for j=1:r

        intrasubject_temporal_ISC(i,j)= corr(squeeze(data(i,j,:,1)), squeeze(data(i,j,:,2))) ;
    
    end
end
%% Inter-subject functional connectivity
% leave one out ISC

loo_ISFC =zeros([n,r,r]);

for i=1:n
    for j=1:r
        for k=1:r

            mean_series_used  = mean_series_session;
            mean_series_used(i,:,:)=[];

            loo_ISFC(i,j,k)= corr(squeeze(mean_series_session(i,j,:)),squeeze(mean(mean_series_used(:,k,:))));
        
        end
    end
end
%% Inter-subject similarity to behavior

behavioral_similarity = zeros([n,n]);

for i=1:n
    for j=1:n
        
        behavioral_similarity(i,j) = sum(behavior(i,:)-behavior(j,:).^2);
    
    end
end

roi_similarity = zeros(r,1);
pvalue = zeros(r,1);

for i=1:r

    [roi_similarity(i),pvalue(i)]=corr(reshape_matrix(behavioral_similarity),reshape_matrix(squeeze(pairwise_temporal_ISC(:,:,i))));

end

save('output.mat','pairwise_temporal_ISC','pairwise_temporal_ISC','loo_dynamic_ISC','loo_spatial_ISC','intrasubject_temporal_ISC','loo_ISFC','behavioral_similarity', 'roi_similarity','pvalue')

%% Answer to the final question:
% The p values I generated are 0.060, 0.739, 0.665, 0.231, 0.054, which
% suggests that ISC is not significantly driven by behavioral similarity.

% I ran the code for the second time, and this time the p values it
% generates are 0.823, 0.755, 0.676, 0.644, 0.017

% Though one of them is lower than .05 which means significant given 0.05
% criterion, the other four is still not significant, so I feel comfortable
% concluding that ISC is not significantly driven by behavioral similarity.
%% define the reshape_matrix function

function [output_vector] = reshape_matrix(square_matrix)
    shape_square = ones(size(square_matrix,1));
    lower_half_mask = tril(shape_square,-1);
    square_vector = square_matrix(:);
    output_vector  = square_vector(lower_half_mask(:)>0);
end