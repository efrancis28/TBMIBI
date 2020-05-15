% % MIBIQuantifySpatialOrganizationMarkerExpression.m
% Author: Leeat Keren, adapted to TB cohort by Erin McCaffrey

% Get enrichment of cells positive for certain proteins to sit together or
% not. Quantify using a defined pixel distance. For each sample get all the
% interactions (defined as pixel distance) between cells positive for protein X
% and protein Y. The use the bootstrapping approach to permute the labels randomly
%
 
% Define paths and read in data
path = '/Users/erinmccaffrey/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/';
pathSegment = '/Users/erinmccaffrey/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/';
% points = [64,65,21,84,42,88,28,89,85,13,35,36,14,15,57,58,19,87,6,...
%           7,33,34,26,27,40,61,47,48,54,55]; %mycobacteria cohort data
points=[67,68,69,70,71,72,73,74,75,76]; %sarcoid cohort data
%points =  [7,13,21,28,33,42,47,85,88,89]; %revised after giant cell correction
dataAll=dataset('File',[pathSegment,'/dataPerCell_3px/granA_cellpheno_CS-asinh-norm_matlab_revised.csv'],'Delimiter',','); %concatenated and annotate matrix of all samples
%converting to matrix by going to cell array and then removing 1st row (headers)
dataAllMat=dataset2cell(dataAll);
dataAllMat=cell2mat(dataAllMat(2:55713,1:46)); %need to excude 47 which is a string

%Go over all pairwise combinations of markers
markerInds = [8,9,11:44]; %indices of the markers in matrix
dataMarkers = dataAllMat(:,markerInds); %subset of data matrix including only markers
markerTitles = dataAll(:,markerInds).Properties.VarNames; %names of the markers
markerNum = length(markerTitles); %number of markers to compare

BootstrapNum = 1000; %number of permutations for bootsrapping to create null
distLim = 100; %cutoff to define the cells that are interacting/close in pixels
closeNum = zeros(markerNum,markerNum); %matrix for number of interactions for each marker pair
closeNumRand = zeros(markerNum,markerNum,BootstrapNum); %matrix for null number of interactions

% define some useful variables for the analysis
patientIdx = 1; %column with patient/point label
cellLabelIdx = 2; %column with cell label

% read in marker thresholds (custom values in same order as matrix columns)
markerThresh=dataset('File',[pathSegment,'/dataPerCell_3px/markerThresholds.csv'],'Delimiter',','); 
thresh_vec = dataset2cell(markerThresh);
thresh_vec = cell2mat(thresh_vec(3:38,2));

for i=1:length(points)
    point=points(i);
    % load relevant data
    disp(['Working on point:',num2str(point)]);
    load([pathSegment,'/Point',num2str(point),'/cellDistances.mat']);
    % get data for current patient
    patientInds = dataAllMat(:,patientIdx) == point; %get row indices of current point
    patientData = dataAllMat(patientInds,:); %get just the data for the current point
    patientDataMarkers = dataMarkers(patientInds,:); %get just the marker data for the current point
    % go over markers
    for j=1:markerNum
        marker1_thresh = thresh_vec(j); %get positiviy threshold for marker
        marker1PosInds = (patientDataMarkers(:,j) > marker1_thresh ); %indices of cells positive for that marker
        marker1PosLabels = patientData(marker1PosInds,cellLabelIdx); %labels of cells positive for that marker
        marker1Num = length(marker1PosLabels); %number of cells positive for that marker
        % iterate over all other markers + curr marker
        parfor k=1:markerNum
            marker2_thresh = thresh_vec(k);
            marker2PosInds = ( patientDataMarkers(:,k) > marker2_thresh ); %indices of cells positive for that marker
            marker2PosLabels = patientData(marker2PosInds,cellLabelIdx); %labels of cells positive for that marker
            marker2Num = length(marker2PosLabels); %number of cells positive for that marker
            truncDistMat = distancesMat(marker1PosLabels,marker2PosLabels); %create a distance matrix of just pos v pos cells
            % turn to binary
            truncDistMatBin = zeros(size(truncDistMat));
            truncDistMatBin(truncDistMat<distLim) = 1; %convert to binary of interacting or not
            % record interaction num
            closeNum(j,k) = sum(sum(truncDistMatBin))  %add number of interactions to matrix of 'real' interactions
            % randomize to get null distribution
            for r=1:BootstrapNum
                % get random labels for marker 1
                marker1LabelsRand = datasample(patientData(:,cellLabelIdx),marker1Num);
                % get random labels for marker 2
                marker2LabelsRand = datasample(patientData(:,cellLabelIdx),marker2Num);
                randTruncDistMat = distancesMat(marker1LabelsRand,marker2LabelsRand); % distance matrix of the random labels
                % turn to binary
                randTruncDistMatBin = zeros(size(randTruncDistMat)); %convert to binary of interacting or not
                randTruncDistMatBin(randTruncDistMat<distLim) = 1; 
                % record interaction num
                closeNumRand(j,k,r) = sum(sum(randTruncDistMatBin)); %add number of interactions to matrix of 'null' interactions
            end
        end
    end

    % Create a vector for the z-scoring (and needed params mean + std)
    z = zeros(markerNum); %z-score
    muhat = zeros(markerNum); %mean
    sigmahat = zeros(markerNum); %std deviation
    
    % Create a vector for the empirical p-value 
    p = zeros(markerNum,markerNum,2);

    % Get the enrichment z-score for each marker
    for j=1:markerNum
        for k=1:markerNum
            tmp= reshape(closeNumRand(j,k,:),BootstrapNum,1);
            [muhat(j,k),sigmahat(j,k)] = normfit(tmp);
            z(j,k) = (closeNum(j,k)-muhat(j,k))/sigmahat(j,k);
            p(j,k,1) = (1+(sum(tmp>=closeNum(j,k))))/(BootstrapNum+1);
            p(j,k,2) = (1+(sum(tmp<=closeNum(j,k))))/(BootstrapNum+1);
        end
    end
    
    % adjust p-values using FDR 0.05 (Inf or NaN will have p value 1)
    p_summary = p(:,:,1);
    for j=1:markerNum
        for k=1:markerNum
            % if interaction is enriched +z grab first p-value
            if (z(j,k) > 0)
                p_summary(j,k) = p(j,k,1);
            % if interaction is avoided -z grab second p-value
            elseif (z(j,k) < 0) 
                p_summary(j,k) = p(j,k,2);
            end
        end
    end
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_summary,.05,'pdep','yes');

    % save results
    save([pathSegment,'/Point',num2str(point),'/spatialAnalysisv2.mat'],'closeNum','closeNumRand','z','muhat','sigmahat','p','adj_p','h');

end


% visualize
resultsDir = [pathSegment,'/dataPerCell_3px/plots/draft_figs/spatial/protein_enrichv2/'];
mkdir(resultsDir);
for i=1:length(points)
    point=points(i);
    disp(['Working on point:',num2str(point)]);
    load([pathSegment,'/Point',num2str(point),'/spatialAnalysisv2.mat'])
    % save results
    save([pathSegment,'/Point',num2str(point),'/spatialAnalysisv2.mat'],'closeNum','closeNumRand','z','muhat','sigmahat','p','adj_p','h');
    zplot = z;
    zplot(isnan(zplot)) = 0;
    zplot(isinf(zplot)) = 0;
    hmap=clustergram(zplot([1:34,36],[1:34,36]),'RowLabels', markerTitles([1:34,36]),...
        'ColumnLabels', markerTitles([1:34,36]),'Colormap', 'redbluecmap','DisplayRange',...
        20, 'DisplayRatio', 0.2);
    addTitle(hmap, ['Point ',num2str(point)])
    fig=plot(hmap);
    saveas(fig,([resultsDir,'Point',num2str(point)]),'epsc')
    save([pathSegment,'/Point',num2str(point),'/spatialAnalysisv2.mat'],'closeNum','closeNumRand','z','muhat','sigmahat','p','adj_p','h','markerTitles');
end

