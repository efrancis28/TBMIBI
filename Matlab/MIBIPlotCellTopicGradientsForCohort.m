% MIBIPlotCellTopicGradientsForCohort.m
% Author: Erin McCaffrey (adapted from original script by Leeat Keren:
% MibiPlotImmunePopulationsForCohort170911.m) 
% Script reads in images for the cohort and then produces a new image of
% the segmentation mask where each cell object is colored by its topic.
% Keeps the colors consistent with those used for plotting in all other
% figures. 

%% Initiate paths to data and define important variables
path = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/';
pathSegment = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/';
points = [33];
imSize=1024;
dataAll=dataset('File',['/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis/all_gran_topic_annotation.csv'],'Delimiter',','); %concatenated and annotate matrix of all samples
%converting to matrix
dataAllMat=dataset2cell(dataAll);
dataAllMat=cell2mat(dataAllMat(2:46980,[1,2,55:65]));

%% Define topic to plot gradient for (ie. will only give values to cells
% within a given topic, but will still plot the probability instead of a
% binary assignment

% topic of interest
topic = 3;

% column containing the topping loadings for the given topic
topic_prob_col = topic + 3;

% column containing the topic assignment
topic_assn_col = 13;


%% Produce topic loading map

for i=1:length(points)
    %load data
    disp(points(i));
    point = points(i);
    load([pathSegment,'/Point',num2str(point),'/segmentationParams3px.mat']);
    %get stats for objects in image
    stats = regionprops(newLmod,'PixelIdxList');
    %get just data for current point
    currInds = (dataAllMat(:,1) == point);
    currCellData = dataAllMat(currInds,:);
    %get labels of objects in image
    currLabelIdentity = labelIdentityNew;
    labelNum = length(currLabelIdentity);
    %create vector of all pixels for assigning color map
    imageL = zeros(size(newLmod,1)*size(newLmod,2),1);
    for j=1:labelNum 
        %get index of current cell
        cellInd = find(ismember(currCellData(:,2),j));
        if ~isempty(cellInd)
            % assign a color based on topic assignment
            if(currCellData(cellInd,topic_assn_col) == topic)
                loading = currCellData(cellInd,topic_prob_col);
                imageL(stats(j).PixelIdxList)=loading;
            end
        end
    end
    %shape back into an image
    imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));
    % plot
    f1=figure;
    colormap('parula');
    imagesc(imageLReshape);
    imwrite(imageLReshape(31:994,31:994),[pathSegment,'/Point',num2str(point),'/topic',num2str(topic),'_gradient.tif'],'tif');
end