% MIBIPlotCellTopicsForCohort.m
% Author: Leeat Keren (original script:
% MibiPlotImmunePopulationsForCohort170911.m) 
% Script reads in images for the cohort and then produces a new image of
% the segmentation mask where each cell object is colored by its topic.
% Keeps the colors consistent with those used for plotting in all other
% figures. 

%% Initiate paths to data and define important variables
path = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/';
pathSegment = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/';
% points = [14,64,65,21,84,42,88,28,89,85,13,35,36,15,57,58,19,87,6,...
%     7,33,34,26,27,40,61,47,48,54,55]; %cohort data
points = [61];
imSize=1024;
dataAll=dataset('File',['/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis/all_TB_topic_annotation.csv'],'Delimiter',','); %concatenated and annotate matrix of all samples
%converting to matrix. This is not a clean way to do it but oh well for
%now..
dataAllMat=dataset2cell(dataAll);
dataAllMat=cell2mat(dataAllMat(2:33693,[1:46,52:63]));
dataAllMat(:,58) = dataAllMat(:,58) +1;

% column containing the numerically encoded FlowSOM phenotype
phenoCol=58;
numPheno=length(unique(dataAllMat(:,phenoCol)));

%load pheno key and define the num/name of phenos
% code myeloid=1, lymph=2, granulocyte=3, imm_other=4, non_immune=5
phenoKey=dataset('File',[pathSegment,'/dataPerCell_3px/topic_numkey.csv'],'Delimiter',','); 

%% Define color code for topic

%load color key with hexcodes and RGB values
colorKey=dataset('File',[pathSegment,'/dataPerCell_3px/custom_color_keys/colorkey_topics.csv'],'Delimiter',','); 

% go through all cell types and produce a key with numerical codes for
% pheno and color
phenoCodes=[];
phenoColors=[];
for i=1:numPheno
    phenoCodes(i,1)=i;
    pheno = string(phenoKey(phenoKey.Code == i, :).x___Pheno);
    phenoR = colorKey(colorKey.x___topic == pheno, :).R;
    phenoG = colorKey(colorKey.x___topic == pheno, :).G;
    phenoB = colorKey(colorKey.x___topic == pheno, :).B;
    phenoRGB = [phenoR, phenoG, phenoB];
    phenoColors(i,1:3) = phenoRGB;
end

%combine into single matrix
phenoRGBkey = [phenoCodes,phenoColors];

% colorcode:
% [ 0 bg - white, cells - rainbow ]
cmap = phenoColors;
cmap01 = cmap/255;

%% Recolor cell masks by topic assignment

% large matrix containing all images
images=zeros(imSize,imSize,3,length(points),'uint8');
images_lin=zeros(imSize,imSize,3,length(points),'uint8');

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
    %go through all labels to determine color assignment
    for j=1:labelNum
        %get index of current cell
        cellInd = find(ismember(currCellData(:,2),j));
        if ~isempty(cellInd)
            % assign a color based on pheno
            cellVal = currCellData(cellInd,phenoCol);
            imageL(stats(j).PixelIdxList)=cellVal;
        end
    end
    %shape back into an image and store in image stack
    imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));
    images(:,:,:,i) = label2rgb(imageLReshape,cmap01,'w');
    % plot single
    f1=figure;
    colormap('parula');
    imagesc(label2rgb(imageLReshape,cmap01,'w'));
    title(['Point ', num2str(point)])
    resultsDir = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Paper/Figure2_LDA';
    imwrite(label2rgb(imageLReshape(31:994,31:994),cmap01,'w'),[resultsDir,'/Point',num2str(point),'_topic_overlay.tif'],'tif');
end