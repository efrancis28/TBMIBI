%% MIBIannotateCellsinGranulomaZones.m
% This script takes in a mask of the myeloid-rich and peripheral zones. For
% each point and cell it creates a binary matrix of whether or not that
% cell label is in the myeloid zone or the peripheral zone. It can be found
% in both or neither.

%define path and points
pathMask = '/Users/erinmccaffrey/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/summed_data';
path = '/Users/erinmccaffrey/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/';
points = [6,7,13,21,28,33,34,35,36,40,47,48,54,55,57,61,84,88,89]; %cohort data
dataAll=dataset('File',[path,'/dataPerCell_3px/granA_cellpheno_CS-asinh-norm_matlab_revised.csv'],'Delimiter',','); %concatenated and annotate matrix of all samples
%converting to matrix. This is not a clean way to do it but oh well for
%now..
dataAllMat=dataset2cell(dataAll);
dataAllMat=cell2mat(dataAllMat(2:55713,1:46)); %need to excude 47 which is a string for some reason -___-

%create matrix with just patient IDs and cell labels
dataAllPatientAndCells=dataAllMat(:,1:2);
dataAllPatientAndCells = dataAllPatientAndCells(ismember(dataAllPatientAndCells(:,1), points), :);

%define column with patient Idx and cell label for filetering
patientIdx = 1; %column with patient label
cellLabelIdx = 2; %column with cell label

%create output matrix
cell_regions = zeros(size(dataAllPatientAndCells,1),size(dataAllPatientAndCells,2)+2);
cell_regions(:,patientIdx)=dataAllPatientAndCells(:,patientIdx);
cell_regions(:,cellLabelIdx)=dataAllPatientAndCells(:,cellLabelIdx);
myeloidIdx = 3;
periphIdx = 4; 

%counter for indexing summary matrix
count=0;

for p=1:length(points)
    point=points(p);
    disp(['point',num2str(point)]);
    
    % load data
    load([pathMask,'/Point',num2str(point),'/regionmask.mat']);
    load([path,'/Point',num2str(point),'/segmentationParams3px.mat']);
    
    % combine masks
    mask(mask>0)=1;
    periph_mask(periph_mask>0)=2;
    combined_mask=mask+periph_mask;
    
    % get region stats for mask
    stats_mask = regionprops(combined_mask,'PixelIdxList');
    myeloid_px = stats_mask(1).PixelIdxList;
    periph_px = stats_mask(2).PixelIdxList;
    
    % get data just for current point
    patientInds = dataAllPatientAndCells(:,patientIdx) == point; %get row indices of current point
    patientData = dataAllPatientAndCells(patientInds,:); %get just the data for the current point
    stats_cells = regionprops(newLmod,'PixelIdxList');
    
    % iterate through all cells and get their location based on pixel idxs
    cells = patientData(:,2); 
    for j = 1:length(cells)
        % define cell object including the boundary
        cell = cells(j);
        cell_px = stats_cells(cell).PixelIdxList;
        %check myeloid zone
        myeloid_overlap = intersect(cell_px,myeloid_px);
        if (~isempty(myeloid_overlap))
            cell_regions(j+count,myeloidIdx) = 1;
        end
        %check periph zone
        periph_overlap = intersect(cell_px,periph_px);
        if (~isempty(periph_overlap))
            cell_regions(j+count,periphIdx) = 1;
        end
    end
    %update counter
    count=count+length(cells);
end

channelLabels = {'SampleID','cellLabelInImage','Myeloid_Zone','Periph_Zone'};
TEXT.PnS = channelLabels;
csvwrite_with_headers([pathMask,'/granA_cellpheno_CS-asinh-norm_matlab_revised.csv'],cell_regions,TEXT.PnS)