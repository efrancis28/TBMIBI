% % MIBICreateDistanceMatrixBetweenCells
% Author: Leeat Keren, modified by Erin McCaffrey
% For each point, reads in the massDS and the segmentation path containing
% the relevant segmentation params. It next produces a distance matrix of
% the center of all cells to all other cells for downstream spatial
% analysis. 

% path to the segmentation params and single cell data
path = '/Users/erinmccaffrey/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised';

% points to create distance matrix for
points = [64,65,21,84,42,88,28,89,85,13,35,36,14,15,57,58,19,87,6,...
    7,33,34,26,27,40,61,47,48,54,55,67,68,69,70,71,72,73,74,75,76];

for i=1:length(points)
    point=points(i);
    disp(['point',num2str(point)]);
    % load the segmentation and the cell data
    load([path,'/Point',num2str(point),'/segmentationParams3px.mat']);
    load([path,'/Point',num2str(point),'/cellData3px.mat']);
    % get centers of cells
    stats = regionprops(newLmod,'centroid');
    % get x,y coordinates in a mat
    centroidCoord = zeros(length(stats),2);
    for c=1:length(stats)
        centroidCoord(c,1) = stats(c).Centroid(1);
        centroidCoord(c,2) = stats(c).Centroid(2);
    end
    % create a matrix of distances between all cells
    distancesMat = pdist2(centroidCoord,centroidCoord);
    save([path,'/Point',num2str(point),'/cellDistances.mat'],'distancesMat');
end
