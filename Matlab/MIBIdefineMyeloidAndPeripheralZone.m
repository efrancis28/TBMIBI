%% MIBIdefineMyeloidAndPeripheralZone.m
% This script takes in a mask of the myeloid-rich region of selected
% samples. Next it, defines a peripheral zone based on the defined pixel
% expansion.

%define path and points
path = '/Users/erinmccaffrey/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/summed_data';
points = [6,7,13,21,28,33,34,35,36,40,47,48,54,55,57,61,84,88,89]; %cohort data

%define expansion for peripheral zone
periph_distance = 100;
se=strel('disk', periph_distance);

for p=1:length(points)
    % load data
    point=points(p);
    disp(['point',num2str(point)]);
    mask=imread([path,'/Point',num2str(point),'/s',num2str(point),'_mask','.tif']);
    mask_outline=boundarymask(mask);
    imwrite(uint16(mask_outline),[path,'/Point',num2str(point),'/regionmask_outline.tif'],'tif');
    mask_combined=imdilate(mask,se);
    periph_mask=mask_combined-mask;
    save ([path,'/Point',num2str(point),'/regionmask.mat'],'mask','periph_mask');
end
    
    
mask_outline=boundary_mask(periph_mask);