%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                
%  Author: Dr. Yasmin M. Kassim  
%  Corresponding Author: Dr. Stefan Jaeger
%  Lister Hill National Center for Biomedical Communications,
%  National Library of Medicine - National Institutes of Health
%
%  For more information, contact:
%
%      Dr. Yasmin M. Kassim 
%      Location:  38A / B1N-28N
%      Phone Number: (301) 827-4730
%      E-mail: yasmin.kassim@nih.gov
% or
%      Dr. Stefan Jaeger
%      Location:  38A / 10N1003O
%      Phone Number: (301) 435-3198
%      E-mail: stefan.jaeger@nih.gov
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Program Name: RBCNet_RunMe.m
%  Inputs:  Malaria images should be placed in the Data folder (patient input folder)
%           Ex: .\Data\234C92P53ThinF\Img\IMG_20150821_150457.png
%           
%  Procedure: 1. Load learning models (U-Net + Faster R-CNN)
%             2. Produce segmentation mask using U-Net
%             3. Detect the cells for each cluster of blobs of the
%             raw image corresponding to the segmented mask from U-Net using Faster R-CNN
%             4. Filter out WBC
%        
%  Outputs:   Store bounding boxes, scores, and visualization in the
%             Output folder (patient output folder)
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath '.\Scripts\';
clear all;
%------------------------- Flags --------------------------------------------------------
watch_blob_detection=0;  % Set this flag if you want to watch the detection of each blob
                         % Be aware it will take much more extra time 
%----------------------------------------------------------------------------------------
opath='..\Data\';
folders = dir(opath);
Output_path='..\Output\';
% Define the two colormaps.
cellcolors={'blue', 'green', 'cyan', 'red', 'magenta', 'yellow','white'};
cc_no=0;
fprintf('Upload Learning models: U-Net for cell segmentation + Faster R-CNN to detect clusters of blobs')
%%-------- Load learning models (U-Net + Faster R-CNN)
load('UNet_segmentation.mat');
load('FRCNN_detector.mat');
SE=strel('disk',10);
SEC=strel('disk',12);
color_co=1;
for i_i = 1:numel(folders)
    i_i
inputr = strcat(opath,filesep,folders(i_i).name);    

if folders(i_i).isdir && ~strcmp(folders(i_i).name,'.') && ~strcmp(folders(i_i).name,'..')
    
folderPath =strcat(opath,filesep,folders(i_i).name,filesep,'Img',filesep);
rawimglist=dir(fullfile(folderPath,'*.jpg'));
n = size(rawimglist,1);

predictPatchSize = [1536 1536];
for i=1:n
val_data=imread(strcat(folderPath,rawimglist(i).name));

val_data=imresize(val_data,0.3);
val_dataf=val_data;
val_data(:,:,4)=ones(897,1594);
%%-------- Produce segmentation mask using U-Net
tic
segmentedImage = segmentImage(val_data,net,predictPatchSize);
%toc
segmentedImage = uint8(val_data(:,:,4)~=0) .* segmentedImage;
% figure
% imshow(segmentedImage,[])

temp=segmentedImage;
segmentedImage(temp==2)=255;
segmentedImage(temp==1)=0;

idxrect=1;
imOS=imread(strcat(folderPath,rawimglist(i).name));
imOS2=imOS;
binaryImage=logical(segmentedImage);

im2=imresize(imOS,0.3);

binaryImage=check_CP(binaryImage,im2);

im2rgb=rgb2gray(im2);
im2rgb_B=uint8(binaryImage).*im2rgb;
thresholdWBC=mean(mean(im2rgb_B(im2rgb_B>0)));

temp=rgb2gray(im2);
circularmask=zeros(size(im2,1),size(im2,2));
circularmask(temp>20)=1;
circularmask(temp<20)=0;
circularmask2=imerode(circularmask,SE);

binaryImage=binaryImage.*circularmask2;
%
binaryImage = bwareaopen(binaryImage, 250);

R=im2(:,:,1);
G=imadjust(im2(:,:,2));
B=im2(:,:,3);

im2(:,:,1)=R;im2(:,:,2)=G;im2(:,:,3)=B;


[labeledImage, numberOfBlobs] = bwlabel(binaryImage);

measurements = regionprops(labeledImage, 'BoundingBox', 'area');
Allareas=[measurements.Area];
blobBoundingBox = regionprops(labeledImage, 'BoundingBox');
[sortedAreas, sortIndexes] = sort(Allareas, 'ascend');

%%------ Detect the cells for each cluster of blobs of the
%%------ raw image corresponding to the segmented mask from U-Net using Faster R-CNN

for k =1: numberOfBlobs           % Loop through all blobs.
    mask_with_requiredBlob = ismember(labeledImage, k);
    mask_with_requiredBlob=imfill(mask_with_requiredBlob,'holes');
    mask_with_requiredBlob=imdilate(mask_with_requiredBlob,SEC);
    im3(:,:,1)=im2(:,:,1).*uint8(mask_with_requiredBlob);
    im3(:,:,2)=im2(:,:,2).*uint8(mask_with_requiredBlob);
    im3(:,:,3)=im2(:,:,3).*uint8(mask_with_requiredBlob);
    Rmax=max(max(im2(:,:,1)));    Gmax=max(max(im2(:,:,2)));    Bmax=max(max(im2(:,:,3)));
    
    im2_connected_comp=imcrop(im3,blobBoundingBox(k).BoundingBox);
    R_comp=im2_connected_comp(:,:,1);
    G_comp=im2_connected_comp(:,:,2);
    B_comp=im2_connected_comp(:,:,3);
    R_comp(R_comp==0)=Rmax;
    G_comp(G_comp==0)=Gmax;
    B_comp(B_comp==0)=Bmax;
    im2_connected_comp(:,:,1)=R_comp;
    im2_connected_comp(:,:,2)=G_comp;
    im2_connected_comp(:,:,3)=B_comp;

    
    if (size(im2_connected_comp,1)>20 && size(im2_connected_comp,1)<33) 
    im2_connected_comp=imresize(im2_connected_comp,[33 size(im2_connected_comp,2)]);    
    end
    if (size(im2_connected_comp,2)>20 && size(im2_connected_comp,2)<33) 
    im2_connected_comp=imresize(im2_connected_comp,[size(im2_connected_comp,1) 33]);    
    end
    if size(im2_connected_comp,1)>32 && size(im2_connected_comp,2)>32

    [bboxes,scores] = detect(detector,im2_connected_comp);
    bboxes=bboxes/0.3;
    n=(size(bboxes,1));
    
    for ib=1:n
    bboxes(ib,1)=bboxes(ib,1)+(blobBoundingBox(k).BoundingBox(1)/0.3);
    bboxes(ib,2)=bboxes(ib,2)+(blobBoundingBox(k).BoundingBox(2)/0.3);
    %%---- Checking if WBC
    RBC=imresize(imcrop(imOS,bboxes(ib,:)),0.3);    
    RBCM=imresize(imcrop(binaryImage,bboxes(ib,:)*0.3),[size(RBC,1) size(RBC,2)]);
    RBCM=imfill(RBCM,'holes');
    WBCchecking=rgb2gray(RBC).*uint8(RBCM);
    WBCchecking=WBCchecking(WBCchecking~=0);
    dark_region=find(WBCchecking<thresholdWBC-50);
    %%---------------------
    
    if numel(dark_region)<200 && (numel(dark_region)<(1/3)*size(WBCchecking,1))
    allrect(idxrect,:)=bboxes(ib,:);allscores(idxrect)=scores(ib);idxrect=idxrect+1;
    else
    bboxes(ib,:)=[0 0 0 0];
    end
    end
    
    if ~isempty(bboxes) 
    if watch_blob_detection==1
    cc_no=cc_no+1;
    cellcolor=cellcolors(mod(color_co,7)+1);color_co=color_co+1;
    imOS = insertObjectAnnotation(imOS,'rectangle',bboxes,scores,'LineWidth',5, 'FontSize',20,'Color',cellcolor);
    imshow(imOS,[])
    end
    else
    bboxes=[1 1 40 40];scores=1;
    bboxes=bboxes/0.3;
    n=(size(bboxes,1));
    
    for ib=1:n
    bboxes(ib,1)=bboxes(ib,1)+(blobBoundingBox(k).BoundingBox(1)/0.3);
    bboxes(ib,2)=bboxes(ib,2)+(blobBoundingBox(k).BoundingBox(2)/0.3);
    %%---- Checking if WBC
    RBC=imresize(imcrop(imOS,bboxes(ib,:)),0.3);    
    RBCM=imresize(imcrop(binaryImage,bboxes(ib,:)*0.3),[size(RBC,1) size(RBC,2)]);
    RBCM=imfill(RBCM,'holes');
    WBCchecking=rgb2gray(RBC).*uint8(RBCM);
    WBCchecking=WBCchecking(WBCchecking~=0);
    dark_region=find(WBCchecking<thresholdWBC-50);
    %%---------------------
    if numel(dark_region)<200 && (numel(dark_region)<(1/3)*size(WBCchecking,1))
    allrect(idxrect,:)=bboxes(ib,:);allscores(idxrect)=scores(ib);idxrect=idxrect+1;
    else
    bboxes(ib,:)=[0 0 0 0];
    end
    end
    end

end
end
toc
%%------ Outputs:   Store the bounding boxes, scores, and visulaization in the
%%------ Output folder in the corresponding patient
allscores=allscores';
imOS2 = insertObjectAnnotation(imOS2,'rectangle',allrect,allscores,'LineWidth',5, 'FontSize',20);
close all
imshow(imOS2,[]);

allrect_Path =strcat(Output_path,filesep,folders(i_i).name,filesep,rawimglist(i).name(1:end-4),filesep,'Bounding Boxes',filesep);
allscores_path =strcat(Output_path,filesep,folders(i_i).name,filesep,rawimglist(i).name(1:end-4),filesep,'Scores',filesep);
vis_path =strcat(Output_path,filesep,folders(i_i).name,filesep,rawimglist(i).name(1:end-4),filesep,'Visualization',filesep);
    if ~exist(allrect_Path, 'dir')
       mkdir(allrect_Path)
    end
    if ~exist(allscores_path, 'dir')
       mkdir(allscores_path)
    end
    if ~exist(vis_path, 'dir')
       mkdir(vis_path)
    end
rectpath=strcat(allrect_Path,rawimglist(i).name(1:end-4),'.mat');
save(rectpath,'allrect');
scorepath=strcat(allscores_path,rawimglist(i).name(1:end-4),'.mat');
save(scorepath,'allscores');
wpathFRCNN=strcat(vis_path,rawimglist(i).name(1:end-4),'.png');
imwrite(imOS2,wpathFRCNN);  
pause(5);
clear allrect
clear allscores
end
end
end
