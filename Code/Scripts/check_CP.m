function [newmask] = check_CP(im2,im1)

binarymask=im2;

[labeledImage, numberOfBlobs] = bwlabel(im2);
measurements = regionprops(labeledImage, 'BoundingBox', 'area');
Allareas=[measurements.Area];
blobBoundingBox = regionprops(labeledImage, 'BoundingBox');
[sortedAreas, sortIndexes] = sort(Allareas, 'ascend');

if max(sortedAreas)>70000

SE=strel('disk',100);
temp=rgb2gray(im1);
circularmask=zeros(size(im2,1),size(im2,2));
circularmask(temp>20)=1;
circularmask(temp<=20)=0;
circularmask2=imerode(circularmask,SE);

im1=im1(:,:,2);

temp=im1(im1>10);

level=median(temp)-3;

tmp=im1;

im1((tmp<level))=0;

im1((im1~=0))=1;

im2(im1==1)=0;

SE=strel('disk',10);

mask=imopen(im2,SE);
gg=im2-mask;

ggclean = bwareaopen(gg, 120);

newmask=logical(ggclean+mask+logical((~circularmask2.*binarymask)));

else
   newmask= binarymask;
end
end

