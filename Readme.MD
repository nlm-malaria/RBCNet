[![img-29.png](https://i.postimg.cc/LsCMjpPs/img-29.png)](https://postimg.cc/B8KyGRMW)
https://user-images.githubusercontent.com/16782894/92193000-42d2a580-ee35-11ea-911b-06a9337613a9.png
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
