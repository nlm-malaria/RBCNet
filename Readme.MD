# RBCNet for cell detecton using  a  dual deep learning architecture
Computer-assisted algorithms have become a mainstay of biomedical applications to improve accuracy and reproducibility of repetitive tasks like manual segmentation and annotation. We propose a novel pipeline for red blood cell detection and counting in thin blood smear microscopy images, named RBCNet, using a dual deep learning architecture. RBCNet consists of a  U-Net first stage for cell-cluster segmentation, followed by a second stage Faster R-CNN for detecting small cell objects within clusters, identified as connected components from the U-Net stage. RBCNet uses cell clustering instead of region proposals, which is robust to cell fragmentation, is highly scalable for detecting small objects or fine scale morphological structures in very large images, can be trained using non-overlapping tiles, and during inference is adaptive to the scale of cell-clusters with a low memory footprint.   
We tested our method on an archived collection of human malaria smears with nearly 200,000 labeled cells across 965 images from 193 patients, acquired in Bangladesh, with each patient contributing five images. Cell detection accuracy using RBCNet was higher than 97%.
The novel dual cascade RBCNet architecture provides more accurate cell detections because the foreground cell-cluster masks from U-Net adaptively guide the detection stage, resulting in a notably higher true positive and lower false alarm rates, compared to traditional and other deep learning methods. The RBCNet pipeline implements a crucial step towards automated malaria diagnosis. 

[![pipeline.png](https://i.postimg.cc/SxZcvCZV/pipeline.png)](https://postimg.cc/xXzkkJVb)

# How to use CNV

There are Three main folders in our repository:

**Code:** contains the algorithm scripts.   
**Data:** contains input images.   
**Output:** contains bounding box detections, scores, and visualization.   
**Examples:** contains some images and video examples.   


## How it works :  

1. Malaria images should be placed in the Data folder (patient input folder):

``` Ex: .\Data\234C92P53ThinF\Img\IMG_20150821_150457.png```

2. Run ```RBCNet_RunMe.m```

This script will:
1. Load learning models (U-Net + Faster R-CNN).  
2. Produce segmentation mask using U-Net.  
3. Detect the cells for each cluster of blobs of the raw image corresponding to the segmented mask from U-Net using Faster R-CNN.  
4. Filter out WBC. 
5. Store bounding boxes, scores, and visualization in the Output folder (patient output folder).  

[![img-29.png](https://i.postimg.cc/LsCMjpPs/img-29.png)](https://postimg.cc/B8KyGRMW)
[![img-29.png](https://i.postimg.cc/LsCMjpPs/img-29.png)](https://youtu.be/--mDaNat03Y)

### Project Collaborators and Contact

**Author:** Dr. Yasmin M. Kassim    
**Corresponding Author:** Dr. Stefan Jaeger


**Created by:** Dr. Yasmin Kassim  
Lister Hill National Center for Biomedical Communications,  
National Library of Medicine - National Institutes of Health

For more information, contact:

* **Dr. Yasmin Kassim**  
Location:  38A / B1N-28N  
Phone Number: (301) 827-4730  
E-mail: yasmin.kassim@nih.gov  


* **Dr. Stefan Jaeger**  
Location:  38A / 10N1003O  
Phone Number: (301) 435-3198  
E-mail: stefan.jaeger@nih.gov


