# Unsupervised_ML_for_Hyperspectral_imaging



For cell segmentation on hyperspectral imagimg , use the function cell_segm(data , hne, imaging, dbg)
with 
- data being a struct object containing: rgb ( rgb vizualisation of the tissue) , spec ( the 3D matrix WidthxHeightxnbr_of_wavelength_channles containing the wavelength spectrum of each pixel in the image) , lambda ( the frequency of each wavelength channels respectively - in nm.).
- hne  a struct object containing  H  and E  : distance matrix after svd of H&E img   - H and E can be calcultated using [H,E] = get_HnE(wl_mat,lambda_mat);
- imaging = 'x20' or 'x40' according to the resolution of the microscope/objective
- dbg = 0  - debugging tool.


for preprocessing the data before applying cell segmentation , use preprocess_tissue(data, whitening) , whitening being a boolean .





the repository doesn’t contain the tissue I worked on, as well as the manually labelled ground truth, but they are available on a USB key, whenever requested to me .



![image](https://user-images.githubusercontent.com/30235346/209973718-7c186d8c-1ab6-4020-8c18-d63ddd1c1fc6.png)
