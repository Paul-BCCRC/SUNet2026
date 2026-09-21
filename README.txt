SuNet:

1) Open Segment_dir.py

2) Set "imagefilepath" to the directory with TIFF images to segment.

3) Set "outpath" to where you want the output. Can be same as imagefilepath.

4) Open "UNET_Segmentation_Sequence_Function.py" Change the Run Parameters to the desired values.

5) Change the weight paths to .pth files in folder SegmentationWeights. 
	PathA = weights for predicting nuclei centers.
	PathB = weights for predicting nuclei boundaries.
	PathC = weights for predicting cytoplasm boundaries. 

6) Run Segment_dir.py



CCG Editor:

1) Open/Run CCG_Editor_App.mlapp in MATLAB

2) Click Source button and select cmg/ccg you want to open and wait for it to load.

3) Hover mouse over Edit View to see what the buttons are. Click it to open up the main display showing the image and mask data within the cmg/ccg.

4) "Channels" button lets you change how the Image data is viewed.
   "Masks" button changes how the mask data is viewed.
   "Feature View" and "Threshold Feature" will be available if the cmg/ccg has a corresponding feature file (.cb4)

5) After done editing select output directory with the output button and type in the file name with extension .cmg or .ccg.

6) Pressing "Save CMG" button will save a new cmg/ccg in the output directory with the output name.

