# NuclearEdge
Custom-written MATLAB (tested in MATLAB R2021b) scripts to find nuclear edge, nuclear area and FISH spots, then calculate the the relative distance from the FISH spots to the nuclear edge (i.e., the distance from a FISH spot to the nuclear edge divided the square root of the nuclear area)

This package includes:
ReadTifFindedges_batch_forCy5_autopick.m
nuclei_counter.m
point_to_line_distance.m
Example input single cell images (cell001.tif, cell002.tif, cell003.tif)

Input: single cell images in TIF format. Each image should contain only a single cell. The image should have three channels: channel 1 is sgGOLDFISH signals; channel 2 is MUC4-R GOLDFISH signals; channel 3 is Hoechst signals staining nucleus.

Output:  
*_EdgeDlist.dat: contains the relative distances from the FISH spots to the nuclear edge for each cell (i.e., the distance from a FISH spot to the nuclear edge divided the square root of the nuclear area).

How to use:
1. Put single cell images (e.g., cell001.tif, cell002.tif, cell003.tif) in the same folder with ReadTifFindedges_batch_forCy5_autopick.m.
2. Run ReadTifFindedges_batch_forCy5_autopick.m. in MATLAB.
3. The program will analyze all TIF image files in the folder one by one. The program will display its first try of identifing nuclear edge/FISH spots. User will need to visually inspect if the program correctly identify nuclear edge and FISH spots. There will be prompts in the Command Window to instruct user to either confirm that the program correctly identified nuclear edge/FISH spots or follow the prompt information to adjust parameters. For later case, the program will use the user-adjusted parameters and try to identifing nuclear edge/FISH spots again and display the results to user (user can adjust the parameter again if necessary). 
4. The program will make a folder called "AfterAnalysis". Analyzed images and output *_EdgeDlist.dat files are in the folder.

Expected output: 
cell001: 0.207689463203132
cell002: 0.157824831313985
cell003: 0.226008597725398
Expected run time: 2 seconds per cell.
