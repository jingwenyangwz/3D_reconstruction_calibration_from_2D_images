# 3D_reconstruction_calibration_from_2D_images
Obtention of the intrinsic parameters of a camera; Finding local matches between several views of an object; 3D reconstruction and calibration.

- section 1: obtaining intrinsic parameters of the camera
  - run Obtention_intrinsic_paras_of_a_camera/jingwen_yang_1.m
  - we calibrate the camera with a checkerboard
  - With :
    - Size, in millimeters, of the checkerboard in your screen.
    - The set of images of the screen checkerboard that you have used for calibration.
    - The resolution of these captured images (in pixels).
    
we can calculate the internal matrix of a camera.


- section 2: Finding local matches between several views of an object
  - object/scene capture
  - Detection, description and matching of feature points: 
  - Qualitative and quantitative evaluation:
    a. Include the estimated Fundamental matrix.
    b. Qualitative evaluate the quality of the estimated fundamental matrix trough the
vgg_gui_F.m function of the vgg-mvg toolbox. Discuss on the obtained results in the report based on the views’ capture conditions. 
    c. Quantitative evaluate the quality of the estimated fundamental matrix by accounting for the number of inliers matchings (i.e. those that agree with the estimated fundamental
matrix) that are returned by estimateFundamentalMatrix. Include illustrative Tables (i.e. referencing views) in your exam’s report to aggregate and present the obtained
results.

  -Selection:
  
  Try out different combination of detectors and descriptors, choose the best combination. 
 
- section 3:  3D reconstruction and calibration

  - Compute consistent point matches among N views. 
  - Compute the Fundamental matrix and an **initial projective reconstruction** from 2 of the cameras. Visualize the mean re-projection error and the reprojection error histogram.
  - Improve this initial reconstruction by means of a Projective **Bundle Adjustment**, using a higher number (maybe all) of your images. In the report:  Visualize the mean re-projection error and the reprojection error histogram.
  - Re-compute the Fundamental matrix between two of the cameras, using the projection matrices obtained after the Projective Bundle Adjustment step
  - Use the properties of the Essential matrix (between two cameras) to obtain a Euclidean reconstruction of the scene (use the re-projected points obtained after the Projective Bundle Adjustment step). 
  Visualize the mean re-projection error and the reprojection error histogram, illustrative results from several viewpoints and the 3D Matlab figure.
  
  
  
  
  
  
