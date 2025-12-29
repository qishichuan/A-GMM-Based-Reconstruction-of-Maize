# 3D Reconstruction Method for Organ-Level Morphology of Field Maize Populations
## Project Description
This project is based on the thesis "A GMM Based Organ-Level Reconstruction of Maize in Field-Scale by Fixed-Point Sliding Scanning Strategy of Mobile Robots" and implements a point cloud-based 3D reconstruction method for maize plants in real field environments. The proposed method employs a two-stage pipeline:  
### First-stage: 
Individual Plant Registration. This stage utilizes a batch registration algorithm based on Gaussian Mixture Models (GMM). Building upon the Batch Expectation-Maximization (EM) framework proposed by Evangelidis et al. in 'Joint Alignment of Multiple Point Sets with Batch and Incremental Expectation-Maximization'(https://ieeexplore.ieee.org/document/7968037), the algorithm performs joint alignment of multi-frame point clouds to effectively mitigate cumulative errors.  
### Second-stage: 
Large-Scale Population Reconstruction. This phase achieves global registration by integrating pose data derived from GNSS/IMU systems.  
Designed specifically for agricultural phenotyping, this approach emphasizes the robustness inherent in non-deep learning architectures.
Note: This method relies on traditional geometric processing algorithms (non-deep learning), thus it does not involve trainable model weights. Algorithm parameters are derived through statistical estimation (e.g., expectation-maximization).
## Environment Dependencies
### Operating System: Linux (Ubuntu 18.04)  
### Core Libraries:  
PCL (Point Cloud Library) 1.8+  
ROS (Robot Operating System) Melodic  
RVIZ  
### Hardware:
16-line LiDAR (e.g., Velodyne VLP-16)  
GNSS/IMU modules (e.g., CGI-410 + Xsens MTi-300)  
Industrial computer (i7 processor, 16GB RAM) Usage
## Data Collection
Use the self-built field mobile platform to collect maize point cloud data.  
Ensure GNSS/IMU data is synchronized and stored with point clouds.
## Model Weights Explanation
This method does not involve deep learning models; hence, no model weight files are included. The algorithm depends on the following key parameters (optimized through experiments in the thesis):
### GMM Registration
Number of Gaussian components K: Set based on point cloud density (default 10-20)  
Iteration count Q: 50-100 (based on convergence criteria)
## Model Testing Instructions
Testing is based on real field experiments (Hefei Dawei Experimental Field) and covers the following aspects:  
### Test Dataset
Environment: Maize at V6-V10 leaf stages, row spacing 0.65m, average plant height 0.85m.  
Data Volume: 30 single-plant point clouds + 3 population point cloud blocks (22m × 7.5m).  
### Evaluation Metrics
Registration Accuracy:
Precision: GMM algorithm achieves 0.97 (indoor) / 0.95 (field)
Target reconstruction error: Area ratio R1=1.12, aspect ratio R2=0.99
Population Reconstruction Validation:
Phenotypic measurements vs. manual ground truth:
Plant height: R2=0.99, RMSE = 0.015 m
Stem width: R2=0.79, RMSE = 0.0039 m (after correction K1=0.67)
Leaf width: R2=0.79, RMSE = 0.0093 m (after correction K2=0.71)
## Result Examples
Single-plant registration results: Point clouds retain 3D leaf morphology   
Result Examples
Single-plant registration results: Point clouds retain 3D leaf morphology
<img width="868" height="331" alt="图片" src="https://github.com/user-attachments/assets/b528e3fb-a736-4da7-a703-f32dceb54de9" />
Population reconstruction: Extracted single-plant point clouds show clear hierarchy
<img width="1010" height="897" alt="图片" src="https://github.com/user-attachments/assets/f6a4d05f-bb1c-4709-8f41-8868be81d75a" />
## *Corresponding author
Lu liu, Associate professor
Anhui Agricultural University
email:vliulu@ahau.edu.cn
