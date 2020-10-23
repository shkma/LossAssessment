# LossAssessment
This repository contains seismic loss assessment MATLAB code.

Prepared by: Shoma Kitayama (shomakit@buffalo.edu)
Last modified: 5.Oct.2020
Version: 1.0

The folder contains MATLAB codes that computes the following values that are the results of seismic loss assessment of seismically isolated buildings with SCBF designed by RI=2 and TFP-1 based on Conditional Spectra approach (notations are explained in the manuscript):

- Loss vulnerability functions
- Expected annual loss (EAL)
- Expected Loss (=EL) Over Time (=t)

To implement this example code, you need to download the data of results of nonlinear seismic response analysis (i.e., multiple stripe analysis; Jalayer, 2003) from the following link:

https://drive.google.com/file/d/1pNAoIicndCVjONx4P9fMqTAlJawJbsfy/view?usp=sharing

The downloaded data (i.e., a folder, after extraction of zip file, "Response_analysis_data_CS_SCBF_RI2.0_Iso1_ASCE16") of results of seismic respone analysis should be located to the folder in the same place as the matlab files are located.

The zip file should be extracted to get a folder that contains data.
The data is the results of analysis of seimically isolated building with SCBF (superstructure is designed per RI=2.0) and with the smallest required isolators (TFP-1) based on ASCE7-16 (2017) seismic design standard.

Please note that one of the MATLAB file "fn_mle_pc.m" used in this folder was obtained from the data set that is associated with the following journal article:
Baker JW. (2015). “Efficient analytical fragility function fitting using dynamic structural analysis.” Earthquake Spectra, 31(1), 579-599.

---
Cited article that does not appear in the manuscript:
Jalayer F. (2003). "Direct probabilitic seismic analysis: Implementing non-linear dynamic assessments." PhD. Thesis. Stanford University.
