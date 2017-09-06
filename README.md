# TMI-2017-Ultrasound-Toews
Code for Tracked Ultrasound Auto-calibration, see article in IEEE Transactions on Medical Imaging 2017


Copyright (c) 2017, Matthew Toews
All rights reserved.

By downloading, copying, installing or using the software you agree to this license.
If you do not agree to this license, do not download, install, copy or use the software.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

This archive contains C++ code for the algorithm described in the article:
   "Phantomless Auto-Calibration and Online Calibration Assessment for a Tracked Freehand 2D Ultrasound Probe",
  Matthew Toews, William Wells, IEEE TMI 2017 

Contact information:
   matt.toews@gmail.com
   matthew.toews@etsmtl.ca
   sw@bwh.harvard.edu

Files:
  feature_extraction_matching.cpp : feature extraction and correspondence, links to library OpenCV2.2
  autocalibration.cpp             : autocalibration from feature correspondences, links to library nlopt-2.2.4
  sweep_2c.txt                    : sample ultrasound sequence input file
  unrelated.US.txt                : sample unrelated ultrasound sequence file
  unrelated.US.surf.dsc           : sample unrelated ultrasound descriptor file
  unrelated.US.surf.key           : sample unrelated ultrasound keypoint file

Additional Data:
  example MNI BITE data set available at www.matthewtoews.com/projects/MNI-BITE.group1_02.zip

Description:
  The code in this archive performs offline processing with the following steps, described below:
     1) feature extraction
     2) feature correspondence
     3) autocalibration

1) Feature Extraction: extract features from one ultrasound sequence.
File: feature_extraction_matching.cpp, compile with flag #FEATURE_EXTRACT set
Inputs:
   1. flag for keypoint type (OpenCV flag)
   2. flag for descriptor type (OpenCV flag)
   3. ultrasound sequence name (*.txt, text file listing images and uncalibrated probe coordinate matrices)
   4. image mask (*.png)
Outputs:
   1. keypoint file (*.surf.key)
   2. descriptor file (*.surf.dsc)
      
Command line example: 
  ./feature_extraction.exe SURF SURF sweep_2c.txt mask.png
 Produces files
     sweep_2c.key : list of extracted keypoint coordinates
     sweep_2c.dsc : list of keypoint descriptors

Note: sweep_2c.txt is a text file describing an ultrasound sequence, tab separated, where each line lists
 <US image filename> <frame index> <4x4 probe coordinate matrix>

MNI\bite-us\group1\02\pre\sweep_2c\2d\2c.2dus.00001sm.png     1       0.122571591958704       0.165615503867502       -0.00555814942305405    164.86931874045607      -0.0408304890414928     0.0234881769734075      -0.200644221855862      -60.794919764371159     -0.160581528742823      0.120448203638904       0.0467753872252562      -111.09545269274065     0       0       0       1
MNI\bite-us\group1\02\pre\sweep_2c\2d\2c.2dus.00002sm.png     2       0.122719820963848       0.16551148268963        -0.00538234853132442    164.824508910182        -0.0405278595156183     0.023526388336327       -0.200701092975177      -61.113228064347624     -0.160544976594485      0.120583653487478       0.0465514257792677      -111.21025587215171     0       0       0       1
...


2) Feature Matching: identify correspondences between two US sequences of the same tissue/scene
File: feature_extraction_matching.cpp, compile with flag #FEATURE_MATCHING set
Inputs:
   1. training ultrasound sequence 1 (*.txt *.surf.key *.surf.dsc in same directory)
   2. unrelated/background ultrasound sequence (*.txt *.surf.key *.surf.dsc in same directory)
   4. output directory name
   3. fitting  ultrasound sequence 2 (*.mha *.surf.key *.surf.dsc in same directory)

Outputs:
   Feature correspondence coordinates and associated probe coordinate matrices, in header(*.txt)/data(*.bin) format
   1. coordinates in sequence 1 (*.txt.non-linear.matches.X.txt, *.txt.non-linear.matches.X.bin)
   2. coordinates in sequence 2 (*.txt.non-linear.matches.Y.txt, *.txt.non-linear.matches.Y.bin)

Example: identify correspondences between sweep_2c and sweep_2a
  > ./feature_matching sweep_2c.txt unrelated.US.txt output sweep_2a.txt 
 Produces files:
     sweep_2c.txt.non-linear.matches.X.txt, sweep_2c.txt.non-linear.matches.X.bin
     sweep_2c.txt.non-linear.matches.Y.txt, sweep_2c.txt.non-linear.matches.Y.bin

3) Auto-calibration: estimate 4x4 calibration matrix from correspondences
File: autocalibration.cpp
Inputs:
   1. correspondence / probe coordinate file (*.txt) 

Example: estimate calibration from correspondences between sequences sweep_2c and sweep_2a
  > ./autocalibration sweep_2c.txt 
 Reads in 4 files sweep_2c.txt.non-linear.matches.(X/Y).(.txt/.bin), outputs calibration matrix.


Known Issues
 Recent versions of OpenCV package SURF feature extraction/matching in a separate library 'nonfree'.  
 The following code can be added for compilation:
    #include <opencv2/nonfree/nonfree.hpp> 
    ... 
    cv::initModule_nonfree();

Libraries/Dependencies
 http://ab-initio.mit.edu/wiki/index.php/NLopt
 http://opencv.org/

Sample MNI Data:
 www.matthewtoews.com/projects/MNI-BITE.group1_02.zip

