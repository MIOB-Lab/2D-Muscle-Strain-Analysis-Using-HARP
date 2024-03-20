function plane3D = Slice3d_info_dicom(SliceInfo)

% 03/03/10
% This code is originaly taken from XL's harp3d
% Modified by SS to be used in HARP4

% generate the plan infomation from dicom file
% SliceInfo is the information of the slice from
% dicominfo function

% Sliceinfo = dicominfo(dicomfile);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Software: 
% 		   Harmonic Phase (HARP)
% 
% 	  Image Analysis and Communications Lab (IACL)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 
% 
% This HARP software is intended for research and educational purposes
% only. You agree to use the HARP software according to the terms and
% conditions outlined in the HARP Disclaimer and End User License
% Agreement (HARP-EULA). In brief, this software has NOT been cleared by
% the Food and Drug Administration, and it may NOT be used for clinical
% or diagnostic purposes.  This software implements only the most basic
% functions of the HARP method and may not be useful for more advanced
% acquisition modes, unusual tag patterns, or advanced tag processing
% goals.  Users are permitted to modify this software for research and
% educational purposes provided that the research and core software is
% properly cited or acknowledged as provided herein and in the
% HARP-EULA.
% 
% JHU owns the entire right, title, and interest in and to the "HARP"
% trademark and the patents related to the HARP software. JHU has
% exclusively licensed the HARP software to Diagnosoft, Inc.  A
% commercial, FDA cleared product is available from Diagnosoft,
% Inc. Permission for IACL and Johns Hopkins University to distribute
% this software for research educational purposes has been granted by
% Diagnosoft, Inc. However, the software may not be redistributed and no
% commercial use can be based on this software.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disclaimer: 
% 
% THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
% APPLICABLE LAW. THE COPYRIGHT HOLDERS AS WELL AS JOHNS HOPKINS
% UNIVERSITY AND DIAGNOSOFT, INC. PROVIDE THE PROGRAM AS IS WITHOUT
% WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED. THE ENTIRE RISK AS
% TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE
% PROGRAM PROVE DEFECTIVE, YOU ASSUME ALL THE COSTS AND RESPONSIBILITIES
% THAT MIGHT RESULT FROM ITS USE.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contributors:
% 
% Fangxu Xing, fxing1@jhu.edu
% Sahar Soleimanifard
% Khaled Z. Abd-Elmoniem
% Harsh K. Agarwal
% Yuanming Suo
% Xiaofeng Liu
% Smita Sampath
% Vijay Parthasarathy
% Jonghye Woo
% Aaron Carass
% Nael F. Osman
% Jerry L. Prince, prince@jhu.edu
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acknowledgement: 
% 
% If you use this HARP software or methods derived from this software to
% process or produce any data for publication, please refer to the
% following publications:
% 
% N. F. Osman, W. S. Kerwin, E. R. McVeigh, and J. L. Prince, "Cardiac
% Motion Tracking Using CINE Harmonic Phase (HARP) Magnetic Resonance
% Imaging", Mag. Res. Med., vol. 42, pp. 1048-1060, 1999.
% 
% N. F. Osman, E. R. McVeigh, and J. L. Prince, "Imaging Heart Motion
% Using Harmonic Phase MRI", IEEE Trans. on Medical Imaging, vol. 19,
% No. 3, pp. 186-202, March 2000.
% 
% X. Liu and J.L. Prince, "Shortest path refinement for motion
% estimation from tagged MR images," IEEE Trans Med Imaging, vol. 29,
% no. 8, pp.1560-1572, March 2010.
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More information: 
% 
% IACL HARP page - http://www.iacl.ece.jhu.edu/static/harp/
% 
% Wikipedia HARP page - http://en.wikipedia.org/wiki/HARP_algorithm
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: 
% 
% Copyright (c) 1999-2012 Johns Hopkins University, Image Analysis
% and Communications Lab (IACL)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plane3D.e1 = SliceInfo.ImageOrientationPatient(4:6); % the column direction
plane3D.e1 = plane3D.e1/norm(plane3D.e1);
plane3D.e2 = SliceInfo.ImageOrientationPatient(1:3); % the row direction
plane3D.e2 = plane3D.e2/norm(plane3D.e2);

plane3D.pixelSpacing = SliceInfo.PixelSpacing;       % mm/pixel
plane3D.rows = double(SliceInfo.Rows);               % #of rows
plane3D.cols = double(SliceInfo.Columns);            % #of columns

plane3D.origin = SliceInfo.ImagePositionPatient - ...
                 plane3D.pixelSpacing(1)*plane3D.e1-plane3D.pixelSpacing(2)*plane3D.e2;
             
plane3D.pnt    = plane3D.origin;  % one point on the plane3D
plane3D.norm   = cross(plane3D.e1, plane3D.e2);  % the normal direction
plane3D.center = plane3D.origin + (plane3D.cols/2+0.5)*plane3D.pixelSpacing(1)*plane3D.e1 ...
                            + (plane3D.rows/2+0.5)*plane3D.pixelSpacing(2)*plane3D.e2;

% (1,1), left top corner in mm 
plane3D.lt = plane3D.origin + 1*plane3D.pixelSpacing(1)*plane3D.e1 ...
                        + 1*plane3D.pixelSpacing(2)*plane3D.e2;
% (128, 1), right top corner in mm 
plane3D.lb = plane3D.origin + plane3D.rows*plane3D.pixelSpacing(1)*plane3D.e1 ...
                        + 1*plane3D.pixelSpacing(2)*plane3D.e2;
% (1, 128), left bottom in mm
plane3D.rt = plane3D.origin + 1*plane3D.pixelSpacing(1)*plane3D.e1 ...
                        + plane3D.cols*plane3D.pixelSpacing(2)*plane3D.e2;
% (128 128), right bottom corner in mm
plane3D.rb = plane3D.origin + plane3D.rows*plane3D.pixelSpacing(1)*plane3D.e1 ...
                        + plane3D.cols*plane3D.pixelSpacing(2)*plane3D.e2;

% image width  in mm
plane3D.width = plane3D.cols*plane3D.pixelSpacing(1);
% image height in mm
plane3D.height = plane3D.rows*plane3D.pixelSpacing(2);

% image plane3D parameters, =[a b c d] such that image function is:
% ax+by+cz+d=0
plane3D.para = [plane3D.norm', -dot(plane3D.pnt, plane3D.norm)];