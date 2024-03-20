function [va,grad] = GetMeshLValue2(HImage,VImage,RealPoints,CalcGrad, pixelSize,Mask_recon)

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

[Ny,Nx]=size(HImage);
NPoints=size(RealPoints,2);
va=zeros(2,NPoints);

grad=zeros(2,2,NPoints);
FloorPoint=fix(RealPoints);
%FloorPoint = round(RealPoints);
dPoint=RealPoints-FloorPoint;


FloorPoint(1,FloorPoint(1,:)>Ny-1)=Ny-2;
FloorPoint(2,FloorPoint(2,:)>Nx-1)=Nx-2;

FloorPoint(1,FloorPoint(1,:)<2)=2;
FloorPoint(2,FloorPoint(2,:)<2)=2;


Index=FloorPoint(1,:)+(FloorPoint(2,:)-1)*Ny;
IndexPlusX=FloorPoint(1,:)+(FloorPoint(2,:))*Ny;
IndexPlusY=FloorPoint(1,:)+1+(FloorPoint(2,:)-1)*Ny;
IndexPlusXY=FloorPoint(1,:)+1+(FloorPoint(2,:))*Ny;


at1=-wrap(HImage(IndexPlusY )-HImage(Index))...
    -wrap(HImage(IndexPlusX )-HImage(Index))...
    +wrap(HImage(IndexPlusXY)-HImage(Index));
bt1=wrap(HImage(IndexPlusX)-HImage(Index));
ct1=wrap(HImage(IndexPlusY)-HImage(Index));

va(1,:)=wrap(at1.*dPoint(1,:).*dPoint(2,:)+...
    bt1.*dPoint(2,:)+...
    ct1.*dPoint(1,:)+...
    HImage(Index));

at2=-wrap(VImage(IndexPlusY)-VImage(Index))...
    -wrap(VImage(IndexPlusX)-VImage(Index))...
    +wrap(VImage(IndexPlusXY)-VImage(Index));
bt2=wrap(VImage(IndexPlusX)-VImage(Index));
ct2=wrap(VImage(IndexPlusY)-VImage(Index));

va(2,:)=wrap(at2.*dPoint(1,:).*dPoint(2,:)+...
    bt2.*dPoint(2,:)+...
    ct2.*dPoint(1,:)+...
    VImage(Index));

if CalcGrad
    
    grad(1,1,:)=at1.*dPoint(1,:)+bt1;
    grad(1,2,:)=at1.*dPoint(2,:)+ct1;
    
    grad(2,1,:)=at2.*dPoint(1,:)+bt2;
    grad(2,2,:)=at2.*dPoint(2,:)+ct2;

    % probably this gradient needs to be symmteric..., but we need the bilinear
    % interpolation anyways. 
    
else
    grad=[];
end
