function Planes3D=Slice3d_info(SliceInfo)

%------------------------------------------------------------------
% -SS-
% From Philips Document:
% orientation = SliceList{i}.Orientation
% orientation: 1,2,3
% preparation direction = SliceList{1}.Planes3D.PrpDir
% preparation direction: RL,AP,FH
%
% Patient coordintae system is fixed, right handed. 
% MPS slice system is right/left handed. The chart shows relation
% of these two systems. Origin is fixed in both systems. They are
% just rotated. 
%
% PAT vs XYZ(magnet system righthanded): RL-->x AP-->y FH-->z 
% 
% orientation prep.dir  M   P   S  MPS
% ----------- -------- --- --- --- ---
% transversal    AP    -RL -AP  FH  R
%                RL    -AP -RL  FH  L
% sagittal       AP     FH -AP -RL  L
%                FH    -AP  FH -FL  R
% coronal        RL     FH -RL  AP  L
%                FH    -RL  FH  AP  R
% 
% Angulation(ap,fh,rl)[degr]:
%
% Xmps = Mori.Mrot.Xpat
%
% Mrot = Rfh.Rap.Rrl is responsible for rotation
% Rrl  = [1,        0,       0]
%        [0,  cos(rl), sin(rl)]
%        [0, -sin(rl), cos(rl)]
%
% Rap  = [ cos(ap), 0, sin(ap)]
%        [       0, 1,       0]
%        [-sin(ap), 0, cos(ap)]
%
% Rfh  = [ cos(fh), sin(fh), 0]
%        [-sin(fh), cos(fh), 0]
%        [       0,       0, 1]
% Mori is responsible for orientation of the two systems
% for example: Orientation = 1 & Prep.dir = Ap
%              Mori = [-1 0 0;0 -1 0;0 0 1]
%
% rot_3d: converts MPS to PAT
% Xpat = inv(Mrot).inv(Mori).Xmps
% mat_rot = inv(Mrot) = inv(Rrl).inv(Rap).inv(Rfh)
% mat_ori = inv(Mori)
%
% inv(Rrl) = [1,       0,       0]
%            [0, cos(rl),-sin(rl)]
%            [0, sin(rl), cos(rl)]
%
%mat_rot = 
%[                         cos(ap)*cos(fh),                        -cos(ap)*sin(fh),        -sin(ap)]
%[-sin(rl)*sin(ap)*cos(fh)+cos(rl)*sin(fh), sin(rl)*sin(ap)*sin(fh)+cos(rl)*cos(fh),-sin(rl)*cos(ap)]
%[ cos(rl)*sin(ap)*cos(fh)+sin(rl)*sin(fh),-cos(rl)*sin(ap)*sin(fh)+sin(rl)*cos(fh), cos(rl)*cos(ap)]
%
% -SS-
%------------------------------------------------------------------

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

tempA=(SliceInfo.PxlSpacing).*(SliceInfo.ReconResolution);

r = -tempA(1)/2;t=tempA(1)/2;b=-tempA(1)/2;l=tempA(1)/2;

ini_norm = [0 ,0 ,1]';
ini_e1   = [-1,0 ,0]';
ini_e2   = [0 ,-1,0]';
origin   = [tempA/2-SliceInfo.PxlSpacing 0]'; % Origing = [0,0]

lt       = [l t 0]'; % lt = [1,1]
rt       = [r t 0]'; % rt = [1,128]
rb       = [r b 0]'; % rb = [128,1]
lb       = [l b 0]'; % lb = [128,128]

ang.ap = SliceInfo.Angulation(1); 
ang.fh = SliceInfo.Angulation(2);
ang.rl = SliceInfo.Angulation(3);

center.ap = SliceInfo.OffCentre(1); 
center.fh = SliceInfo.OffCentre(2); 
center.rl = SliceInfo.OffCentre(3);

Planes3D.center = [center.rl;center.ap;center.fh]; 

orientation     = SliceInfo.Orientation; 
prpdir          = SliceInfo.PrpDir; %preparation dir

Planes3D.norm   = rot_3d(ang,orientation,prpdir,ini_norm);
Planes3D.e1     = rot_3d(ang,orientation,prpdir,ini_e1);
Planes3D.e2     = rot_3d(ang,orientation,prpdir,ini_e2);

Planes3D.origin = rot_3d(ang,orientation,prpdir,origin) + Planes3D.center;

Planes3D.lt     = rot_3d(ang,orientation,prpdir,lt) + Planes3D.center;
Planes3D.rt     = rot_3d(ang,orientation,prpdir,rt) + Planes3D.center;
Planes3D.rb     = rot_3d(ang,orientation,prpdir,rb) + Planes3D.center;
Planes3D.lb     = rot_3d(ang,orientation,prpdir,lb) + Planes3D.center;

Planes3D.pnt    = Planes3D.origin;
Planes3D.wide   = norm(Planes3D.origin - Planes3D.rt);
Planes3D.length = norm(Planes3D.origin - Planes3D.lb);
%Planes3D.e1     = (Planes3D.rt - Planes3D.lt)/norm(Planes3D.rt - Planes3D.lt);
%Planes3D.e2     = (Planes3D.lb - Planes3D.lt)/norm(Planes3D.lb - Planes3D.lt);

return;

% --------------------------------------------------------------------
%                     rot_3d
% --------------------------------------------------------------------

function rot_v = rot_3d(ang,ori,dir,v)
    rad = pi/180;

% -SS-
% Old code:
%     sx = sin( -ang.rl * rad);
%     sy = sin( -ang.ap * rad);
%     sz = sin( -ang.fh * rad);
%     cx = cos( -ang.rl * rad);
%     cy = cos( -ang.ap * rad);
%     cz = cos( -ang.fh * rad);
%     
%     matrix(1,1) = cy * cz;
%     matrix(2,1) = -sz * cx + sx * sy * cz;
%     matrix(3,1) = sx * sz + sy * cx * cz;
% 
%     matrix(1,2) = sz * cy;
%     matrix(2,2) = cx * cz + sx * sy * sz;
%     matrix(3,2) = -sx * cz + sy * sz * cx;
% 
%     matrix(1,3) = -sy;
%     matrix(2,3) = sx * cy;
%     matrix(3,3) = cx * cy;
% this way of finding matrix uses -ang.ap I don't know why.
% gives us another view at the image

% New code:
    sx = sin( ang.rl * rad);
    sy = sin( ang.ap * rad);
    sz = sin( ang.fh * rad);
    cx = cos( ang.rl * rad);
    cy = cos( ang.ap * rad);
    cz = cos( ang.fh * rad);
    
    matrix(1,1) =  cy*cz;
    matrix(2,1) = -sx*sy*cz + cx*sz;
    matrix(3,1) =  cx*sy*cz + sx*sz;
    
    matrix(1,2) = -cy*sz;
    matrix(2,2) =  sx*sy*sz + cx*cz;
    matrix(3,2) = -cx*sy*sz + sx*cz;
    
    matrix(1,3) = -sy;
    matrix(2,3) = -sx*cy;
    matrix(3,3) = cx*cy;

% -SS- 
    
    mat_rot = matrix;


    switch ori
        case 1
            switch dir
                case 'Anterior-Posterior'
                    matrix = [-1 0 0;0 -1 0;0 0 1];
                case 'Right-Left'
                    matrix = [0 -1 0;-1 0 0;0 0 1];
                otherwise
                    msg('ms_geometry_slice_to_apat:illegal view_axis')
            end
        case 2
            switch dir
                case 'Anterior-Posterior'
                    matrix = [0 0 -1;0 -1 0;1 0 0];
                case 'Feet-Head'
                    matrix = [0 0 -1;-1 0 0;0 1 0];
                otherwise
                    msg('ms_geometry_slice_to_apat:illegal view_axis')
            end
        case 3
            switch dir
                case 'Right-Left'
                    matrix = [0 -1 0;0 0 1;1 0 0];
                case 'Feet-Head'
                    matrix = [-1 0 0;0 0 1;0 1 0];
                otherwise
                    msg('ms_geometry_slice_to_apat:illegal view_axis')
            end
        otherwise
            msg('ms_geometry_slice_to_apat:illegal view_axis')
    end
    
    mat_ori = matrix;
    rot_v = mat_rot * mat_ori * v;
return;

