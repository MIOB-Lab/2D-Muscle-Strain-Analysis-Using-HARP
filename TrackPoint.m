function pt_track = TrackPoint(pt, CurrentSlice, CurrentPhase, dim, RefineOrNot)

% This script is written by Xiaofeng Liu.
% dim = 1: 1D HARP tracking. 
% dim = 2: 2D HARP tracking.  This is the default value

% RefineOrNot = 1: use refinement
% RefineOrNot = 0: use traditional HARP
% if this is not specified when calling the function, the function will
% look at SliceList.AdvancedSettings.TrackUseRefinement;

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

SliceList = getappdata(0, 'SliceList');

numPhase = size(SliceList{CurrentSlice}.dynamic{1}.Data, 1);

pt_track = zeros(2,numPhase);
pt_track(1, CurrentPhase) = pt(1);
pt_track(2, CurrentPhase) = pt(2);

if ~exist('dim','var')
    dim = 2;
end

if ~exist('RefineOrNot','var')
    RefineOrNot = SliceList{CurrentSlice}.AdvancedSettings.TrackUseRefinement;
end

if RefineOrNot      % is use refinement
    
    %mx = interp2(squeeze(SliceList{CurrentSlice}.MotionBackwardX(CurrentPhase, :,:)), pt(2), pt(1));
    %my = interp2(squeeze(SliceList{CurrentSlice}.MotionBackwardY(CurrentPhase, :,:)), pt(2), pt(1));
    %pt_ref = pt_track(:, CurrentPhase) + [mx; my];
    
    for phase = CurrentPhase:numPhase-1
        pt_cur = pt_track(:, phase);
        mx = interp2(squeeze(SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardX(phase, :,:)), pt_cur(2), pt_cur(1));
        my = interp2(squeeze(SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardY(phase, :,:)), pt_cur(2), pt_cur(1));
        
        pt_track(:,phase+1) = pt_cur + [mx; my];
    end
    for phase = CurrentPhase:-1:2
        pt_cur = pt_track(:, phase);
        mx = interp2(squeeze(SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardX(phase, :,:)), pt_cur(2), pt_cur(1));
        my = interp2(squeeze(SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardY(phase, :,:)), pt_cur(2), pt_cur(1));
        
        pt_track(:,phase-1) = pt_cur + [mx; my];
    end
      
else % use traditional HARP
    % Track forward
    curPt = pt;
    for phase = (CurrentPhase+1) : numPhase

        hx1 = squeeze(SliceList{CurrentSlice}.HarpPhase1(phase-1, :,:));
        hx2 = squeeze(SliceList{CurrentSlice}.HarpPhase1(phase, :,:));

        if dim == 1
            hy1 = SliceList{CurrentSlice}.RefPhase2;
            hy2 = hy1;
        else
            hy1 = squeeze(SliceList{CurrentSlice}.HarpPhase2(phase-1, :,:));
            hy2 = squeeze(SliceList{CurrentSlice}.HarpPhase2(phase, :,:));
        end
        [new_pt, motion] = TrackPointHarp(curPt, hx1, hx2, hy1, hy2);
        pt_track(:, phase) = new_pt;
        curPt = new_pt;
    end

    % Track backward
    curPt = pt;
    for phase = (CurrentPhase-1) : -1 : 1

        hx1 = squeeze(SliceList{CurrentSlice}.HarpPhase1(phase+1, :,:));
        hx2 = squeeze(SliceList{CurrentSlice}.HarpPhase1(phase, :,:));

        if dim == 1
            hy1 = SliceList{CurrentSlice}.RefPhase2;
            hy2 = hy1;
        else
            hy1 = squeeze(SliceList{CurrentSlice}.HarpPhase2(phase+1, :,:));
            hy2 = squeeze(SliceList{CurrentSlice}.HarpPhase2(phase, :,:));
        end
        [new_pt, motion] = TrackPointHarp(curPt, hx1, hx2, hy1, hy2);
        pt_track(:, phase) = new_pt;
        curPt = new_pt;
    end

end


% ----------------------------------------------------------------
% Function: Track one point using traditional HARP tracking over one frame
function [new_pt, motion] = TrackPointHarp(pt, hx1, hx2, hy1, hy2)

iter = 1;
y = pt;
[ny,nx] = size(hx1);

% the phase value of pt at the first time frame
current_A(1) = getLValue(hx1,pt);
current_A(2) = getLValue(hy1,pt);
h2 = [hx2, hy2];

ephi1 = ones(2,1);
while ((norm(ephi1(:,iter))>.0001)&&(iter<15))
    iter=iter+1;

    m=fix(y);
    dm=y-m;
    if (y(1)>=ny)
        m(1)=ny-1;dm(1)=.99;
    end;

    if (y(1)<1)
        m(1)=1; dm(1)=0;
    end;

    if (y(2)>=nx)
        m(2)=nx-1;dm(2)=.99;
    end;

    if (y(2)<1)
        m(2)=1; dm(2)=0;
    end;

    for k=1:2,
        at=-wrap(h2(m(1)+1,m(2)+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx))-...
            wrap(h2(m(1),m(2)+1+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx))+...
            wrap(h2(m(1)+1,m(2)+1+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx));
        bt=wrap(h2(m(1),m(2)+1+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx));
        ct=wrap(h2(m(1)+1,m(2)+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx));
        va=wrap(at*dm(1)*dm(2)+bt*dm(2)+ct*dm(1)+h2(m(1),m(2)+(k-1)*nx));
        ephi1(k,iter)=wrap(va-current_A(k));
        grad(k,1)=(at*dm(1)+bt);
        grad(k,2)=(at*dm(2)+ct);
    end;

    vt=-inv(grad'*grad)*grad'*ephi1(1:2,iter);

    if (abs(vt(1))>1),
        vt(1)=sign(vt(1));
    end;

    if (abs(vt(2))>1),
        vt(2)=sign(vt(2));
    end;
    y=y+[vt(2);vt(1)];

end % while

new_pt = y;
motion = new_pt - pt;

