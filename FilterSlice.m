function FilterSlice(SliceIndex)

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

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

% If scan type is not supported.
if ~sum(strcmp(SliceList{SliceIndex}.ScanType,{'CSPAMM-LINES','SPAMM-LINES','MICSR'}))
    errordlg('No filter available for this scan type!');
    return;
end

nx = SliceList{SliceIndex}.ReconResolution(1);
ny = SliceList{SliceIndex}.ReconResolution(2);
NumDyn = length(SliceList{SliceIndex}.dynamic);

for dyn_indx = 1:NumDyn
    
    omega    = SliceList{SliceIndex}.dynamic{dyn_indx}.Omega;
    fltrshft = SliceList{SliceIndex}.dynamic{dyn_indx}.Filter_Shift;
    ratios   = SliceList{SliceIndex}.dynamic{dyn_indx}.Ratios;
    rotation = SliceList{SliceIndex}.dynamic{dyn_indx}.Rotation*pi/180;
    decay    = SliceList{SliceIndex}.dynamic{dyn_indx}.Decay;
    
    FILTER_T = fftshift(design_filter(omega+fltrshft,ratios,decay,rotation,nx,ny));
    FILTER_C = fftshift(design_filter(-(omega+fltrshft),ratios,decay,rotation,nx,ny));
    
    frames   = size(SliceList{SliceIndex}.dynamic{dyn_indx}.Data,1);
    TLC      = SliceList{SliceIndex}.dynamic{dyn_indx}.TLC;
    NXY      = SliceList{SliceIndex}.NXY;
    
    HrtData_T  = zeros ([frames NXY(1) NXY(2)]);
    HrtData_C  = zeros ([frames NXY(1) NXY(2)]);
    
    for phs_indx = 1:frames
        A = ifft2(FILTER_T.*(fft2(squeeze(SliceList{SliceIndex}.dynamic{dyn_indx}.Data(phs_indx,:,:)))));
        B = ifft2(FILTER_C.*(fft2(squeeze(SliceList{SliceIndex}.dynamic{dyn_indx}.Data(phs_indx,:,:)))));
        HrtData_T(phs_indx,:,:) = A(floor(TLC(1)+(1:NXY(1))),floor(TLC(2)+(1:NXY(2))));
        HrtData_C(phs_indx,:,:) = B(floor(TLC(1)+(1:NXY(1))),floor(TLC(2)+(1:NXY(2))));
    end
    
    %-SS- 01/28/10
    % Ignore scalings and assume f(x) is real and w>0
    % f(x) --FT--> F(v)
    % f(x)cos(wx) --FT--> F(v-w)+F(v+w)
    % (F(v-w)+F(v+w)).*FILTER_T --FTinverse--> f(x)e^{-jwx}
    % FilteredPeaks{1} = f(x)e^{-jwx}
    % FilteredPeaks{2} = f(x)e^{+jwx}
    % Phase = -jwx
    % Magnitude = f(x)
    
    % When we are done with all dynamics, wait until a few other things
    % are completed.
    
    SliceList{SliceIndex}.dynamic{dyn_indx}.FilteredPeaks{1} = HrtData_T;
    SliceList{SliceIndex}.dynamic{dyn_indx}.FilteredPeaks{2} = HrtData_C;
    if((dyn_indx == 1))
        SliceList{SliceIndex}.HarpPhase1 = angle((HrtData_T + conj(HrtData_C)));
        %SliceList{SliceIndex}.HarpPhase1      = angle(sqrt(HrtData_T.*conj(HrtData_C)));
        %SliceList{SliceIndex}.HarpPhase1      = wrap(0.5*(wrap(angle(HrtData_T)-angle(HrtData_C))));
        %SliceList{SliceIndex}.HarpPhase1      = wrap(0.5*((angle(HrtData_T)-angle(HrtData_C))));
        SliceList{SliceIndex}.HarpMagnitude = 0.5*(abs(HrtData_T)+abs(HrtData_C));
        SliceList{SliceIndex}.AvailResults = ...
            AddToAvailableResults(SliceList{SliceIndex}.AvailResults,[{'HarpMagnitude'};{'HarpPhase1'}]);
        
    elseif((dyn_indx == 2))
        SliceList{SliceIndex}.HarpPhase2 = angle((HrtData_T + conj(HrtData_C)));
        %SliceList{SliceIndex}.HarpPhase2= angle(sqrt(HrtData_T.*conj(HrtData_C)));
        %SliceList{SliceIndex}.HarpPhase2= wrap(0.5*(wrap(angle(HrtData_T)-angle(HrtData_C))));
        SliceList{SliceIndex}.HarpMagnitude = 0.5*(0.5*(abs(HrtData_T)+abs(HrtData_C))+SliceList{SliceIndex}.HarpMagnitude);
        SliceList{SliceIndex}.AvailResults = ...
            AddToAvailableResults(SliceList{SliceIndex}.AvailResults,{'HarpPhase2'});
    end
end

SliceList{SliceIndex}.ProcessInfo.FilterStatus.flag = 'clean';
setappdata(0,'SliceList',SliceList);


function MyList = AddToAvailableResults(MyList, New2Add)

New2AddLength = length(New2Add);

for i = 1:New2AddLength,
    OldLength = length(MyList);
    found = false;
    for j = 1:OldLength,
        if strcmp(MyList(j),New2Add(i))
            found = true;
            break;
        end
    end
    if ~found
        MyList(OldLength+1) = New2Add(i);
    end
end
