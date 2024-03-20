function [output_data,Info_data] = siemens_Dicom_Load_Data(pathname,filein,ScanType)

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

% --------------------------------------------
% output_data.slice is SliceList. Info_data is original dicom header.
% Initialization

% First File of First Dir needed for initialization
info = dicominfo([pathname{1} filein{1,1}]);

% NumberOfK == 2 is SPAMM data
% NumberOfK == 4 is CSPAMM data
[NumberOfK,NumberOfFiles] = size(filein);

% find number of phases
NumberOfPhases = 1;
if isfield(info,'CardiacNumberOfImages')
    NumberOfPhases = info.CardiacNumberOfImages;
end
NumberOfSlices = NumberOfFiles/NumberOfPhases;

% initialization
Info_data = cell(NumberOfSlices,1); % original dicom header
slice = cell(NumberOfSlices,1); % SliceList
SliceOrigin = [NaN;NaN;NaN];
PhsIndex = ones(NumberOfSlices,1);
wb = waitbar(0,'Loading DICOM data...');
clear info;

for FileIndex = 1:NumberOfFiles
    
    waitbar(FileIndex/NumberOfFiles,wb);
    
    Image = [];
    for KIndex = 1:NumberOfK
        info{KIndex} = dicominfo([pathname{KIndex} [filein{KIndex,FileIndex}]]);
        Image(:,:,KIndex) = double(dicomread([pathname{KIndex} [filein{KIndex,FileIndex}]])); %#ok<AGROW>
    end
    
    % found is to find an image with the same slice location (aka another time frame)
    found = false;
    SlcIndex = 1;
    for NoSlc = 1:size(SliceOrigin,3)
        if info{1}.ImagePositionPatient == SliceOrigin(:,:,NoSlc)
            SlcIndex = NoSlc - 1;
            found = true;
            break;
        end
    end
    
    if ~found % then this is a new slice
        SlcIndex = size(SliceOrigin,3);
        SliceOrigin(:,:,end+1) = info{1}.ImagePositionPatient;     %#ok<AGROW>       
        slice{SlcIndex}.ReconResolution(1,1) = double(info{1}.Rows);
        slice{SlcIndex}.ReconResolution(1,2) = double(info{1}.Columns);        
        slice{SlcIndex}.Thickness = double(info{1}.SliceThickness);
        slice{SlcIndex}.PxlSpacing = double(info{1}.PixelSpacing');
        
        Info_data{SlcIndex} = info{1};
    end
    
    if NumberOfK == 1 % High-Res Single Phase with no RR interval
        
        slice{SlcIndex}.dynamic{1}.phase{PhsIndex(SlcIndex)}.magnitude = Image(:,:,1);
        if isfield(info{1},'TriggerTime')
            slice{SlcIndex}.dynamic{1}.ttime(PhsIndex(SlcIndex)) = info{1}.TriggerTime;
        else
            slice{SlcIndex}.dynamic{1}.ttime(PhsIndex(SlcIndex)) = 0; %joker
        end
        slice{SlcIndex}.dynamic{1}.freq(PhsIndex(SlcIndex))  = 1; %joker
        
    else
        if strcmpi(ScanType,'CSPAMM-LINES') % CSPAMM image            
            if NumberOfK == 8 % 2 dynamics
                cmpIm1a = Image(:,:,1).*exp(sqrt(-1)*(Image(:,:,2)-2048)*pi/2048); % reconstruct complex image
                cmpIm1b = Image(:,:,3).*exp(sqrt(-1)*(Image(:,:,4)-2048)*pi/2048); % reconstruct complex image
                DynImage{1}  = real(cmpIm1a - cmpIm1b);
                cmpIm2a = Image(:,:,5).*exp(sqrt(-1)*(Image(:,:,6)-2048)*pi/2048); % reconstruct complex image
                cmpIm2b = Image(:,:,7).*exp(sqrt(-1)*(Image(:,:,8)-2048)*pi/2048); % reconstruct complex image
                DynImage{2}  = real(cmpIm2a - cmpIm2b);
            elseif NumberOfK == 4 % 1 dynamic
                cmpIm1a = Image(:,:,1).*exp(sqrt(-1)*(Image(:,:,2)-2048)*pi/2048); % reconstruct complex image
                cmpIm1b = Image(:,:,3).*exp(sqrt(-1)*(Image(:,:,4)-2048)*pi/2048); % reconstruct complex image
                DynImage{1}  = real(cmpIm1a - cmpIm1b);
            end
            % In philips there are two min and max RR interval but in siemens
            % it is only one nominal interval
            slice{SlcIndex}.dynamic{1}.minRRintrvl(PhsIndex(SlcIndex)) = min(info{1}.NominalInterval,info{3}.NominalInterval);
            slice{SlcIndex}.dynamic{1}.maxRRintrvl(PhsIndex(SlcIndex)) = max(info{1}.NominalInterval,info{3}.NominalInterval);
            
            if NumberOfK == 8 % 2 dynamics
                slice{SlcIndex}.dynamic{2}.minRRintrvl(PhsIndex(SlcIndex)) = min(info{5}.NominalInterval,info{7}.NominalInterval);
                slice{SlcIndex}.dynamic{2}.maxRRintrvl(PhsIndex(SlcIndex)) = max(info{5}.NominalInterval,info{7}.NominalInterval);
            end
            
        elseif strcmpi(ScanType,'MICSR') % MICSR image
            if NumberOfK == 4 % 2 dynamics
                DynImage{1}  = (Image(:,:,1) - Image(:,:,2)).*(Image(:,:,1) + Image(:,:,2));
                DynImage{2}  = (Image(:,:,3) - Image(:,:,4)).*(Image(:,:,3) + Image(:,:,4));
            elseif NumberOfK == 2 % 1 dynamic
                DynImage{1}  = (Image(:,:,1) - Image(:,:,2)).*(Image(:,:,1) + Image(:,:,2));
            end
            % In philips there are two min and max RR interval but in siemens
            % it is only one nominal interval
            slice{SlcIndex}.dynamic{1}.minRRintrvl(PhsIndex(SlcIndex)) = min(info{1}.NominalInterval,info{2}.NominalInterval);
            slice{SlcIndex}.dynamic{1}.maxRRintrvl(PhsIndex(SlcIndex)) = max(info{1}.NominalInterval,info{2}.NominalInterval);
            
            if NumberOfK == 4 % 2 dynamics
                slice{SlcIndex}.dynamic{2}.minRRintrvl(PhsIndex(SlcIndex)) = min(info{3}.NominalInterval,info{4}.NominalInterval);
                slice{SlcIndex}.dynamic{2}.maxRRintrvl(PhsIndex(SlcIndex)) = max(info{3}.NominalInterval,info{4}.NominalInterval);
            end
            
        elseif strcmpi(ScanType,'SPAMM-LINES') % SPAMM image
            DynImage{1}  = Image(:,:,1);
            DynImage{2}  = Image(:,:,2);
            
            slice{SlcIndex}.dynamic{1}.minRRintrvl(PhsIndex(SlcIndex)) = info{1}.NominalInterval;
            slice{SlcIndex}.dynamic{1}.maxRRintrvl(PhsIndex(SlcIndex)) = info{1}.NominalInterval;
            
            slice{SlcIndex}.dynamic{2}.minRRintrvl(PhsIndex(SlcIndex)) = info{2}.NominalInterval;
            slice{SlcIndex}.dynamic{2}.maxRRintrvl(PhsIndex(SlcIndex)) = info{2}.NominalInterval;
        end
        slice{SlcIndex}.dynamic{1}.phase{PhsIndex(SlcIndex)}.magnitude = DynImage{1};
        if length(DynImage) == 2
            slice{SlcIndex}.dynamic{2}.phase{PhsIndex(SlcIndex)}.magnitude = DynImage{2};
        end
        
        if isfield(info{1},'TriggerTime') % need check down here
            slice{SlcIndex}.dynamic{1}.ttime(PhsIndex(SlcIndex)) = info{1}.TriggerTime;
            if length(DynImage) == 2
                slice{SlcIndex}.dynamic{2}.ttime(PhsIndex(SlcIndex)) = info{2}.TriggerTime;
            end
        else
            slice{SlcIndex}.dynamic{1}.ttime(PhsIndex(SlcIndex)) = 0;
            if length(DynImage) == 2
                slice{SlcIndex}.dynamic{2}.ttime(PhsIndex(SlcIndex)) = 0;
            end
        end
        slice{SlcIndex}.dynamic{1}.freq(PhsIndex(SlcIndex))  = 1; %joker
        if length(DynImage) == 2
            slice{SlcIndex}.dynamic{2}.freq(PhsIndex(SlcIndex))  = 1; %joker
        end
    end
    
    PhsIndex(SlcIndex) = PhsIndex(SlcIndex) + 1;
end

close(wb);
output_data.slice = slice;