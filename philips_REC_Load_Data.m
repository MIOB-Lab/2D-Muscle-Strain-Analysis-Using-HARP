function output_data = philips_REC_Load_Data(idir,filein,par,out_slice_no,out_scan_dynamic,out_cardiac_phase,out_data_type)

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

if(~isfield(par,'info'))
    disp('Data file is empty! Please select another file.');
    return;
end
Xsize=par.info(1,10);
Ysize=par.info(1,11);
PHASES=par.CrdPhs;
NODYN=par.NumDyn;
NOSLICES=par.NumSlc;
coils=1;

NumDataVariables=size(par.info,1)/PHASES/NODYN/NOSLICES;

fidin = fopen([idir filein],'r','ieee-le');

INDEX      = squeeze(par.info(:,[1 3:5]));
INDSCN_TYP = squeeze(par.info(:,6));
INDSCN_TYP=zeros(size(INDSCN_TYP));%Scan_type can hold interesting information
                                   %Rightnow, it is better to just ignore
                                   %it!
SCAN_TYPES = unique(INDSCN_TYP);
NOSCANS    = length(unique(INDSCN_TYP));

par.info_reor=zeros(size(par.info));
found=false;
data_info=0;
%This is the order data is saved
if(nargin>3)
    if(out_scan_dynamic>NODYN)
        disp(sprintf('Requested dynamic not found!. dynamic must be <=%d',NODYN));
        disp(sprintf('Dynamic No. %d will be used instead!',NODYN));
        out_scan_dynamic=NODYN;
    end;    
    if(out_slice_no>NOSLICES)
        disp(sprintf('Requested slice not found!. slice must be <=%d',NOSLICES));
        disp(sprintf('Slice No. %d will be used instead!',NOSLICES));
        out_slice_no=NOSLICES;
    end;    
    if(out_cardiac_phase>PHASES)
        disp(sprintf('Requested cardiac phase not found!. phase must be <=%d',PHASES));
        disp(sprintf('Phase No. %d will be used instead!',PHASES));
        out_cardiac_phase=PHASES;
    end;    
    for slc=1:NOSLICES % First slice then the second then, ....
        for dyn=1:NODYN % First dynamic then the second then, ....
            for data_type=0:4 % Magnitude --> real --> imag --> phase --> recon
                for phase=1:PHASES % Cardiac phase 1 --> 2 --> 3 --> ....
                    LST=(INDEX==ones(NumDataVariables*PHASES*NODYN*NOSLICES,1)*[out_slice_no out_scan_dynamic out_cardiac_phase out_data_type]);  
                    LST=LST(:,1).*LST(:,2).*LST(:,3).*LST(:,4);
                    I=find(LST);
                    if(~isempty(I))
                        fseek(fidin,2*(I-1)*Xsize*Ysize,'bof');
                        output_data = fread(fidin,[Xsize Ysize],'short'); %read the partial k-space
                        data_info=par.info(I,:);
                        output_data=(output_data*data_info(13)+data_info(12))/(data_info(13)*data_info(14));
                        found=true;
                    end
                    if(found)break;end;
                end
                if(found)break;end;
            end
            if(found)break;end;
        end
        if(found)break;end;
    end;
else % load all data in the file
    output_data=cell(0);
    for data_image=1:size(INDEX,1) % First dynamic then the second then, ....
        fseek(fidin,2*(data_image-1)*Xsize*Ysize,'bof');
        out_data=double(fread(fidin,[Xsize Ysize],'short'));
        data_info=par.info(data_image,:);
        out_data=(out_data*data_info(13)+data_info(12))/(data_info(13)*data_info(14));
        switch INDEX(data_image,4)
            case 0
                output_data.slice{INDEX(data_image,1)}.dynamic{INDEX(data_image,2)+NODYN*(find(SCAN_TYPES==INDSCN_TYP(data_image))-1)}.phase{INDEX(data_image,3)}.magnitude=out_data;
            case 1
                output_data.slice{INDEX(data_image,1)}.dynamic{INDEX(data_image,2)+NODYN*(find(SCAN_TYPES==INDSCN_TYP(data_image))-1)}.phase{INDEX(data_image,3)}.real     =out_data;
            case 2
                output_data.slice{INDEX(data_image,1)}.dynamic{INDEX(data_image,2)+NODYN*(find(SCAN_TYPES==INDSCN_TYP(data_image))-1)}.phase{INDEX(data_image,3)}.imag     =out_data;
            case 3
                output_data.slice{INDEX(data_image,1)}.dynamic{INDEX(data_image,2)+NODYN*(find(SCAN_TYPES==INDSCN_TYP(data_image))-1)}.phase{INDEX(data_image,3)}.angle    =out_data;
            case 4
                output_data.slice{INDEX(data_image,1)}.dynamic{INDEX(data_image,2)+NODYN*(find(SCAN_TYPES==INDSCN_TYP(data_image))-1)}.phase{INDEX(data_image,3)}.recon    =out_data;
        end
        output_data.slice{INDEX(data_image,1)}.dynamic{INDEX(data_image,2)}.ttime(INDEX(data_image,3))=data_info(33);
        output_data.slice{INDEX(data_image,1)}.dynamic{INDEX(data_image,2)}.freq(INDEX(data_image,3))=data_info(37);
        output_data.slice{INDEX(data_image,1)}.dynamic{INDEX(data_image,2)}.minRRintrvl(INDEX(data_image,3))=data_info(38);
        output_data.slice{INDEX(data_image,1)}.dynamic{INDEX(data_image,2)}.maxRRintrvl(INDEX(data_image,3))=data_info(39);
        output_data.slice{INDEX(data_image,1)}.ReconResolution  = data_info(10:11);
        output_data.slice{INDEX(data_image,1)}.Angulation       = data_info(17:19);
        output_data.slice{INDEX(data_image,1)}.OffCentre        = data_info(20:22);
        output_data.slice{INDEX(data_image,1)}.Thickness        = data_info(23);
        output_data.slice{INDEX(data_image,1)}.Orientation      = data_info(26);
        output_data.slice{INDEX(data_image,1)}.PxlSpacing       = data_info(29:30);
    end
end

fclose(fidin);