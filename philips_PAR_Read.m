function PAR= philips_PAR_Read(idir,filein)

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

PatNamStr='.    Patient name                       :';%   ROT_PHNT
ExmNamStr='.    Examination name                   :';%   ROT_PHNT
PrtNamStr='.    Protocol name                      :';%   EPI_CSPM_WRK CLEAR
ExmDatStr='.    Examination date/time              :';%   2005.07.08 / 18:41:59
SrsTypStr='.    Series Type                        :';%   Image   MRSERIES
AcqNumStr='.    Acquisition nr                     :';%   30
RcnNumStr='.    Reconstruction nr                  :';%   1
ScnDurStr='.    Scan Duration [sec]                :';%   25.2
CrdPhsStr='.    Max. number of cardiac phases      :';%   20
NumEcoStr='.    Max. number of echoes              :';%   1
NumSlcStr='.    Max. number of slices/locations    :';%   1
NumDynStr='.    Max. number of dynamics            :';%   1
NumMixStr='.    Max. number of mixes               :';%   1
PntPosStr='.    Patient position                   :';%   Head First Supine
PrpDirStr='.    Preparation direction              :';%   Anterior-Posterior
TchniqStr='.    Technique                          :';%   FFE
ScnResStr='.    Scan resolution  (x, y)            :';%   192  10
RecResStr='.    Recon resolution (x, y)            :';%   256 256
ScnMdeStr='.    Scan mode                          :';%   MS
RepTmeStr='.    Repetition time [ms]               :';%   30.00  
FOVAFRStr='.    FOV (ap,fh,rl) [mm]                :';%   300.00  8.00  300.00
WtrFtSStr='.    Water Fat shift [pixels]           :';%   0.73
AngAFRStr='.    Angulation midslice(ap,fh,rl)[degr]:';%   0.00  0.00  0.00
OfCntrStr='.    Off Centre midslice(ap,fh,rl) [mm] :';%   -84.55  17.60  -18.18
FlwCmpStr='.    Flow compensation <0=no 1=yes> ?   :';%   0
PreSatStr='.    Presaturation     <0=no 1=yes> ?   :';%   0
PhsEncStr='.    Phase encoding velocity [cm/sec]   :';%   0.00  0.00  0.00
MTCYesStr='.    MTC               <0=no 1=yes> ?   :';%   0
SPRYesStr='.    SPIR              <0=no 1=yes> ?   :';%   0
EPI_NoStr='.    EPI factor        <0,1=no EPI>     :';%   1
DynScnStr='.    Dynamic scan      <0=no 1=yes> ?   :';%   0
DffYesStr='.    Diffusion         <0=no 1=yes> ?   :';%   0
DffEcoStr='.    Diffusion echo time [ms]           :';%   0.00


fid=fopen([idir filein]);
PAR=cell(0);


while 1
	tline = fgetl(fid);
	if ~ischar(tline), break, end

    if(findstr(PatNamStr,tline))
    	 PAR.PatNam = sscanf(tline,strcat(PatNamStr,'%s'));
     elseif(findstr(ExmNamStr,tline))
    	 PAR.ExmNam = sscanf(tline,strcat(ExmNamStr,'%s'));
     elseif(findstr(PrtNamStr,tline))
    	 PAR.PrtNam = sscanf(tline,strcat(PrtNamStr,'%s'));
     elseif(findstr(ExmNamStr,tline))
    	 PAR.ExmNam = sscanf(tline,strcat(ExmNamStr,'%s'));
     elseif(findstr(ExmDatStr,tline))
    	 PAR.ExmDat = sscanf(tline,strcat(ExmDatStr,'%d.%d.%d / %d:%d:%d'))';
     elseif(findstr(SrsTypStr,tline))
    	 PAR.SrsTyp = sscanf(tline,strcat(SrsTypStr,'%s %s %s %s'));
     elseif(findstr(AcqNumStr,tline))
    	 PAR.AcqNum = sscanf(tline,strcat(AcqNumStr,'%d'));
     elseif(findstr(RcnNumStr,tline))
    	 PAR.RcnNum = sscanf(tline,strcat(RcnNumStr,'%d'));
     elseif(findstr(ScnDurStr,tline))
    	 PAR.ScnDur = sscanf(tline,strcat(ScnDurStr,'%f'));
     elseif(findstr(CrdPhsStr,tline))
    	 PAR.CrdPhs = sscanf(tline,strcat(CrdPhsStr,'%d'));
     elseif(findstr(NumEcoStr,tline))
    	 PAR.NumEco = sscanf(tline,strcat(NumEcoStr,'%d'));
     elseif(findstr(NumSlcStr,tline))
    	 PAR.NumSlc = sscanf(tline,strcat(NumSlcStr,'%d'));
     elseif(findstr(NumDynStr,tline))
    	 PAR.NumDyn = sscanf(tline,strcat(NumDynStr,'%d'));
     elseif(findstr(NumMixStr,tline))
    	 PAR.NumMix = sscanf(tline,strcat(NumMixStr,'%d'));
     elseif(findstr(PntPosStr,tline))
    	 PAR.PntPos = sscanf(tline,strcat(PntPosStr,'%s %s %s %s'));
     elseif(findstr(PrpDirStr,tline))
    	 PAR.PrpDir = sscanf(tline,strcat(PrpDirStr,'%s %s %s %s'));
     elseif(findstr(TchniqStr,tline))
    	 PAR.Tchniq = sscanf(tline,strcat(TchniqStr,'%s %s %s %s'));
     elseif(findstr(ScnResStr,tline))
    	 PAR.ScnRes = sscanf(tline,strcat(ScnResStr,'%d','%d'))';
     elseif(findstr(RecResStr,tline))
    	 PAR.RecRes = sscanf(tline,strcat(RecResStr,'%d','%d'))';
     elseif(findstr(ScnMdeStr,tline))
    	 PAR.ScnMde = sscanf(tline,strcat(ScnMdeStr,'%s %s %s %s'));
     elseif(findstr(RepTmeStr,tline))
    	 PAR.RepTme = sscanf(tline,strcat(RepTmeStr,'%f'))';
     elseif(findstr(FOVAFRStr,tline))
    	 PAR.FOVAFR = sscanf(tline,strcat(FOVAFRStr,'%f','%f','%f'))';
     elseif(findstr(WtrFtSStr,tline))
    	 PAR.WtrFtS = sscanf(tline,strcat(WtrFtSStr,'%f'))';
     elseif(findstr(AngAFRStr,tline))
    	 PAR.AngAFR = sscanf(tline,strcat(AngAFRStr,'%f','%f','%f'))';
     elseif(findstr(OfCntrStr,tline))
    	 PAR.OfCntr = sscanf(tline,strcat(OfCntrStr,'%f','%f','%f'))';
     elseif(findstr(FlwCmpStr,tline))
    	 PAR.FlwCmp = sscanf(tline,strcat(FlwCmpStr,'%d'))';
     elseif(findstr(PreSatStr,tline))
    	 PAR.PreSat = sscanf(tline,strcat(PreSatStr,'%d'))';
     elseif(findstr(PhsEncStr,tline))
    	 PAR.PhsEnc = sscanf(tline,strcat(PhsEncStr,'%f','%f','%f'))';
     elseif(findstr(MTCYesStr,tline))
    	 PAR.MTCYes = sscanf(tline,strcat(MTCYesStr,'%d'))';
     elseif(findstr(SPRYesStr,tline))
    	 PAR.SPRYes = sscanf(tline,strcat(SPRYesStr,'%d'))';
     elseif(findstr(EPI_NoStr,tline))
    	 PAR.EPI_No = sscanf(tline,strcat(EPI_NoStr,'%d'))';
     elseif(findstr(DynScnStr,tline))
    	 PAR.DynScn = sscanf(tline,strcat(DynScnStr,'%d'))';
     elseif(findstr(DffYesStr,tline))
    	 PAR.DffYes = sscanf(tline,strcat(DffYesStr,'%d'))';
     elseif(findstr(DffEcoStr,tline))
    	 PAR.DffEco = sscanf(tline,strcat(DffEcoStr,'%f'))';
    end

    if(findstr('=== IMAGE INFORMATION DEFINITION ===',tline))
        index=0;
        Exp1='\(([0-9]*)\**integer\)';
        Exp2='\(([0-9]*)\**float\)';
		while 1
			tline = fgetl(fid);
		    if(findstr('=== IMAGE INFORMATION ===',tline)), break, end;
            [s1,f1,t1]=regexp(tline,Exp1);s=s1;f=f1;data_type='(integer)';
            [s2,f2,t2]=regexp(tline,Exp2);if(isempty(s))s=s2;f=f2;data_type='(float)';end;
            if(isempty(s))continue;end;
            data_length=1;
            if((f(length(f))-s(length(s))+1)>length(data_type))
                data_length=str2num(tline((s(length(s))+1):(f(length(f))-length(data_type))));    
            end;
            index=index+1;
            PAR.InfoNameAndSize{index}.Name=deblank(tline(1:(s(length(s))-1)));
            PAR.InfoNameAndSize{index}.Type=data_type;
            PAR.InfoNameAndSize{index}.Width=data_length;
            PAR.InfoNameAndSize{index}.StartCol=1;
            if(index>1)
                PAR.InfoNameAndSize{index}.StartCol=PAR.InfoNameAndSize{index-1}.StartCol+PAR.InfoNameAndSize{index-1}.Width;
            end;
        end
    end;

    if(findstr('=== IMAGE INFORMATION ===',tline))
        %disp('=== IMAGE INFORMATION ===');
        index=0;
		while 1
			tline = fgetl(fid);
		    if(findstr('END OF DATA DESCRIPTION FILE',tline)), break, end;
            [A,ok]=str2num(tline);
            if(ok==1)
                index=index+1;
                PAR.info(index,:)=A';
            end
        end
    end;
end
fclose(fid);

