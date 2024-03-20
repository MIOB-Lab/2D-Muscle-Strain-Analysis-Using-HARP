function ReorganizeDICOMFiles(prefileString)

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

initialPATH = pwd;

display('Please choose the directory to import data from.');
dicomDIR = uigetdir('./', 'Choose the directory to IMPORT data from ..');
if dicomDIR==0
    display('Canceled. Program now exit.');
    return;
end

display('Please choose the directory to export data to.');
saveDIR = uigetdir('./', 'Choose the directory to EXPORT data to ...');
if saveDIR==0
    display('Canceled. Program now exit.');
    return;
end

% dicomDIR = '/Users/xiaofengliu/Desktop/source';
% saveDIR = '/Users/xiaofengliu/Desktop/target';

display(' I am working. Pleasse do NOT exit... ....');
tic;
if nargin == 0
    prefileString = 'IM';
end

parseFileName(dicomDIR, saveDIR, prefileString);
toc

cd(initialPATH);

display('Done.');

% ------------------------------------------------
function parseFileName(dicomDIR, saveDIR, prefileString)

content = dir(dicomDIR);

for n = 1:length(content)
	if strcmp(content(n).name, '.') || strcmp(content(n).name, '..')
        continue;
    end

	if content(n).isdir
    	% if this is a folder, check contents of the folder
        folder = [dicomDIR, filesep, content(n).name];
        parseFileName(folder, saveDIR, prefileString);
    else
        % if it is a file, check whether it is a dicom file.
        % If not, do nothing
        % if yes, save it according to the header information.
        try
        	filename = [dicomDIR, filesep, content(n).name];
            saveFile(filename, saveDIR, prefileString);

        catch
            warning(['The file ', filename, ' CANNOT be saved!']);
        end
    end
end

% --------------------------------------------------------
function saveFile(dicomFile, saveDIR, prefileString)

initialPATH = pwd;
info = dicominfo(dicomFile);

patientName    = info.PatientName.FamilyName;
studyName      = info.StudyDescription;
seriesName     = info.SeriesDescription;
seriesNumber   = info.SeriesNumber;
instanceNumber = info.InstanceNumber;

cd(saveDIR)

% 1 enter first level folder
content  = dir;
found = false;
studyName = RemoveSlashFromString(studyName);
for n = 1:length(content)
    if strcmp(content(n).name, studyName)
        found = true;
        break;
    end
end

if ~found
   mkdir(studyName);
end
cd(studyName);

% 2, enter second level folder

folder_name = strcat(sprintf('%04d', seriesNumber), '_', seriesName);

content  = dir;
found = false;
folder_name = RemoveSlashFromString(folder_name);
for n = 1:length(content)
    if strcmp(content(n).name, folder_name)
        found = true;
        break;
    end
end

if ~found
   mkdir(folder_name);
end
cd(folder_name);

if ~isempty(patientName)
    % If more than one patient
    % 3, enter third level folder
    content     = dir;
    found       = false;
    patientName = RemoveSlashFromString(patientName);
    for n = 1:length(content)
        if strcmp(content(n).name, patientName)
            found = true;
            break;
        end
    end

    if ~found
        mkdir(patientName);
    end
    cd(patientName);
end

% 4, create the file name

targetFile = [prefileString, '-', sprintf('%04d', seriesNumber), '-', sprintf('%04d', instanceNumber), '.dcm'];

copyfile(dicomFile, targetFile);

% 5, done, go back to the old path

cd(initialPATH);

% remove '/' from string because it may create trouble when used for folder
% name in windows
function str = RemoveSlashFromString(str)

for n = 1:length(str)
    if str(n) == '/'  || str(n) == '\' || str(n) == ':'
        str(n) = '-';
    end
end

