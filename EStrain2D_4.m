function EStrain2D_4(SliceIndex)

% This function calculates 2D strain from filtered tagged images.
% FP1, FP2: Filtered Peaks - FP1: first dynamic, FP2: second dynamic
% HM: Harp Magnitude
% ReconeRes: Number of Pixels in rows and cols of the images 
% PxlSpacing: Pixel Spacing
% Center: Center of ROI 
% Omega1, Omega2: Omega of each tag direction (position of peaks) - Omega1:
% first dynamic, Omega2: second dynamic
% Scan Type: type of scan. Strain calculation only implemented for
% zHarp,d-zharp and TrueHarp for Harp4 software.
% RemoveOffset: if 1 considers first phase as reference phase and subtracts
% if from all following phases. In case of zHarp this flag needs to be one
% otherwise calculations are wrong.
% NLDiffSmooth: enables using NLD Smoothing
% HP3: HarpPhase3. If available, Ecz,Erz,Exz,Eyz strain will be computed.
% zEnc: Z_Encoding, needed if HP3 is available. 

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
    return
end;

ScanType = SliceList{SliceIndex}.ScanType;
if strcmp(ScanType,'TrueHARP (5 DNY)')
    FP1{1} = SliceList{SliceIndex}.dynamic{3}.FilteredPeaks{1};
    FP1{2} = SliceList{SliceIndex}.dynamic{5}.FilteredPeaks{1};
    Omega1 = SliceList{SliceIndex}.dynamic{3}.Omega;
    FP2{1} = SliceList{SliceIndex}.dynamic{2}.FilteredPeaks{1};
    FP2{2} = SliceList{SliceIndex}.dynamic{4}.FilteredPeaks{1};
    Omega2 = SliceList{SliceIndex}.dynamic{4}.Omega;
elseif strcmp(ScanType,'TrueHARP (6 DNY)')
    FP1{1} = SliceList{SliceIndex}.dynamic{3}.FilteredPeaks{1};
    FP1{2} = SliceList{SliceIndex}.dynamic{5}.FilteredPeaks{1};
    Omega1 = SliceList{SliceIndex}.dynamic{3}.Omega;
    FP2{1} = SliceList{SliceIndex}.dynamic{4}.FilteredPeaks{1};
    FP2{2} = SliceList{SliceIndex}.dynamic{6}.FilteredPeaks{1};
    Omega2 = SliceList{SliceIndex}.dynamic{4}.Omega;
else % for zHarp and d-zHarp

    FP1    = SliceList{SliceIndex}.dynamic{1}.FilteredPeaks; % from first dynamic
    FP2    = SliceList{SliceIndex}.dynamic{2}.FilteredPeaks; % from second dynamic
    Omega1 = SliceList{SliceIndex}.dynamic{1}.Omega;
    Omega2 = SliceList{SliceIndex}.dynamic{2}.Omega;
    % FilteredPeaks{1} is at +omega and {2} is at -omega
end
HM         = SliceList{SliceIndex}.HarpMagnitude;
ReconRes   = SliceList{SliceIndex}.ReconResolution;
PxlSpacing = SliceList{SliceIndex}.PxlSpacing; 
CENTER     = SliceList{SliceIndex}.ROI_CNTR; %[row col]

RemoveOffset = 1;
NLDiffSmooth = 0;

%(FP1,FP2,HM,ReconRes,PxlSpacing,CENTER,Omega1,Omega2,RemoveOffset,NLDiffSmooth,HP3,zEnc)

FOV=PxlSpacing(1)*ReconRes(1);
% assumes same pixel spacing and resolution in both rows and cols
OMEGA_S = [Omega1*2*pi/FOV Omega2*2*pi/FOV]; %2x2 matrix

% OMEGA_S = [a 0] % Horizontal
%         = [0 a] % Vertical

% CALCULATION OF GRADIENTS
% @: this sign is used as \rho
% z12: @dyn1/@cols
% z11: @dyn1/@rows
% z22: @dyn2/@cols
% z21: @dyn2/@rows
[z12,z11] = gradHARP_4(FP1,OMEGA_S(:,1),PxlSpacing,RemoveOffset); %[col row]
[z22,z21] = gradHARP_4(FP2,OMEGA_S(:,2),PxlSpacing,RemoveOffset); %[col row]

% m32 & m31: @HarpPhase3
if isfield(SliceList{SliceIndex},'HarpPhase3')
    if ~isfield(SliceList{SliceIndex},'z_encoding')
        disp('Choose z-encoding for the slice');
        return;
    end
    HP3 = SliceList{SliceIndex}.HarpPhase3;
    zENC = SliceList{SliceIndex}.z_encoding;
    [m32,m31] = gradZHARP_4(HP3,zENC,PxlSpacing,RemoveOffset);
end

% Modified by SS 12/20/09
% Never checked. Not sure if it works
No_Frames = size(z12,1);

if(NLDiffSmooth)
    for frame=1:No_Frames
        if exist('m31','var')
        [z12_t,z11_t,z22_t,z21_t,m32_t,m31_t]=NLDiffSmooth(squeeze(z12(frame,:,:)),squeeze(z11(frame,:,:)),...
                                               squeeze(z22(frame,:,:)),squeeze(z21(frame,:,:)),...
                                               squeeze(m32(frame,:,:)),squeeze(m31(frame,:,:)),...
                                               squeeze(HM(frame,:,:)));
        m32(frame,:,:)=m32_t;
        m31(frame,:,:)=m31_t;
        else
        [z12_t,z11_t,z22_t,z21_t]=NLDiffSmooth(squeeze(z12(frame,:,:)),squeeze(z11(frame,:,:)),...
                                               squeeze(z22(frame,:,:)),squeeze(z21(frame,:,:)),...
                                               squeeze(HM(frame,:,:)));
        end 
        z11(frame,:,:)=z11_t;
        z21(frame,:,:)=z21_t;
        z12(frame,:,:)=z12_t;
        z22(frame,:,:)=z22_t;
        pause(1);
    end;
end;

% STRAIN CALCULATION
% A\b = inv(A)*b
omega=(OMEGA_S*OMEGA_S')\OMEGA_S;

NX        = size(z12,2);
NY        = size(z12,3);

m11=ones(size(z12));
m12=zeros(size(z12));
m21=zeros(size(z12));
m22=ones(size(z12));

ymat=reshape(rem((1:(NX*NY))'-1,NY)+1,NY,NX);%rows
%xmat=ymat';                                  %cols

[xmat,ymat] = meshgrid(1:NY,1:NX);
rvec_1=zeros((size(z12)));
rvec_2=zeros((size(z12)));
cvec_1=zeros((size(z12)));
cvec_2=zeros((size(z12)));
%--SS-
%|----------------- x
%|\  
%| \
%|  \
%|   \       ec
%|    \      /
%|theta\ 90 /
%|      \  /
%y     er /
%If Theta is angle with axis y(rows of image)
% x-y and ec-er : right handed systems
%[er]=[cos(theta)  sin(theta)] [ey]=[rows]
%[ec] [-sin(theta) cos(theta)] [ex] [cols]

for frame=1:No_Frames
    rvec_t1=(ymat-CENTER(1,frame))*PxlSpacing(1);    
    rvec_t2=(xmat-CENTER(2,frame))*PxlSpacing(2);    
    rvec_norm=sqrt(rvec_t1.^2+rvec_t2.^2)+1e-10;
    rvec_1(frame,:,:)=rvec_t1./rvec_norm;   % rvec_1 = cos(theta)
    rvec_2(frame,:,:)=rvec_t2./rvec_norm;   % rvec_2 = sin(theta)
    cvec_1(frame,:,:)=-rvec_2(frame,:,:);   % cvec_1 = -sin(theta)
    cvec_2(frame,:,:)=rvec_1(frame,:,:);    % cvec_2 = cos(theta)
end;

%--SS-: june 17-09
%M matrix is Finv = I - du/dx where dx is in deformed coordinate system
%Here what we calculate and call z is actually -du/dx. Refer to HA's
%Strain document for more information.

% -SS- 12/29/09
% Refer to Khaled's Paper in Med Imag Analysis 08
% Phase stores x-ux and by filtering HARP gives -u


m11=m11+omega(1,1)*z11;%Converting angle into distance by * TAG_SEP/(2*pi)
m12=m12+omega(1,1)*z12;
m21=m21+omega(2,1)*z11;
m22=m22+omega(2,1)*z12;

m11=m11+omega(1,2)*z21;%Converting angle into distance by * TAG_SEP/(2*pi)
m12=m12+omega(1,2)*z22;
m21=m21+omega(2,2)*z21;
m22=m22+omega(2,2)*z22;

EC_STRAIN= 1- sqrt((m11.*cvec_1+m12.*cvec_2).^2+(m21.*cvec_1+m22.*cvec_2).^2);
ER_STRAIN= 1- sqrt((m11.*rvec_1+m12.*rvec_2).^2+(m21.*rvec_1+m22.*rvec_2).^2);
EX_STRAIN= 1- sqrt(m11.^2+m21.^2);
EY_STRAIN= 1- sqrt(m12.^2+m22.^2);

if exist('m31','var')
    
    ax= PxlSpacing(1)*ones(size(m31)); % Segment x-length
    cx=-PxlSpacing(1)*m31;             % Segment z-length || The (-)sign is because 
    %the direction of the gradient was opposite to the direction of the
    %element segment.
    %Segment y-length = 0
    ay= PxlSpacing(1)*ones(size(m32));
    cy=-PxlSpacing(1)*m32;
    ECz_STRAIN= 1- sqrt(((m11.*ax.*cvec_1+m12.*ay.*cvec_2).^2+(m21.*ax.*cvec_1+m22.*ay.*cvec_2).^2+(m31.*ax.*cvec_1+m32.*ay.*cvec_2+cx.*cvec_1+cy.*cvec_2).^2)./((ax.*cvec_1).^2+(ay.*cvec_2).^2+(cx.*cvec_1).^2+(cy.*cvec_2).^2+eps));
    ERz_STRAIN= 1- sqrt(((m11.*ax.*rvec_1+m12.*ay.*rvec_2).^2+(m21.*ax.*rvec_1+m22.*ay.*rvec_2).^2+(m31.*ax.*rvec_1+m32.*ay.*rvec_2+cx.*rvec_1+cy.*rvec_2).^2)./((ax.*rvec_1).^2+(ay.*rvec_2).^2+(cx.*rvec_1).^2+(cy.*rvec_2).^2+eps));

    EXz_STRAIN= 1- sqrt(((m11.*ax).^2+(m21.*ax).^2+(m31.*ax+cx).^2)./(ax.^2+cx.^2+eps));
    EYz_STRAIN= 1- sqrt(((m12.*ay).^2+(m22.*ay).^2+(m32.*ay+cy).^2)./(ay.^2+cy.^2+eps));
end

% comput m transpose m: 
% m = Finv
% c = Finv^{T}Finv = Cauchy Deformation Tensor
% e = 1/2(I-c) = Eulerian Strain Tensor
mt11=m11.*m11+m21.*m21;
mt12=m11.*m12+m21.*m22;
%mt21=m12.*m11+m22.*m21;
mt22=m12.*m12+m22.*m22;

% compute the eigen-values and -vectors
ga=sqrt((mt11-mt22).^2+4*mt12.^2)/2;
eig1=(mt11+mt22)/2+ga;
eig2=(mt11+mt22)/2-ga;

vect1_1 = mt12;
vect1_2 = (mt22-mt11)/2+ga;
vect1_norm = sqrt(vect1_1.^2+vect1_2.^2)+1e-10;
vect1_1 = vect1_1./vect1_norm;
vect1_2 = vect1_2./vect1_norm;

%--SS-
%Comment by SS: June 17-09
%lambda1 and lambda2 are the singular values of 
%F = I + du/dX. here we first find Finv = I - du/dx, then we find
%Eigen values of Finv'*Finv and that's why lambda1 = 1./sqrt(). 
%Principle values are sorted from large to small.

lambda1 = eig2+eps;  % min eigen values of Finv (eig2)
lambda(:,:,:,1)=1./(sqrt(lambda1))-1; %larger Sig-Val of F
%lambda(:,:,:,1)=(1-lambda1)/2; %larger eig-Val of e
lambda2 = eig1+eps;  % max eigen values of Finv (eig1)
lambda(:,:,:,2)=1./(sqrt(lambda2))-1; %smaller Sig-Val of F
%lambda(:,:,:,2)=(1-lambda2)/2; %smaller eig-Val of e

% lambda1 and lambda2 are singular values of F minus "1". 
% Eigen values of Eulerian strain tensor are 
% (1-eig2)/2 and (1-eig1)/2

% [vect1_1,vect1_2] corresponds to eig1 corresponds to lambda 2
% [vect2_1,vect2_2] corresponds to eig2 corresponds to lambda 1

%Storing Eigen Vectors
%First eigen vector is calculated ([vect1_1 vect1_2])
%To find second eigen vector, firest eigen vector is rotated by 90 degrees
eigvec(:,:,:,1,2) = vect1_1; %2nd Eig-Vect corresponding to smaller Eig-Val
eigvec(:,:,:,2,2) = vect1_2; %[row col]
eigvec(:,:,:,1,1) = -vect1_2;%1st Eig-Vect corresponding to larger Eig-Val
eigvec(:,:,:,2,1) = vect1_1; %[row col]

% first eigen vector: [vect1_1 vect1_2] = [row col]
%angle of first eigen vector and circumferential direction;
angle(:,:,:,1) = asin(sin(acos((-vect1_2).*cvec_1+vect1_1.*cvec_2))); 
%angle of second eigen vector and circumferential direction;
angle(:,:,:,2) = asin(sin(acos(vect1_1.*cvec_1+vect1_2.*cvec_2))); 

Strain2DInfoSet.EC_STRAIN  = EC_STRAIN;
Strain2DInfoSet.ER_STRAIN  = ER_STRAIN;
Strain2DInfoSet.EX_STRAIN  = EX_STRAIN;
Strain2DInfoSet.EY_STRAIN  = EY_STRAIN;
Strain2DInfoSet.lambda     = lambda;
Strain2DInfoSet.eigvec     = eigvec;
Strain2DInfoSet.angle      = angle;  % angle with circum strain

if exist('ECz_STRAIN','var')
    Strain2DInfoSet.ECz_STRAIN = ECz_STRAIN;
    Strain2DInfoSet.ERz_STRAIN = ERz_STRAIN;
    Strain2DInfoSet.EXz_STRAIN = EXz_STRAIN;
    Strain2DInfoSet.EYz_STRAIN = EYz_STRAIN;
end

SliceList{SliceIndex}.Strain2DInfoSet = Strain2DInfoSet;
setappdata(0,'SliceList',SliceList);
return;
