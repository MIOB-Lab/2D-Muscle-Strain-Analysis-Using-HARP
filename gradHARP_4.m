function [gradx,grady]=gradHARP_4(CINE, OMEGA_S, PXL_WDTH,RMVOFST,NLD_SMOTH,FILTER_T)

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

    No_Peaks  = length(CINE);
    No_Frames = size(CINE{1},1);
    NX        = size(CINE{1},2);
    NY        = size(CINE{1},3);

    RemoveOffsetFirst = RMVOFST;
    
    gradient_factor   = 1;
    %ymat              = reshape(rem((1:(NX*NY))'-1,NY)+1,NY,NX)
    [xmat,ymat] = meshgrid(1:NY,1:NX);
    % ymat: rows

    harp0  = 0;
    harp01 = 0;
    harp02 = 0;

    gradx  = zeros(size(CINE{1}));
    grady  = zeros(size(CINE{1}));

    
    %-SS- 01/28/10
    %CINE{1} = f(X)e^{-jwX} %true
    %CINE{2} = f(X)e^{+jwX} %conj
    
    %angle(CINE{1}) = wX = wx + w(X-x) = wx -ux
    %x: deformed coordinate X: undeformed coordinate
    %shifting = -wx
    %harp+shifting = -ux
    % --> same thing for the other peak: harp2-shifting
    % --> gradient farctor = 2 (results from two peaks)
    
    shifting=-ymat*PXL_WDTH(1)*OMEGA_S(1)-xmat*PXL_WDTH(2)*OMEGA_S(2);

    if(RemoveOffsetFirst)
        if (No_Peaks>=2)
            harp01=squeeze(angle(CINE{1}(1,:,:)))+shifting;
            harp02=squeeze(angle(CINE{2}(1,:,:)))+flipud(fliplr(shifting));
        else
            harp0 =squeeze(angle(CINE{1}(1,:,:)))+shifting;
        end;
    end;

    for frame=1:No_Frames,
       
        if (No_Peaks>=2)
            harp1       = squeeze(angle(CINE{1}(frame,:,:)))-harp01;
            harp2       = squeeze(angle(CINE{2}(frame,:,:)))-harp02;
            Mag         = 0.5*squeeze(abs(CINE{1}(frame,:,:)) + abs(CINE{2}(frame,:,:)));

            shifted_harp= wrap((harp1+shifting)-(harp2+flipud(fliplr(shifting)))); % shifted harp
            % Above line is similar to
            % shifted_harp = wrap(wrap(harp1+shifting)-wrap(harp2-shifting));
            
            % Added by -HA-
            if(exist('FILTER_T','var'))
                shifted_harp = angle(ifft2(FILTER_T.*fft2(Mag.*exp(sqrt(-1)*shifted_harp))));
            end
            gradient_factor=2;
        else
            harp=squeeze(angle(CINE{1}(frame,:,:)))-harp0;

            shifted_harp=wrap(harp+shifting); % shifted harp
            Mag         = 0.5*squeeze(abs(CINE{1}(frame,:,:)));
            shifted_harp = angle(ifft2(FILTER_T.*fft2(Mag.*shifted_harp)));
        end;
        % Added by -HA-
        usePolynomialDerivativeEstimation = 0;
        if (usePolynomialDerivativeEstimation && (frame ==7))
            dx_shifted_harp(1:NY,1:(NX-1))=wrap(shifted_harp(:,2:(NX))-shifted_harp(:,1:(NX-1)));
            dx_shifted_harp(1:NY,NX)=dx_shifted_harp(1:NY,NX-1);
            dy_shifted_harp(1:(NY-1),1:NX)=wrap(shifted_harp(2:(NY),:)-shifted_harp(1:(NY-1),:));
            dy_shifted_harp(NY,1:NX)=dx_shifted_harp(NY-1,1:NX);
            
            numOfPoints = 5; % odd number
            orderOfPolynomial = 3;
            [X Y] = meshgrid(0:orderOfPolynomial,0:orderOfPolynomial);
            PowerMatrix = [X(:) Y(:)];
            Nexponent = size(PowerMatrix,1);
            clear X Y

            [X Y] = meshgrid(1:numOfPoints,1:numOfPoints);
            X = X - (numOfPoints+1)/2;X = X(:).';
            Y = Y - (numOfPoints+1)/2;Y = Y(:).';
            J = (repmat(X,[Nexponent 1]).^repmat(PowerMatrix(:,1),[1 numOfPoints*numOfPoints])).*...
                (repmat(Y,[Nexponent 1]).^repmat(PowerMatrix(:,2),[1 numOfPoints*numOfPoints]));
            for iy = ((numOfPoints+1)/2 + 1):(NY-(numOfPoints+1)/2)
                for ix = ((numOfPoints+1)/2 + 1):(NX-(numOfPoints+1)/2)
                    IX = X(1:numOfPoints:end)+ix;
                    IY = Y(1:numOfPoints)+iy;
                    mag = squeeze(abs(CINE{1}(frame,IY,IX)))*0+1;
                    W = diag(exp((X.^2 + Y.^2)/(2*400)).*(mag(:)'));
                    A = inv(J*W*J')*J;
                    A = reshape(A(find(sum(PowerMatrix,2) == 1),:),[2 numOfPoints numOfPoints]);
                    
                    dx_shifted_harp(ix,iy) = sum(sum(shifted_harp(IY,IX).*squeeze(A(1,:,:))));
                    dy_shifted_harp(ix,iy) = sum(sum(shifted_harp(IY,IX).*squeeze(A(2,:,:))));
                end
            end
        else
            %-Suo- the NX and NY in following section is interchanged.
            dx_shifted_harp(1:NX,1:(NY-1))=wrap(shifted_harp(:,2:(NY))-shifted_harp(:,1:(NY-1)));
            dx_shifted_harp(1:NX,NY)=dx_shifted_harp(1:NX,NY-1);
            dy_shifted_harp(1:(NX-1),1:NY)=wrap(shifted_harp(2:(NX),:)-shifted_harp(1:(NX-1),:));
            dy_shifted_harp(NX,1:NY)=dx_shifted_harp(NX-1,1:NY);
        end
        
        %Added by -HA-
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform 1D smoothing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        oneDSmoothing =0;
        if (oneDSmoothing)
            filterN = 3;
%             oneDSmoothingFilter = ones(2*filterN+1,1);
            oneDSmoothingFilter=exp(-(-filterN:filterN).^2/(0.4163*filterN)^2)';
%            oneDSmoothingFilter=cos((-filterN:filterN)'/filterN*pi/4);

            dx_shifted_harp = dx_shifted_harp.*squeeze(abs(CINE{1}(frame,:,:)));
            dx_shifted_harp = conv2(dx_shifted_harp,oneDSmoothingFilter.');
            dx_shifted_harp = dx_shifted_harp(:,filterN +(1:NY));
            normaliseFilter = squeeze(abs(CINE{1}(frame,:,:)));
            normaliseFilter = conv2(normaliseFilter,oneDSmoothingFilter.');
            normaliseFilter = normaliseFilter(:,filterN +(1:NY));
            dx_shifted_harp = dx_shifted_harp./(eps+normaliseFilter);
            
            dy_shifted_harp = dy_shifted_harp.*squeeze(abs(CINE{2}(frame,:,:)));
            dy_shifted_harp = conv2(dy_shifted_harp,oneDSmoothingFilter);
            dy_shifted_harp = dy_shifted_harp(filterN +(1:NX),:);            
            normaliseFilter = squeeze(abs(CINE{2}(frame,:,:)));
            normaliseFilter = conv2(normaliseFilter,oneDSmoothingFilter);
            normaliseFilter = normaliseFilter(filterN +(1:NX),:);
            dy_shifted_harp = dy_shifted_harp./(eps+normaliseFilter);
        end
        
        gradx(frame,:,:)=dx_shifted_harp/PXL_WDTH(1)/gradient_factor;
        grady(frame,:,:)=dy_shifted_harp/PXL_WDTH(2)/gradient_factor;
    end;
return;
