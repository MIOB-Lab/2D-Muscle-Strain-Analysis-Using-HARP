function trackGrid = VectorTrack(sliceListForThisSlice, currentFrame, xRange, yRange, frames, ...
    maskForReferenceFrame, NXY, RefineOrNot)

% sliceList: Data Structure for this slice. 
% currentFrame: reference frame
% xRange: ROI in x direction
% yRange: ROI in y direction
% frames: # of frames
% maskForReferenceFrame: Mask for the reference frame
% NX: Number of pixels in the x-direction
% NY: Number of pixels in the y-direction
% RefineOrNot: use refinement or not, added by Xiaofeng Liu
% RefineOrNot = 0: use original vector tracking; other: use refinement. 
% RefineOrNot default value is 0. 

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

% Do the tracking using vectorized tracking (Khaled
% and Ayman) technique. It has been adapted to suit this program.

% Changed by Xiaofeng Liu, to incorporate it with refinement

maxIterations = 10;

% Form the grid on which tracking needs to be done. Only a
% few points on the grid will be tracked because of
% masking.

trackGrid = 0 * ones(length(yRange), length(xRange), frames,2);

trackingDetails = cell(1);
% trackingDetails stores all the information that is
% applicable to the tracking routine, like number of
% tracked points, the co-ordinates of tracked points etc.

clear Y vt trackingDetails trackThisPoint nPoints currentA grad A trackedPointsFinal
clear numTrackedPointsFinal IND1 IND2

% Initialize the number of tracked points
trackingDetails.numTrackedPoints = 0;

% Initialize 2 index  variables
x = 0;
y = 0;

for yROIcoord = yRange
    x=0; y = y+1;
    
    for xROIcoord = xRange
        x = x+1;
        currentY = [yROIcoord;xROIcoord];
        
        % getValue from mask Undeformed
        % because, we are tracking forward
        trackThisPoint = GetValue(maskForReferenceFrame,currentY);
        
        if trackThisPoint ~=0
            % The actual starting location of the point being tracked.
            trackingDetails.pointsRefFrame(:,trackingDetails.numTrackedPoints+1)=currentY;
            
            % Above are the actual coordinates. Do the
            % same thing for the index as well.
            
            trackingDetails.index(:,trackingDetails.numTrackedPoints+1) = [y;x];
            
            % The number of points being tracked
            trackingDetails.numTrackedPoints=trackingDetails.numTrackedPoints+1;
        end
    end
end


% Start the vectorized tracking.

% First initialize the matrix where all the tracked
% points will be stored.

trackingDetails.trackedPoints=zeros([size(trackingDetails.pointsRefFrame) frames]);

% Store the current location of points for the
% reference frame
trackingDetails.trackedPoints(:,:,currentFrame)=trackingDetails.pointsRefFrame;

frameList=[currentFrame+1:frames (currentFrame-1):-1:1];
prevFrame=[currentFrame:(frames - 1) currentFrame:-1:2];

% Set and initialize few more variables that will be
% useful later.
nPoints = trackingDetails.numTrackedPoints;

% currentA = zeros(nPoints,2);
% A = zeros(1,2);
% grad = zeros(2,2);

% added by Xiaofeng
% Now judge whether use the refinement to track or the original vector
% tracking code
if ~exist('RefineOrNot','var')
    RefineOrNot = 0; % default value: use original vector tracking
end

if RefineOrNot  % use refinement
    
    a = sliceListForThisSlice.TLC(1);
    b = sliceListForThisSlice.TLC(2);
    
%     pts_cur = trackingDetails.trackedPoints(:,:,currentFrame);
    %mx = interp2(squeeze(sliceListForThisSlice.MotionBackwardX(currentFrame, :,:)), pts_cur(2,:)+b, pts_cur(1,:)+a);
    %my = interp2(squeeze(sliceListForThisSlice.MotionBackwardY(currentFrame, :,:)), pts_cur(2,:)+b, pts_cur(1,:)+a);
    
    %pts_ref = pts_cur + [mx; my];
    
    for phase = currentFrame : frames-1
        pts_cur =  trackingDetails.trackedPoints(:,:,phase);
        mx = interp2(squeeze(sliceListForThisSlice.MotionForwardX(phase, :,:)), pts_cur(2,:)+b, pts_cur(1,:)+a);
        my = interp2(squeeze(sliceListForThisSlice.MotionForwardY(phase, :,:)), pts_cur(2,:)+b, pts_cur(1,:)+a);
        
        trackingDetails.trackedPoints(:,:,phase+1) = pts_cur + [mx; my];
    end
    for phase = currentFrame : -1 : 2
        pts_cur =  trackingDetails.trackedPoints(:,:,phase);
        mx = interp2(squeeze(sliceListForThisSlice.MotionBackwardX(phase, :,:)), pts_cur(2,:)+b, pts_cur(1,:)+a);
        my = interp2(squeeze(sliceListForThisSlice.MotionBackwardY(phase, :,:)), pts_cur(2,:)+b, pts_cur(1,:)+a);
        
        trackingDetails.trackedPoints(:,:,phase-1) = pts_cur + [mx; my];
    end
    

else    % use the original vector tracking

    K = 2;
    phaseImages= zeros(K,NXY(1), NXY(2), frames);



%     a   = sliceListForThisSlice.TLC(1);
%     b   = sliceListForThisSlice.TLC(2);

    % Fill the HARP images for all the K's and all time frames
    for frm_indx = 1:frames

        %     scope1=(1:NY)+fix((frm_indx-1)/COL)*NY;
        %     scope2=(1:NX)+rem(frm_indx-1,COL)*NX;

        %     for k=1:K
        %         phaseImages(k,:,:,frm_indx)=angle(CINE_I(scope1,scope2+(k-1)*COL*NX));
        %     end;

        phaseImages(1,:,:,frm_indx) = sliceListForThisSlice.HarpPhase1(frm_indx,(1:NXY(1)),((1:NXY(2))));
        phaseImages(2,:,:,frm_indx) = sliceListForThisSlice.HarpPhase2(frm_indx,(1:NXY(1)),((1:NXY(2))));

    end;

    [currentA,grad]=GetMeshLValue2(squeeze(phaseImages(2,:,:,currentFrame)),...
        squeeze(phaseImages(1,:,:,currentFrame)),trackingDetails.trackedPoints(:,:,currentFrame),0);


    for frmIndx = 1:frames-1
        Iteration=1;
        ephi1=ones(2,nPoints);
        Y=squeeze(trackingDetails.trackedPoints(:,:,prevFrame(frmIndx)));

        % while ((norm(sum(ephi1,2)) > .0001) & (Iteration <
        % maxIterations)),

        % Changed above by Vijay so that every tracked point satifies the
        % condition, and not just the overall norm.

        while (~isempty((((sqrt(ephi1(1,:).^2 + ephi1(2,:).^2)) > 0.0001)))) & (Iteration < maxIterations)
            Iteration=Iteration+1;
            [A,grad]=GetMeshLValue2(squeeze(phaseImages(2,:,:,frameList(frmIndx))),...
                squeeze(phaseImages(1,:,:,frameList(frmIndx))),Y,1);
            ephi1=wrap(A-currentA);

            % Compute the step direction
            for n=1:nPoints,
                vt(:,n)=-inv(grad(:,:,n)'*grad(:,:,n))*grad(:,:,n)'*ephi1(:,n);
            end
            vt(vt>1)=1;
            vt(vt<-1)=-1;
            Y=Y+vt(2:-1:1,:);

        end

        trackingDetails.trackedPoints(:,:,frameList(frmIndx))=Y;

    end

end % end for vector tracking or refinement vector tracking

trackedPointsFinal=zeros(trackingDetails.numTrackedPoints*2,frames);
numTrackedPointsFinal=trackingDetails.numTrackedPoints;

% Now all the tracked points are in place
trackedPointsFinal(:,1:frames)=reshape(squeeze(trackingDetails.trackedPoints(:,:,:)),...
    [numTrackedPointsFinal*2 frames]);


% Now the trackedPointsFinal is form,
% -                                -T
% |x1 y1 x2 y2....xN yN for frame 1|
% |x1 y1 x2 y2....xN yN for frame 2|
% -                                -
% Now we need to rearrange it to the original grid.

a1 = trackingDetails.index(1,:);
a2 = trackingDetails.index(2,:);
% a3=ones(size(a1));
% a4=ones(size(a1));

for myFrame = 1:frames
%     
%     scope1=(1:NY)+fix((frm_indx-1)/COL)*NY;
%     scope2=(1:NX)+rem(frm_indx-1,COL)*NX;
%     
    
    IND1=(sub2ind(size(trackGrid),a1,a2,myFrame*ones(size(a1)),1*ones(size(a1))));
    trackGrid(IND1)=reshape(trackingDetails.trackedPoints(1,:,myFrame),size(IND1));
    
    IND2=(sub2ind(size(trackGrid),a1,a2,myFrame*ones(size(a1)),2*ones(size(a1))));
    trackGrid(IND2)=reshape(trackingDetails.trackedPoints(2,:,myFrame),size(IND2));
    
end
return;
