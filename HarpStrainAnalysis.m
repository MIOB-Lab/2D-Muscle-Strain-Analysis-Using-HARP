%% The purpose of this code is to extract a ROI from the first principal strain map acquired from the HARP software. You first need to do the analysis using the HARP technique, save the results, load it into the workspace and then execute this code. This code analyzes the first principal strain but you can modify it to analyze any strain you want
%% Code should be executed block by block, specially the block of segementation. You can run upto line 39 at one go, segment the muscle of interest and save as mentioned, then run the rest of the code
%% The last block is for overlaying the colormap onto the initial scan. This is not necessary if you are happy to just see the colormap alone and not superimposed onto the scan
%% This code was written by Divya Pradip Roy (divya-pradip.roy@uvm.edu)

%% Inputs: Add values to the variables as needed 
n= 12 ; % Time frame of interest
I= dicomread("C:\Users\droy6\OneDrive - University of Vermont\Desktop\Semesters\Spring 2024\Research\Session 14 Feb 2024\Tagging\Axial\A1 Magnitude\IM_0133.dcm"); % Selecting the directory of the image where analysis will be done
Nxy=[150 106]; % Plug in the size of the ROI from the software
ROI_Location= [48 169]; % Plug in the location of the topmost corner of the ROI from the software

%% Calculations using the inputs
rows_to_add_before= ROI_Location(2)-1;  %For Zero-padding
cols_to_add_before= ROI_Location(1)-1;    % For Zero-padding
rows_to_add_after= 480- Nxy(1)- rows_to_add_before  % For Zero-padding
cols_to_add_after= 480 - Nxy(2)- cols_to_add_before;  % For Zero-padding
xcords= [ROI_Location(1),ROI_Location(1)+Nxy(2), ROI_Location(1) + Nxy(2), ROI_Location(1)];            % X co-ordinates for the mask  
ycords= [ROI_Location(2), ROI_Location(2), ROI_Location(2)+ Nxy(1),ROI_Location(2)+ Nxy(1)];            % Y co-ordinates for the mask


%% Visualizing the image 
figure
imshow(I,[])

%% Extracting the strain value of that frame in a 2D matrix and zero-padding it to make it the same size as the image
P1= SliceList{1,1}.Strain2DInfoSet.lambda(n,:,:,1); % If second principal strain is wanted, replace 1 with 2 at the end
P1=squeeze(P1);
P1= padarray(P1,[rows_to_add_before,cols_to_add_before],'pre');
P1= padarray(P1,[rows_to_add_after,cols_to_add_after],'post');

%% Creating a mask similar to the ROI chosen in the software and visualizing the cropped image
mask= roipoly(I,xcords,ycords);
ROI_I=I;
ROI_I(~mask)=0;
figure
imshow(ROI_I,[])

%% Segmenting the cropped image to get the muscle of interest location. Save the logical mask as BW and the image as VL_Segmented
imageSegmenter(ROI_I)
I(~BW)=0; % Zero out everything in the image except the muscle of interest and show it
figure
imshow(I,[])

%% Showcase P1 as it is and then crop everything from it other that the strain values of the muscle of interest
figure
imagesc(P1,[-0.5,0.5])
colorbar
P1(~BW)=0;
figure
imagesc(P1,[-0.5,0.5])
colorbar

%% Getting numbers: Calculating strain statistics
[r,c]= find(BW==1);
Strain= zeros(size(r));
for i=1:length(r)
    Strain(i,1)= P1(r(i,1),c(i,1));
end

Max_Strain_Before_Normalization= max(Strain(:));
Min_Strain_Before_Normalization= min(Strain(:));
Mean_Strain_Before_Normalization= mean(Strain(:));
STD_Strain_Before_Normalization= std(Strain(:));

%% Overlaying colormap on the scan
I_O= dicomread("C:\Users\droy6\OneDrive - University of Vermont\Desktop\Semesters\Spring 2024\Research\Session 14 Feb 2024\Tagging\Axial\A1 Magnitude\IM_0133.dcm");
figure
ax1=axes;
imagesc(I_O);
colormap(ax1,'gray')  %Upto this part we just visualized the scan onto which we are going to overlay the colormap

ax2=axes;
imagesc(ax2,P1,"AlphaData",BW>0); %This makes only the mask area visible in the plot and everything else transparent
colormap(ax2,'jet') % Just plotting a jet colormap initially which will be changed to RWB manually after plotting
caxis(ax2,[-0.5 0.5]);
ax2.Visible= 'off'; % This switces off everthing in the colormap except the area defined by the mask, which is our area of interest 
linkprop([ax1 ax2], 'Position');
colorbar