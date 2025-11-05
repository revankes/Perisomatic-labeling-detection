%%%%%%%%%%%%%%%%%%%%%%%%%%%%INTRO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% with the ijm script, we made sure to select one z-stack image, where we
% take the ROI out of the pink channel. These are the cells. Information is then measured in other channels, 
% which cells are PV +/- (blue), mCherry +/- (red) and fos +/- (green). These measurements are added into a dataframe. 
% 
% To measure the PV input around different types of cells we make an
% enlarged the circumference around the Nissl cells, and measure the PV
% signal data.

% different methods of PV measurement are added, such as binarized and
% non-binarized, averaged to the area of the cell and not.


% In this script, we add all the information together from the different
% channels into one data frame. We then sort the cells into mCherry +/- and
% PV +/-, fos+/- whereafter we are able to say something about the PV input on
% these cells. 

% it is important that you only select one cell with your fiji script,
% since the cutoffs are based on single cells, not double cells

%
clear all
close all

addpath(genpath('')) ; %add file location here

%% load files and import data
filenames = {'_Blue', '_BlueLarge', '_BlueThreshold', '_BlueThresholdLarge', '_GREEN','_RED','_PINK'};
fnames = dir('*.txt');
myfiles  = length(fnames); % how many files you have

Finaldata(1,:) = cellstr({'filename', '#Nissl', '#PV+', '#mCherry+',....
        '#mCherry-', '#mCherry- VIC', '#fos+', '#mCh+PV-', 'exp mCh+PV- ', 'mean exp mCh+PV-', 'BINexp mCh+PV-', 'meanBINexp mCh+PV-', ...
        '#mCh+PV-Fos+', 'exp mCh+PV-Fos+', 'mean exp mCh+PV-Fos+', 'BINexp mCh+PV-Fos+', 'meanBINexp mCh+PV-Fos+', ...
        '#mCh+PV-Fos-', 'exp mCh+PV-Fos-', 'mean exp mCh+PV-Fos-' , 'BINexp mCh+PV-Fos-', 'meanBINexp mCh+PV-Fos-', ...
        '#mChVIC-PV-', 'exp mChVIC-PV-', 'mean exp mChVIC-PV-', 'BINexp mChVIC-PV-', 'meanBINexp mChVIC-PV-',  ...
        '#mChVIC-PV-Fos+', 'exp mChVIC-PV-Fos+', 'mean exp mChVIC-PV-Fos+', 'BINexp mChVIC-PV-Fos+', 'meanBINexp mChVIC-PV-Fos+',  ...
        '#mChVIC-PV-Fos-', 'exp mChVIC-PV-Fos-', 'mean exp mChVIC-PV-Fos-', 'BINexp mChVIC-PV-Fos-', 'meanBINexp mChVIC-PV-Fos-', ...
        '#Fos+PV-', 'expFos+PV-', 'mean exp Fos+PV-', 'BINexp Fos+PV-', 'meanBINexp Fos+PV-', ...
        '#Fos-PV-', 'expFos-PV-', 'mean exp Fos-PV-', 'BINexp Fos-PV', 'meanBINexp Fos-PV' }) ; %final data frame

% fill in a 0 for no mCherry present, fill a 1 in if you do have mCherry cells in your image.
mCherry_present = 1;


%% start of loop
rownr = 1 ;

for i = 1:7:myfiles % for 1 to the number of myfiles, with steps of 7

    clear cell_data
    clear DataPVpos
    clear DataPVneg
    clear DatamCherrypos
    clear DatamCherryneg
    clear DatamCherryposPVneg
    clear DatamCherrynegVinPVneg 
    clear datamCherry_neg_vicinity 
    clear dataFos_neg_vicinity
    clear DataFosneg
    clear DataFosnegPVneg
    clear DataFospos
    clear DataFosposPVneg
    clear DatamCherrynegVinPVnegFosneg
    clear DatamCherrynegVinPVnegFospos
    clear DatamCherryposPVnegFosneg
    clear DatamCherryposPVnegFospos


    B = importdata((fnames(i,1).name)); %Blue, blue cells
    BL = importdata((fnames(i+1,1).name)); %BlueLarge,  blue enlarged cells
    BBin = importdata((fnames(i+2,1).name)); %Blue binarized cells
    BLBin = importdata((fnames(i+3,1).name)) ; %BlueLarge, blue enlarged binarized cells
    G = importdata((fnames(i+4,1).name)); %Green, Fos/green cells
    P = importdata((fnames(i+5,1).name)); %Pink, Nissl/pink cells
    R = importdata((fnames(i+6,1).name)); %Red, mCherry/red cells
   
    currFilename=struct2cell(fnames(i,1)); %current file which you are working on, but then in cell form
    currName = currFilename(1,1); % name of the file you're working on 

    rownr = rownr + 1 ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % concatinate data into one dataset
    %%%%%%%%%%%%%%% %%%%%%%%%%%%%%%  

    %general
    cellnr = P.data(:,1) ; %cell number
    areacell = P.data(:,2) ; %area of cell
    xcoordcell = P.data(:,7) ; %x-coordinate cell
    ycoordcell = P.data(:,8) ; %y-coordatine cell
    
    %pink/Nissl
    meanP = P.data(:,3) ; %pink mean 
    expP = meanP .* areacell ; %total exp Pink
    
    %blue/PV
    meanB = B.data(:,3) ; %blue mean 
    expB = meanB .* areacell ; %total exp blue
    percentileB = meanB/(max(meanB)-min(meanB)) ; % percentile cell 
    
    %binarized blue
    meanBin = BBin.data(:,3) ; %blue binarized data 
    expBin = meanBin .* areacell ; %total binarized exp blue
    
    %green/fos
    meanG = G.data(:,3) ; %green mean 
    expG = meanG .* areacell ; %total exp green
    percentileG = meanG/(max(meanG)-min(meanG)) ; % percentile cell 
    
    %red/mCherry
    meanR = R.data (:,3) ; %red mean 
    expR = meanR .* areacell ; %total exp red
    percentileR = meanR/(max(meanR)-min(meanR)); % percentile cell 
    
    %blue large/enlarged PV 
    areacellLarge = BL.data(:,2); %area of enlarged cell
    meanBL = BL.data(:,3); %blue mean enlarged cell
    expBL = meanBL .* areacellLarge ; %total exp blue enlarged
    expBdonut = expBL-expB ; % expression BL - B (so only the part/donut around the cell)
    areadonut = areacellLarge-areacell ;
    meanBdonut = expBdonut./areadonut ; %average expression in donut
    
    %binarized blue large/enlarged PV
    meanBLBin= BLBin.data(:,3); %blue large binarized data
    expBLBin = meanBLBin .* areacellLarge ; % total binarized exp blue
    expBINdonut = expBLBin - expBin ; % expression binarized BL-B (so only the part around the cell)
    meanexpBINdonut = expBINdonut ./areadonut ; % expression binarized and averaged per area
   
    % concatinate data into table
    cell_data = table (cellnr, areacell, xcoordcell, ycoordcell, meanP, expP, meanB, ...
       expB, percentileB, meanBin, expBin, meanR, expR, percentileR, meanG, expG, percentileG, areacellLarge, meanBL, expBL, expBdonut, areadonut, meanBdonut, meanBLBin, expBLBin, expBINdonut,meanexpBINdonut);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sorting cells: mCherry pos/neg, PV pos/neg, Fos pos/neg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PV 
    % because we saw that there was quite a difference in range of
    % brightness, we decided to do a cutoff at the 75% in the 'brightest
    % cell-least brightest' cell range. 
    % !Test this for your image set!
    
    %This cutoff does not work perfectly if you do not have any PV cells!

    thresholdB = 0.75 ; %positive PV cell threshold
    thresholdB_low = 0.50 ; %negative PV cell threshold

    % PV pos cells 
    idxPVpos = cell_data.percentileB >= thresholdB;
    DataPVpos= cell_data(idxPVpos,:);
    
    % PV neg cells
    idxPVneg = cell_data.percentileB < thresholdB_low;
    DataPVneg= cell_data(idxPVneg,:);

    %mCherry 

    % mCherry thresholding is also based on percentiles, 
    % !Test this for your image set!
    
    thresholdRhigh = 0.35 ; % positive mCherry cell threshold
    thresholdRlow  = 0.20 ; % negative mCherry cell threshold

    % if you have no mCherry cells, it won't work that well. Therefore, if
    % mCherry_present = 0 then a mean threshold will be used. If the
    % mCherry_present = 1, then the normal thresholding will be used.

     thresholdR_noMch = 192500.000000000 ; %test this

    if mCherry_present == 1 % in the case of having mCherry cells
         % all mCherry pos cells
        idxmCherrypos = cell_data.percentileR>=thresholdRhigh;
        DatamCherrypos = cell_data(idxmCherrypos,:);
        % all mCherry neg cells
        idxmCherryneg = cell_data.percentileR<=thresholdRlow; 
        DatamCherryneg = cell_data(idxmCherryneg,:);

    elseif mCherry_present == 0 % in the case of having no mCherry cells present
        % all mCherry pos cells
        idxmCherrypos = cell_data.expR>=thresholdR_noMch;
        DatamCherrypos = cell_data(idxmCherrypos,:);
         % all mCherry neg cells
        idxmCherryneg = cell_data.expR<thresholdR_noMch;
        DatamCherryneg = cell_data(idxmCherryneg,:);
    end 

    %Fos

    %again, Fos cells are based on percentiles. 
    % !Test this for your image set!

    thresholdGhigh =0.4 ; % positive Fos cell threshold
    thresholdGlow  =0.25 ; % negative Fos cell threshold

    %all fos pos cells
    idxFospos = cell_data.percentileG >=thresholdGhigh;
    DataFospos = cell_data(idxFospos,:);

    %all fos neg cells
    idxFosneg = cell_data.percentileG <=thresholdGlow;
    DataFosneg = cell_data(idxFosneg,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mCherry positive cells, excl PV pos cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % mCherry cells that are PV negative 
    idxmCherryposPVneg = DatamCherrypos.percentileB<thresholdB;
    DatamCherryposPVneg = DatamCherrypos(idxmCherryposPVneg,:);
    
    totexpmCherryposPVneg= mean(DatamCherryposPVneg.expBdonut); % average total exp blue around mCherry pos cells (but PV neg)
    meanexpmCherryposPVneg = mean (DatamCherryposPVneg.meanBdonut); % average mean exp blue around mChherry pos cells (but PV neg)
    
    binexpmCherryposPVneg = mean (DatamCherryposPVneg.expBINdonut) ; %binarized 
    meanbinexpmCherryposPVneg = mean (DatamCherryposPVneg.meanexpBINdonut) ;% average per area binarized
    
    countmCherryposPVneg= height(DatamCherryposPVneg.cellnr);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mCherry negative cells - only in the vicinity of the mCherry positive cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mCherryneg_vectors = [DatamCherryneg.xcoordcell, DatamCherryneg.ycoordcell] ; %x and y coordinates of mCherry neg cells
    mCherryposPVneg_vectors = [DatamCherryposPVneg.xcoordcell, DatamCherryposPVneg.ycoordcell] ;%x and y coordinates of mCherry pos cells

    vicinityN = 3; %the number of mCherryneg cells you want returned lying in the vicinity of the mCherry pos cells
    
    % Return an array with Distances from 'mCherryneg_vectors' to points in the 'mCherryposPVneg_vectors' sorted 
    % small to large, and keep only the 3 smallest/closest

    [DisttoPos,IndextoPos] = pdist2(mCherryneg_vectors,mCherryposPVneg_vectors,'euclidean','Smallest',vicinityN);
    
    datamCherry_neg_vicinity = DatamCherryneg(IndextoPos, :); %makes a new array with only the mCherry neg cells
    % that are in the vicinity of the mCherry pos cells (this should be 3x the amount of mCherryposPVneg cells)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mCherry neg vicinity cells, excl PV pos cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % all mCherry neg vicinity cells that are PV negative 
    idxmCherrynegPVnegv = datamCherry_neg_vicinity.percentileB<thresholdB;
    DatamCherrynegVinPVneg = datamCherry_neg_vicinity(idxmCherrynegPVnegv,:);

    totexpmCherrynegVinPVneg= mean(DatamCherrynegVinPVneg.expBdonut); % average total exp blue around mCherry negative cells (and PV neg)
    meanexpmCherrynegVinPVneg = mean(DatamCherrynegVinPVneg.meanBdonut); % average mean exp blue around mCherry neg cells (and PV neg)
    
    binexpmCherrynegVinPVneg = mean (DatamCherrynegVinPVneg.expBINdonut) ; %binarized 
    meanbinexpmCherrynegVinPVneg = mean (DatamCherrynegVinPVneg.meanexpBINdonut) ;% average per area binarized
    
    countmCherrynegVinPVneg= height(DatamCherrynegVinPVneg.cellnr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Overlap of cells (PV/mCherry/Fos)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % all mCherry cells that are PV negative, and Fos negative
    idxmCherryposPVnegFosneg = DatamCherryposPVneg.percentileG<=thresholdGlow ;
    DatamCherryposPVnegFosneg = DatamCherryposPVneg(idxmCherryposPVnegFosneg,:);

    totexpmCherryposPVneg_Fosneg=mean(DatamCherryposPVnegFosneg.expBdonut); % average total exp blue around mCherry pos, PV neg and Fos NEG cells
    meanexpmCherryposPVneg_Fosneg = mean(DatamCherryposPVnegFosneg.meanBdonut); %average mean exp blue around mCherry pos, PV neg and Fos NEG cells
   
    binexpmCherryposPVneg_Fosneg = mean (DatamCherryposPVnegFosneg.expBINdonut) ; %binarized 
    meanbinexpmCherryposPVneg_Fosneg = mean (DatamCherryposPVnegFosneg.meanexpBINdonut) ;% average per area binarized

    countmCherryposPVneg_Fosneg=height(DatamCherryposPVnegFosneg.cellnr);

    % all mCherry cells that are PV negative, but Fos positive
    idxmCherryposPVnegFospos = DatamCherryposPVneg.percentileG>=thresholdGhigh ;
    DatamCherryposPVnegFospos = DatamCherryposPVneg(idxmCherryposPVnegFospos,:);

    totexpmCherryposPVneg_Fospos=mean(DatamCherryposPVnegFospos.expBdonut); %average exp blue around mCherry pos, PV neg and Fos POS cells
    meanexpmCherryposPVneg_Fospos=mean(DatamCherryposPVnegFospos.meanBdonut); % average mean exp blue around mCherry pos, PV neg and Fos POS cells
    
    binexpmCherryposPVneg_Fospos = mean (DatamCherryposPVnegFospos.expBINdonut) ; %binarized 
    meanbinexpmCherryposPVneg_Fospos = mean (DatamCherryposPVnegFospos.meanexpBINdonut) ;% average per area binarized

    countmCherryposPVneg_Fospos=height(DatamCherryposPVnegFospos.cellnr);

    % all mCherry neg vicinity cell that are PV negative, and Fos negative
    idxmCherrynegVinPVnegFosneg = DatamCherrynegVinPVneg.percentileG<=thresholdGlow ;
    DatamCherrynegVinPVnegFosneg = DatamCherrynegVinPVneg(idxmCherrynegVinPVnegFosneg,:);
   
    totexpmCherrynegVinPVneg_Fosneg=mean(DatamCherrynegVinPVnegFosneg.expBdonut); %average total exp blue around mCherry neg, PV neg and Fos NEG cells
    meanexpmCherrynegVinPVneg_Fosneg=mean(DatamCherrynegVinPVnegFosneg.meanBdonut); %average mean exp blue around mCherry neg, PV neg and Fos NEG cells
    
    binexpmCherrynegVinPVneg_Fosneg = mean (DatamCherrynegVinPVnegFosneg.expBINdonut) ; %binarized 
    meanbinexpmCherrynegVinPVneg_Fosneg = mean (DatamCherrynegVinPVnegFosneg.meanexpBINdonut) ;% average per area binarized
    
    countmCherrynegVinPVneg_Fosneg=height(DatamCherrynegVinPVnegFosneg.cellnr);

    % all mCherry neg cells vicinity that are PV negative, but Fos positive
    idxmCherrynegVinPVnegFospos = DatamCherrynegVinPVneg.percentileG>=thresholdGhigh ;
    DatamCherrynegVinPVnegFospos = DatamCherrynegVinPVneg(idxmCherrynegVinPVnegFospos,:);
    
    totexpmCherrynegVinPVneg_Fospos=mean(DatamCherrynegVinPVnegFospos.expBdonut); %average total exp blue around mCherry neg, PV neg and Fos POS cells
    meanexpmCherrynegVinPVneg_Fospos=mean(DatamCherrynegVinPVnegFospos.meanBdonut); %average mean exp blue around mCherry neg, PV neg and Fos POS cells
    
    binexpmCherrynegVinPVneg_Fospos = mean (DatamCherrynegVinPVnegFospos.expBINdonut) ; %binarized 
    meanbinexpmCherrynegVinPVneg_Fospos = mean (DatamCherrynegVinPVnegFospos.meanexpBINdonut) ;% average per area binarized
    
    countmCherrynegVinPVneg_Fospos=height(DatamCherrynegVinPVnegFospos.cellnr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fos positive and negative cells, excluding PV pos cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Fos cells that are PV negative 
    idxFosposPVneg = DataFospos.percentileB<thresholdB;
    DataFosposPVneg = DataFospos(idxFosposPVneg,:);
    
    totexpFosposPVneg= mean(DataFosposPVneg.expBdonut); % average total exp blue around Fos pos cells (but PV neg)
    meanexpFosposPVneg = mean(DataFosposPVneg.meanBdonut); % average mean exp blue around Fos pos cells (but PV neg)
    
    binexpFosposPVneg = mean (DataFosposPVneg.expBINdonut) ; %binarized 
    meanbinexpFosposPVneg = mean (DataFosposPVneg.meanexpBINdonut) ;% average per area binarized
    
    countFosposPVneg= height(DataFosposPVneg.cellnr);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fos negative cells - only in the vicinity of the Fos positive cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Fosneg_vectors = [DataFosneg.xcoordcell, DataFosneg.ycoordcell] ; %x and y coordinates of Fos neg cells
    FosposPVneg_vectors = [DataFosposPVneg.xcoordcell, DataFosposPVneg.ycoordcell] ;%x and y coordinates of Fos pos cells

    vicinityN = 3; %the number of Fos cells you want returned lying in the vicinity of the Fos pos cells
    
    % Return an array with Distances from 'mCherryneg_vectors' to points in the 'mCherryposPVneg_vectors' sorted 
    % small to large, and keep only the 3 smallest/closest

    [DisttoPos,IndextoPos] = pdist2(Fosneg_vectors,FosposPVneg_vectors,'euclidean','Smallest',vicinityN);
    
    dataFos_neg_vicinity = DataFosneg(IndextoPos, :); %make a new array with only the Fos neg cells
    % that are in the vicinity of the Fos pos cells (this should be 3x the amount of FosposPVneg cells)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fos neg vicinity cells, excl PV pos cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    idxFosnegPVneg = dataFos_neg_vicinity.percentileB<thresholdB;
    DataFosnegPVneg = dataFos_neg_vicinity(idxFosnegPVneg,:);
    
    totexpFosnegPVneg= mean(DataFosnegPVneg.expBdonut); % average total exp blue around Fos pos cells (but PV neg)
    meanexpFosnegPVneg = mean(DataFosnegPVneg.meanBdonut); %average mean exp blue around Fos pos cells (but PV neg)
    
    binexpFosnegPVneg = mean (DataFosnegPVneg.expBINdonut) ; %binarized 
    meanbinexpFosnegPVneg = mean (DataFosnegPVneg.meanexpBINdonut) ;% average per area binarized
    
    countFosnegPVneg= height(DataFosnegPVneg.cellnr);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % General info of the image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    countNisslcells= height(cell_data.cellnr); %amount of total Nissl cells
    countPvPoscells= height(DataPVpos.cellnr); %amount of total PV cells
    countmCherryPoscells= height(DatamCherrypos.cellnr); %amount of total mCherry cells
    countmCherryNegcells= height(DatamCherryneg.cellnr); %amount of total mCherry neg cells
    countmCherryNegviccells = height(datamCherry_neg_vicinity); %amount of total mCherry neg cells in vicinity of mCherry cells
    countFosCells = height(DataFospos.cellnr); %amount of total Fos cells 

    % add data to a final data frame 
   
    Finaldata(rownr, 1)  = currName ; %file name
    Finaldata(rownr, 2)  = num2cell(countNisslcells) ; %amount of nissl cells
    Finaldata(rownr, 3)  = num2cell(countPvPoscells) ; %amount of PV cells
    Finaldata(rownr, 4)  = num2cell(countmCherryPoscells) ; %amount of mCherry cells
    Finaldata(rownr, 5)  = num2cell(countmCherryNegcells) ; %amount of mCherry neg cells
    Finaldata(rownr, 6)  = num2cell(countmCherryNegviccells) ;%amount of mCherry neg cells in vicinity to mCherry pos cells
    Finaldata(rownr, 7)  = num2cell(countFosCells) ; %amount of Fos cells 

    %mCherry pos
    Finaldata(rownr, 8)  = num2cell(countmCherryposPVneg) ; %count of mCherry Pos cells, no PV cell
    Finaldata(rownr, 9)  = num2cell(totexpmCherryposPVneg) ; %tot exp 
    Finaldata(rownr, 10) = num2cell(meanexpmCherryposPVneg) ; %mean exp
    Finaldata(rownr, 11) = num2cell(binexpmCherryposPVneg) ; % tot bin exp
    Finaldata(rownr, 12) = num2cell(meanbinexpmCherryposPVneg) ; %av bin exp

    Finaldata(rownr, 13) = num2cell(countmCherryposPVneg_Fospos) ; %count mc+ PV - Fos pos
    Finaldata(rownr, 14) = num2cell(totexpmCherryposPVneg_Fospos); %tot exp
    Finaldata(rownr, 15) = num2cell(meanexpmCherryposPVneg_Fospos) ; %mean exp
    Finaldata(rownr, 16) = num2cell(binexpmCherryposPVneg_Fospos) ; %tot bin exp
    Finaldata(rownr, 17) = num2cell(meanbinexpmCherryposPVneg_Fospos) ; %av bin exp
    
    Finaldata(rownr, 18) = num2cell(countmCherryposPVneg_Fosneg) ; %count mC+ PV- Fos neg
    Finaldata(rownr, 19) = num2cell(totexpmCherryposPVneg_Fosneg) ; %tot exp
    Finaldata(rownr, 20) = num2cell(meanexpmCherryposPVneg_Fosneg) ; %mean exp
    Finaldata(rownr, 21) = num2cell(binexpmCherryposPVneg_Fosneg) ;%tot bin exp
    Finaldata(rownr, 22) = num2cell(meanbinexpmCherryposPVneg_Fosneg) ;%av bin exp

    %mCherry neg, vicinity cells 
    Finaldata(rownr, 23) = num2cell(countmCherrynegVinPVneg) ;%count of mCherry neg cells, no PV cell, vicinity
    Finaldata(rownr, 24) = num2cell(totexpmCherrynegVinPVneg) ; %tot exp 
    Finaldata(rownr, 25) = num2cell(meanexpmCherrynegVinPVneg) ; %mean exp
    Finaldata(rownr, 26) = num2cell(binexpmCherrynegVinPVneg) ;%tot bin exp
    Finaldata(rownr, 27) = num2cell(meanbinexpmCherrynegVinPVneg) ;%av bin exp

    Finaldata(rownr, 28) = num2cell(countmCherrynegVinPVneg_Fospos) ; % count mC- PV- Fos pos
    Finaldata(rownr, 29) = num2cell(totexpmCherrynegVinPVneg_Fospos) ; %tot exp
    Finaldata(rownr, 30) = num2cell(meanexpmCherrynegVinPVneg_Fospos) ; %mean exp
    Finaldata(rownr, 31) = num2cell(binexpmCherrynegVinPVneg_Fospos) ; %tot bin exp
    Finaldata(rownr, 32) = num2cell(meanbinexpmCherrynegVinPVneg_Fospos) ;%av bin exp

    Finaldata(rownr, 33) = num2cell(countmCherrynegVinPVneg_Fosneg) ; %count of mC- PV- Fos neg
    Finaldata(rownr, 34) = num2cell(totexpmCherrynegVinPVneg_Fosneg) ; %tot exp
    Finaldata(rownr, 35) = num2cell(meanexpmCherrynegVinPVneg_Fosneg); %mean exp
    Finaldata(rownr, 36) = num2cell(binexpmCherrynegVinPVneg_Fosneg) ;%tot bin exp
    Finaldata(rownr, 37) = num2cell(meanbinexpmCherrynegVinPVneg_Fosneg) ;%av bin exp

    %Fos pos
    Finaldata(rownr, 38) = num2cell(countFosposPVneg) ; % count Fos pos, no PV cell
    Finaldata(rownr, 39) = num2cell(totexpFosposPVneg) ; %tot exp
    Finaldata(rownr, 40) = num2cell(meanexpFosposPVneg) ; %mean exp
    Finaldata(rownr, 41) = num2cell(binexpFosposPVneg) ;%tot bin exp
    Finaldata(rownr, 42) = num2cell(meanbinexpFosposPVneg) ;%av bin exp

    %Fos neg, vicinity
    Finaldata(rownr, 43) = num2cell(countFosnegPVneg) ; % count Fos neg, no PV cell
    Finaldata(rownr, 44) = num2cell(totexpFosnegPVneg) ; %tot exp
    Finaldata(rownr, 45) = num2cell(meanexpFosnegPVneg); %mean exp
    Finaldata(rownr, 46) = num2cell(binexpFosnegPVneg) ;%tot bin exp
    Finaldata(rownr, 47) = num2cell(meanbinexpFosnegPVneg) ;%av bin exp

end  

 %removes NaN and replaces with 0's
 for k = 1:size(Finaldata,1)
  for j = 1:size(Finaldata,2)
    Finaldata{k,j}(isnan(Finaldata{k,j}))=0;
  end
end   



