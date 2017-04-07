%% Count the masked objects as tags 
function [ind_measures] = counting_areas_GFAP (masks_path)
   
    %% Number of subsets %%
    
    imagefiles = dir([masks_path, 'nucglio*.jpg']);
    nfiles = length(imagefiles);
    proteinfiles = dir([masks_path, 'health_protein*.jpg']);
    glioprotfiles = dir([masks_path, 'glio_protein*.jpg']);
    control_nuc = dir([masks_path, 'nucleus*.jpg']);
    astro = dir([masks_path, 'astro_protein*.jpg']);
    
    
    nucleus_number = zeros(nfiles,1);
    protein_area = zeros(nfiles, 1);
    sick_exp_area = zeros(nfiles, 1);
    control_nucleus_number = zeros(nfiles, 1);
    astro_area = zeros(nfiles, 1);
    
    %% Image treatement %%
 
    % treating each binary image
    for j=1:nfiles
        %Chargin nucleus glioma image
        currentnucleusfile = [masks_path imagefiles(j).name];
        I = imread(currentnucleusfile);
        
        % Noise removal 
        seD = strel('disk',1);
        denoise = imerode(I,seD);
        
        % Filling holes
        filled = imfill(denoise, 'holes');
        
        % Steve method from here %

        % Open area
        bw2 = ~bwareaopen(~filled, 10);
        
        % Distance transform
        D = -bwdist(~filled);

        % Compute the watersheed transform of D
        Ld = watershed(D);
        
% The watershed ridge lines, in white, correspond to Ld == 0. Let's use these ridge 
%lines to segment the binary image by changing the corresponding pixels into background.

        bw2 = filled;
        bw2(Ld == 0) = 0;

% The "raw" watershed transform is known for its tendency to "oversegment" 
% an image. The reason is something I mentioned above: each local minimum, no 
% matter how small, becomes a catchment basin. A common trick, then, in 
% watershed-based segmentation methods is to filter out tiny local minima 
% using imextendedmin and then modify the distance transform so that no 
% minima occur at the filtered-out locations. This is called "minima imposition" 
% and is implemented via the function imimposemin.

% The following call to imextendedmin should ideally just produce small 
% spots that are roughly in the middle of the cells to be segmented. I'll use 
% imshowpair to superimpose the mask on the original image.

        mask = imextendedmin(D,2);

% Modify the distance transform so it only has minima at the desired locations, 
% and then repeat the watershed steps above.

        D2 = imimposemin(D,mask);
        Ld2 = watershed(D2);
        bw3 = filled;
        bw3(Ld2 == 0) = 0;

        % Continuing with my work %
        
        %Noise removal to separate connected objects
        seD = strel('disk',6);
        denoise2 = imerode(bw3,seD);
        
        % Compute the watershed transform and display the resulting label matrix as an RGB image.
        [L2, num] = bwlabel(denoise2);
      

        % Writing in a vector
        nucleus_number(j) = num*10000;% * 10000 for significative number when comparing the density
        
    end
%%    
    for j=1:nfiles
        %Chargin control nucleus image (healthy nucleus)
        currentcontrolfile = [masks_path control_nuc(j).name];
        I = imread(currentcontrolfile);
        

        % Noise removal 
        seD = strel('disk',1);
        denoise = imerode(I,seD);
        
        % Filling holes
        filled = imfill(denoise, 'holes');
        
        % Steve method from here %

        % Open area
        bw2 = ~bwareaopen(~filled, 10);
        
        % Distance transform
        D = -bwdist(~filled);

        % Compute the watersheed transform of D
        Ld = watershed(D);
        
% The watershed ridge lines, in white, correspond to Ld == 0. Let's use these ridge 
%lines to segment the binary image by changing the corresponding pixels into background.

        bw2 = filled;
        bw2(Ld == 0) = 0;

% The "raw" watershed transform is known for its tendency to "oversegment" 
% an image. The reason is something I mentioned above: each local minimum, no 
% matter how small, becomes a catchment basin. A common trick, then, in 
% watershed-based segmentation methods is to filter out tiny local minima 
% using imextendedmin and then modify the distance transform so that no 
% minima occur at the filtered-out locations. This is called "minima imposition" 
% and is implemented via the function imimposemin.

% The following call to imextendedmin should ideally just produce small 
% spots that are roughly in the middle of the cells to be segmented. I'll use 
% imshowpair to superimpose the mask on the original image.

        mask = imextendedmin(D,2);

% Modify the distance transform so it only has minima at the desired locations, 
% and then repeat the watershed steps above.

        D2 = imimposemin(D,mask);
        Ld2 = watershed(D2);
        bw3 = filled;
        bw3(Ld2 == 0) = 0;

        % Continuing with my work %
        
        %Noise removal to separate connected objects
        seD = strel('disk',6);
        denoise2 = imerode(bw3,seD);
        
        % Compute the watershed transform and display the resulting label matrix as an RGB image.
        [L2, num] = bwlabel(denoise2);
      
         
         control_nucleus_number(j) = num*10000;% *10000 for significative number when calculating the density
        
    end

%%    
    for j=1:nfiles
        %Protein
        currentproteinfile = [masks_path proteinfiles(j).name]; %Health protein
        currentsick_protfile = [masks_path glioprotfiles(j).name]; %Glioma protein
        currentastro = [masks_path astro(j).name]; %Astrocyte protein
        I = imread(currentproteinfile);
        S = imread(currentsick_protfile);
        A = imread(currentastro);

        
        % Smooth 
        seL = strel('line',3,3);
        %final = imerode(I, seL);
        finalS = imerode(S, seL);
        
        %Quantification of the area in pixels
        % white = bwarea(final);
        white = sum(I(:,1));
        whiteS = sum(finalS(:,1));
        whiteA = sum(A(:,1));
        % Area of reference for this subset
        [rows, columns] = size (I);
        ref = rows * columns; 
        
        %Counting the number of pixels (in comments hoe to count the area)
        %area = white * 0.01768 * 0.01768; %1 pixel is 0.01768 um
        area = white;
        %areaS = whiteS * 0.01768 * 0.01768;
        areaS = whiteS;
        areaA = whiteA;
        %areaREF = ref * 0.01768 * 0.01768;% Area of reference for this subset
        areaREF = ref;
        
        %Writing it in a vector
        protein_area(j) = area;% Area of the protein in the normal tissue
        sick_exp_area(j) = areaS;% Area of the protein in the glioma
        astro_area(j) = areaA;
        reference(j) = areaREF;% Area of reference for this subset
        
    end

%%     
    ind_measures = zeros(nfiles, 5);
    ind_measures (:,1) = protein_area;
    ind_measures (:,2) = nucleus_number;
    ind_measures (:,3) = sick_exp_area;
    ind_measures (:,4) = control_nucleus_number;
    ind_measures (:,5) = astro_area;
    ind_measures (:,6) = reference;
   
end
