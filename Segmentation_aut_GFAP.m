function [masks_path]= Segmentation_aut_GFAP (ind_tissue_path, project_path) %(path, tissue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GFAP segmentation %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Reading the images %%
    imagefiles = dir([ind_tissue_path, '/figure*.jpg']);
    nfiles = length(imagefiles);
    
    
  
    %% Segmenting 5 tissues %%
    
    % NUCLEUS: nuc_control.csv contains the pixels for healthy nucleus. The
    % file produced by nuc_control is called nucleus*.jpg. The image is
    % segmented considering 1.5std for RED, 1.3std for GREEN and 1std for
    % BLUE. NOTE: Not worth, deleted.
    
    % NUCLEUS GLIOMA: nuc_high.csv contains the pixels for nucleus in
    % glioma. The file produced is called nuc_glio*.jpg. The image is
    % segmented considering 2.5std
    
    % PROTEIN GLIOMA: protein_high.csv contains the pixels for protein in high
    % degree of the glioma with particular aglomeration. The file produced is 
    % called glio_protein*.jpg. The image is segmented considering 1.5std
    % for RED, 1.4 for GREEN and 1.3 for BLUE. The consensus is 1.5std.
    
    % PROTEIN ASTROCYTE: astrocyte.csv contains the pixes for protein
    % expressed in normal tissue, based on the control images. The file 
    % generated is called strocyte*.jpg. The image is segmented 
    % considering 1std for RED, from 0 to 28 pixels value for GREEN and from  
    % 0 to 15 pixels value for BLUE. The consensus is 1std.
    
    % PROTEIN NEUTROPIL: protein_control.csv contains the pixes for protein
    % expressed in normal tissue, based on the control images. The file 
    % generated is called control_protein*.jpg. The image is segmented 
    % considering 1.5std for RED, 1.3std for GREEN and 1.5std for BLUE. The
    % consensus is 1.5std.

    mkdir(ind_tissue_path, 'Masks');
    masks_path = [ind_tissue_path, 'Masks/'];
    for j=1:nfiles
        image=imagefiles(j).name;
        I = imread([ind_tissue_path,image]);
        %This makes the order of the files to change maintaining the number
        %they got in the meander sampling.  Therefore it is necessary to
        %write them all again and to collect same number for figure,
        %epithelium and mucous. That won't happen whenever we use directlu
        %the imagefile imb_ab from the meander sampling without writing the
        %figures*.jpg.
        
%         imagefile = img_ab(:,:,:,j);
%         I = imread (imagefile);
        
        
           %% Segmenting nucleus healthy
        filename = 'nuc_control_prot.csv';  % Pixels library for Nucleus with normal response
        pixels = csvread([project_path, filename]);

        % Look for the standard deviation to choose the values of the threshold.
        color = mean(pixels); %It deletes the distribution of the selected pixels.
        threshold = 1.5*std(pixels);

        redBand = I(:,:, 1);
        greenBand = I(:,:, 2);
        blueBand = I(:,:, 3);
 
        % Threshold each color band

        thresholdmax = color + threshold;
        redthresholdmax = thresholdmax (1);
        greenThresholdmax = thresholdmax (2);
        blueThresholdmax = thresholdmax (3);

        thresholdmin = color - threshold;
        redthresholdmin = thresholdmin (1);
        greenThresholdmin = thresholdmin (2);
        blueThresholdmin = thresholdmin (3);

        redMaskmax = (redBand < redthresholdmax);
        greenMaskmax = (greenBand < greenThresholdmax);
        blueMaskmax = (blueBand < blueThresholdmax);

        redMaskmin = (redBand > redthresholdmin);
        greenMaskmin = (greenBand > greenThresholdmin);
        blueMaskmin = (blueBand > blueThresholdmin);

        redMask = uint8( redMaskmax & redMaskmin);
        greenMask = uint8(greenMaskmax & greenMaskmin);
        blueMask = uint8(blueMaskmax & blueMaskmin);

        % Combine the masks to find where all 3 are "true."
        resultObjectsMask = uint8(redMask & greenMask & blueMask);
        [a, b] = size(resultObjectsMask);
        white = double(ones(a , b));
        finalMask = double (resultObjectsMask) .* double(white);
        
        imwrite(finalMask,[masks_path,'nucleus',(num2str(j)),'.jpg']);
    
        
%         %% Evaluating if the image is empty area
%         [evaluation] = image_eval (finalMask);
%         if strcmp(evaluation, 'negative')
%             imwrite(finalMask,[masks_path,'epithelium',(num2str(j)),'.jpg']);
%         else
%         end
       

          %% Segmenting nucleus with sick expression
          
        filename = 'nuc_high.csv'; % Pixels library for glioma nucleus
        pixels = csvread([project_path, filename]);

        % Look for the standard deviation to choose the values of the threshold.
        color = mean(pixels); %It deletes the distribution of the selected pixels.
        threshold = 2.5*std(pixels);
        
        thresholdmax = color + threshold;
        redthresholdmax = thresholdmax (1);
        greenThresholdmax = thresholdmax (2);
        blueThresholdmax = thresholdmax (3);

        thresholdmin = color - threshold;
        redthresholdmin = thresholdmin (1);
        greenThresholdmin = thresholdmin (2);
        blueThresholdmin = thresholdmin (3);

        redMaskmax = (redBand < redthresholdmax);
        greenMaskmax = (greenBand < greenThresholdmax);
        blueMaskmax = (blueBand < blueThresholdmax);

        redMaskmin = (redBand > redthresholdmin);
        greenMaskmin = (greenBand > greenThresholdmin);
        blueMaskmin = (blueBand > blueThresholdmin);

        redMask = uint8( redMaskmax & redMaskmin);
        greenMask = uint8(greenMaskmax & greenMaskmin);
        blueMask = uint8(blueMaskmax & blueMaskmin);

        % Combine the masks to find where all 3 are "true."
        resultObjectsMask = uint8(redMask & greenMask & blueMask);
        [a, b] = size(resultObjectsMask);
        white = double(ones(a , b));
        finalMask = double (resultObjectsMask) .* double(white);
        imwrite(finalMask,[masks_path,'nucglio',(num2str(j)),'.jpg']);
        
        
        %% Segmenting extracellular protein from normal response (neutropil)
        filename = 'protein_control.csv'; % Pixels library for protein from neutropil
        pixels = csvread([project_path, filename]);

        % Look for the standard deviation to choose the values of the threshold.
        color = mean(pixels); %It deletes the distribution of the selected pixels.
        threshold = 1.5*std(pixels);
        
        thresholdmax = color + threshold;
        redthresholdmax = thresholdmax (1);
        greenThresholdmax = thresholdmax (2);
        blueThresholdmax = thresholdmax (3);

        thresholdmin = color - threshold;
        redthresholdmin = thresholdmin (1);
        greenThresholdmin = thresholdmin (2);
        blueThresholdmin = thresholdmin (3);

        redMaskmax = (redBand < redthresholdmax);
        greenMaskmax = (greenBand < greenThresholdmax);
        blueMaskmax = (blueBand < blueThresholdmax);

        redMaskmin = (redBand > redthresholdmin);
        greenMaskmin = (greenBand > greenThresholdmin);
        blueMaskmin = (blueBand > blueThresholdmin);

        redMask = uint8( redMaskmax & redMaskmin);
        greenMask = uint8(greenMaskmax & greenMaskmin);
        blueMask = uint8(blueMaskmax & blueMaskmin);

        % Combine the masks to find where all 3 are "true."
        resultObjectsMask = uint8(redMask & greenMask & blueMask);
        [a, b] = size(resultObjectsMask);
        white = double(ones(a , b));
        finalMask = double (resultObjectsMask) .* double(white);
        imwrite(finalMask,[masks_path,'health_protein',(num2str(j)),'.jpg']);

        %% Segmenting protein expressed in glioma high degree
        filename = 'protein_high.csv'; % Pixels library for protein
        pixels = csvread([project_path, filename]);

        % Look for the standard deviation to choose the values of the threshold.
        color = mean(pixels); %It deletes the distribution of the selected pixels.
        threshold = 1.5*std(pixels);
        
        thresholdmax = color + threshold;
        redthresholdmax = thresholdmax (1);
        greenThresholdmax = thresholdmax (2);
        blueThresholdmax = thresholdmax (3);

        thresholdmin = color - threshold;
        redthresholdmin = thresholdmin (1);
        greenThresholdmin = thresholdmin (2);
        blueThresholdmin = thresholdmin (3);

        redMaskmax = (redBand < redthresholdmax);
        greenMaskmax = (greenBand < greenThresholdmax);
        blueMaskmax = (blueBand < blueThresholdmax);

        redMaskmin = (redBand > redthresholdmin);
        greenMaskmin = (greenBand > greenThresholdmin);
        blueMaskmin = (blueBand > blueThresholdmin);

        redMask = uint8( redMaskmax & redMaskmin);
        greenMask = uint8(greenMaskmax & greenMaskmin);
        blueMask = uint8(blueMaskmax & blueMaskmin);

        % Combine the masks to find where all 3 are "true."
        resultObjectsMask = uint8(redMask & greenMask & blueMask);
        [a, b] = size(resultObjectsMask);
        white = double(ones(a , b));
        finalMask = double (resultObjectsMask) .* double(white);
        imwrite(finalMask,[masks_path,'glio_protein',(num2str(j)),'.jpg']);

%% Segmenting GFAP protein expressed in the astrocytes
        filename = 'astrocyte.csv'; % Pixels library for protein
        pixels = csvread([project_path, filename]);

        % Look for the standard deviation to choose the values of the threshold.
        color = mean(pixels); %It deletes the distribution of the selected pixels.
        threshold = 1*std(pixels);
        
        thresholdmax = color + threshold;
        redthresholdmax = thresholdmax (1);
        greenThresholdmax = thresholdmax (2);
        blueThresholdmax = thresholdmax (3);

        thresholdmin = color - threshold;
        redthresholdmin = thresholdmin (1);
        greenThresholdmin = thresholdmin (2);
        blueThresholdmin = thresholdmin (3);

        redMaskmax = (redBand < redthresholdmax);
        greenMaskmax = (greenBand < greenThresholdmax);
        blueMaskmax = (blueBand < blueThresholdmax);

        redMaskmin = (redBand > redthresholdmin);
        greenMaskmin = (greenBand > greenThresholdmin);
        blueMaskmin = (blueBand > blueThresholdmin);

        redMask = uint8( redMaskmax & redMaskmin);
        greenMask = uint8(greenMaskmax & greenMaskmin);
        blueMask = uint8(blueMaskmax & blueMaskmin);

        % Combine the masks to find where all 3 are "true."
        resultObjectsMask = uint8(redMask & greenMask & blueMask);
        [a, b] = size(resultObjectsMask);
        white = double(ones(a , b));
        finalMask = double (resultObjectsMask) .* double(white);
        imwrite(finalMask,[masks_path,'astro_protein',(num2str(j)),'.jpg']);
        
    end

end

