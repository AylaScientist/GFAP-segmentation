function [ind_measures] = image_processor (imagefiles, nfiles, ind_tissue_path, project_path, path, i, j, sample_id, sample_tissue, C)

    %Subsetting the image for this id in this tissue
    Periodic_sampling_aut(imagefiles, nfiles, ind_tissue_path, project_path);
    %[writing] = Meander_sampling_aut(imagefiles, nfiles, ind_tissue_path, project_path);
    %writing
    %Segmenting the images
    masks_path = Segmentation_aut_GFAP (ind_tissue_path, path);
    %Contabilizing histological areas for this id in this tissue
    [ind_measures] = counting_areas_GFAP (masks_path);
    rmdir([ind_tissue_path],'s');
end

