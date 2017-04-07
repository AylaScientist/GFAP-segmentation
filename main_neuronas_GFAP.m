%%%
% This program segments the GFAP images from a particular experimental
% project and extracts:
%   @number of cancer nuclei
%   @GFAP Neuropile (pixels)
%   @GFAP Glioma cells (pixels)
%   @GFAP Astrocytes (pixels)
% 
% The reference volume is the size of the subset previously set @ 200x250 pixels


function main_neuronas_GFAP ()
% Initializing the project:
    project_name = input ('Please, insert the project name. ', 's');
    %path = 'C:/Users/Ayla/Dropbox/Matlab/';
    %path = what('Matlab');
    %path='/media/Data Storage/Dropbox/Dropbox/Matlab/'
    %path = '/Users/user/Dropbox/Matlab/';
    path = 'D:\Dropbox\Dropbox\Matlab\';
    project_path = [path,project_name,'/'];
    project_description = readtable ([project_path, project_name,'.xls']);
    project_variables = table2cell(project_description);
    [individuals, columns] = size (project_description);
    tissues = columns - 1;
    mkdir('temp');
    images_path =what('temp');%Temp folder must be created already
    C = [];
    %%
% Processing the data    
row_counter = 1;
    for i = 1:tissues
        for j = 1:individuals
            %Reading the image files for this individual in this tissue
            sample_id = project_variables(j,1);
            sample_tissue = project_variables(j,(i+1));
            sample = strcat(num2str(j), sample_tissue);
            sample_name = sample{1};
            %Creating temporary dir
            mkdir([images_path.path, '/'],sample_name);
            ind_tissue_path = [images_path.path,'/',sample_name,'/'];
            %Finding the original image to segment
            filename = [sample_name,'*.jpg'];
            imagefiles = dir([project_path, filename]);  
            [nfiles, ncols] = size(imagefiles); 
            [ind_measures] = image_processor (imagefiles, nfiles, ind_tissue_path, project_path, path, i, j, C);
            [final_nfiles, length] = size(ind_measures);
            for ii = 1:final_nfiles
                sample_id_str = sample_id{1};
                sample_tissue_str = sample_tissue{1};
                loc_number = num2str(ii);
                protein_area = num2str(ind_measures(ii,1));
                nucleus_number = num2str(ind_measures(ii,2));
                no_exp_area = num2str(ind_measures(ii,3));
                healthy_nucleus = num2str(ind_measures(ii,4));
                astro_area = num2str(ind_measures(ii,5));
                reference = num2str(ind_measures(ii,6));
                C{row_counter, 1} = sample_id_str;
                C{row_counter, 2} = sample_tissue_str;
                C{row_counter, 3} = loc_number;
                C{row_counter, 4} = protein_area;
                C{row_counter, 5} = nucleus_number;
                C{row_counter, 6} = no_exp_area;
                C{row_counter, 7} = healthy_nucleus;
                C{row_counter, 8} = astro_area;
                C{row_counter, 9} = reference;
                row_counter = row_counter + 1;
            end
        end
    end
    %% Writting data into a table
    T = cell2table(C);
    T.Properties.VariableNames = {'individual', 'tissue', 'location', 'protein', 'nucleus', 'glioma_protein','healthy_nucleus','astro_protein','area_ref'};
    %writetable (T, [project_path,project_name,'results.xls']); %Only if
    %you have excel installed.
    writetable (T, [project_path,project_name,'results.csv']); 
    %,'Delimiter',',','QuoteStrings',true
end