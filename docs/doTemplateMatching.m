% Script is used for template matching. 
% KS 0708 -- Cleaned Nirav's script.

%Define your participants folder
supercomputer_ic_folder = '/fs/plum3_share2/Users/elena/Projects/ICA_img_clustering/AD/ICs/';
num_of_subjects = 38;

%Define the output folder where the results of this script will be
%stored.
out_folder = '/fs/plum3_share2/Users/elena/Projects/ICA_img_clustering/AD/template_matching_results';
final_table=[];

%Define template that you wish to use
template_file = '/fs/plum3_share3/templates/image_templates_2004PNAS/3T_ICAcontrol_radio_001_001.img'; 
template_name = 'DMN';

% Load the template ROI
V = spm_vol( template_file );
[Y,XYZ] = spm_read_vols(V,0);
roi_intensity = 200; % The template image is a binary file -- activation indicated by 200.

% For each subject compute the similarity score of all the components
% belonging to that subject with the template.
for ii = 1:num_of_subjects

disp(['Processing Subject#' num2str(ii)]);

subj_folder = sprintf('%ss%02d',supercomputer_ic_folder, ii);
cd(subj_folder);
comp_count = dir('vol*.nii');

    for jj = 1:length(comp_count)
        num_component = jj - 1;
        disp(['Processing Component#' num2str(num_component+1)]);

        activation_in_roi = 0;voxel_in_roi = 0;
        activation_out_roi = 0;voxel_out_roi = 0;

        if num_component <= 9
          component_file_name= sprintf('%s/vol000%d.nii', subj_folder, num_component);
        elseif num_component <= 99 && num_component > 9
          component_file_name = sprintf('%s/vol00%d.nii', subj_folder, num_component);
        else 
          component_file_name = sprintf('%s/vol0%d.nii', subj_folder, num_component);
        end

        V1 = spm_vol(component_file_name);
        [Y1, XYZ1] = spm_read_vols(V1,0);

        indexes = find( Y == roi_intensity );
        for i=1:numel(indexes)
            if Y1(indexes(i)) ~= 0
                activation_in_roi = activation_in_roi + Y1(indexes(i));
                voxel_in_roi = voxel_in_roi+1;
            end
        end
        Y1(indexes) = 0;
        for i=1:numel(Y)
            if Y1(i) ~= 0
                activation_out_roi = activation_out_roi+Y1(i);
                voxel_out_roi = voxel_out_roi+1;
            end
        end

       activation_difference(jj) = (activation_in_roi/voxel_in_roi) - (activation_out_roi/voxel_out_roi);

    end 
    
    % Compute the Top 3 best similarity score components
    [comp_sort_descend idx] = sort(abs(activation_difference),'descend');

    high_freq = 0; % assumption: all the components are good/non-high frequency components. Run the ICA_high frequency script if you are skeptical about this assumption
    table_components_final = [ii, idx(1),comp_sort_descend(1), high_freq, idx(2),comp_sort_descend(2),high_freq, idx(3),comp_sort_descend(3),high_freq]; 
    
    final_table = [final_table; table_components_final];
    clear activation_difference
    
end

save([out_folder filesep 'best_' template_name '_components.txt'], 'final_table','-ascii','-tabs');




