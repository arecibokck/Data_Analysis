function Analyze_Insect(data_folder)

    treatment_subdirs = subDirs(data_folder);
    for i = 1:length(treatment_subdirs)
        treatment_dir = treatment_subdirs(i).name;
        motion_type_data_files = dir(strcat(data_folder,'\',...
                                            treatment_dir,'\*.csv'));
        for j = 1:length(motion_type_data_files)
            
            motion_type_data_file = motion_type_data_files (j).name;
            
            if isempty(strfind(motion_type_data_file,'angle'))
                
                current_data_file = strcat(data_folder,'\',...
                                           treatment_dir,'\',...
                                           motion_type_data_file);
                current_data = Arrange_Points(current_data_file);
                angle_data = calculateAngle(current_data);
                angle_data_file = strrep(current_data_file,...
                                         '.csv','-angle.csv');
                csvwrite(angle_data_file, angle_data);
                
            end
            
        end
        
        motion_type_data_files = dir(strcat(data_folder,'\',...
                                    treatment_dir,'\*.csv'));
                                
        for j = 1:length(motion_type_data_files)
            
            motion_type_data_file = motion_type_data_files (j).name;
            
            if ~isempty(strfind(motion_type_data_file,'angle'))                    
                if ~isempty(strfind(motion_type_data_file,'Hz'))
                    [stimulus, response] = ReadCSV(data_folder, treatment_dir, motion_type_data_file)
                    AnalyzeSineMotion(stimulus, response, strcat(data_folder,'\', treatment_dir, '\'));
                    delete(strcat(data_folder,'\', treatment_dir, '\',motion_type_data_file));
                else
                    [stimulus, response] = ReadCSV(data_folder, treatment_dir, motion_type_data_file)
                    AnalyzeChirpMotion(stimulus, response, strcat(data_folder,'\', treatment_dir, '\'));
                    delete(strcat(data_folder,'\', treatment_dir, '\', motion_type_data_file));
                end
            end
        end
    end
    
end 

function dirListing = subDirs(directory)

    dirListing = dir(directory);
    dirListing(~[dirListing.isdir]) = []; %remove items which are not directories
    tf = ismember( {dirListing.name}, {'.', '..'});
    dirListing(tf) = [];  %remove current and parent directory
    
end

function angle_data = calculateAngle(data)

    Head_1 = [data(:, 1) data(:, 2)]; % specifying coordinates
    Head_2 = [data(:, 3) data(:, 4)];
    Thorax_1 = [data(:, 5) data(:, 6)];
    Thorax_2 = [data(:, 7) data(:, 8)];
    Head = Head_2 - Head_1;  % specifying lines
    Thorax = Thorax_2 - Thorax_1;
    theta1 = (atan2(Head(:, 2), Head(:, 1))); % datangle between Head and X-axis
    theta2 = (atan2(Thorax(:, 2), Thorax(:, 1))); % datangle between Thorax and X-axis
    thetahead = 180 * (unwrap(theta1)) / pi; % to make angles continous
    thetathorax = 180 * (unwrap(theta2)) / pi;
    theta = thetahead - thetathorax; % angle between head and thorax
    angle_data = [thetahead thetathorax theta];
    
end

   
function AnalyzeSineMotion(s, r, folder)
    
    fast_FT(s, r, 600, folder, 'sine');
    Correlation(s, r, 600, folder, 'sine');
    
end
function AnalyzeChirpMotion(s, r, folder)
    
    %fast_FT(s, r, 600, folder, 'chirp');
    Correlation(s, r, 600, folder, 'chirp');
    %Spec_Analysis((s, r, 600, 'chirp_bode');
    
end

function [stimulus, response] = ReadCSV(data_folder, treatment_dir, motion_type_data_file)

    fileID = fopen(strcat(data_folder,'\',...
                                           treatment_dir,'\',...
                                           motion_type_data_file));
    stimulus = cell2mat(textscan(fileID,'%*s %f %*[^\n]', 'Delimiter' , ','));
    fclose(fileID);
    fileID = fopen(strcat(data_folder,'\',...
                           treatment_dir,'\',...
                           motion_type_data_file));
    response = cell2mat(textscan(fileID,'%*s %*s %f', 'Delimiter' , ','));
    fclose(fileID);
    
end

       