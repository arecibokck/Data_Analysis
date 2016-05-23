function [dirListing] = Analyze_Insect(insect)
    folder = 'E:\NCBS_Head_Stabilization_Project\Data Analysis\Data\';
    dirListing = dir(strcat(folder, insect));
    dirListing(~[dirListing.isdir]) = [];
    tf = ismember( {dirListing.name}, {'.', '..'});
    dirListing(tf) = [];  %remove current and parent directory.
    
    for i = 1:length(dirListing)
        if strcmp(dirListing(i).name,'Pre-Anaesthesia')
            currentdir = strcat(strcat(folder,insect),'\Pre-Anaesthesia');
            subList_Pre = dir(currentdir);
            tf = ismember( {subList_Pre.name}, {'.', '..'});
            subList_Pre(tf) = [];  %remove current and parent directory.
            
            %SINE
            data_sf = Arrange_Points(subList_Pre(1).name);
            Plot_SandR(data_sf, currentdir, 'stimulus_s.csv', 'response_s.csv');
            subList_Pre = dir(currentdir);
            tf = ismember( {subList_Pre.name}, {'.', '..'});
            subList_Pre(tf) = [];  %remove current and parent directory.
            s_index = 0;
            r_index = 0;
            for j = 1:length(subList_Pre)
                if strcmp(subList_Pre(j).name,'stimulus_s.csv')
                    s_index = j;
                elseif strcmp(subList_Pre(j).name, 'response_s.csv')
                    r_index = j;
                end
            end
            fast_FT(subList_Pre(s_index).name, subList_Pre(r_index).name, 600, currentdir, 'sine');
            %Correlation
            
            %CHIRP
            data_ch = Arrange_Points(subList_Pre(2).name);
            Plot_SandR(data_ch, currentdir, 'stimulus_ch.csv', 'response_ch.csv');
            subList_Pre = dir(currentdir);
            tf = ismember( {subList_Pre.name}, {'.', '..'});
            subList_Pre(tf) = [];  %remove current and parent directory.
            s_index = 0;
            r_index = 0;
            for j = 1:length(subList_Pre)
                if strcmp(subList_Pre(j).name,'stimulus_ch.csv')
                    s_index = j;
                elseif strcmp(subList_Pre(j).name, 'response_ch.csv')
                    r_index = j;
                end
            end
            fast_FT(subList_Pre(s_index).name, subList_Pre(r_index).name, 600, currentdir, 'chirp');
            %Correlation
            %Spec_Analysis(subList_Pre(s_index).name, subList_Pre(r_index).name, currentdir, 'chirp_bode');
            %delete(strcat(currentdir,'\stimulus_s.csv'), strcat(currentdir,'\response_s.csv'), strcat(currentdir,'\stimulus_ch.csv'), strcat(currentdir,'\response_ch.csv'));
        
        elseif strcmp(dirListing(i).name,'Post-Anaesthesia')
            currentdir = strcat(strcat(folder,insect),'\Post-Anaesthesia');
            subList_Pre = dir(currentdir);
            tf = ismember( {subList_Pre.name}, {'.', '..'});
            subList_Pre(tf) = [];  %remove current and parent directory.
            
            %SINE
            data_sf = Arrange_Points(subList_Pre(1).name);
            Plot_SandR(data_sf, currentdir, 'stimulus_s.csv', 'response_s.csv');
            subList_Pre = dir(currentdir);
            tf = ismember( {subList_Pre.name}, {'.', '..'});
            subList_Pre(tf) = [];  %remove current and parent directory.
            s_index = 0;
            r_index = 0;
            for j = 1:length(subList_Pre)
                if strcmp(subList_Pre(j).name,'stimulus_s.csv')
                    s_index = j;
                elseif strcmp(subList_Pre(j).name, 'response_s.csv')
                    r_index = j;
                end
            end
            fast_FT(subList_Pre(s_index).name, subList_Pre(r_index).name, 600, currentdir, 'sine');
            %Correlation
            
            %CHIRP
            data_ch = Arrange_Points(subList_Pre(2).name);
            Plot_SandR(data_ch, currentdir, 'stimulus_ch.csv', 'response_ch.csv');
            subList_Pre = dir(currentdir);
            tf = ismember( {subList_Pre.name}, {'.', '..'});
            subList_Pre(tf) = [];  %remove current and parent directory.
            s_index = 0;
            r_index = 0;
            for j = 1:length(subList_Pre)
                if strcmp(subList_Pre(j).name,'stimulus_ch.csv')
                    s_index = j;
                elseif strcmp(subList_Pre(j).name, 'response_ch.csv')
                    r_index = j;
                end
            end
            fast_FT(subList_Pre(s_index).name, subList_Pre(r_index).name, 600, currentdir, 'chirp');
            %Correlation
            %Spec_Analysis(subList_Pre(s_index).name, subList_Pre(r_index).name, currentdir, 'chirp_bode');
            %delete(strcat(currentdir,'\stimulus_s.csv'), strcat(currentdir,'\response_s.csv'), strcat(currentdir,'\stimulus_ch.csv'), strcat(currentdir,'\response_ch.csv'));
        end
    end
end