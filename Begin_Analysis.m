function Begin_Analysis()
    path = 'E:\NCBS_Head_Stabilization_Project\Data Analysis';
    dirL = dir(path);
    dirL(~[dirL.isdir]) = [];
    r = ismember( {dirL.name}, {'.', '..'});
    dirL(r) = [];  %remove current and parent directory.

    i = 0;
    count = 0;
    while 1
        insect = strcat('.\Data\A_', num2str(i))
        fPath = strcat(strcat(path,'\'), insect)
        if isequal(exist(fPath, 'dir'),7)
            Analyze_Insect(insect);
            count = count + 1;
            if count == length(dirL)
                break;
            end
        end
        i = i + 1;
    end

end