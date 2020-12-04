function CellOut = convertstruct2cell(StructIn)
    % CellOut = Convertstruct2cell(StructIn)
    % converts a struct into a cell-matrix where the first column contains
    % the fieldnames and the second the contents
    CellOut = [fieldnames(StructIn),struct2cell(StructIn)]'; 
end