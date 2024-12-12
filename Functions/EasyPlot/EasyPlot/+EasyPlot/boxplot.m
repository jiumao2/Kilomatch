function boxplot(varargin)
    fields_to_save = {'UserData', 'Box', 'Color'};
    data_fields = cell(1, length(fields_to_save));
    for k = 1:length(fields_to_save)
        data_fields{k} = varargin{1}.(fields_to_save{k});
    end

    boxplot(varargin{:});

    for k = 1:length(fields_to_save)
        varargin{1}.(fields_to_save{k}) = data_fields{k};
    end
end