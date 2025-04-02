function options = SetDefaultOptions(options, default_field_value_pairs)
%SETDEFAULTOPTIONS
for i = 1 : length(default_field_value_pairs)
    field = default_field_value_pairs{i}{1};
    if(~isfield(options, field))
        value = default_field_value_pairs{i}{2};
        options.(field) = value;
    end
end

end