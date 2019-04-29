function settings=add_default(settings, field, entry)
    if ~isfield(settings, field)
        settings.(field)=entry;
    end
end

