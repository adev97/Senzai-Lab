function c = get_cell_type_color(ct)
    if strcmpi(ct, 'RGC')
        c = [0.85 0.15 0.15];  % red
    elseif strcmpi(ct, 'SC')
        c = [0.10 0.45 0.80];  % blue
    else
        c = [0.4 0.4 0.4];     % gray for unknown
    end
