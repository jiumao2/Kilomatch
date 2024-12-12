function ytickangle(axes_all, angle)

if ~iscell(axes_all)
    ytickangle(axes_all, angle);
    return
end

for k = 1:size(axes_all,1)
    for j = 1:size(axes_all,2)
        ytickangle(axes_all{k,j}, angle);
    end
end

end