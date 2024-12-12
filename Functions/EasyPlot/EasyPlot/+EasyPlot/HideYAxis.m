function HideYAxis(axes_all)
    if iscell(axes_all)
        for k = 1:size(axes_all,1)
            for j = 1:size(axes_all,2)
                axes_all{k,j}.YAxis.Visible = 'off';
            end
        end
    else
        axes_all.YAxis.Visible = 'off';
    end
end