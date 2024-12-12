function setMargin(axes_all, position, value)
    % value: [MarginTop MarginLeft MarginBottom MarginRight] or 1x1 double
    if length(value) == 4
        if iscell(axes_all)
            for k = 1:length(axes_all)
                axes_all{k}.UserData.MarginTop = value(1);
                axes_all{k}.UserData.MarginLeft = value(2);
                axes_all{k}.UserData.MarginBottom = value(3);
                axes_all{k}.UserData.MarginRight = value(4);
            end
        else
            axes_all.UserData.MarginTop = value(1);
            axes_all.UserData.MarginLeft = value(2);
            axes_all.UserData.MarginBottom = value(3);
            axes_all.UserData.MarginRight = value(4);  
        end
        return
    end
    
    if strcmpi(position, 'left')
        fieldName = 'MarginLeft';
    elseif strcmpi(position, 'right')
        fieldName = 'MarginRight';
    elseif strcmpi(position, 'top')
        fieldName = 'MarginTop';
    elseif strcmpi(position, 'bottom')
        fieldName = 'MarginBottom';
    end

    if iscell(axes_all)
        for k = 1:length(axes_all)
            axes_all{k}.UserData.(fieldName) = value;
        end
    else
        axes_all.UserData.(fieldName) = value; 
    end    

end
