function place(axes_all, reference, location)
    if iscell(reference)
        top_ref = -Inf;
        bottom_ref = Inf;
        left_ref = Inf;
        right_ref = -Inf;
        for k = 1:size(reference, 1)
            for j = 1:size(reference, 2)
                top_this = reference{k,j}.Position(2)+reference{k,j}.Position(4)+reference{k,j}.UserData.MarginTop;
                bottom_this = reference{k,j}.Position(2)-reference{k,j}.UserData.MarginBottom;
                left_this = reference{k,j}.Position(1)-reference{k,j}.UserData.MarginLeft;
                right_this = reference{k,j}.Position(1)+reference{k,j}.Position(3)+reference{k,j}.UserData.MarginRight;
                
                top_ref = max(top_ref,top_this);
                bottom_ref = min(bottom_ref,bottom_this);
                left_ref = min(left_ref,left_this);
                right_ref = max(right_ref,right_this);
            end
        end         
    else
        top_ref = reference.Position(2)+reference.Position(4)+reference.UserData.MarginTop;
        bottom_ref = reference.Position(2)-reference.UserData.MarginBottom;
        left_ref = reference.Position(1)-reference.UserData.MarginLeft;
        right_ref = reference.Position(1)+reference.Position(3)+reference.UserData.MarginRight;
    end
    
    if ~iscell(axes_all)
        axes_all = {axes_all};
    end

    top = -Inf;
    bottom = Inf;
    left = Inf;
    right = -Inf;
    for k = 1:size(axes_all, 1)
        for j = 1:size(axes_all, 2)
            top_this = axes_all{k,j}.Position(2)+axes_all{k,j}.Position(4)+axes_all{k,j}.UserData.MarginTop;
            bottom_this = axes_all{k,j}.Position(2)-axes_all{k,j}.UserData.MarginBottom;
            left_this = axes_all{k,j}.Position(1)-axes_all{k,j}.UserData.MarginLeft;
            right_this = axes_all{k,j}.Position(1)+axes_all{k,j}.Position(3)+axes_all{k,j}.UserData.MarginRight;
            
            top = max(top,top_this);
            bottom = min(bottom,bottom_this);
            left = min(left,left_this);
            right = max(right,right_this);
        end
    end

    if strcmpi(location, 'left')
        dx = left_ref-right;
        dy = 0;
    elseif strcmpi(location, 'right')
        dx = right_ref-left;
        dy = 0;
    elseif strcmpi(location, 'top')
        dx = 0;
        dy = top_ref-bottom;
    elseif strcmpi(location, 'bottom')
        dx = 0;
        dy = bottom_ref-top;
    end
    
    EasyPlot.move(axes_all, 'dx', dx, 'dy', dy);
end