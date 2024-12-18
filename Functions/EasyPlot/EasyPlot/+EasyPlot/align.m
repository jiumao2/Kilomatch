function align(axes_all, reference, position)
    if iscell(reference)
        if any(strcmpi({'left', 'right', 'top', 'bottom'}, position))
            idx = EasyPlot.Utils.getSideAxes(reference, position);
            reference = reference{idx};
        else
            top_ref = -Inf;
            bottom_ref = Inf;
            left_ref = Inf;
            right_ref = -Inf;
            for k = 1:size(reference, 1)
                for j = 1:size(reference, 2)
                    top_this = reference{k,j}.Position(2)+reference{k,j}.Position(4);
                    bottom_this = reference{k,j}.Position(2);
                    left_this = reference{k,j}.Position(1);
                    right_this = reference{k,j}.Position(1)+reference{k,j}.Position(3);
                    
                    top_ref = max(top_ref,top_this);
                    bottom_ref = min(bottom_ref,bottom_this);
                    left_ref = min(left_ref,left_this);
                    right_ref = max(right_ref,right_this);
                end
            end         
        end
    else
        top_ref = reference.Position(2)+reference.Position(4);
        bottom_ref = reference.Position(2);
        left_ref = reference.Position(1);
        right_ref = reference.Position(1)+reference.Position(3);
    end

    if any(strcmpi({'horizontalCenter', 'verticalCenter'}, position))
        if iscell(axes_all)
            top = -Inf;
            bottom = Inf;
            left = Inf;
            right = -Inf;
            for k = 1:size(axes_all, 1)
                for j = 1:size(axes_all, 2)
                    top_this = axes_all{k,j}.Position(2)+axes_all{k,j}.Position(4);
                    bottom_this = axes_all{k,j}.Position(2);
                    left_this = axes_all{k,j}.Position(1);
                    right_this = axes_all{k,j}.Position(1)+axes_all{k,j}.Position(3);
                    
                    top = max(top,top_this);
                    bottom = min(bottom,bottom_this);
                    left = min(left,left_this);
                    right = max(right,right_this);
                end
            end
        else
            top = axes_all.Position(2)+axes_all.Position(4);
            bottom = axes_all.Position(2);
            left = axes_all.Position(1);
            right = axes_all.Position(1)+axes_all.Position(3);
        end
    end

    if strcmpi(position, 'left')
        out = reference.Position(1);
        if ~iscell(axes_all)
            axes_all.Position(1) = out;
            return
        end
    
        idx_side = EasyPlot.Utils.getSideAxes(axes_all, position);
        dx = out-axes_all{idx_side}.Position(1);
        dy = 0;
    elseif strcmpi(position, 'right')
        out = reference.Position(1)+reference.Position(3);
        if ~iscell(axes_all)
            axes_all.Position(1) = out-axes_all.Position(3);
            return
        end
    
        idx_side = EasyPlot.Utils.getSideAxes(axes_all, position);
        dx = out-axes_all{idx_side}.Position(1)-axes_all{idx_side}.Position(3);
        dy = 0;
    elseif strcmpi(position, 'top')
        out = reference.Position(2)+reference.Position(4);
        if ~iscell(axes_all)
            axes_all.Position(2) = out-axes_all.Position(4);
            return
        end
    
        idx_side = EasyPlot.Utils.getSideAxes(axes_all, position);
        dx = 0;
        dy = out-axes_all{idx_side}.Position(2)-axes_all{idx_side}.Position(4);
    elseif strcmpi(position, 'bottom')
        out = reference.Position(2);
        if ~iscell(axes_all)
            axes_all.Position(2) = out;
            return
        end
    
        idx_side = EasyPlot.Utils.getSideAxes(axes_all, position);
        dx = 0;
        dy = out-axes_all{idx_side}.Position(2);
    elseif strcmpi(position, 'horizontalCenter')
        out = 0.5*(left_ref+right_ref);
        dx = out-0.5*(left+right);
        dy = 0;
    elseif strcmpi(position, 'verticalCenter')
        out = 0.5*(top_ref+bottom_ref);
        dx = 0;
        dy = out-0.5*(top+bottom);
    end
    
    EasyPlot.move(axes_all, 'dx', dx, 'dy', dy);
end