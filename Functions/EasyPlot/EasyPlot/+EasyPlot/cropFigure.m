function cropFigure(fig)
    drawnow;
    top = -Inf;
    bottom = Inf;
    left = Inf;
    right = -Inf;
    for k = 1:length(fig.Children)
        top_child = fig.Children(k).Position(2)+fig.Children(k).Position(4)+fig.Children(k).UserData.MarginTop;
        bottom_child = fig.Children(k).Position(2)-fig.Children(k).UserData.MarginBottom;
        left_child = fig.Children(k).Position(1)-fig.Children(k).UserData.MarginLeft;
        right_child = fig.Children(k).Position(1)+fig.Children(k).Position(3)+fig.Children(k).UserData.MarginRight;
        
        top = max(top,top_child);
        bottom = min(bottom,bottom_child);
        left = min(left,left_child);
        right = max(right,right_child);
    end
    
    if isfield(fig.UserData,'OtherHandles')
        for k = 1:length(fig.UserData.OtherHandles)
            top_child = fig.UserData.OtherHandles{k}.Position(2)+fig.UserData.OtherHandles{k}.Position(4)+fig.UserData.OtherHandles{k}.UserData.MarginTop;
            bottom_child = fig.UserData.OtherHandles{k}.Position(2)-fig.UserData.OtherHandles{k}.UserData.MarginBottom;
            left_child = fig.UserData.OtherHandles{k}.Position(1)-fig.UserData.OtherHandles{k}.UserData.MarginLeft;
            right_child = fig.UserData.OtherHandles{k}.Position(1)+fig.UserData.OtherHandles{k}.Position(3)+fig.UserData.OtherHandles{k}.UserData.MarginRight;
            
            top = max(top,top_child);
            bottom = min(bottom,bottom_child);
            left = min(left,left_child);
            right = max(right,right_child);
        end
    end

    new_width = right-left+fig.UserData.MarginLeft+fig.UserData.MarginRight;
    new_height = top-bottom+fig.UserData.MarginTop+fig.UserData.MarginBottom;
    fig.Position(3) = new_width;
    fig.Position(4) = new_height;

    for k = 1:length(fig.Children)
        fig.Children(k).Position(1) = fig.Children(k).Position(1)-(left-fig.UserData.MarginLeft);
        fig.Children(k).Position(2) = fig.Children(k).Position(2)-(bottom-fig.UserData.MarginBottom);
    end
    
    if isfield(fig.UserData,'OtherHandles')
        for k = 1:length(fig.UserData.OtherHandles)
            fig.UserData.OtherHandles{k}.Position(1) = fig.UserData.OtherHandles{k}.Position(1)-(left-fig.UserData.MarginLeft);
            fig.UserData.OtherHandles{k}.Position(2) = fig.UserData.OtherHandles{k}.Position(2)-(bottom-fig.UserData.MarginBottom);    
        end
    end

    fig.Position(1) = EasyPlot.DefaultValue.PositionX;
    fig.Position(2) = EasyPlot.DefaultValue.PositionY;
end