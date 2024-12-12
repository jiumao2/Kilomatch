function addHandle(fig, handle)
    if isfield(fig.UserData, 'OtherHandles')
        fig.UserData.OtherHandles{end+1} = handle;
    else
        fig.UserData.OtherHandles{1} = handle;
    end
end