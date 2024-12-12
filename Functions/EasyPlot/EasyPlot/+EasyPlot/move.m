function move(axes_all, varargin)
    dx = 0;
    dy = 0;
    if nargin > 1
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'dx')
                dx = varargin{k+1};
            elseif strcmpi(varargin{k},'dy')
                dy = varargin{k+1};
            else
                error('Wrong argument!');
            end
        end
    end

    if iscell(axes_all)
        for k = 1:size(axes_all,1)
            for j = 1:size(axes_all,2)
                pos = axes_all{k,j}.Position;
                pos(1) = pos(1)+dx;
                pos(2) = pos(2)+dy;
                axes_all{k,j}.Position = pos;
            end
        end
    else
        pos = axes_all.Position;
        pos(1) = pos(1)+dx;
        pos(2) = pos(2)+dy;
        axes_all.Position = pos;        
    end
end