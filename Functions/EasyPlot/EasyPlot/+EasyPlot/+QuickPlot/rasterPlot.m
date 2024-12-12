function rasterPlot(spike_times, event_times, varargin)
    % spike_times: 1xn_spikes array
    % event_times: 1xn_events array
    tPre = -1000;
    tPost = 1000;
    Width = 3;
    Height = 3;

    if nargin > 2
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'tPre')
                tPre = varargin{k+1};
            elseif strcmpi(varargin{k},'tPost')
                tPost = varargin{k+1};
            elseif strcmpi(varargin{k},'Width')
                Width = varargin{k+1};
            elseif strcmpi(varargin{k},'Height')
                Height = varargin{k+1};
            else
                error('Wrong argument!');
            end
        end
    end

    fig = EasyPlot.figure();
    ax = EasyPlot.axes(fig,...
        "Height", Height,...
        'Width', Width,...
        'MarginBottom', 1,...
        'MarginLeft', 1);    

    x_all = [];
    y_all = [];
    for k = 1:length(event_times)
        spike_times_this = spike_times(spike_times>=event_times(k)+tPre & spike_times<=event_times(k)+tPost)-event_times(k);
        if ~isempty(spike_times_this)
            numspikes=length(spike_times_this);
            x=ones(1,3*numspikes)*nan;
            y=ones(1,3*numspikes)*nan;
    
            y(1:3:3*numspikes)=-0.5+k;
            y(2:3:3*numspikes)=y(1:3:3*numspikes)+1;
            x(1:3:3*numspikes)=spike_times_this;
            x(2:3:3*numspikes)=spike_times_this; 

            x_all = [x_all, x];
            y_all = [y_all, y];
        end
    end
    plot(ax, x_all, y_all, 'k-', 'LineWidth', 1);
    xlabel(ax, 'time');
    ylabel(ax, 'Trials');
    xlim(ax, [tPre, tPost]);
    ylim(ax, [0.5, length(event_times)+0.5]);
    yticks(ax, [1, length(event_times)]);

    EasyPlot.cropFigure(fig);
end