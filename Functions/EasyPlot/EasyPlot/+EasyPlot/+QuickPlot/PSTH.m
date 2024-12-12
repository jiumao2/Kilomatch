function PSTH(spike_times, trigger_times, varargin)
    % spike_times: 1xn double in ms
    % trigger_times: 1xn double in ms
    Width = 5;
    Height = 5;
    binWidth = 20;
    gaussian_kernel = 20;
    tPre = -500;
    tPost = 500;

    if nargin > 1
        for k = 1:2:size(varargin,2)
            if strcmpi(varargin{k},'Width')
                Width = varargin{k+1};
            elseif strcmpi(varargin{k},'Height')
                Height = varargin{k+1};
            elseif strcmpi(varargin{k},'tPre')
                tPre = varargin{k+1};
            elseif strcmpi(varargin{k},'tPost')
                tPost = varargin{k+1};
            elseif strcmpi(varargin{k},'binWidth')
                binWidth = varargin{k+1};
            elseif strcmpi(varargin{k},'gaussian_kernel')
                gaussian_kernel = varargin{k+1};
            else
                error('Wrong argument!');
            end
        end
    end

    fig = EasyPlot.figure("Height",50,'Width',50);
    ax = EasyPlot.axes(fig,...
        "Height",Height,...
        'Width',Width,...
        'MarginBottom',0.8,...
        'MarginLeft',0.8);    
    
    params.pre = -tPre;
    params.post = tPost;
    params.binwidth = binWidth;
    [psth, tpsth] = jpsth(spike_times, trigger_times, params);
    psth = smoothdata(psth,'gaussian',gaussian_kernel*5/binWidth);

    plot(ax,tpsth,psth,'k-');
    xlabel(ax,'time');
    ylabel(ax,'firing rate (Hz)');

    EasyPlot.cropFigure(fig);
    function [psth, tpsth, trialspxmat, tspkmat, trigtimes] = jpsth(spxtimes, trigtimes, params)
        % JY 8.9.2020
        % spktimes in ms
        % trigtimes in ms
        % params.pre, params.post, params.binsize
        
        pre = params.pre;
        post = params.post;
        trigtimes = round(trigtimes);
        if size(trigtimes,1) == 1
            trigtimes = trigtimes';
        end
        
        binwidth = params.binwidth; % also in msec
        
        n_events = length(trigtimes); % time of events
        
        spxtimes                 =      round(spxtimes); 
        tspkall                     =       1:max([spxtimes trigtimes'])+5000;
        spkall                      =      zeros(1, length(tspkall));
        spkall(spxtimes)      =     1; 
        trialspxmat               =     zeros(pre+post+1, n_events);
        
        for i = 1:n_events 
            if trigtimes(i)-pre>0 && trigtimes(i)+post < max(tspkall)
                trialspxmat(:, i) = spkall (trigtimes(i)-pre : trigtimes(i)+post);
            else
               trialspxmat(:, i) = NaN*ones(pre+post+1, 1);
            end
        end
        
        
        tspkmat = (-pre):post;
        
        [~, inan] = find(isnan(trialspxmat));
        
        spkmat = trialspxmat;
        spkmat(:, inan)=[];
        trialspxmat = spkmat;
        trigtimes(inan) = [];
         
        [spkhistos, ts ] = spikehisto(spkmat,1000, (pre+post)/binwidth);
        ts = ts*1000 - pre;
         
        psth = spkhistos;
        tpsth = ts; 
    end
end