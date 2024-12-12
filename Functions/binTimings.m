function [spike_train, t_spike_train] = binTimings(spike_times, binwidth, varargin)
    StartFromZero = 'on';
    tStart = NaN;
    t_edges = NaN;
    if nargin>=3
        for i=1:2:size(varargin,2)
            switch varargin{i}
                case 'StartFromZero'
                    StartFromZero = varargin{i+1}; 
                case 'tStart'
                    tStart = varargin{i+1}; 
                case 't_edges'
                    t_edges = varargin{i+1}; 
                otherwise
                    errordlg('unknown argument')
            end
        end
    end
    
    if isnan(t_edges)
        switch StartFromZero
            case 'on'
                tStart = 0;
            case 'off'
                if isnan(tStart)
                    tStart = spike_times(1);
                end
            otherwise
                errordlg('unknown argument')   
        end
        t_edges = tStart:binwidth:spike_times(end)+binwidth; 
    end
    
    t_spike_train = (t_edges(1:end-1)+t_edges(2:end))/2;
    
    spike_train = zeros(1,length(t_edges)-1);
    
    k = 1; j = 1;
    while j <= length(spike_times) && k < length(t_edges)
        if spike_times(j) < t_edges(1)
            j = j+1;
        elseif spike_times(j) < t_edges(k+1)
            spike_train(k) = spike_train(k)+1;
            j = j+1;
        else
            k = k+1;
        end
    end
        
end