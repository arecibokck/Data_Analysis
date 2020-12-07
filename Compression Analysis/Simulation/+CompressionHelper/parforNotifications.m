% Copyright (c) 2019 Andrea Alberti
%
% All rights reserved.
classdef parforNotifications < handle
    properties
        N;   % number of iterations
        text = 'Please wait ...';   % text to show
        width = 50;
        showWarning = true;
    end
    properties (GetAccess = public, SetAccess = private)
        n;
    end
    properties (Access = private)
        inProgress = false;
        percent;
        DataQueue;
        usePercent;
        Nstr;
        NstrL;
        lastComment;
    end
    methods
        function this = parforNotifications()
            this.DataQueue = parallel.pool.DataQueue;
            afterEach(this.DataQueue, @this.updateStatus);
        end
        % Start progress bar
        function PB_start(this,N,varargin)
            assert(isscalar(N) && isnumeric(N) && N == floor(N) && N>0, 'Error: ''N'' must be a scalar positive integer.');
            
            this.N = N;
            
            p = inputParser;
            addParameter(p,'message','Please wait: ');
            addParameter(p,'usePercentage',true);
            
            parse(p,varargin{:});
            
            this.text = p.Results.message;
            assert(ischar(this.text), 'Error: ''Message'' must be a string.');
            
            this.usePercent = p.Results.usePercentage;
            assert(isscalar(this.usePercent) && islogical(this.usePercent), 'Error: ''usePercentage'' must be a logical scalar.');
            
            this.percent = 0;
            this.n = 0;
            this.lastComment = '';
            if this.usePercent
                fprintf('%s [%s]: %3d%%\n',this.text, char(32*ones(1,this.width)),0);
            else
                this.Nstr = sprintf('%d',this.N);
                this.NstrL = numel(this.Nstr);
                fprintf('%s [%s]: %s/%s\n',this.text, char(32*ones(1,this.width)),[char(32*ones(1,this.NstrL-1)),'0'],this.Nstr);
            end
            
            this.inProgress = true;
        end
        % Iterate progress bar
        function PB_iterate(this,str)
            if nargin == 1
                send(this.DataQueue,'');
            else
                send(this.DataQueue,str);
            end
        end
        function warning(this,warn_id,msg)
            if this.showWarning
                msg = struct('Action','Warning','Id',warn_id,'Message',msg);
                send(this.DataQueue,msg);
            end
        end
        function PB_reprint(this)
            p = round(100*this.n/this.N);
            
            this.percent = p;
            
            cursor_pos=1+round((this.width-1)*p/100);
            
            if p < 100
                sep_char = '|';
            else
                sep_char = '.';
            end
            
            if this.usePercent
                fprintf('%s [%s%s%s]: %3d%%\n', this.text, char(46*ones(1,cursor_pos-1)), sep_char, char(32*ones(1,this.width-cursor_pos)),p);
            else
                nstr=sprintf('%d',this.n);
                fprintf('%s [%s%s%s]: %s/%s\n', this.text, char(46*ones(1,cursor_pos-1)), sep_char, char(32*ones(1,this.width-cursor_pos)),[char(32*ones(1,this.NstrL-numel(nstr))),nstr],this.Nstr);
            end
        end
        function updateStatus(this,data)
            
            if ischar(data)
                
                this.n = this.n + 1;
                
                p = round(100*this.n/this.N);
                
                if p >= this.percent+1 || this.n == this.N
                    this.percent = p;
                    
                    cursor_pos=1+round((this.width-1)*p/100);
                    
                    if p < 100
                        sep_char = '|';
                    else
                        sep_char = '.';
                    end
                    
                    if ~isempty(data)
                        comment = [' (',data,')'];
                    else
                        comment = '';
                    end
                    
                    if this.usePercent
                        fprintf('%s%s%s%s]: %3d%%%s\n',char(8*ones(1,58+numel(this.lastComment))), char(46*ones(1,cursor_pos-1)), sep_char, char(32*ones(1,this.width-cursor_pos)),p,comment);
                    else
                        nstr=sprintf('%d',this.n);
                        fprintf('%s%s%s%s]: %s/%s%s\n',char(8*ones(1,55+2*numel(this.Nstr)+numel(this.lastComment))), char(46*ones(1,cursor_pos-1)), sep_char, char(32*ones(1,this.width-cursor_pos)),[char(32*ones(1,this.NstrL-numel(nstr))),nstr],this.Nstr,comment)
                    end
                    
                    this.lastComment = comment;
                    
                    
                    if p == 100
                        this.inProgress = false;
                    end
                end
                
            else
                switch data.Action
                    case 'Warning'
                        warning(data.Id,[data.Message,newline]);
                        if this.inProgress
                            this.PB_reprint();
                        end
                end
                
            end
            
        end
    end
end


