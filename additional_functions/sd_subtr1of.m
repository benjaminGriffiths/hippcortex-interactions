% This function subtracts 1 over f in time frequency data by fitting a
% linear function and writes saves the slope for each time point;
% Pow is the output of ft_frequanalysis with cfg.output='pow' and
% cfg.keeptrials='yes'; Note that this function will only work if the
% dimord = rpt_chan_freq_time:

% As an input it needs:
% cfg.toi which specifies the time windows which you want to have corrected
% e.g. [-0.5 3.5]
% 
% As an output the script generates the following variables:
% power = 1/f corrected power spectrum
% slp = ERP like structure with the 1/f slope for each time point

% Written by Simon Hanslmayr on 02/02/2017 

function [power,slp]=sh_subtr1of(cfg,pow)

switch pow.dimord
    case 'rpt_chan_freq_time'
        
        power = pow;
        
        t1=nearest(pow.time,cfg.toi(1));
        t2=nearest(pow.time,cfg.toi(2));
        pspc=pow.powspctrm(:,:,:,t1:t2);
        power.time=pow.time(t1:t2);
        
        freq=pow.freq;
        
        logfrq=log10(freq);
        logpow=log10(pspc);
        X = [ones(length(logfrq),1) logfrq'];% we add the ones so that the intercept is also fitted
        
        ntrls=size(pspc,1);
        nchans=size(pspc,2);
        nts=size(pspc,4);
        power.slope=zeros(ntrls,nchans,nts);
        power.powspctrm=NaN(ntrls,nchans,length(freq),nts);
        
        slpdo='rpt_chan_time';
        slp=struct('avg', zeros(nchans,nts), 'time',power.time,'label',{power.label},'trial', ...
            zeros(ntrls,nchans,nts), 'dimord', slpdo); 
        
        for n=1:ntrls
            for ch=1:nchans
                for t=1:nts
                    tmp_ps=squeeze(logpow(n,ch,:,t));
                    idx=find(~isnan(tmp_ps));
                    tmp_ps=tmp_ps(idx);
                    % Fit a linear function to log transformed power
                    b = X(idx,:)\tmp_ps;% Believe it or not, this simple comand does a linear regression; Handy, innit?
                    linft = X(idx,:)*b;% give the linear fit
                    slp.trial(n,ch,t)=b(2);
                    temp=tmp_ps-linft;
                    power.powspctrm(n,ch,idx,t) = 10.^temp;
                    power.slope(n,ch,t) = b(2);
                end
            end
        end
        slp.avg=squeeze(mean(slp.trial,1));

    case 'chan_freq'
        power = pow;
        
        pspc=pow.powspctrm;
        freq=pow.freq;
        
        logfrq=log10(freq);
        logpow=log10(pspc);
        X = [ones(length(logfrq),1) logfrq'];% we add the ones so that the intercept is also fitted
        
        nchans=size(pspc,1);
        power.slope=zeros(nchans,1);
        power.powspctrm=NaN(nchans,length(freq));
        
        slpdo='chan';
        slp=struct('avg', zeros(nchans,1), 'time',[],'label',{power.label},'trial', ...
            zeros(nchans,1), 'dimord', slpdo); 
        
        for ch=1:nchans
            tmp_ps=squeeze(logpow(ch,:));
            idx=find(~isnan(tmp_ps));
            tmp_ps=tmp_ps(idx)';
            % Fit a linear function to log transformed power
            b = X(idx,:)\tmp_ps;% Believe it or not, this simple comand does a linear regression; Handy, innit?
            linft = X(idx,:)*b;% give the linear fit
            slp.trial(ch)=b(2);
            temp=tmp_ps-linft;
            power.powspctrm(ch,idx) = 10.^temp;
            power.slope(ch) = b(2);
        end
        slp.avg=squeeze(mean(slp.trial,1));
        
        
    otherwise
        disp('Error: dimord must be rpt_chan_freq_time');
end