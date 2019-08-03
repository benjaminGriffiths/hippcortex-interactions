function output = sd_subtract_erp(data)

% duplicate input
output = data;

% get timelock of data
tml = ft_timelockanalysis([],data);

% cycle through each trial and subtract ero
for trl = 1 : numel(data.trial)
    output.trial{trl} = output.trial{trl} - tml.avg;
end