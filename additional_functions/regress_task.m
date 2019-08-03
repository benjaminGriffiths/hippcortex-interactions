function res = regress_task(pow,task)

% get dimensions
input_size = size(pow);

% collapse in 2d matrix
pow = pow(:,:);
res = nan(size(pow));

% cycle through each point
for i = 1 : size(pow,2)
    
    % get beta
    B = task\pow(:,i);
    
    % get residual
    res(:,i) = (pow(:,i) - task*B);    
end

% reshape powspctrm
res = reshape(res,input_size);