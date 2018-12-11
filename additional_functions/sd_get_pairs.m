function chancmb = sd_get_pairs(data1,data2)

% --- ORIGINAL FUNCTION --- %
n1 = numel(data1.label);
n2 = numel(data2.label);

chancmb = [];
for i = 1 : n1
    for j = 1 : n2
        chancmb(end+1,:) = [i,j];
    end
end
% 
% % get number of channels
% n1 = numel(data1.label);
% n2 = numel(data2.label);
% 
% % predefine matrix for channel pairs
% chancmb = [];
% 
% % cycle through neocortical channels
% for i = 1 : n1
%     
%     % cycle through hippocampal channels
%     for j = 1 : n2
%         
%         % check that channels do not come from same shaft
%         if ~strcmpi(data1.label{i}(8:9),data2.label{j}(8:9))
%             chancmb(end+1,:) = [i,j];
%         end
%     end
% end

%if isempty(chancmb); error('No pairs found...'); end