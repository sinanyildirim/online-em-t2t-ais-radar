function assert_prerequisites()
%ASSERT_PREREQUISITES Fail early when required MATLAB functions are missing.

required_functions = {
    'assignkbest'
    'exportgraphics'
    'gamrnd'
    'log_sum_exp'
    'mvnrnd'
    'pagemtimes'
    'pagetranspose'
    'randsample'
    'tensorprod'
    };

missing_functions = {};
for i = 1:numel(required_functions)
    function_name = required_functions{i};
    if exist(function_name, 'file') == 0 && exist(function_name, 'builtin') == 0
        missing_functions{end + 1} = function_name; %#ok<AGROW>
    end
end

if ~isempty(missing_functions)
    error('Missing required MATLAB functions: %s', strjoin(missing_functions, ', '));
end
end
