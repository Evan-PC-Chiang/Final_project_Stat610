function [selectedOptionA, selectedOptionB] = validateOptions(group1Options, group2Options, varargin)
    % Initialize selected options
    selectedOptionA = '';
    selectedOptionB = '';
    
    % Check input parameters
    for i = 1:length(varargin)
        % Check for options in Group 1
        if any(strcmp(varargin{i}, group1Options))
            if isempty(selectedOptionA)
                selectedOptionA = varargin{i};
            else
                error('Only one option from Group 1 can be selected.');
            end
        % Check for options in Group 2
        elseif any(strcmp(varargin{i}, group2Options))
            if isempty(selectedOptionB)
                selectedOptionB = varargin{i};
            else
                error('Only one option from Group 2 can be selected.');
            end
        else
            % Display a warning if an invalid option is found
            warning('Invalid option: "%s". This option will be ignored.', varargin{i});
        end
    end
    
    % Ensure that an option from each group is selected
    if isempty(selectedOptionA)
        error('One option from Group 1 must be selected.');
    end
    if isempty(selectedOptionB)
        error('One option from Group 2 must be selected.');
    end
end