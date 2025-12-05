classdef ProgressBar < handle
    %PROGRESSBAR A class to display and update a progress bar with a separate complete button
    %
    %   pb = ProgressBar(n) initializes a progress bar with n steps.
    %   pb.update(i) updates the progress to the i-th step, displaying "Processing i of n".
    %   pb.complete() displays a confirmation button when the progress is done.

    properties (Access = private)
        h          % Handle for the waitbar
        fontSize = 12;  % Fixed font size for the displayed text
        startTime % Timer start time
    end
    
    methods
        % Constructor
        function obj = ProgressBar(title)
            % Start the timer automatically
            obj.startTime = tic;
            % Initialize the progress bar window
            obj.h = waitbar(0, 'Calculation Starts......');
            % Set up initial window title and message
            set(obj.h, 'Name', title);
            % Set the font size for the waitbar text
            textHandle = findall(obj.h, 'Type', 'Text');
            set(textHandle, 'FontSize', obj.fontSize, 'FontName','Helvetica');

            drawnow;
        end
        
        % Update method
        function update(obj, i, n)
            % Calculate progress percentage
            progress = i / n;
            % Calculate elapsed time
            elapsedTime = toc(obj.startTime);
            % Update the waitbar with the current "i of n" message
            waitbar(progress, obj.h, sprintf('Processing %d of %d | Elapsed time: %.2f seconds', i, n, elapsedTime));
            

            % Set the font size for the waitbar text
            textHandle = findall(obj.h, 'Type', 'Text');
            set(textHandle, 'FontSize', obj.fontSize, 'FontName','Helvetica');
        end

        % Update method
        function update_per(obj, i, n)
            % Calculate progress percentage
            progress = i / n;
            % Calculate elapsed time
            elapsedTime = toc(obj.startTime);
            % Update the waitbar with the current "i of n" message
            waitbar(progress, obj.h, sprintf('Processing %.1f%% | Elapsed time: %.2f seconds', progress*100, elapsedTime));
            

            % Set the font size for the waitbar text
            textHandle = findall(obj.h, 'Type', 'Text');
            set(textHandle, 'FontSize', obj.fontSize, 'FontName','Helvetica');
        end
        
        % Complete method to show confirm button
        function complete(obj)
            % Update the waitbar message to indicate completion
            waitbar(1, obj.h, 'Complete!');
            
            % Get the position of the waitbar window
            waitbarPos = get(obj.h, 'Position');  % [left, bottom, width, height]

            % Set the button size
            buttonWidth = 100;
            buttonHeight = 25;

            % Calculate the button's centered position within the waitbar
            buttonX = (waitbarPos(3) - buttonWidth) / 2;
            buttonY = 10;  % Position slightly above the bottom

            % Calculate elapsed time
            elapsedTime = toc(obj.startTime);

            % Add a confirmation button in the center of the waitbar
            uicontrol('Style', 'pushbutton', 'String', sprintf('Run time: %.2f seconds', elapsedTime), ...
                'Position', [buttonX, buttonY, buttonWidth, buttonHeight], ...
                'Parent', obj.h, ...
                'Callback', 'delete(gcbf)');  % Close the waitbar when clicked

            % % remind me
            % fs = 44100;  
            % t = 0:1/fs:1.2;  
            % f1 = 440;  % 1: 440 Hz (A4)
            % f2 = 523;  % 2: 523 Hz (C5)
            % 
            % y1 = sin(2 * pi * f1 * t) ...               
            %     + 0.6 * sin(2 * pi * 2 * f1 * t) ...    
            %     + 0.4 * sin(2 * pi * 3 * f1 * t) ...    
            %     + 0.3 * sin(2 * pi * 4 * f1 * t);       
            % 
            % y2 = sin(2 * pi * f2 * t) ...               
            %     + 0.5 * sin(2 * pi * 2 * f2 * t) ...    
            %     + 0.3 * sin(2 * pi * 3 * f2 * t) ...    
            %     + 0.2 * sin(2 * pi * 4 * f2 * t);       
            % 
            % y = y1 + y2;
            % 
            % y = y / max(abs(y));
            % halfLength = round(length(t) / 2);
            % amplitudeRampUp = linspace(0, 1, halfLength);
            % amplitudeRampDown = linspace(1, 0, length(t) - halfLength);
            % amplitudeRamp = [amplitudeRampUp, amplitudeRampDown];
            % y = y .* amplitudeRamp;
            % 
            % sound(y, fs);
        end
        
        % Destructor
        function delete(obj)
            % Close the waitbar when the object is deleted
            if isvalid(obj.h)
                close(obj.h);
            end
        end
    end
end