classdef Geometry_self_improved
    properties
        Lines
        Simple_lines
        Final_lines
        Fault
        H
    end

    properties (Hidden)
        X
        Y
        Id
        Filename
        AngleThreshold
    end

    methods
        % Constructor
        function obj = Geometry_self(filename,angleThreshold,H)
            % Initialize properties
            obj.H = -H;
            obj.Filename = filename;
            obj.AngleThreshold = angleThreshold;
            [obj.X, obj.Y, obj.Id] = obj.readData(filename);
            obj.Lines = obj.lines_group();
            obj.Simple_lines = obj.simplifyLines();
            [obj.Final_lines, obj.Fault] = obj.editDetach();
            
        end

        % Read data in to [x, y, Id]
        function [x, y, Id] = readData(~, filename)
            % Open the file for reading
            fileID = fopen(filename, 'r');
    
            % Read the header line
            header = fgetl(fileID);
            
            % Initialize arrays to store the data
            x = [];
            y = [];
            Id = [];
            
            % Read the data line by line
            while ~feof(fileID)
                % Read a line of data
                dataLine = fgetl(fileID);
                
                % Split the line into its components
                data = sscanf(dataLine, '%f %f %f');
                
                % Check if data has the expected number of columns (3 numerical values)
                if length(data) == 3
                    % Append the data to the respective arrays
                    x = [x; data(1)];
                    y = [y; data(2)];
                    Id = [Id; data(3)];
                end
            end
            
            % Close the file
            fclose(fileID);
        end

        % Group them by Id
        function lines = lines_group(obj)
            % Read the data from the file

            % Find unique Ids
            uniqueIds = unique(obj.Id);
            
            % Initialize a cell array to hold the lines
            lines = cell(length(uniqueIds), 1);
            
            % Loop over each unique Id
            for i = 1:length(uniqueIds)
                % Find indices corresponding to the current Id
                indices = find(obj.Id == uniqueIds(i));
                
                % Group the points corresponding to the current Id
                line = [obj.X(indices), obj.Y(indices)];
                
                % Store the line in the cell array
                lines{i} = line;
            end
        end

        % calculate angle
        function angle = calculateAngle(~, v1, v2)
            cosTheta = dot(v1, v2) / (norm(v1) * norm(v2));
            cosTheta = max(min(cosTheta, 1), -1);
            angle = acosd(cosTheta);
        end

        % Simplify line
        function simple_lines = simplifyLines(obj)
            % Initialize the cell array for simplified lines
            simple_lines = cell(size(obj.Lines));
            % loop through obj.Lines, skip detachment and the last one
            for i = 1:length(obj.Lines)
                if i == 1
                    simple_lines{i} = obj.Lines{i};
                elseif i == length(obj.Lines)
                    simple_lines{i} = obj.Lines{i};
                else
                % skip fault with only 2 points, sort by Y value
                    if length(obj.Lines{i}) == 2
                        simple_lines{i} = obj.Lines{i};
                    else
                        if obj.Lines{i}(end,2)>obj.Lines{i}(1,2)
                            tmp = flipud(obj.Lines{i});
                        else
                            tmp = obj.Lines{i};
                        end
                        j = 2; % start from second point, while end point remains the same.
                        points_to_remove = [];
                        while j <= length(tmp)-1
                            first_seg = tmp(j-1,:) - tmp(j,:);
                            second_seg = tmp(j,:) - tmp(j+1,:);
                            % calculate the degree
                            angle = obj.calculateAngle(first_seg, second_seg);
                            % remains the same if angle > angleThreshold
                            if angle < obj.AngleThreshold
                                points_to_remove = [points_to_remove; j];
                            else
                                points_to_remove = points_to_remove;
                            end
                            j=j+1;
                        end
                        tmp(points_to_remove,:) = [];
                        simple_lines{i} = tmp;
                    end
                end
            end
        end

        % Snap fault to detachment by either add point on detachment or
        % nearest point.
        function [lines2, fault] = editDetach(obj)
            % find detachment
            loop = 1;
            fault = [];
            lines = obj.Simple_lines;
            lines2 = [];
            % find last point of other faults
            % dist < 5, move fault. otherwise add a point in detachment
            for i = 2:length(obj.Lines)-3
                [nnsIdx, dist] = knnsearch(lines{1}, lines{i}(end,:));
                if min(dist)<5
                    lines{i}(end,:) = lines{1}(nnsIdx,:);
                else
                    lines{1} = [lines{1}(1:nnsIdx-1, :); lines{i}(end,:); lines{1}(nnsIdx:end, :)];
                end
            end
            for i = 1:length(lines)
                for j = 1:length(lines{i})-1
                    if i == length(lines)
                        x1 = lines{i}(j,1);
                        y1 = lines{i}(j,2);
                        x2 = lines{i}(j+1,1);
                        y2 = lines{i}(j+1,2);
                    else
                        x11 = lines{i}(j,1);
                        y11 = lines{i}(j,2);
                        x22 = lines{i}(j+1,1);
                        y22 = lines{i}(j+1,2);
                        [x1, y1, x2, y2] = adjustLineToXAxis(x11, y11, x22, y22);
                        [x1, y1, x2, y2] = adjustLineToWish(x1, y1, x2, y2,obj.H);
                    end
                    if sum(isnan([x1, y1, x2, y2])) == 0
                        lines2{i}(j,1) = x1;
                        lines2{i}(j,2) = y1;
                        lines2{i}(j+1,1) = x2;
                        lines2{i}(j+1,2) = y2;
                        fault(loop,1) = -y1;
                        fault(loop,2) = x1;
                        if y1 >= y2
                            if atand((y1-y2)/(x2-x1)) > 0
                                fault(loop,4) = mod(atand((y1-y2)/(x2-x1)),180);
                                fault(loop,6) = 1;
                            else
                                fault(loop,4) = mod(180+atand((y1-y2)/(x2-x1)),180);
                                fault(loop,6) = -1;
                            end
                        else
                            fault(loop,4) = 180+mod(atand((y1-y2)/(x2-x1)),180);
                        end
                        
                        fault(loop,3) = sqrt((x1-x2)^2+(y1-y2)^2);
                        fault(loop,5) = i;
                        fault(loop,7) = -y2;
                        fault(loop,8) = x2;
                        fault(loop,9) = loop;

                        loop = loop+1;
                    else 
                        continue
                    end
                end
            end
        end

        function plot_geometry(obj)
            fault = obj.Fault;
            figure; hold on
            for k=1:size(fault,1)
                a=[fault(k,2) fault(k,2)+fault(k,3)*cos(fault(k,4)*pi/180)]; 
                b=-[fault(k,1) fault(k,1)+fault(k,3)*sin(fault(k,4)*pi/180)];
                    if k <= sum(fault(:,5) ~= max(fault(:,5)))
                        plot(a,b,'r')
                        text(mean(a),mean(b),num2str(k))
                    else
                        plot(a,b,'c','LineStyle','--')
                        text(mean(a),mean(b)-0.5,num2str(k))
                    end

            end
            ylim([obj.H 0]);
            x1 = min(xlim); x2 = max(xlim);
            plot([x1 x2], [0 0])
            plot([x1 x2], [obj.H obj.H],'b',"LineStyle","--")
            axis equal
            pbaspect([3 1 1]);
            set(gcf, 'Position', [100, 100, 2200, 600]);
        end


    end
end