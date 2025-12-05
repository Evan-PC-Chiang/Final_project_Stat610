classdef build_geometry < handle
    properties
        Lines
        Final_lines
        Fault
        Segment
        H
        extractedString
        Bottom_point
        SegPair
        ShearSegment
        ShearFault
        rowNumbers_de_rate
        ShearNorm
    end
    properties (Hidden, Transient)
        X
        Y
        Id
        Segment_pre
        Filename
        Fault_ori
        move_value
        Segment_ori
        Fault_tmp
    end

    methods
        function obj = build_geometry(filename,H)
            % Initialize properties
            obj.H = -H;
            obj.extractedString = regexp(filename, '(?<=/)[^/]+(?=\.)', 'match', 'once');
            [obj.X, obj.Y, obj.Id] = obj.readData(filename);
            obj.Lines = obj.lines_group();
            [obj.Final_lines, obj.Segment_pre] = obj.editDetach();
            obj.Fault_ori = build_geometry.build_parm(obj.Segment_pre);
            [obj.Fault_tmp, obj.Segment_ori, obj.move_value] = obj.move_it();
            [obj.Segment, obj.Bottom_point, obj.Fault, obj.rowNumbers_de_rate] = obj.add_detach();
        end

        % Read data in to [x, y, Id]
        function [x, y, Id] = readData(~, filename)
            % Open the file for reading
            fileID = fopen(filename, 'r');
    
            % Read the header line
            %header = fgetl(fileID);
            
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

        function [lines2, segments_fault] = editDetach(obj)
            % find detachment
            lines = obj.Lines;
            segments_fault = [];
            lines2 = [];
            % find last point of other faults
            % dist < 5, move fault. otherwise add a point in detachment
            for i = 2:length(obj.Lines)-2
                [~, minIndex] = min(lines{i}(:, 2));
                [nnsIdx, dist] = knnsearch(lines{1}, lines{i}(minIndex,:));
                if min(dist)<0.002
                    lines{i}(minIndex,:) = lines{1}(nnsIdx,:);
                else
                    lines{1} = [lines{1}(1:nnsIdx-1, :); lines{i}(minIndex,:); lines{1}(nnsIdx:end, :)];
                end
                lines{1} = sortrows(lines{1},1);
            end
            
            row_num = 0;
            for i = 1:length(lines)
                for j = 1:length(lines{i})-1
                    row_num = row_num + 1;
                    x1 = lines{i}(j,1);
                    y1 = lines{i}(j,2);
                    x2 = lines{i}(j+1,1);
                    y2 = lines{i}(j+1,2);
                    [x11, y11, x22, y22] = adjustLineToXAxis(x1, y1, x2, y2);
                    [x_1, y_1, x_2, y_2] = adjustLineToWish(x11, y11, x22, y22,obj.H);
                    if sum(isnan([x_1 y_1 x_2 y_2])) > 0
                        row_num = row_num - 1;
                    elseif x_1 == 0 && y_1 == 0
                        row_num = row_num - 1;
                    else
                        segments_fault(row_num,:) = [x_1 y_1 x_2 y_2 i];
                        lines2{i}(j,1) = x_1;
                        lines2{i}(j,2) = y_1;
                        lines2{i}(j+1,1) = x_2;
                        lines2{i}(j+1,2) = y_2;
                    end
                end
                matrix = lines2{i};
                matrix(all(matrix == 0, 2), :) = [];
                lines2{i} = matrix;
            end
        end

        % move Chelungpu to -32
        function [Snapped_Fault, segment, move_value] = move_it(obj)
            data = obj.Fault_ori;
            index = find(data(:,1) == 0, 1);
            move_value = data(index,2) + 32;
            data(:,2) = data(:,2) - move_value;
            Snapped_Fault = data;
            segment_data = obj.Segment_pre;
            segment_data(:,1) = segment_data(:,1) - move_value;
            segment_data(:,3) = segment_data(:,3) - move_value;
            segment = segment_data;
        end

        % add segment on the boundary of h
        function [segment, bottom_point, data_fault, index] = add_detach(obj)
            data_fault = obj.Fault_tmp;
            data = obj.Segment_ori;
            num = max(data(:,5))+1;
            index = data(:,4) == obj.H;
            bottom_point = [data(index,3:4); 100 -20];
            for i = 1:sum(index)
                data(end+1,:) = [bottom_point(i,:) bottom_point(i+1,:) num];
                data_fault(end+1,:) = [abs(bottom_point(i,2)) bottom_point(i,1) abs(bottom_point(i,1)-bottom_point(i+1,1)) 0 num 1];
            end
            segment = data;
        end

        function add_shear(obj, SegPair)
            try
                % Validate that SegPair is numeric and has 2 columns
                validateattributes(SegPair, {'numeric'}, {'2d', 'ncols', 2});
            catch
                error('greensShear:InvalidInput', ...
                      'Input must be a numeric matrix with exactly 2 columns.');
            end
            obj.SegPair = SegPair;
            shear_seg = ones(height(SegPair),5);
            shear_dir = ones(height(SegPair),1);
            shear_norm = ones(height(SegPair),4);
            for i = 1:size(SegPair,1)
                    seg1 = [obj.Segment(SegPair(i,1),:)];
                    seg2 = [obj.Segment(SegPair(i,2),:)];
                    point1 = [seg1(1), seg1(2)];
                    point2 = [seg1(3), seg1(4)];
                    point3 = [seg2(1), seg2(2)];
                    point4 = [seg2(3), seg2(4)]; 
                    % find common point and the other point for both of
                    % them
                    if isequal(round(point1,3), round(point3,3)) || isequal(round(point1,3), round(point4,3))
                        repeatedPoint = point1;
                    elseif isequal(round(point2,3), round(point3,3)) || isequal(round(point2,3), round(point4,3))
                        repeatedPoint = point2;
                    else
                        error(['Segment ', num2str(SegPair(i,1)), ' and ', num2str(SegPair(i,2)), ' are not connected!' ]);
                    end
                    if isequal(point1, repeatedPoint)
                        uniquePoint1 = point2;
                    else
                        uniquePoint1 = point1;
                    end
                    if isequal(point3, repeatedPoint)
                        uniquePoint2 = point4;
                    else
                        uniquePoint2 = point3;
                    end
                    % build the vector
                    vec1 = repeatedPoint - uniquePoint1;
                    vec2 = repeatedPoint - uniquePoint2;
                    vec1_norm = vec1/norm(vec1);
                    vec2_norm = vec2/norm(vec2);
                    tan_1 = abs(vec1_norm(1)/vec1_norm(2));
                    tan_2 = abs(vec2_norm(1)/vec2_norm(2));
                    if tan_1 > tan_2
                        shear_dir(i,:) = -1;
                    else
                        shear_dir(i,:) = 1;
                    end
                    shear_norm(i,:) = [vec1_norm vec2_norm];
                    % find bisector
                    bisector = vec1_norm + vec2_norm;
                    if bisector(2) < 0
                        bisector = -bisector;
                    end
                    % project to surface
                    scaling = (0 - repeatedPoint(2))/bisector(2);
                    surface_point = repeatedPoint + scaling.*bisector;
                    shear_seg(i,:) = [repeatedPoint surface_point(1) 0 i];
            end
            obj.ShearSegment = shear_seg;
            obj.ShearFault = build_geometry.build_parm(shear_seg);
            obj.ShearFault(:,6) = shear_dir;
            obj.ShearNorm = [shear_norm shear_dir];
        end


        function plot_lines(obj)
            figure; hold on
            for i = 1:length(obj.Final_lines)
                plot(obj.Final_lines{i}(:,1),obj.Final_lines{i}(:,2))
                text(mean(obj.Final_lines{i}(:,1)),mean(obj.Final_lines{i}(:,2)),int2str(i))
            end
            axis equal
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
                        plot(a,b,'r')
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

        function plot_Shear(obj)
            fault = obj.ShearFault;
            for k=1:size(fault,1)
                a=[fault(k,2) fault(k,2)+fault(k,3)*cos(fault(k,4)*pi/180)]; 
                b=-[fault(k,1) fault(k,1)+fault(k,3)*sin(fault(k,4)*pi/180)];
                plot(a,b,'k')

            end
        end

        function write_file(obj)
            file_str = [obj.X-obj.move_value obj.Y obj.Id];
            writematrix(file_str,['../nad_data/', obj.extractedString, '_moved.dat'],"Delimiter",'tab')
        end

    end
    methods (Static)
        function Fault = build_parm(data)
            
            % Initialize output arrays
            largerYPoints = zeros(size(data, 1), 2); % Points with larger y-values
            distances = zeros(size(data, 1), 1);    % Distances between points
            angles = zeros(size(data, 1), 1);       % Angles between points
            
            % Loop through each row
            for i = 1:size(data, 1)
                % Extract points
                x1 = data(i, 1);
                y1 = data(i, 2);
                x2 = data(i, 3);
                y2 = data(i, 4);
                
                % Determine the point with the larger y-value
                if y1 > y2
                    largerYPoints(i, :) = [abs(y1), x1];
                else
                    largerYPoints(i, :) = [abs(y2), x2];
                end
                
                % Compute distance between the two points
                distances(i) = sqrt((x2 - x1)^2 + (y2 - y1)^2);
                
                % Compute angle (in degrees) between the two points
                angles(i) = atan2d(y2 - y1, x2 - x1);
            end
            direction = ones(size(angles));
            direction(angles > 0) = -1;
            final_angle = -angles;
            final_angle(direction == -1) = final_angle(direction == -1)+180;
            direction_2 = ones(size(angles));
            direction_2(final_angle > 90) = -1;

            Fault = [largerYPoints distances final_angle data(:,5) direction_2];
        end
    end
end