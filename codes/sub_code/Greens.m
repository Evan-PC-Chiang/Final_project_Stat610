classdef Greens < handle
    %GREENS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xobs
        Zobs
        GreensFunction
        uxe_de
        uze_de
        SegPair
        de_rate_Number
        GFShear
    end
    
    properties (Hidden)
        Xobs
        Fault
        Middle
        H
        Lines
        SplitGroup
        Idx
        Endpoint
        zobs
        Mu = 1
        t = 1
        Range = 300
    end
    
    
    methods
        % Constructor
        function obj = Greens(Geometry, xobs, kind)
            %GREENS Construct an instance of this class
            % Geometry = Output from Geometry
            % slip_rate = Slip rates for each thrust faults
            % H = Elastic thickness
            obj.Fault = Geometry.Fault;
            obj.Lines = Geometry.Final_lines;
            obj.Middle = median(obj.Fault(:,2));
            obj.H = Geometry.H;
            [obj.SplitGroup, obj.Idx, obj.Endpoint] = obj.splipGroup();
            obj.xobs = xobs;
            obj.zobs = 0:abs(obj.H);
            [obj.GreensFunction, obj.Xobs, obj.Zobs] = obj.greensFault(kind);
            [obj.uxe_de, obj.uze_de, obj.de_rate_Number] = obj.deepFlat();
        end

        % update SegPair
        function addSegPair(obj, SegPair)
            try
                % Validate that SegPair is numeric and has 2 columns
                validateattributes(SegPair, {'numeric'}, {'2d', 'ncols', 2});
            catch
                error('greensShear:InvalidInput', ...
                      'Input must be a numeric matrix with exactly 2 columns.');
            end
            obj.SegPair = SegPair;
        end
        
        % Seperate detachment into several groups by fault.
        function [split_group, idx, endpoint] = splipGroup(obj)
            group = obj.Fault(:,5);
            endpoint = zeros(1,length(obj.Lines)-2); 
            for i = 2:length(obj.Lines)-1
                endpoint(i-1) = min(obj.Lines{i}(:,2));
            end
            detachment_line = obj.Lines{1}; % find the break points
            idx0 = [];
            for j = 1:length(endpoint)
                [~, idx1, ~] = intersect(detachment_line(:,2), endpoint(j));
                idx0 = [idx0; idx1];
            end
            idx = [idx0-1; length(detachment_line)-1];
            for i = 1:length(idx)
                group(1:idx(i)) = group(1:idx(i))-30;
            end
            split_group = group(group ~= max(group));
        end
        
        % Assign slip rates
        function seg_slip = assignSlips(obj,slip_rate)
            slips = zeros(size(obj.SplitGroup));
            id_thrust = [0;find(diff(obj.SplitGroup(max(obj.Idx)+1:length(obj.SplitGroup))) ~= 0) + max(obj.Idx); length(obj.SplitGroup)]; % changes
            % Assign slip rates on thrust faults
            for i = 1:length(id_thrust)-1
                slips(id_thrust(i)+1:id_thrust(i+1)) = slip_rate(i);
            end
            % Assign slip rates on detachment
            slips(1:max(obj.Idx)) = sum(abs(slip_rate));
            for i = fliplr(1:length(obj.Idx)-1)
                slips(1:obj.Idx(i)) = slips(1:obj.Idx(i)) - abs(slip_rate(i));
            end
            seg_slip = slips;
            %Assign_slip = [slips(unique(obj.Idx(1:end-1)+1)); slips(id_thrust(2:end-1)+1)];
        end
        % Build Greens for faults
        function [greens_fault, Xobs, Zobs] = greensFault(obj, kind)
            if nargin < 2 || isempty(kind)
                kind = 'firstRow';
            end
        
            validKinds = {'all', 'firstRow'};
            if ~ismember(kind, validKinds)
                error('Invalid value for "kind". Valid options are: %s', strjoin(validKinds, ', '));
            end
        
            [Xobs1, Zobs1] = meshgrid(obj.xobs, obj.zobs);
            switch kind
                case 'all'
                    Xobs = Xobs1;
                    Zobs = Zobs1;
                case 'firstRow'
                    Xobs = Xobs1(1, :);
                    Zobs = Zobs1(1, :);
            end

            pb = ProgressBar('Greens for Faults');
            uniqueVals = unique(obj.SplitGroup);
            ux_Fault = cell(length(obj.SplitGroup),1);
            uz_Fault = cell(length(obj.SplitGroup),1);
            uxe_Fault = cell(length(obj.SplitGroup),1);
            uze_Fault = cell(length(obj.SplitGroup),1);
            for i = 1:size(uniqueVals,1)
                indices = find(obj.SplitGroup == uniqueVals(i));
                for loop = 1:length(indices)
                    [ux1,uz1,ux1e,uz1e,x1] = drive_visc_ramp_body_return_elastic(abs(obj.H),obj.Fault(indices(loop),1),0,obj.Fault(indices(loop),4),obj.Fault(indices(loop),3),obj.Fault(indices(loop),6),obj.Mu,obj.t,obj.zobs); %shift tip of fault to x=0
                    x1 = x1+obj.Fault(indices(loop),2); % shift them back
                    [X1, Z1] = meshgrid(x1,obj.zobs);
                    ux_Fault{indices(loop)} = griddata(X1, Z1, ux1, Xobs, Zobs,'linear');
                    uz_Fault{indices(loop)} = griddata(X1, Z1, uz1, Xobs, Zobs,'linear');
                    uxe_Fault{indices(loop)} = griddata(X1, Z1, ux1e, Xobs, Zobs,'linear');
                    uze_Fault{indices(loop)} = griddata(X1, Z1, uz1e, Xobs, Zobs,'linear');
                    pb.update(indices(loop),length(obj.SplitGroup));
                end
            end
            greens_fault = containers.Map({'ux_Fault','uz_Fault','uxe_Fault','uze_Fault'},{ux_Fault,uz_Fault,uxe_Fault,uze_Fault});
            pb.complete();
        end
        
        % Make Shear zone
        function greensShear(obj)
            pb = ProgressBar('Greens for Shear');
            uxS_Shear = cell(size(obj.SegPair,1),1);
            uzS_Shear = cell(size(obj.SegPair,1),1);
            uxSe_Shear = cell(size(obj.SegPair,1),1);
            uzSe_Shear = cell(size(obj.SegPair,1),1);
            for i = 1:size(obj.SegPair,1)
                seg1 = [obj.Fault(obj.SegPair(i,1),:)];
                seg2 = [obj.Fault(obj.SegPair(i,2),:)];
                if seg1(1) > seg2(1)
                    seg_upper = seg2;
                    seg_lower = seg1;
                else
                    seg_upper = seg1;
                    seg_lower = seg2;
                end
                dip = 90 - (seg_upper(4) + seg_lower(4))/2;
                if dip < 0
                    dip_axial = 180 + dip;
                else
                    dip_axial = dip;
                end
                
                xtip_axial = seg_lower(2) + seg_lower(1)/tan(dip_axial*pi/180);
                L_axial = seg_lower(1)/sin(dip_axial*pi/180);
                [uxS,uzS,uxSe,uzSe,xS] = drive_visc_ramp_body_return_elastic(abs(obj.H),0,0,180-dip_axial,L_axial,-1,obj.Mu,obj.t,obj.zobs); %shift tip of fault to x=0, dip-180 dip and negative sign on slip for surface dipping down to left
                xS = xS+xtip_axial;
                [X1, Z1] = meshgrid(xS,obj.zobs);
                uxS_Shear{i} = griddata(X1, Z1, uxS, obj.Xobs, obj.Zobs,'linear');
                uzS_Shear{i} = griddata(X1, Z1, uzS, obj.Xobs, obj.Zobs,'linear');
                uxSe_Shear{i} = griddata(X1, Z1, uxSe, obj.Xobs, obj.Zobs,'linear');
                uzSe_Shear{i} = griddata(X1, Z1, uzSe, obj.Xobs, obj.Zobs,'linear');
                pb.update(i,size(obj.SegPair,1));
                
            end
            obj.GFShear = containers.Map({'ux_Shear','uz_Shear','uxe_Shear','uze_Shear'},{uxS_Shear,uzS_Shear,uxSe_Shear,uzSe_Shear});
            pb.complete();
        end

        function [uxe_de, uze_de,rowNumbers_de_rate] = deepFlat(obj)
            subsetIndices = 1:sum(obj.Fault(:,5) ~= max(obj.Fault(:,5)));
            rowNumbers_de_rate = subsetIndices(obj.Fault(subsetIndices, 7) == -obj.H);
            BottomDetach = obj.Fault(rowNumbers_de_rate,:);
            uxe_de = cell(height(BottomDetach),1);
            uze_de = cell(height(BottomDetach),1);
            for i = 1:height(BottomDetach())
                m = [BottomDetach(i,8) -abs(obj.H) 0 1*obj.t];
                [uxe_de{i}, uze_de{i}] = EdgeDisp(m,obj.Xobs,-obj.Zobs,0.25);
            end
        end

        function slip_shear = ShearV(obj,seg_slip)
            slip_shear = zeros(size(obj.SegPair,1),1);
            for i = 1:size(obj.SegPair,1)
                seg1 = [obj.Fault(obj.SegPair(i,1),:) seg_slip(obj.SegPair(i,1))];
                seg2 = [obj.Fault(obj.SegPair(i,2),:) seg_slip(obj.SegPair(i,2))];
                if seg1(1) > seg2(1)
                    seg_upper = seg2;
                    seg_lower = seg1;
                else
                    seg_upper = seg1;
                    seg_lower = seg2;
                end
                if seg_upper(4) < seg_lower(4)
                    slip_shear(i) = -1*(abs(seg_upper(end))+abs(seg_lower(end)))*sin(abs(seg_upper(4)-seg_lower(4))/2*pi/180);
                else
                    slip_shear(i) = 1*(abs(seg_upper(end))+abs(seg_lower(end)))*sin(abs(seg_upper(4)-seg_lower(4))/2*pi/180);
                end
                % dip = 90 - (seg_upper(4) + seg_lower(4))/2;
                % if dip < 0
                %     slip_shear{i} = (abs(seg_upper(end))+abs(seg_lower(end)))*sin(abs(seg_upper(4)-seg_lower(4))/2*pi/180);
                % else
                %     slip_shear{i} = -1*(abs(seg_upper(end))+abs(seg_lower(end)))*sin(abs(seg_upper(4)-seg_lower(4))/2*pi/180);
                % end
            end
        end

                % Make Shear zone
       

    end
    methods(Static)
        function plotShear(Geometry, SegPair)
                for i = 1:size(SegPair,1)
                    seg1 = [Geometry.Fault(SegPair(i,1),:)];
                    seg2 = [Geometry.Fault(SegPair(i,2),:)];
                    if seg1(1) > seg2(1)
                        seg_upper = seg2;
                        seg_lower = seg1;
                    else
                        seg_upper = seg1;
                        seg_lower = seg2;
                    end
                    dip = 90 - (seg_upper(4) + seg_lower(4))/2;
                    if dip < 0
                        dip_axial = 180 + dip;
                    else
                        dip_axial = dip;
                    end
                    xtip_axial = seg_lower(2) + seg_lower(1)/tan(dip_axial*pi/180);
                    line_shear = [xtip_axial 0; seg_lower(2) seg_lower(1)];
                    plot(line_shear(:,1),-line_shear(:,2),'k')
                end
        end
    end
end

