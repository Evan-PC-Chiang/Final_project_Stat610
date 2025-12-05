classdef Greens_v2 < handle
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
        ShearFault
        Build_gemo
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
        function obj = Greens_v2(build_geometry, xobs, kind)
            %GREENS Construct an instance of this class
            % Geometry = Output from Geometry
            % slip_rate = Slip rates for each thrust faults
            % H = Elastic thickness
            obj.Build_gemo = build_geometry;
            obj.Fault = obj.Build_gemo.Fault;
            obj.Lines = obj.Build_gemo.Final_lines;
            obj.ShearFault = obj.Build_gemo.ShearFault;
            obj.SegPair = obj.Build_gemo.SegPair;
            obj.de_rate_Number = find(obj.Build_gemo.rowNumbers_de_rate);
            obj.Middle = median(obj.Fault(:,2));
            obj.H = obj.Build_gemo.H;
            [obj.SplitGroup, obj.Idx, obj.Endpoint] = obj.splipGroup();
            obj.xobs = xobs;
            obj.zobs = 0:abs(obj.H);
            [obj.GreensFunction, obj.Xobs, obj.Zobs] = obj.greensFault(kind);
            [obj.uxe_de, obj.uze_de] = obj.deepFlat();
        end
        
        % Seperate detachment into several groups by fault.
        function [split_group, idx, endpoint] = splipGroup(obj)
            group = obj.Fault(:,5);
            endpoint = zeros(1,length(obj.Lines)-1); 
            for i = 2:length(obj.Lines)
                endpoint(i-1) = min([obj.Lines{i}(1,2) obj.Lines{i}(end,2)]);
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
            fault = obj.ShearFault;

            uxS_Shear = cell(size(obj.ShearFault,1),1);
            uzS_Shear = cell(size(obj.ShearFault,1),1);
            uxSe_Shear = cell(size(obj.ShearFault,1),1);
            uzSe_Shear = cell(size(obj.ShearFault,1),1);
            for i = 1:size(fault,1)
                target_shear = fault(i,:);
                xtip_axial = target_shear(2);
                L_axial = target_shear(3);
                [uxS,uzS,uxSe,uzSe,xS] = drive_visc_ramp_body_return_elastic(abs(obj.H),0,0,target_shear(4),L_axial,-1,obj.Mu,obj.t,obj.zobs); %shift tip of fault to x=0, dip-180 dip and negative sign on slip for surface dipping down to left
                xS = xS+xtip_axial;
                [X1, Z1] = meshgrid(xS,obj.zobs);
                uxS_Shear{i} = griddata(X1, Z1, uxS, obj.Xobs, obj.Zobs,'linear');
                uzS_Shear{i} = griddata(X1, Z1, uzS, obj.Xobs, obj.Zobs,'linear');
                uxSe_Shear{i} = griddata(X1, Z1, uxSe, obj.Xobs, obj.Zobs,'linear');
                uzSe_Shear{i} = griddata(X1, Z1, uzSe, obj.Xobs, obj.Zobs,'linear');
                pb.update(i,size(obj.ShearFault,1));
                
            end
            obj.GFShear = containers.Map({'ux_Shear','uz_Shear','uxe_Shear','uze_Shear'},{uxS_Shear,uzS_Shear,uxSe_Shear,uzSe_Shear});
            pb.complete();
        end

        function [uxe_de, uze_de] = deepFlat(obj)
            BottomDetach = obj.Build_gemo.Bottom_point(1:end-1,:);
            uxe_de = cell(height(BottomDetach),1);
            uze_de = cell(height(BottomDetach),1);
            for i = 1:height(BottomDetach())
                m = [BottomDetach(i,1) -abs(obj.H) 0 1*obj.t];
                [uxe_de{i}, uze_de{i}] = EdgeDisp(m,obj.Xobs,-obj.Zobs,0.25);
            end
        end

        function slip_shear = ShearV(obj,seg_slip)
            slip_shear = zeros(size(seg_slip,1),1);
            shear_norm = obj.Build_gemo.ShearNorm;
            for i = 1:size(obj.SegPair,1)
                v1 = shear_norm(i,1:2).*seg_slip(obj.SegPair(i,1));
                v2 = shear_norm(i,3:4).*seg_slip(obj.SegPair(i,2));
                v3 = norm(v1+v2);
                if obj.ShearFault(i,4)<90
                    slip_shear(i) = v3*shear_norm(i,5)*-1;
                else
                    slip_shear(i) = v3*shear_norm(i,5);
                end
            end
        end
    end
end

