function [ux,uz, uz_lt, operator] = get_model_result_elastic(GF, seg_slip, locking, varargin)
    % define valid model by str
    group1Options = {'Plate', 'Elastic'};
    group2Options = {'Shear', 'NoShear'};
    [selectedOptionA, selectedOptionB] = validateOptions(group1Options, group2Options, varargin{:});
    
    if strcmp(selectedOptionB, 'Shear')
        operator = 1;
    else
        operator = 0;
    end
    slip_shear = GF.ShearV(seg_slip);
    Get_rates4inv_st(GF.GreensFunction, seg_slip, locking);
    Get_rates4inv(GF.GreensFunction, seg_slip);
    Get_rates4inv(GF.GFShear, slip_shear);

    if strcmp(selectedOptionA, 'Plate')
        ux = ux_Fault_st + operator*ux_Shear;
        uz = uz_Fault_st + operator*uz_Shear;
        uz_lt = uz_Fault + operator*uz_Shear;
    else
        ux_de = zeros(size(GF.uxe_de{1}));
        uz_de = zeros(size(GF.uze_de{1}));
        for i=1:length(GF.de_rate_Number)
            ux_de = ux_de + GF.uxe_de{i}*seg_slip(GF.de_rate_Number(i));
            uz_de = uz_de + GF.uze_de{i}*seg_slip(GF.de_rate_Number(i));
        end
        ux = uxe_Fault_st + ux_de + operator*uxe_Shear;
        uz = uze_Fault_st + uz_de + operator*uze_Shear;
        uz_lt = uze_Fault + uz_de + operator*uze_Shear;
    end
end