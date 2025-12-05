% this function is to get the likelihood of the posterior model compare to the original data
% input: ux, uz, data_matrix, locking
% output: logrho, st_vert, st_horz, lt_vert

function [logrho, model_st_vert, model_st_horz, model_lt_vert] = get_logrho_2D(ux, uz, uz_lt, data_matrix)
    % get data_type
    data_type = data_matrix(:,1);
    % find all the nan-nan value in the third and fifth column of data_matrix
    % data_type distance_from_start V_parrallel V_perpendicular V_vertical sigma_parrallel sigma_perpendicular sigma_vertical
    hor_data = data_matrix(:,3);
    ver_data = data_matrix(:,5);
    sigma_horz = ones(size(data_matrix(:,6)));
    sigma_vert = ones(size(data_matrix(:,8)));
    % % set floor
    threshold_hor = 2;
    threshold_ver = 2;
    sigma_horz(sigma_horz < threshold_hor) = threshold_hor;
    sigma_vert(sigma_vert < threshold_ver) = threshold_ver;
    sigma_horz(6) = 0.01;

    % get the data for st_vert and st_horz, and their sigma
    data_st_vert = ver_data(data_type == 1);
    data_st_horz = hor_data(data_type == 1);
    sigma_st_vert = sigma_vert(data_type == 1)*2;
    sigma_st_horz = sigma_horz(data_type == 1)*0.8;
    % same for lt data
    data_lt_vert = ver_data(data_type == 2);
    sigma_lt_vert = sigma_vert(data_type == 2)*0.6;

    % find the model data
    model_st_vert = uz(1,data_type == 1);% - uz(1);
    model_st_horz = ux(1,data_type == 1) - ux(1);
    model_lt_vert = uz_lt(1,data_type == 2);% - uz_lt(1);

    % calculate the total likelihood by vertical stacking the data as well as the model
    data_all = [data_st_vert; data_st_horz; data_lt_vert];
    model_all = [model_st_vert'; model_st_horz'; model_lt_vert'];
    sigma_all = [sigma_st_vert; sigma_st_horz; sigma_lt_vert];
    % calculate the log likelihood
    logrho = -.5*(data_all./sigma_all-model_all./sigma_all)'*(data_all./sigma_all-model_all./sigma_all);
end