clear
% mcmc inversion for slip rates and locling ratios
addpath ../tools/
addpath sub_code/

model = 'Elastic'; % Elastic or Plate
kind = 'Shear';  % NoShear or Shear

% setup bounds
upper_slip_bounds = [40 40 40 40];
lower_slip_bounds = [0 0 0 0];

% How many runs
niter = 20000;

% Load the Green's functions and data
load('Greens/2RampGeometries.mat')
GF.SegPair=GF.Build_gemo.SegPair;
GF_all.SegPair=GF_all.Build_gemo.SegPair;
%% end of input

data = GF.Build_gemo;
% locking
locking = 0.5*ones(sum(data.Fault(:,5) ~= max(data.Fault(:,5))),1);

% set up bounds for slip and locking separately, first column is lower bound, second column is upper bound
bounds = [[lower_slip_bounds'; 0*ones(length(locking),1)] [upper_slip_bounds'; ones(length(locking),1)]];

% set up the initial values
% setup slip-rate, same order with fault input, no input for detachment
slip_rate = (upper_slip_bounds + lower_slip_bounds)/2;


% get fault slip-rate
seg_slip = GF.assignSlips(slip_rate);
seg_slip(29) = seg_slip(14);
% get shear slip-rate
slip_shear = GF.ShearV(seg_slip);
% get velocity
[ux, uz, uz_lt, operator] = get_model_result_elastic(GF, seg_slip, locking, model,kind);

% get likelihood compare to original data
[logrho, model_st_vert, model_st_horz, model_lt_vert] = get_logrho_2D(ux, uz, uz_lt, data_matrix);

%% MCMC inversion
% set up the parameters
npara = length(slip_rate) + length(locking);
% set up the matrix to store the results
para = zeros(niter+1, npara);
% put initial values in it
para(1,:) = [slip_rate locking'];
% set up the matrix to store the likelihood
logrho_all = zeros(niter+1, 1);
logrho_all(1) = logrho;
% set up the matrix to store the model results
model_st_vert_all = zeros(niter+1, length(model_st_vert));
model_st_horz_all = zeros(niter+1, length(model_st_horz));
model_lt_vert_all = zeros(niter+1, length(model_lt_vert));
% ser up the step size as 10% of the bounds
step_size = 0.1*(bounds(:,2)-bounds(:,1))';
pb = ProgressBar('MCMC Inversion');
randomInt = randi([3500, 5400]);
% start the loop
for i = 1:niter
    % pick one parameter to update
    para_new = para(i,:);
    % randomly pick one parameter to update
    para_new(randi(npara)) = para(i,randi(npara)) + step_size(randi(npara))*randn;
    % check if the new parameter is within the bounds
    if all(para_new >= bounds(:,1)') && all(para_new <= bounds(:,2)')
        % get the new slip rate and locking
        slip_rate_new = para_new(1:length(slip_rate));
        locking_new = para_new(length(slip_rate)+1:end);
        % get the new slip rate
        seg_slip_new = GF.assignSlips(slip_rate_new);
        seg_slip_new(29) = seg_slip_new(14);
        % get the new shear slip rate
        slip_shear_new = GF.ShearV(seg_slip_new);
        % get the new velocity
        [ux_new, uz_new, uz_lt_new, operator_new] = get_model_result_elastic(GF, seg_slip_new, locking_new, model,kind);
        % get the new likelihood
        [logrho_new, model_st_vert_new, model_st_horz_new, model_lt_vert_new] = get_logrho_2D(ux_new, uz_new, uz_lt_new, data_matrix);
        % calculate the acceptance probability
        alpha = exp(logrho_new - logrho_all(i));
        % accept the new parameter with probability alpha
        if alpha>1
            accept = 1;
        else
            if rand < alpha
                accept = 1;
            else
                accept = 0;
            end
        end
    else
        accept = 0;
    end
    % update the parameter, likelihood, and model results
    if accept == 1
        para(i+1,:) = para_new;
        logrho_all(i+1) = logrho_new;
        model_st_vert_all(i+1,:) = model_st_vert_new;
        model_st_horz_all(i+1,:) = model_st_horz_new;
        model_lt_vert_all(i+1,:) = model_lt_vert_new;
    else
        para(i+1,:) = para(i,:);
        logrho_all(i+1) = logrho_all(i);
        model_st_vert_all(i+1,:) = model_st_vert_all(i,:);
        model_st_horz_all(i+1,:) = model_st_horz_all(i,:);
        model_lt_vert_all(i+1,:) = model_lt_vert_all(i,:);
    end
    if mod(i,randomInt) == 0
        pb.update_per(i,niter);
    else 
        continue
    end
end
pb.complete();

%% plot result (mean of last 100 run)
para_mean = mean(para(end-100:end,:));
slip_rate_mean = para_mean(1:length(slip_rate));
locking_mean = para_mean(length(slip_rate)+1:end);
seg_slip_mean = GF.assignSlips(slip_rate_mean);
seg_slip_mean(29) = seg_slip_mean(14);

[model_st_horz_mean, model_st_vert_mean, model_lt_vert_mean, ~] = get_model_result_elastic(GF_all, seg_slip_mean, locking_mean, model,kind);
model_st_horz_mean = model_st_horz_mean - model_st_horz_mean(1);

beautiful_plot_cross_result
