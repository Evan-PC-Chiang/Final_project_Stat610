clear
close all

addpath ../tools/
addpath sub_code/

filename = '../data/2RampGeometries.dat';
h = 20;

data = build_geometry(filename,h);

SegPair = [15 16; 16 1; 17 18; 20 3; 22 4; 4 5; 5 6; 6 7; 8 9; 9 10; 10 11; 11 12; 12 13; 23 13; 14 29; 21 22; 26 27];


data.add_shear(SegPair)

data.plot_geometry
data.plot_Shear
%%
% load data 
load('../data/2D_profile_data.mat','data_matrix_profile_output')
data_matrix_profile_output = data_matrix_profile_output(data_matrix_profile_output(:,8)~=9999,:);

ref_data = [1 -80 0 0 0 1 1 1];
data_matrix = [ref_data;data_matrix_profile_output];

data_type = data_matrix(:,1);

xobs_all = [linspace(-1*80,80,81)]; % for plot
xobs = data_matrix(:,2); % for inversion

GF_all = Greens_v2(data,xobs_all,'all');
GF_all.greensShear();
GF = Greens_v2(data,xobs,'firstRow');
GF.greensShear();

matFilename = strcat('Greens/', data.extractedString, '.mat');
save(matFilename,'GF','GF_all','data_matrix');