% This script performs basic housekeeping tasks and path management. It is
% called by many of the scripts used to construct results.
% The results/figures folders can be changed according to user preference.

%% BASIC TIDYING UP
clear variables;
close all;
clc;

%% PATH MANAGEMENT
restoredefaultpath;
addpath(genpath('MAPSlite'));
addpath(genpath('Functions'));
addpath(genpath('Toolkit'));
addpath('Models');

%% SPECIFY FOLDERS FOR FIGURES, RESULTS ETC
figureFolder = ['Figures' filesep];
resultsFolder = ['Results' filesep];
inputsFolder = [resultsFolder 'Inputs' filesep];
outputsFolder = [resultsFolder 'Outputs' filesep];
dataFolder = ['Data' filesep];