function RunData = create_simulation_journey_base(...
    Model,H,T,Ytilde0,tt)
% This macro creates a LSS model simulation steady state base.
% It creates a complete MAPS linear state space (LSS) model forecast run 
% dataset consistent with the model input and with the assumptions
% underpinning LSS models.
%
% INPUTS:
%   -> Model: MAPS model structure
%   -> H: numeric scalar describing the forecast horizon
%   -> T (optional): number of periods to include in the back data
%   -> Ytilde0 (optional): initial conditions for raw observables
%   -> tt (optional/model dependent): complete set of time-varying trends 
%      over past & forecast
%
% OUTPUTS:
%   -> RunData: MAPS forecast run dataset structure
%
% DETAILS:
%   -> This macro creates a steady state simulation data base for
%      simulation in MAPS & EASE simulation journeys.
%   -> The dataset produced is a complete MAPS forecast dataset structure 
%      with a complete set of forecast data and forecast constraints, as 
%      well as T periods of past data (with a default of 3 if not input, 
%      chosen to allow for annual growth rates to be computed in every 
%      period of the forecast horizon).
%   -> Model variables & shocks are set equal to zero, consistent with the
%      LSS model being at steady state.
%   -> If the model input contains measurement equations, model observables 
%      are set equal model variables with the addition of the constants 
%      from the measurement equations.
%   -> If the model input contains data transformation equations and if 
%      the optional initial conditions input is not provided, then raw 
%      observables are given initial conditions that are set using a 
%      heuristic based on the sort of behaviour implicit in the data 
%      transformation equations.
%   -> Finally, if the model input includes time-varying trends in the data 
%      transformation equations and they are not provided on input, then 
%      they are set using another heuristic based on their effect in the 
%      data transformation equations. If the time-varying trends are
%      provided on input, then both must be provided (i.e. both the past &
%      forecast bits) in one matrix.
%   -> Consistent with standard assumptions about the effect of data
%      transformation in LSS model forecast runs (in "data" journeys), both 
%      the raw observable and time-varying trend simulation bases span one
%      additional period over the past.
%
% NOTES:
%   -> Note that the heuristics to generate the time-varying deterministic
%      trends do depend on some assumptions about the rules for data 
%      transformation equations. Changes to those rules may invalidate this
%      code.
%
% This version: 27/09/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~is_positive_real_integer(H)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin>2 && ~is_positive_real_integer(T)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);   
elseif nargin>3 && ~is_finite_real_numeric_column_vector(Ytilde0)
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin>4 && ~is_finite_real_two_dimensional_numeric_matrix(tt)
    errId = ['MAPS:',mfilename,':BadInput5'];
    generate_and_throw_MAPS_exception(errId);      
end

%% HANDLE OPTIONAL INPUTS
% The default value is for T=3, which guarantees that users will be able to
% compute annual growth rates of data over the whole forecast horizon. Note
% that this number applies to model variables and observables. Any raw
% observables and time-varying trends will have one additional data point 
% (reflecting that data transformation means 1 period of data is lost).
if nargin < 3
    T = 3;    
end
if nargin < 4
    initialConditionsForROprovided = false;
else
    initialConditionsForROprovided = true;
end
if nargin < 5
    timeVaryingTrendsProvided = false;
else
    timeVaryingTrendsProvided = true;
end

%% UNPACK MODEL CLASS INFO
modelIsLinearStateSpace = unpack_model(Model,{'modelIsLinearStateSpace'});
if ~modelIsLinearStateSpace
    errId = ['MAPS:',mfilename,':BadModelClass'];
    generate_and_throw_MAPS_exception(errId);
end

%% UNPACK MODEL CHARACTERISTICS INFO
[modelIsForwardLooking,modelHasMeasurementEqs,modelHasMeasurementErrors,...
    modelHasDataTransformationEqs,modelHasTimeVaryingTrends] = ...
    unpack_model(Model,{'modelIsForwardLooking',...
    'modelHasMeasurementEqs','modelHasMeasurementErrors',...
    'modelHasDataTransformationEqs','modelHasTimeVaryingTrends'});

%% CHECK CONSISTENCY OF MODEL CHARACTERISTICS WITH OPTIONAL INPUTS
if initialConditionsForROprovided && ~modelHasDataTransformationEqs
    errId = ['MAPS:',mfilename,':ROinitCondsInputInvalidWithModel'];
    generate_and_throw_MAPS_exception(errId);     
elseif timeVaryingTrendsProvided && ~modelHasTimeVaryingTrends
    errId = ['MAPS:',mfilename,':ttInputInvalidWithModel'];
    generate_and_throw_MAPS_exception(errId);        
end

%% CHECK OPTIONAL INPUTS ARE VALID GIVEN MODEL & OTHER INPUTS
if initialConditionsForROprovided
    YtildeMnems = unpack_model(Model,'YtildeMnems');
    if size(Ytilde0,1)~=size(YtildeMnems,1)
        errId = ['MAPS:',mfilename,...
            ':ROinitCondsInputIncompatibleWithModel'];
        generate_and_throw_MAPS_exception(errId);
    end
end
if timeVaryingTrendsProvided
    etatMnems = unpack_model(Model,'etatMnems');
    if size(tt,1) ~= size(etatMnems,1)
        errId = ['MAPS:',mfilename,':ttInputIncompatibleWithModel'];
        generate_and_throw_MAPS_exception(errId);
    elseif size(tt,2) ~= T+1+H
        errId = ['MAPS:',mfilename,':ttInputIncompatibleWithPeriodsInput'];
        generate_and_throw_MAPS_exception(errId);
    end
end

%% MODEL VARIABLES
% Compute the dimensions of the model variables and set the data over the
% past (the equivalent to the smoothed estimates of the model variables
% from the Kalman smoother) to zeros (and the forecast constraints to a 
% vector of zeros. Set the data over the forecast to a matrix of zeros.
PHI = unpack_model(Model,{'PHI'});
nx = size(PHI,1);
xs = zeros(nx,T);
xT = zeros(nx,1);
xf = zeros(nx,H);
RunData.Past.modelVariables = xs;
RunData.Constraint.modelVariables = xT;
RunData.Forecast.modelVariables = xf;

%% SHOCKS
% Compute the dimensions of the shocks and set the unanticipated shocks
% over the past to a vector of zeros (consistent with the assumption that
% model variables are at steady state made above). Set the unanticipated 
% forecast shocks to zero as a base. If the model is forward looking do the
% same for the anticipated shocks.
nz = size(PHI,2);
zs = zeros(nz,T);
u = zeros(nz,H);
RunData.Past.Shocks.unanticipated = zs;
RunData.Forecast.Shocks.unanticipated = u;
if modelIsForwardLooking
    a = zeros(nz,H);
    RunData.Forecast.Shocks.anticipated = a;
end

%% MODEL OBSERVABLES
% If the model contains model observables & measurement equations, set the
% model observables over the past and forecast to be equal to the constants 
% from the measurement equations (which is consistent with the base for 
% shocks and model variables described above). 
if modelHasMeasurementEqs
    D = unpack_model(Model,{'D'});
    Y = D*ones(1,T);
    Yf = D*ones(1,H);
    RunData.Past.modelObservables = Y;
    RunData.Forecast.modelObservables = Yf;
end

%% MEASUREMENT ERRORS
% If the model contains measurement errors, set their values over the past
% & forecast to zero.
if modelHasMeasurementErrors 
    V = unpack_model(Model,{'V'});
    nw = size(V,2);
    ws = zeros(nw,T);
    wf = zeros(nw,H);
    RunData.Past.measurementErrors = ws;
    RunData.Forecast.measurementErrors = wf;
end

%% RAW OBSERVABLES & TIME TRENDS
% If the model contains raw observables and data transformation equations
% set the raw observables over the past and forecast to be consistent with
% the assumptions made above. The level for the raw observables in the
% period prior to period 1 for the model observables above is arbitrary for 
% all raw observables with data transformation equations that include the 
% difference operator (since the level for the raw observables cannot be 
% determined from the one period of the model observables alone). The 
% heuristic for these are set in a sub-function below. If the model 
% contains deterministic time-varying trends another heuristic is used to 
% setup some dummy data in another sub-function below, which is again based
% on the way in which the time-varying trends enter the data transformation
% equations. See the sub-functions below for more details.
if modelHasDataTransformationEqs
    [DTfunHandle,RTfunHandle,YtildeEqStrs] = unpack_model(...
        Model,{'DTfunHandle','RTfunHandle','YtildeEqStrs'});
    if modelHasTimeVaryingTrends
        etatMnems = unpack_model(Model,{'etatMnems'});
        if ~timeVaryingTrendsProvided
            Ymnems = unpack_model(Model,{'Ymnems'});
            [etat0,etat,etatf] = set_time_varying_trends(...
                DTfunHandle,etatMnems,Ymnems,YtildeEqStrs,T,H);
        else
            etat0 = tt(:,1);
            etat = tt(:,2:T+1);
            etatf = tt(:,T+2:T+1+H);
        end
        etatT = etat(:,T);
        RunData.Past.timeVaryingTrends = [etat0 etat];
        RunData.Constraint.timeVaryingTrends = etatT;
        RunData.Forecast.timeVaryingTrends = etatf;
        if ~initialConditionsForROprovided
            Ytilde0 = set_raw_observables_initial_conditions(...
                RTfunHandle,Y,etat0,etat);
        end
        Ytilde = transform_observables_from_model_to_raw_space(...
            RTfunHandle,Y,Ytilde0,etat,etat0);
        YtildeT = Ytilde(:,T);
        Ytildef = transform_observables_from_model_to_raw_space(...
            RTfunHandle,Yf,YtildeT,etatf,etatT);
    else
        if ~initialConditionsForROprovided
            Ytilde0 = set_raw_observables_initial_conditions(...
                RTfunHandle,Y);
        end
        Ytilde = transform_observables_from_model_to_raw_space(...
            RTfunHandle,Y,Ytilde0);
        YtildeT = Ytilde(:,T);
        Ytildef = transform_observables_from_model_to_raw_space(...
            RTfunHandle,Yf,YtildeT);
    end 
    RunData.Past.rawObservables = [Ytilde0 Ytilde];
    RunData.Constraint.rawObservables = Ytilde(:,T);
    RunData.Forecast.rawObservables = Ytildef;
end

end

%% FUNCTION TO SET RAW OBSERVABLES INITIAL CONDITIONS
function Ytilde0 = set_raw_observables_initial_conditions(...
    RTfunHandle,Y,etat0,etat)
% This helper sets the initial conditions for the raw observables. 
% It uses a heuristic to set the initial conditions based on the output 
% from retransformations of the observables.
%
% INPUTS:
%   -> RTfunHandle: data retransformation function handle
%   -> Y: simulation base for model observables over the past
%   -> etat0 (model dependent): initial condition for time-varying trends
%   -> etat (model dependent): past simulation base for time trends
%
% OUTPUTS:
%   -> Ytilde0: initial condistions for the raw observables

%% INITIALSE THE RAW OBSERVABLES
% Set the raw observables in the initial period to an index of 100. 
nYtilde = size(Y,1);
Ytilde0 = 100*ones(nYtilde,1);

%% TRANSFORM THE MODEL OBSERVBALES
% Use the initial conditions for the raw observables and the past data for
% the model observables to compute a consistent series for the raw
% observables using the data retransformation symbolics.
if nargin > 2
    Ytilde = transform_observables_from_model_to_raw_space(...
        RTfunHandle,Y,Ytilde0,etat,etat0);
else
    Ytilde = transform_observables_from_model_to_raw_space(...
        RTfunHandle,Y,Ytilde0);    
end

%% CORRECT THE "LEVELS" RAW OBSERVBALES
% For raw obswervables that are related to the level of the corresponding
% model observable (ie those that exist in log or level operators), set the
% initial condition equal to the value found in the first period of the
% past on retransformation above (which will be comppletely independent of
% the level of the raw observables in the initial period.
DTlevelLogicals = ~any(diff(Ytilde')',2);
Ytilde0(DTlevelLogicals,1) = Ytilde(DTlevelLogicals,1);

end

%% FUNCTION TO SET THE TIME-VARYING TREND SIMULATION DATA
function [etat0,etat,etatf] = set_time_varying_trends(...
    DTfunHandle,etatMnems,Ymnems,YtildeEqStrs,T,H)
% This helper sets the simulation base for determinstic time trends. 
% It uses a heuristic to set the values for the time trends to constants
% (over the past and the forecast) depending on the effect they have on the
% observables.
%
% INPUTS:
%   -> DTfunHandle: data transformation function handle
%   -> Y: simulation base for model observables over the past
%   -> etatMnems: mnemonics for the time-varying trends
%   -> Ymnems: model observable mnemonics
%   -> YtildeEqStrs: data transformation equations
%   -> T: periods over the past
%   -> H: forecast horizon
%
% OUTPUTS:
%   -> etat0: simulation base initial conditions for the time trends
%   -> etat: simulation base past data for the time trends
%   -> etatf: simulation base forecast data for the time trends

%% COMPUTE DIMENSIONS OF INPUTS
% Compute the dimensions of the raw observables and time-varying trends
% input.
nYtilde = size(YtildeEqStrs,1);
netat = size(etatMnems,1);

%% INITIALISE TIME TRENDS
% Initialise the time trends for the initial period, the past and the
% forecast. Set these equal to zero, which is a netural assumption for
% additive time trends (eg a time-varying deterministic trend for 
% inflation).
etat0 = zeros(netat,1);
etat = zeros(netat,T);
etatf = zeros(netat,H);

%% COMPUTE MODEL OBSERVABLES FOR THE INITIALISED TIME TRENDS
% Compute the model observables for the initialised time trend given the 
% dummy data for raw observables input and the data transformation 
% equations. 
Ytilde0 = ones(nYtilde,1);
Ytildef = ones(nYtilde,1);
Yf = transform_observables_from_raw_to_model_space(...
    DTfunHandle,Ytildef,Ytilde0,etat(:,T),etat0);

%% FIND TRANSFORMED RAW OBSERVABLES EQUAL TO INFINITY
% Find the logical indices of any model observables that equal NaN or inf
% on transformation. These should be associatd with time trends that enter
% as part of a log operator or as part of as scalars in a diff operator 
% (ie Y = diff(Ytilde/etat) or Y = log(Ytilde/etat) or some combination).
YfInfLogics = any(~isfinite(Yf),2);
if ~any(YfInfLogics)
    return
end

%% FIND THE TIME TRENDS ASSOCIATED WITH THE BAD TRANSFORMED RAW OBSERVABLES
% Find the mnemonics of the time trends from the data transformation 
% equations which are associated with the bad model observables.
YetatIncCellStr = compute_equations_incidence_matrix_in_string_format(...
    YtildeEqStrs,etatMnems,Ymnems);
YfInfIncLogics = ismember(YetatIncCellStr(:,1),Ymnems(YfInfLogics));
etatToChange = unique(YetatIncCellStr(YfInfIncLogics,2));

%% CHANGE THE INVALID TIME TRENDS
% Lookup the indices of the time trends associated with the bad raw
% observables and change the defaults for these to 1 which will solve the 
% problem as long as they were the original cause of the bad
% transformations.
etatToChangeInd = lookup_model_index_numbers(etatMnems,etatToChange);
etat0(etatToChangeInd) = 1;
etat(etatToChangeInd,:) = 1;
etatf(etatToChangeInd,:) = 1;

end