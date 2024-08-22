function UpdatedModel = solve_LSS_model(Model,theta)
% This macro solves LSS models in a standard way or using a bespoke solver.
% It updates the MAPS model structure for the solution to the model (and 
% the model's numeric structural matrices) consistent with a new parameter
% set or with the parameter set in the model object passed in. The default
% behaviour is for the model to be solved using a standard RE algorithm 
% (AIM). Alternatively, a custom solver can be added to the model object
% and the standard solver will be by-passed to a solver of the user's
% choice.

% INPUTS:   
%   -> Model: MAPS model structure
%   -> theta (optional): column vector of parameter values
%
% OUTPUTS:  
%   -> UpdatedModel: Updated MAPS model structure
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> is_finite_real_numeric_column_vector
%   -> unpack_model
%   -> <custom_solve>
%   -> generate_MAPS_exception_add_cause_and_throw
%   -> solve_LSS_model_steady_state
%   -> create_LSS_model_equation_numerics
%   -> solve_LSS_model_RE_dynamics
%   -> create_measurement_equation_numerics
%   -> compute_model_variable_covariance_matrix
%   -> pack_model
%
% DETAILS: 
%   -> If used in a standard way (without a custom solver), this macro
%      updates the model structure with a new RE solution consistent with
%      the parameter set passed in. It also saves all intermediate outputs
%      like the numeric structural matrices into the updated model object.
%   -> It is also possible to by-pass the standard solution implemented in
%      MAPS using a handle to a function that implements a bespoke solution
%      (eg flexible-price model conditional on the sticky-price state).
%   -> In either case, the aim is to produce the following numeric 
%      matrices: HB, HC, HF & PSI from the model equations, HB*x{t-1}+
%      HC*x{t}+HF*x{t+1} = PSI*z{t}; D, G & V from the measurement
%      equations, Y{t} = D+G*x{t}+V*w{t}; ss, the steady state solution;
%      B, PHI & F from the dynamic solution of the model, x{t} =
%      B*x{t-1}+PHI*u{t} with F used to discount future anticipated shocks;
%      P, the variance-covariance matrix of the model variables in steady 
%      state.
%      
% NOTES:
%   -> See <> for a description of the standard LSS solution algorithm in 
%      MAPS. 
%   -> See <> for a description of what custom solvers are and how to use 
%      them.
%
% This version: 16/05/2011
% Author(s): Alex Haberis and Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of inputs is as expected.
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin']; 
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin>1 && ~is_finite_real_numeric_column_vector(theta)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% UNPACK DEFAULT PARAMETERS
% If only 1 input was passed in, unpack the parameters from the model
% object as default.
if nargin < 2
    theta = unpack_model(Model,{'theta'}); 
end
                            
%% UNPACK MODEL CHARACTERISTICS INFO
% Unpack information about the model's characteristics to be used below to
% determine which parts of the code to execute.
[modelHasMeasurementEqs,modelHasMeasurementErrors,...
    modelHasSteadyStateEqs] = ...
    unpack_model(Model,{'modelHasMeasurementEqs',...
    'modelHasMeasurementErrors','modelHasSteadyStateEqs'});

%% SOLVE LSS MODEL
% Solve the model using either a custom solver (which can be used to solve 
% for steady state & RE dynamics in a bespoke manner) or using the standard 
% steady state procedure and rational expectations solver (AIM algorithm)
% in MAPS.
if isfield(Model.Numerics,'customSolver')
    custom_solve = Model.Numerics.customSolver;
    if ~strcmp(class(custom_solve),'function_handle')
        errId = ['MAPS:',mfilename,':BadCustomSolverModelField'];
        generate_and_throw_MAPS_exception(errId);
    else
        try
            [B,PHI,F,ss]  = custom_solve(Model,theta);
        catch CustomSolveE
            errId = ['MAPS:',mfilename,':CustomSolveE'];
            generate_MAPS_exception_add_cause_and_throw(...
                CustomSolveE,errId,{func2str(custom_solve)});
        end
    end
else
    if modelHasSteadyStateEqs
        SSfunHandle = unpack_model(Model,{'SSfunHandle'});
        ss = solve_LSS_model_steady_state(SSfunHandle,theta);
        thetaAndSS = {theta,ss};
    else
        thetaAndSS = {theta};
    end
    [HBfunHandle,HCfunHandle,HFfunHandle,PSIfunHandle] = unpack_model(...
        Model,{'HBfunHandle','HCfunHandle','HFfunHandle','PSIfunHandle'});
    [HB,HC,HF,PSI] = create_LSS_model_equation_numerics(thetaAndSS,...
            HBfunHandle,HCfunHandle,HFfunHandle,PSIfunHandle);
    [B,PHI,F] = solve_LSS_model_RE_dynamics(HB,HC,HF,PSI);  
end

%% CREATE MEASUREMENT EQUATION NUMERICS
% Create the measurement equation numerics in Y{t} = D + G*x{t}+V*w{t} with
% or without loadings on the measurement errors depending on whether the 
% model contains them or not. 
if modelHasMeasurementEqs
    [DfunHandle,GfunHandle] = unpack_model(...
        Model,{'DfunHandle','GfunHandle'});
    if modelHasMeasurementErrors
        VfunHandle = unpack_model(Model,{'VfunHandle'});
        [D,G,V] = create_measurement_equation_numerics(...
            thetaAndSS,DfunHandle,GfunHandle,VfunHandle);
    else
        [D,G] = create_measurement_equation_numerics(...
            thetaAndSS,DfunHandle,GfunHandle);
    end
end

%% CREATE COMPUTE MODEL VARIABLE COVARIANCE MATRIX
% If the model has data transformation equations, compute the covariance
% matrix of the model variables (in steady state) to initialise the Kalman
% filter. If, for whatever reason, the computation fails and the model is
% forward looking, throw an exception. If the model is backward-looking and
% the computation fails, set the covariance matrix equal to a scaled
% identity matrix and warn the user that the computation failed (possibly
% because the model has a unit root).
if modelHasMeasurementEqs
    try
        P = compute_model_variable_covariance_matrix(PHI,B);
    catch CovarComputeE
        modelIsForwardLooking = unpack_model(...
            Model,{'modelIsForwardLooking'});
        if ~modelIsForwardLooking
            warning(['MAPS:',mfilename,':CovarComputeOverwrite'],...
                ['It was not possible to compute the covariance ',...
                'matrix of the model variables. The routine to do so ',...
                'returned the following message: ',...
                CovarComputeE.message,' The matrix is instead being ',...
                'initialised as a scaled identity matrix - if this is ',...
                'not as expected, please carefully review the model ',...
                'for possible errors']);
            P = 0.01*eye(size(B));
        else
            rethrow(CovarComputeE);
        end
    end
end

%% PACK MODEL
% Update the model structure input with the structural matrices,
% parameters, steady states and solution computed above.
UpdatedModel = pack_model(Model,...
    {'B','PHI','F','HB','HC','HF','PSI','theta'},...
    {B,PHI,F,HB,HC,HF,PSI,theta});
if modelHasSteadyStateEqs
    UpdatedModel = pack_model(UpdatedModel,{'ss'},{ss});
end
if modelHasMeasurementEqs
    UpdatedModel = pack_model(UpdatedModel,{'D','G','P'},{D,G,P});
    if modelHasMeasurementErrors
        UpdatedModel = pack_model(UpdatedModel,{'V'},{V});
    end
end
 


end