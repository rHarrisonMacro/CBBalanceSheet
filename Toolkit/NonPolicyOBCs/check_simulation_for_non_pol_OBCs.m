function [simIsValid,OBCsAreActiveSim,diagnostics] = ...
    check_simulation_for_non_pol_OBCs(xSim,nonPolOBCinfo,...
    OBCsAreActive,constraintTolerance)
% This function checks whether a simulation satisfies occasionally binding
% constraints (OBCs), in the form of contemporary slackness conditions.
% See Section 6 of the paper for definitions of "active" and "inactive"
% constraints.
% 
% INPUTS
% ->    xSim: nx*H matrix of simulation paths
% ->    nonPolOBCinfo: structure containing the information about
% "inactive" and "active".
% ->    OBCsAreActive: N*H matrix of logicals indicating whether 
% each of the N OBCs is "active" (=true) or "inactive" (=false) in each 
% simulation period t=1,...,H.
% ->	constraintTolerance: (optional) positive scalar specifying
% numerical tolerance for checking contemporary slackness conditions.
%
% OUTPUTS
% ->    simIsValid: logical scalar indicating whether the simulation
% satisfies the contemporary slackness conditions implied by
% nonPolOBCstatus
% ->    OBCsAreActiveSim: N*H matrix of logicals indicating whether 
% each of the N OBCs *should be* "active" (=true) or "inactive" (=false) in 
% each simulation period t=1,...,H.
% ->    diagnostics: structure containing information about the specific
% contemporary slackness conditions, with fields:
%       -- inactiveConstraintsSatisfied: N*H matrix of logicals indicating
%       whether the inequality constraint of the "inactive" state of the
%       OBC is valid (i.e., if C_{i}xtilde_{i,t}>=0 for i=1,..,N and
%       t=1,...,H.
%       -- constraintsAssumptionMatchesSim: scalar logical indicating
%       whether the assumptions regarding the active OBCs are valid in the
%       simulation.
%       -- inactiveConstraintsValid: scalar logical indicating whether the
%       assumptions regarding the inactive OBCs are valid in the
%       simulation.
%
% Author(s): Richard Harrison
%
% This version: 16/12/20

%% SET DEFAULT CONSTRAINT TOLERANCE IF NECESSARY
if nargin < 3 || nargin>4
    error('Function requires 3 or 4 inputs.');
elseif nargin == 3
    constraintTolerance = 1e-10;
end
if constraintTolerance<0
    error('Constraint tolerance must be positive.');
end

%% EXTRACT SIMULATION PATHS FOR XTILDE
xtildeSim = xSim(nonPolOBCinfo.xtildeInds,:);

%% GATHER INFORMATION ABOUT THE SIMULATION AND NON-POLICY OBCs
H = size(xtildeSim,2);
N = size(nonPolOBCinfo.G,1);
kMat = repmat(nonPolOBCinfo.k,1,H) + repmat(constraintTolerance,N,H);
diagnostics = struct;
simIsValid = true;

%% BUILD INTERMEDIATE VARIABLES FOR LOGICAL CHECKS OF OBCS
Cx = xtildeSim(nonPolOBCinfo.inactive.VarInds,:);
Gx = xtildeSim(nonPolOBCinfo.active.VarInds,:);
% Inactive constraints first
CxEqualsZero = ~OBCsAreActive;
GxMinuskValues = Gx-kMat;
GxMinuskIsNonNegative = (GxMinuskValues>-constraintTolerance);
% Active constraints
GxEqualsk = OBCsAreActive;
CxIsNonNegative = (Cx>-constraintTolerance);

%% CHECK THE VALIDITY OF THE CONTEMPORARY SLACKNESS CONDITIONS INDIVIDUALLY
validInactiveStates = (CxEqualsZero & GxMinuskIsNonNegative);
invalidInactiveStates = (CxEqualsZero & ~GxMinuskIsNonNegative);
validActiveStates = (CxIsNonNegative & GxEqualsk);
invalidActiveStates = (~CxIsNonNegative & GxEqualsk);

%% RECORD THE STATUS OF THE CONSTRAINTS ACCORDING TO THE SIM
% Note that the logical returns false in the inactive state and true if the
% OBC is active: hence the logical indicating whether the constraint is
% inactive must be reversed.
OBCsAreActiveSim = ~validInactiveStates;

%% CHECKS FOR VALIDITY OF THE SIMULATION
% Check if assumptions about inactive constraints under which the simulation
% was produced are consistent with the simulation output
if any(any(OBCsAreActiveSim~=OBCsAreActive))
    simIsValid = false;
    constraintsAssumptionMatchesSim = false;
    [row,col] = find(ismember(GxMinuskValues, min(GxMinuskValues(:))));
    largestViolation.row = row;
    largestViolation.col = col;
else
    constraintsAssumptionMatchesSim = true;
end
% Check for any periods in which the contemporary slackness conditions are
% violated
if any(any(invalidActiveStates)) || any(any(invalidInactiveStates))
    simIsValid = false;
    contemporarySlacknessConditionsHold = false;
else
    contemporarySlacknessConditionsHold = true;
end 

%% PACK DIAGNOSTICS STRUCTURE WITH ADDITIONAL RESULTS
diagnostics.validActiveStates = validActiveStates;
diagnostics.invalidActiveStates = invalidActiveStates;
diagnostics.validInactiveStates = validInactiveStates;
diagnostics.invalidInactiveStates = invalidInactiveStates;
diagnostics.constraintsAssumptionMatchesSim = ...
    constraintsAssumptionMatchesSim;
diagnostics.contemporarySlacknessConditionsHold = ...
    contemporarySlacknessConditionsHold;
if ~simIsValid && ~constraintsAssumptionMatchesSim
    diagnostics.largestViolation = largestViolation;
end

end