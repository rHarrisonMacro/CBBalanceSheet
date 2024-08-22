function OBCsAreActive = update_non_policy_OBC_guess(...
    OBCsAreActiveGuess,OBCsAreActiveSim,diagnostics,method,guessUpdateHandle)
% This function updates the guess about the periods in which occasionally
% binding constraints are "active".
%
% Function to update logicals that dictate whether non-policy OBCs are
% active (=true) or not (=false).
%
% INPUTS
% -> OBCsAreActiveGuess: N*H matrix of logicals indicating whether 
% each of the N OBCs is "active" (=true) or "inactive" (=false) in each 
% simulation period t=1,...,H. This represents the "guess" (for the status
% of the non-policy OBCs in each period of the simulation) on which the
% simulation was constructed.
%
% -> OBCsAreActiveSim: N*H matrix of logicals indicating whether 
% each of the N OBCs is "active" (=true) or "inactive" (=false) in each 
% simulation period t=1,...,H. This is based on a test of the OBCs on the
% simulation output (using "check_simulation_for_non_pol_OBCs.m").
%
% -> diagnostics: structure with diagnostics of OBC checks, from 
% "check_simulation_for_non_pol_OBCs.m"
% 
% -> method: (optional) string specifying the method by which the
% OBCsAreActive guess is updated.
%
% -> guessUpdateHandle: (optional) function handle for user-defined 
% approach to updating OBCsAreActive guess.
%
% OUTPUTS
% -> OBCsAreActive: updated N*H matrix of logicals indicating 
% whether each of the N OBCs is "active" (=true) or "inactive" (=false) in 
% each simulation period t=1,...,H. The updated matrix constitutes a better
% guess for the solution.
%
% NOTES
% -> The function is designed to allow for additional methods for updating
% the guess for OBCs to be incorporated.
% 
% Author(s): Richard Harrison
% 
% This version: 15/12/20

%% CHECK INPUTS & SET DEFAULT METHOD IF REQUIRED
if nargin<3 || nargin>5
    error('Function requires 3, 4 or 5 inputs.');
elseif nargin==3
    method = 'fullUpdateInactiveFirst';
end

%% UPDATE THE GUESS ACCORDING TO SPECIFIED UPDATE SCHEME
switch method
    case {'fullUpdateActiveFirst','fullUpdateInactiveFirst',...
            'singleUpdateInactiveFirst'}
        if ~diagnostics.constraintsAssumptionMatchesSim
            switch method
                 case 'singleUpdateInactiveFirst'
                     OBCsAreActive = OBCsAreActiveGuess;
                     OBCsAreActive(diagnostics.largestViolation.row,...
                         diagnostics.largestViolation.col) = true;
                otherwise
                    OBCsAreActive = OBCsAreActiveSim;
                    % This updates the guess to impose OBCs in periods for 
                    % which the "inactive" contemporary slackness 
                    % conditions are violated based on the simulation 
                    % output.
            end
        else
            % If the constraints assumptions match the simulation, then
            % there must be invalid assumptions about active or inactive
            % states.
            switch method
                case 'fullUpdateInactiveFirst'
                    if any(any(diagnostics.invalidInactiveStates))
                        OBCsAreActive = flip_invalid_states(...
                            OBCsAreActiveGuess,diagnostics,...
                            'inactive');
                    elseif any(any(diagnostics.invalidActiveStates))
                        OBCsAreActive = flip_invalid_states(...
                            OBCsAreActiveGuess,diagnostics,...
                            'active');
                    else
                        error(['The OBC status in the simulation appears ',...
                            'to be inconsistent.']);
                    end
                case 'fullUpdateActiveFirst'
                    if any(any(diagnostics.invalidInactiveStates))
                        OBCsAreActive = flip_invalid_states(...
                            OBCsAreActiveGuess,diagnostics,...
                            'inactive');
                    elseif any(any(diagnostics.invalidActiveStates))
                        OBCsAreActive = flip_invalid_states(...
                            OBCsAreActiveGuess,diagnostics,...
                            'active');
                    else
                        error(['The OBC status in the simulation appears ',...
                            'to be inconsistent.']);
                    end
                case 'singleUpdateInactiveFirst'
                    if any(any(diagnostics.invalidActiveStates))
                        OBCsAreActive = switch_first_invalid_state(...
                            OBCsAreActiveGuess,diagnostics,...
                            'active');
                    elseif any(any(diagnostics.invalidInactiveStates))
                        OBCsAreActive = switch_first_invalid_state(...
                            OBCsAreActiveGuess,diagnostics,...
                            'inactive');
                    else
                        error(['The OBC status in the simulation appears ',...
                            'to be inconsistent.']);
                    end
                otherwise
                    error(['Update method "' method '" is not supported.']);
            end
        end
    case 'userDefined'
        if nargin<5
            error(['Function handle for "userDefined" update method ',...
                'was not supplied.']);
        end
        try
            OBCsAreActive = guessUpdateHandle(OBCsAreActiveGuess,...
                OBCsAreActiveSim,diagnostics);
        catch OBCerror
            disp(['Attempt to update guess for OBCs using user defined ',...
                'function failed with the following error message:']);
            throw(OBCerror);
        end
    otherwise
        error(['Update method "' method '" is not supported.']);
end

end

%% SUB-FUNCTION TO CHANGE INVALID ASSUMPTIONS STATE OF OBCS
function constraintsAreActive = ...
    flip_invalid_states(constraintsAreActiveGuess,diagnostics,type)

switch type
    case 'inactive'
        % Cases in which the guess incorrectly assumed that OBCs were
        % inactive are switched to assume that the OBCs are active.
        constraintsAreActive = constraintsAreActiveGuess;
        constraintsAreActive(diagnostics.invalidInactiveStates) = true;
    case 'active'
        % Cases in which the guess incorrectly assumed that OBCs were
        % active are switched to assume that the OBCs are inactive.
        constraintsAreActive = constraintsAreActiveGuess;
        constraintsAreActive(diagnostics.invalidActiveStates) = false;
    otherwise
        error(['The string argument "type" passed to this function must ',...
            'take the value of "active" or "inactive" only.']);
end

end


%% SUB-FUNCTION TO CHANGE INVALID ASSUMPTIONS STATE OF OBCS
function constraintsAreActive = ...
    switch_first_invalid_state(constraintsAreActiveGuess,diagnostics,type)

switch type
    case 'inactive'
        % Cases in which the guess incorrectly assumed that OBCs were
        % inactive are switched to assume that the OBCs are active.
        constraintsAreActive = constraintsAreActiveGuess;
        invalidInactiveInds = find(diagnostics.invalidInactiveStates);
        constraintsAreActive(invalidInactiveInds(1)) = true;
    case 'active'
        % Cases in which the guess incorrectly assumed that OBCs were
        % active are switched to assume that the OBCs are inactive.
        constraintsAreActive = constraintsAreActiveGuess;
        invalidActiveInds = find(diagnostics.invalidActiveStates);
        constraintsAreActive(invalidActiveInds(1)) = false;
    otherwise
        error(['The string argument "type" passed to this function must ',...
            'take the value of "active" or "inactive" only.']);
end

end

