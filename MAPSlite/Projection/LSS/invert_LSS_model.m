function NewShocks = invert_LSS_model(B,PHI,F,xT,Shocks,Instructions)
% This forecast module implements a LSS model inversion.
% It inverts the solution of a linear state space (LSS) model to find the
% the values for the anticipated, and/or unanticipated shocks, that deliver
% judgemental forecasts for model variables.
%
% INPUTS:
%   -> B: nx*nx matrix of backward loadings on the model variables
%   -> PHI: nx*nz matrix of loadings on the shocks
%   -> F: nx*nx matrix of forward loadings on the anticipated shocks
%   -> xT: nx*1 vector of initial conditions for the model variables
%   -> Shocks: structure containing an existing set of shocks
%       - anticipated (model dependent/optional): nz*H matrix of
%         anticipated shocks
%       - unanticipated (optional): nz*H matrix of unanticipated shocks
%   -> Instructions: structure detailing the inversion required:
%       - algorithmIsMinimumDamage: true or false
%       - AnticipatedFixes (optional): row cell array with columns denoting
%         periods (e.g. AnticipatedFixes.modelVarIndices{1} corresponds to
%         period 1 etc):
%           - modelVarIndices: indices of the model variables to be fixed
%           - modelVarValues: values for the fixes
%           - shockUsageIndices: indices of shocks to implement the
%             anticipated fixes
%       - UnanticipatedFixes (optional): as for anticipated fixes:
%           - modelVarIndices: indices of the model variables to be fixed
%           - modelVarValues: values for the fixes
%           - shockUsageIndices: indices of shocks to implement the
%             unanticipated fixes
%
% OUTPUTS:
%   -> NewShocks: structure containing a new set of shocks
%       - anticipated (model/input dependent): nz*H matrix of shocks
%       - unanticipated (input dependent): nz*H matrix of shocks
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> validate_LSS_model_solution_matrices
%   -> interrogate_LSS_model_variable_initial_conditions
%   -> interrogate_forecast_shocks_structure
%   -> validate_inversion_instructions (sub-function)
%   -> create_inversion_info (sub-function)
%   -> check_validity_of_inversion (sub-function)
%   -> compute_reorder_operators (sub-function)
%   -> reorder_model_solution (sub-function)
%   -> check_feasibility_of_inversion (sub-function)
%   -> stack_fixes_and_shocks (sub-function)
%   -> compute_indexes_of_reordered_variables_and_shocks (sub-function)
%   -> partition_and_stack_model_solution (sub-function)
%   -> compute_matrices_in_simultaneous_system (sub-function)
%   -> combine_and_reorder_shocks (sub-function)
%
% DETAILS:
%   -> This forecast module implements a LSS model inversion. It inverts a
%      LSS model to find the set of shocks that deliver a set of exactly or
%      over-identified anticipated and/or unanticipated fixes.
%   -> It contains two algoirthms that can be used to deliver over-
%      identified fixes (i.e. using more instruments (shocks) than targets
%      (fixes)): a minimum damage algorithm which seeks to minimise the sum
%      squared deviation of the shocks used to implement the inversion from
%      their existing values (algorithmIsMinimumDamage = true); a minimum
%      variance algorithm which seeks to minimise the sum squared values of
%      the shocks used to implement the inversion regardless of their
%      previous values (algorithmIsMinimumDamage = false).
%   -> In either case, the algorithm works by reordering the model
%      solution, partitioning it into vectors that have values that are
%      known and then inverting to find the unknowns (which are the values
%      for the shocks that deliver the inversion).
%   -> This module will throw an exception if: the inversion instructions
%      contain anticipated fixes and the input does not include anticipated
%      shocks or the model used is backward-looking (the F matrix is
%      zeros); the inversion requested is unanticipated and the input does
%      not include unanticipated shocks; the inversion violates nay of the
%      basic inversion rules (e.g. must have at least as many instruments 
%      as targets etc), or if the inversion is infeasible (based on rank 
%      conditions - eg trying to fix a forcing process using any of the 
%      shocks that do not drive that process).
%
% NOTES:
%   -> See <> for more details on MAPS linear state space forecast
%      functionality. See also <> for a detailed description of LSS model
%      inversion.
%   -> Note that the algorithm is coded for the general case, but it
%      implicitly allows for specific cases. For example, if the inversion
%      is an unanticipated inversion, then the vector of anticipated fixes
%      made over model variables is set to zeros(0,1) in all S periods of
%      the inversion. Similarly, if no anticipated shocks are input
%      (perhaps because the model is backward-looking), then the matrix of
%      anticipated shocks is set to zeros over the whole forecast horizon.
%      The matrix algebra follows through correctly without having to
%      litter the code with a bunch of "if" "else" or "switch" "case"
%      statements.
%
% This version: 14/03/2011
% Author(s): Matt Waldron

%% CHECK NUMBER OF INPUTS
% Check that the number of inputs passed in is as expected.
if nargin < 6
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
end

%% VALIDATE SOLUTION MATRICES
% Call a MAPS module helper to validate that the solution matrices passed
% in have the expected shape and are consistent with each other.
validate_LSS_model_solution_matrices(B,PHI,F);

%% VALIDATE THE MODEL VARIABLE INITIAL CONDITIONS
% Call a MAPS module helper to validate that the model variable initial
% conditions. If successful, the function returns a numeric scalar
% denoting the the dimensionality of the model variables.
nx = interrogate_LSS_model_variable_initial_conditions(xT,B);

%% VALIDATE THE SHOCKS INPUT
% Call a MAPS forecast module helper function to validate the existing
% forecast shocks input. If successful, the function returns a structure
% describing the content of the shocks structure and numeric scalars
% denoting the forecast horizon and the dimensionality of the shocks.
[ShocksInclude,H,nz] = interrogate_forecast_shocks_structure(Shocks,PHI,F);

%% SETUP MODEL INDICES FOR VARIABLES & SHOCKS
% Setup a vector of model indices for the model variables and shocks.
xInd = (1:nx)';
zInd = (1:nz)';

%% VALIDATE THE INVERSION INSTRUCTIONS
% Call a sub-function to check the inversion instructions input structure.
% If sucessful, the function returns another structure, describing the
% content of the input.
InstructionsInclude = validate_inversion_instructions(...
    Instructions,H,xInd,zInd,ShocksInclude,F);

%% CREATE INVERSION INFO
% Call a sub-function to create all the information necessary to complete
% the inversion. The function unpacks the inputs that exist in the
% instructions. Note that it also creates "empty" inputs, where the
% particular type of input did not exist in the inversion instructions.
% For example, if "InstructionsInclude.unanticipatedFixes" = false, then xu
% & xuInd are set to 1*S cell arrays of zeros(0,1).
[xa,xaInd,xu,xuInd,an,anInd,un,unInd,abari,aiInd,ubari,uiInd] = ...
    create_inversion_info(...
    Instructions,InstructionsInclude,Shocks,ShocksInclude,H,zInd);

%% DETERMINE DIMENSIONALITY OF INPUTS
% Compute the dimensionality of each of the partitioned variables in each
% of the periods of the inversion, S. For example, nxa is a 1*S row vector
% with the number of model variables to be fixed using anticipated shocks
% in period, s, given by nxa(s).
nxa = cellfun('size',xa,1);
nxu = cellfun('size',xu,1);
nxn = nx-nxa-nxu;
nan = cellfun('size',an,1);
nun = cellfun('size',un,1);
nai = nz-nan;
nui = nz-nun;

%% CHECK THAT THE INVERSION IS VALID
% Check that the information unpacked above constitutes a valid inversion
% in the sense that it meets all of the inversion rules.
check_validity_of_inversion(nxa,xaInd,nxu,xuInd,nai,aiInd,nui,uiInd,H);

%% COMPUTE THE OPERATOR MATRICES REQUIRED TO REORDER THE SYSTEM
% Call a sub-function to compute the matrices required to reorder the
% model solution appopriately in each of the periods of the inversion. For
% example, Tx is a 1*S cell array with Tx{s} containing an nx*nx matrix to
% reorder the model variables in period s so that:
% Tx{s}*x = xT{s} = [xa{s};xu{s};xn{s}].
[Tx,Ta,Tu] = compute_reorder_operators(nx,nxa,nxu,xInd,xaInd,xuInd,...
    nz,zInd,nai,aiInd,nan,anInd,nui,uiInd,nun,unInd,H);

%% COMPUTE THE REORDERED MODEL SOLUTION MATRICES
% Call a sub-function to reorder the model solution appopriately in each of
% the periods of the inversion. For example, BT is a 1*S cell with BT{s}
% containing B reordered appropriately in period s so that:
% BT{s} = [Ba{s}a{s-1} Ba{s}u{s-1} Ba{s}n{s-1};
%          Bu{s}a{s-1} Bu{s}u{s-1} Bu{s}n{s-1};
%          Bn{s}a{s-1} Bn{s}u{s-1} Bn{s}n{s-1}].
[BT,FTkTimesPHITa,PHITu] = reorder_model_solution(...
    B,F,PHI,Tx,Ta,Tu,H,nx);

%% CHECK FEASIBILITY OF THE REQUIRED INVERSIONS
% Check that the anticipated & unanticipated inversions are feasible. They
% are not feasible if the model variables to be inverted are not linearly
% independent or if the shocks chosen to implement the inversions cannot
% simultaneously affect all the variables.
check_feasibility_of_inversion(FTkTimesPHITa,PHITu,nxa,nxu,nai,nui,H);

%% STACK VALUES FOR KNOWN ENDOGENOUS VARIABLES & SHOCKS
% Call a sub-function to stack fixes and shocks not being used in the
% inversion across the whole horizon S in individual vectors. For example,
% xaVec = [xa{1};xa{2};...;xa{S-1};xa{s}].
[xaVec,xuVec,xalagVec,xulagVec,anVec,unVec,zbariVec] = ...
    stack_fixes_and_shocks(xa,xu,nxa,nxu,an,un,abari,ubari,H);

%% COMPUTE INDICES OF COMPONENTS OF REORDERED SYSTEM
% Call a sub-function to compute the indexes of the reordered model
% variables & shocks and the reordered variables & shocks in the stacked
% partitioned system. For example, xaInd & xalagInd are 1*S cells such that
% xaInd{s} & xalagInd{s} refer to the indexes of xa{s} & xalag{s-1} in
% xT{s} & xT{s-1} respectively. While xaVecInd & xalagVecInd are 1*S cells
% such that xaInd{s} & xalagInd{s} refer to the indexes of xa{s} & xalag{s}
% in the stacked vectors xaVec & xalagVec.
[xTaInd,xTuInd,xTnInd,xTalagInd,xTulagInd,xTnlagInd,xTaVecInd,...
    xTuVecInd,xTnVecInd,xTalagVecInd,xTulagVecInd,xTnlagVecInd,aTiInd,...
    aTnInd,aTiVecInd,aTnVecInd,uTiInd,uTnInd,uTiVecInd,uTnVecInd] = ...
    compute_indexes_of_reordered_variables_and_shocks(...
    nx,nxa,nxu,nxn,nai,nan,nui,nun,H);

%% PARTITION AND STACK MODEL SOLUTION
% Call a sub-function to partition the reordered model solution & stacks it
% together across the whole inversion horizon. For example, BaaMat is a
% sum(nxa)*sum(nxa(1:S-1)) matrix containing the individual matrices Baa{s}
% across its diagonal such that:
% BaaMat(xaVecInd{s},xalagVecInd{s}) = Ba{s}a{s-1}.
[BaaMat,BauMat,BanMat,BuaMat,BuuMat,BunMat,BnaMat,BnuMat,BnnMat,...
    AaiLoadMat,AanLoadMat,AuiLoadMat,AunLoadMat,AniLoadMat,AnnLoadMat,...
    UaiLoadMat,UanLoadMat,UuiLoadMat,UunLoadMat,UniLoadMat,UnnLoadMat] ...
    = partition_and_stack_model_solution(BT,FTkTimesPHITa,PHITu,nx,nxa,...
    nxu,nxn,nai,nan,nui,nun,xTaInd,xTuInd,xTnInd,xTalagInd,xTulagInd,...
    xTnlagInd,xTaVecInd,xTuVecInd,xTnVecInd,xTalagVecInd,xTulagVecInd,...
    xTnlagVecInd,aTiInd,aTnInd,aTiVecInd,aTnVecInd,uTiInd,uTnInd,...
    uTiVecInd,uTnVecInd,H);

%% COMPUTE THE LOADING & CONSTANT MATRIX IN THE STACKED SYSTEM
% Call a sub-function to compute the constant & loading in the siumltanous
% system of equations sych that: Wzi = C, where zi = [aiVec;uiVec].
[C,W] = compute_matrices_in_simultaneous_system(nx,nxa,nxu,nxn,nai,nui,...
    xT,xaVec,xalagVec,xuVec,xulagVec,anVec,unVec,xTaVecInd,xTuVecInd,...
    xTnVecInd,xTalagVecInd,xTulagVecInd,xTnlagVecInd,BaaMat,BauMat,...
    BanMat,BuaMat,BuuMat,BunMat,BnaMat,BnuMat,BnnMat,AaiLoadMat,...
    AanLoadMat,AuiLoadMat,AunLoadMat,AniLoadMat,AnnLoadMat,UaiLoadMat,...
    UanLoadMat,UuiLoadMat,UunLoadMat,UniLoadMat,UnnLoadMat,H);

%% IMPLEMENT THE INVERSION
% Invert the system to compute the shocks, zi = [aiVec;uiVec], depending on
% the choice of inversion algorithm in the inputs. If the algorithm chosen
% is minimum variance (i.e. not minimum damage), define the variance
% covariance matrix of the shocks (because of the assumption implemented in
% the construction of MAPS models that all shocks are standard normally
% distributed, this is just the identity).
if Instructions.algorithmIsMinimumDamage
    zi = zbariVec+(W'/(W*W'))*(C-W*zbariVec);
else
    SIGMAzi = eye(sum(nai)+sum(nui),sum(nai)+sum(nui));
    zi = (SIGMAzi\(W'/(W/SIGMAzi*W')))*C;
end

%% CREATE A COMPLETE MATRIX OF SHOCKS ACROSS ALL TIME PERIODS
% Call a sub-function to combine the shocks and reorder them such that a &
% u are nz*S matrices with, for example, a(:,1) containing the original,
% model ordering of the shocks for period s = 1 of the inversion.
[a,u] = combine_and_reorder_shocks(zi,nz,nai,nui,aTiVecInd,uTiVecInd,...
    anVec,unVec,aTnVecInd,uTnVecInd,Ta,Tu,H);

%% PACK SHOCKS OUTPUT
% Pack the shocks output into a structure, ready for use in other MAPS
% functions.
if ShocksInclude.anticipated
    NewShocks.anticipated = a;
end
if ShocksInclude.unanticipated
    NewShocks.unanticipated = u;
end

end

%% FUNCTION TO COMBINE SHOCKS USED & NOT USED AND THEN REORDER THEM
function [a,u] = combine_and_reorder_shocks(zi,nz,nai,nui,aiVecInd,...
    uiVecInd,anVec,unVec,anVecInd,unVecInd,Ta,Tu,S)
% This helper creates a complete model ordered set of shocks.
% It separates out the anticipated and unanticipated shocks computed in the
% inversion. It then combines them with the shocks that were not used to
% achieve the inversion before reordering the vector using the reordering
% operators.
%
% INPUTS:
%   -> zi: combined set of shocks computed in the inversion
%   -> nz: number of shocks in the model
%   -> nai: number of anticipated shocks in the inversion in each period
%   -> nui: number of unanticipated shocks in the inversion in each period
%   -> aiVecInd: model index numbers of stacked anticipated shocks used
%   -> uiVecInd: model index numbers of stacked unanticipated shocks used
%   -> anVec: stacked values for anticipated shocks not used
%   -> unVec: stacked values for unanticipated shocks not used
%   -> anVecInd: model indexes of stacked anticipated shocks not used
%   -> unVecInd: model indexes of stacked unanticipated shocks not used
%   -> Ta: reorder operators for the anticipated shocks
%   -> Tu: reorder operators for the unanticipated shocks
%   -> S: horizon of the inversion
%
% OUTPUTS:
%   -> a: complete matrix of anticipated shocks
%   -> u: complete matrix of unanticipated shocks
%
% CALLS:
%   -> none

%% EXTRACT THE VALUES FOR THE SHOCKS IN THE INVERSION
% Extract the anticipated shocks and the unanticipated shocks computed in
% the inversion using the fact that inversion orders the anticipated shocks
% first.
aiVec = zi(1:sum(nai));
uiVec = zi(sum(nai)+1:sum(nai)+sum(nui));

%% SETUP MATRICES FOR THE SHOCKS
% Setup matrices for the shocks to be filled in below. Setup a matrix of
% inversion ordered and model ordered shocks for each type.
aT = zeros(nz,S);
uT = zeros(nz,S);
a = zeros(nz,S);
u = zeros(nz,S);

%% COMPUTE THE SHOCKS
% For each period of the inversion from 1 to S compute the inversion
% ordered shock vector by stacking the inverted shocks on top of the
% shocks being held fixed. Use the reorder operator for that period to
% compute the model ordered shocks.
for s = 1:S
    aT(:,s) = [aiVec(aiVecInd{s});anVec(anVecInd{s})];
    uT(:,s) = [uiVec(uiVecInd{s});unVec(unVecInd{s})];
    a(:,s) = Ta{s}'*aT(:,s);
    u(:,s) = Tu{s}'*uT(:,s);
end

end

%% FUNCTION TO COMPUTE THE CONSTANT & LOADING IN THE SIMULTANEOUS SYSTEM
function [C,W] = compute_matrices_in_simultaneous_system(nx,nxa,nxu,nxn,...
    nai,nui,xT,xaVec,xalagVec,xuVec,xulagVec,anVec,unVec,xTaVecInd,...
    xTuVecInd,xTnVecInd,xTalagVecInd,xTulagVecInd,xTnlagVecInd,BaaMat,...
    BauMat,BanMat,BuaMat,BuuMat,BunMat,BnaMat,BnuMat,BnnMat,AaiLoadMat,...
    AanLoadMat,AuiLoadMat,AunLoadMat,AniLoadMat,AnnLoadMat,UaiLoadMat,...
    UanLoadMat,UuiLoadMat,UunLoadMat,UniLoadMat,UnnLoadMat,S)
% This helper creates the constant and loadings in the simultaneous system.
% It uses the reordered, partitioned model solution to run through the
% matrix algebra necessary to compute the constants and loadings required
% for the inversion.
%
% INPUTS:
%   -> nx: number of model variables in the model
%   -> nxa: number of anticipated model variable fixes in each period
%   -> nxu: number of unanticipated model variable fixes in each period
%   -> nxn: number of model variable fixes not being fixed in each period
%   -> nai: number of anticipated shocks in the inversion in each period
%   -> nui: number of unanticipated shocks in the inversion in each period
%   -> xT: vector of initial conditions for the model variables
%   -> xaVec: stacked values for the anticipated fixes
%   -> xalagVec: stacked lagged values for the anticipated fixes
%   -> xuVec: stacked values for the unanticipated fixes
%   -> anVec: stacked values for anticipated shocks not used
%   -> unVec: stacked values for unanticipated shocks not used
%   -> xTaVecInd: reordered model indexes of stacked anticipated fixes
%   -> xTuVecInd: reordered model indexes of stacked unanticipated fixes
%   -> xTnVecInd: reordered model indexes of model variables not fixed
%   -> xTalagVecInd: reordered inds of stacked lagged anticipated fixes
%   -> xTulagVecInd: reordered inds of stacked lagged unanticipated fixes
%   -> xTnlagVecInd: reordered inds of stacked lagged variables not fixed
%   -> BaaMat: anticipated fixes loadings on lagged anticipated fixes
%   -> BauMat: anticipated fixes loadings on lagged unanticipated fixes
%   -> BanMat: anticipated fixes loadings on lagged variables not fixed
%   -> BuaMat: unanticipated fixes loadings on lagged anticipated fixes
%   -> BuuMat: unanticipated fixes loadings on lagged unanticipated fixes
%   -> BunMat: unanticipated fixes loadings on lagged variables not fixed
%   -> BnaMat: variables not fixed loadings on lagged anticipated fixes
%   -> BnuMat: variables not fixed loadings on lagged unanticipated fixes
%   -> BnnMat: variables not fixed loadings on lagged variables not fixed
%   -> AaiLoadMat: anticipated fixes loads on anticipated shocks used
%   -> AanLoadMat: anticipated fixes loads on anticipated shocks not used
%   -> AuiLoadMat: unanticipated fixes loads on anticipated shocks used
%   -> AunLoadMat: unanticipated fixes loads on anticipated shocks not used
%   -> AniLoadMat: variables not fixed loads on anticipated shocks used
%   -> AnnLoadMat: variables not fixed loads on anticipated shocks not used
%   -> UaiLoadMat: anticipated fixes loads on unanticipated shocks used
%   -> UanLoadMat: anticipated fixes loads on unanticipated shocks not used
%   -> UuiLoadMat: unanticipated fixes loads on unanticipated shocks used
%   -> UunLoadMat: unanticipated fixes loads on unanticipated shocks not ..
%   -> UniLoadMat: variables not fixed loads on unanticipated shocks used
%   -> UnnLoadMat: variables not fixed loads on unanticipated shocks not ..
%   -> S: horizon of the inversion
%
% OUTPUTS:
%   -> C: constant in the simultaneous system
%   -> W: loadings in the simultaneous system
%
% CALLS:
%   -> none

%% SETUP COME COMPONENTS OF THE SIMULTANEOUS SYSTEM
% Setup some of the componenents of teh simultaneous system and some
% intermediate variables required for the algebra.
W1aa = AaiLoadMat;
W1ua = AuiLoadMat;
W1au = UaiLoadMat;
W1uu = UuiLoadMat;
W2aa = zeros(sum(nxa),sum(nai));
W2ua = zeros(sum(nxu),sum(nai));
W2au = zeros(sum(nxa),sum(nui));
W2uu = zeros(sum(nxu),sum(nui));
Aan = AanLoadMat*anVec;
Aun = AunLoadMat*anVec;
Ann = AnnLoadMat*anVec;
Uan = UanLoadMat*unVec;
Uun = UunLoadMat*unVec;
Unn = UnnLoadMat*unVec;
xnpartsumlag = zeros(nx+sum(nxn(1:S-1)),1);
AnilagLoadMat = zeros(nx+sum(nxn(1:S-1)),sum(nai));
AnilagLoadMat(nx+1:nx+sum(nxn(1:S-1)),:) = AniLoadMat(1:sum(nxn(1:S-1)),:);
AsumnilagLoad = zeros(nx,sum(nai));
UnilagLoadMat = zeros(nx+sum(nxn(1:S-1)),sum(nui));
UnilagLoadMat(nx+1:nx+sum(nxn(1:S-1)),:) = UniLoadMat(1:sum(nxn(1:S-1)),:);
UsumnilagLoad = zeros(nx,sum(nui));

%% COMPUTE THE COMPONENTS OF THE SYSTEM
% Compute all the componenents of the system starting from the first period
% of the inversion 1 and ending in the last period, S.
for s = 1:S
    if s == 1
        xnpartsumlag(xTnlagVecInd{s},1) = xT;
    else
        xnpartsumlag(xTnlagVecInd{s},1) = ...
            BnnMat(xTnVecInd{s-1},xTnlagVecInd{s-1})...
            *xnpartsumlag(xTnlagVecInd{s-1},1)...
            +BnaMat(xTnVecInd{s-1},xTalagVecInd{s-1})...
            *xalagVec(xTalagVecInd{s-1})...
            +BnuMat(xTnVecInd{s-1},xTulagVecInd{s-1})...
            *xulagVec(xTulagVecInd{s-1})...
            +Ann(xTnVecInd{s-1},1)+Unn(xTnVecInd{s-1},1);
        AsumnilagLoad = BnnMat(xTnVecInd{s-1},xTnlagVecInd{s-1})...
            *AsumnilagLoad+AnilagLoadMat(xTnlagVecInd{s},:);
        UsumnilagLoad = BnnMat(xTnVecInd{s-1},xTnlagVecInd{s-1})...
            *UsumnilagLoad+UnilagLoadMat(xTnlagVecInd{s},:);
    end
    W2aa(xTaVecInd{s},:) = ...
        BanMat(xTaVecInd{s},xTnlagVecInd{s})*AsumnilagLoad;
    W2ua(xTuVecInd{s},:) = ...
        BunMat(xTuVecInd{s},xTnlagVecInd{s})*AsumnilagLoad;
    W2au(xTaVecInd{s},:) = ...
        BanMat(xTaVecInd{s},xTnlagVecInd{s})*UsumnilagLoad;
    W2uu(xTuVecInd{s},:) = ...
        BunMat(xTuVecInd{s},xTnlagVecInd{s})*UsumnilagLoad;
end

%% COMPUTE THE MATRICES IN THE SIMULTANEOUS SYSTEM
% Compute the anticipated and unanticipated constants and the loadings on
% the shocks to invert.
Ca = xaVec-BaaMat*xalagVec-BauMat*xulagVec-Aan-Uan-BanMat*xnpartsumlag;
Cu = xuVec-BuaMat*xalagVec-BuuMat*xulagVec-Aun-Uun-BunMat*xnpartsumlag;
C = [Ca;Cu];
W = [W1aa+W2aa W1au+W2au; W1ua+W2ua W1uu+W2uu];

end

%% FUNCTION TO PARTITION AND VECTORIZE MODEL SOLUTION MATRICES
function [BaaMat,BauMat,BanMat,BuaMat,BuuMat,BunMat,BnaMat,BnuMat,...
    BnnMat,AaiLoadMat,AanLoadMat,AuiLoadMat,AunLoadMat,AniLoadMat,...
    AnnLoadMat,UaiLoadMat,UanLoadMat,UuiLoadMat,UunLoadMat,UniLoadMat,...
    UnnLoadMat] = partition_and_stack_model_solution(BT,FTkTimesPHITa,...
    PHITu,nx,nxa,nxu,nxn,nai,nan,nui,nun,xTaInd,xTuInd,xTnInd,xTalagInd,...
    xTulagInd,xTnlagInd,xTaVecInd,xTuVecInd,xTnVecInd,xTalagVecInd,...
    xTulagVecInd,xTnlagVecInd,aTiInd,aTnInd,aTiVecInd,aTnVecInd,uTiInd,...
    uTnInd,uTiVecInd,uTnVecInd,S)
% This helper partitions and stacks the components of the model solution.
% It uses the model index information to partition and stack the reordered
% model solution across all the inversion forecast quarters, 1 to S.
%
% INPUTS:
%   -> BT: reordered solution backward loadings for each period
%   -> FTkTimesPHITa: reordered anticipated shock loadings for each period
%   -> PHITu: reordered unanticipated shock loadings for each period
%   -> nx: number of model variables in the model
%   -> nxa: number of anticipated model variable fixes in each period
%   -> nxu: number of unanticipated model variable fixes in each period
%   -> nxn: number of model variable fixes not being fixed in each period
%   -> nai: number of anticipated shocks in the inversion in each period
%   -> nan: number of anticipated shocks not used in each period
%   -> nui: number of unanticipated shocks in the inversion in each period
%   -> nun: number of unanticipated shocks not used in each period
%   -> xTaInd: reordered model indexes of anticipated fixes
%   -> xTuInd: reordered model indexes of unanticipated fixes
%   -> xTnInd: reordered model indexes of model variables not fixed
%   -> xTalagInd: reordered model indexes of lagged anticipated fixes
%   -> xTulagInd: reordered model indexes of lagged unanticipated fixes
%   -> xTnlagInd: reordered model indexes of lagged model variables not ...
%   -> xTaVecInd: reordered model indexes of stacked anticipated fixes
%   -> xTuVecInd: reordered model indexes of stacked unanticipated fixes
%   -> xTnVecInd: reordered model indexes of model variables not fixed
%   -> xTalagVecInd: reordered inds of stacked lagged anticipated fixes
%   -> xTulagVecInd: reordered inds of stacked lagged unanticipated fixes
%   -> xTnlagVecInd: reordered inds of stacked lagged variables not fixed
%   -> aTiInd: reordered model indexes of inversion anticipated shocks
%   -> aTnInd: reordered model indexes of non-inversion anticipated shocks
%   -> aTiVecInd: reordered model inds of stacked inversion anticipated ...
%   -> aTnVecInd: reordered model inds of stacked non-inversion anticipated
%   -> uTiInd: reordered model inds of inversion unanticipated shocks
%   -> uTnInd: reordered model inds of non-inversion unanticipated shocks
%   -> uTiVecInd: reordered model inds of stacked inversion unanticipated..
%   -> uTnVecInd: reordered model inds of stacked inversion non-inversion..
%   -> S: horizon of the inversion
%
% OUTPUTS:
%   -> BaaMat: anticipated fixes loadings on lagged anticipated fixes
%   -> BauMat: anticipated fixes loadings on lagged unanticipated fixes
%   -> BanMat: anticipated fixes loadings on lagged variables not fixed
%   -> BuaMat: unanticipated fixes loadings on lagged anticipated fixes
%   -> BuuMat: unanticipated fixes loadings on lagged unanticipated fixes
%   -> BunMat: unanticipated fixes loadings on lagged variables not fixed
%   -> BnaMat: variables not fixed loadings on lagged anticipated fixes
%   -> BnuMat: variables not fixed loadings on lagged unanticipated fixes
%   -> BnnMat: variables not fixed loadings on lagged variables not fixed
%   -> AaiLoadMat: anticipated fixes loads on anticipated shocks used
%   -> AanLoadMat: anticipated fixes loads on anticipated shocks not used
%   -> AuiLoadMat: unanticipated fixes loads on anticipated shocks used
%   -> AunLoadMat: unanticipated fixes loads on anticipated shocks not used
%   -> AniLoadMat: variables not fixed loads on anticipated shocks used
%   -> AnnLoadMat: variables not fixed loads on anticipated shocks not used
%   -> UaiLoadMat: anticipated fixes loads on unanticipated shocks used
%   -> UanLoadMat: anticipated fixes loads on unanticipated shocks not used
%   -> UuiLoadMat: unanticipated fixes loads on unanticipated shocks used
%   -> UunLoadMat: unanticipated fixes loads on unanticipated shocks not ..
%   -> UniLoadMat: variables not fixed loads on unanticipated shocks used
%   -> UnnLoadMat: variables not fixed loads on unanticipated shocks not ..
%
% CALLS:
%   -> none

%% SETUP THE OUTPUTS
% Setup the outputs using the information about the numbers of each fix and
% shock type passed in.
BaaMat = zeros(sum(nxa),sum(nxa(1:S-1)));
BauMat = zeros(sum(nxa),sum(nxu(1:S-1)));
BanMat = zeros(sum(nxa),nx+sum(nxn(1:S-1)));
BuaMat = zeros(sum(nxu),sum(nxa(1:S-1)));
BuuMat = zeros(sum(nxu),sum(nxu(1:S-1)));
BunMat = zeros(sum(nxu),nx+sum(nxn(1:S-1)));
BnaMat = zeros(sum(nxn),sum(nxa(1:S-1)));
BnuMat = zeros(sum(nxn),sum(nxu(1:S-1)));
BnnMat = zeros(sum(nxn),nx+sum(nxn(1:S-1)));
AaiLoadMat = zeros(sum(nxa),sum(nai));
AanLoadMat = zeros(sum(nxa),sum(nan));
AuiLoadMat = zeros(sum(nxu),sum(nai));
AunLoadMat = zeros(sum(nxu),sum(nan));
AniLoadMat = zeros(sum(nxn),sum(nai));
AnnLoadMat = zeros(sum(nxn),sum(nan));
UaiLoadMat = zeros(sum(nxa),sum(nui));
UanLoadMat = zeros(sum(nxa),sum(nun));
UuiLoadMat = zeros(sum(nxu),sum(nui));
UunLoadMat = zeros(sum(nxu),sum(nun));
UniLoadMat = zeros(sum(nxn),sum(nui));
UnnLoadMat = zeros(sum(nxn),sum(nun));

%% EXTRACT THE PARTITIONED MODEL SOLUTION INFO
% For each period in the inversion used the stacked and non-stacked
% reordered model index information to create the stacked, partitioned
% components of the model solution.
for s = 1:S
    BaaMat(xTaVecInd{s},xTalagVecInd{s}) = BT{s}(xTaInd{s},xTalagInd{s});
    BauMat(xTaVecInd{s},xTulagVecInd{s}) = BT{s}(xTaInd{s},xTulagInd{s});
    BanMat(xTaVecInd{s},xTnlagVecInd{s}) = BT{s}(xTaInd{s},xTnlagInd{s});
    BuaMat(xTuVecInd{s},xTalagVecInd{s}) = BT{s}(xTuInd{s},xTalagInd{s});
    BuuMat(xTuVecInd{s},xTulagVecInd{s}) = BT{s}(xTuInd{s},xTulagInd{s});
    BunMat(xTuVecInd{s},xTnlagVecInd{s}) = BT{s}(xTuInd{s},xTnlagInd{s});
    BnaMat(xTnVecInd{s},xTalagVecInd{s}) = BT{s}(xTnInd{s},xTalagInd{s});
    BnuMat(xTnVecInd{s},xTulagVecInd{s}) = BT{s}(xTnInd{s},xTulagInd{s});
    BnnMat(xTnVecInd{s},xTnlagVecInd{s}) = BT{s}(xTnInd{s},xTnlagInd{s});
    UaiLoadMat(xTaVecInd{s},uTiVecInd{s}) = PHITu{s}(xTaInd{s},uTiInd{s});
    UanLoadMat(xTaVecInd{s},uTnVecInd{s}) = PHITu{s}(xTaInd{s},uTnInd{s});
    UuiLoadMat(xTuVecInd{s},uTiVecInd{s}) = PHITu{s}(xTuInd{s},uTiInd{s});
    UunLoadMat(xTuVecInd{s},uTnVecInd{s}) = PHITu{s}(xTuInd{s},uTnInd{s});
    UniLoadMat(xTnVecInd{s},uTiVecInd{s}) = PHITu{s}(xTnInd{s},uTiInd{s});
    UnnLoadMat(xTnVecInd{s},uTnVecInd{s}) = PHITu{s}(xTnInd{s},uTnInd{s});
    for k = 0:S-s
        AaiLoadMat(xTaVecInd{s},aTiVecInd{s+k}) = ...
            FTkTimesPHITa{s,s+k}(xTaInd{s},aTiInd{s+k});
        AanLoadMat(xTaVecInd{s},aTnVecInd{s+k}) = ...
            FTkTimesPHITa{s,s+k}(xTaInd{s},aTnInd{s+k});
        AuiLoadMat(xTuVecInd{s},aTiVecInd{s+k}) = ...
            FTkTimesPHITa{s,s+k}(xTuInd{s},aTiInd{s+k});
        AunLoadMat(xTuVecInd{s},aTnVecInd{s+k}) = ...
            FTkTimesPHITa{s,s+k}(xTuInd{s},aTnInd{s+k});
        AniLoadMat(xTnVecInd{s},aTiVecInd{s+k}) = ...
            FTkTimesPHITa{s,s+k}(xTnInd{s},aTiInd{s+k});
        AnnLoadMat(xTnVecInd{s},aTnVecInd{s+k}) = ...
            FTkTimesPHITa{s,s+k}(xTnInd{s},aTnInd{s+k});
    end
end

end

%% FUNCTION TO COMPUTE INDEXES OF REORDERED VARIABLES & SHOCKS
function [xTaInd,xTuInd,xTnInd,xTalagInd,xTulagInd,xTnlagInd,...
    xTaVecInd,xTuVecInd,xTnVecInd,xTalagVecInd,xTulagVecInd,...
    xTnlagVecInd,aTiInd,aTnInd,aTiVecInd,aTnVecInd,uTiInd,uTnInd,...
    uTiVecInd,uTnVecInd] = ...
    compute_indexes_of_reordered_variables_and_shocks(...
    nx,nxa,nxu,nxn,nai,nan,nui,nun,S)
% This helper computes the indices of the reordered model solution.
% It uses information about the number of each type of fix and the
% convention that model variables are reordered with anticipated fixes
% first, then unanticipated fixes, then variables not fixed and that both
% anticipated & unanticipated shocks are separately reordered so that the
% shocks being used to achieve the inversion are ordered first. This
% function also computes the model variable indices in the lagged system
% and both the contemporaenous model variable / shock indices and lagged
% model variable indices in the time-stacked system.
%
% INPUTS:
%   -> nx: number of model variables in the model
%   -> nxa: number of anticipated model variable fixes in each period
%   -> nxu: number of unanticipated model variable fixes in each period
%   -> nxn: number of model variable fixes not being fixed in each period
%   -> nai: number of anticipated shocks in the inversion in each period
%   -> nan: number of anticipated shocks not used in each period
%   -> nui: number of unanticipated shocks in the inversion in each period
%   -> nun: number of unanticipated shocks not used in each period
%   -> S: horizon of the inversion
%
% OUTPUTS:
%   -> xTaInd: reordered model indexes of anticipated fixes
%   -> xTuInd: reordered model indexes of unanticipated fixes
%   -> xTnInd: reordered model indexes of model variables not fixed
%   -> xTalagInd: reordered model indexes of lagged anticipated fixes
%   -> xTulagInd: reordered model indexes of lagged unanticipated fixes
%   -> xTnlagInd: reordered model indexes of lagged model variables not ...
%   -> xTaVecInd: reordered model indexes of stacked anticipated fixes
%   -> xTuVecInd: reordered model indexes of stacked unanticipated fixes
%   -> xTnVecInd: reordered model indexes of model variables not fixed
%   -> xTalagVecInd: reordered inds of stacked lagged anticipated fixes
%   -> xTulagVecInd: reordered inds of stacked lagged unanticipated fixes
%   -> xTnlagVecInd: reordered inds of stacked lagged variables not fixed
%   -> aTiInd: reordered model indexes of inversion anticipated shocks
%   -> aTnInd: reordered model indexes of non-inversion anticipated shocks
%   -> aTiVecInd: reordered model inds of stacked inversion anticipated ...
%   -> aTnVecInd: reordered model inds of stacked non-inversion anticipated
%   -> uTiInd: reordered model inds of inversion unanticipated shocks
%   -> uTnInd: reordered model inds of non-inversion unanticipated shocks
%   -> uTiVecInd: reordered model inds of stacked inversion unanticipated..
%   -> uTnVecInd: reordered model inds of stacked inversion non-inversion..
%
% CALLS:
%   -> none

%% SETUP THE OUTPUTS
% All outputs are setup as row cell vectors with the same number of columns
% as periods in the inversion.
xTaInd = cell(1,S);
xTuInd = cell(1,S);
xTnInd = cell(1,S);
xTalagInd = cell(1,S);
xTulagInd = cell(1,S);
xTnlagInd = cell(1,S);
xTaVecInd = cell(1,S);
xTuVecInd = cell(1,S);
xTnVecInd = cell(1,S);
xTalagVecInd = cell(1,S);
xTulagVecInd = cell(1,S);
xTnlagVecInd = cell(1,S);
aTiInd = cell(1,S);
aTnInd = cell(1,S);
aTiVecInd = cell(1,S);
aTnVecInd = cell(1,S);
uTiInd = cell(1,S);
uTnInd = cell(1,S);
uTiVecInd = cell(1,S);
uTnVecInd = cell(1,S);

%% COMPUTE THE INDICES
% Run through each time period of the inversion and compute the model
% indices in the partitioned system using the convention about how the
% model is reordered.
for s = 1:S
    xTaInd{s} = (1:nxa(s))';
    xTuInd{s} = nxa(s)+(1:nxu(s))';
    xTnInd{s} = nxa(s)+nxu(s)+(1:nxn(s))';
    aTiInd{s} = (1:nai(s))';
    aTnInd{s} = nai(s)+(1:nan(s))';
    uTiInd{s} = (1:nui(s))';
    uTnInd{s} = nui(s)+(1:nun(s))';
    if s == 1
        xTalagInd{s} = (1:0)';
        xTulagInd{s} = (1:0)';
        xTnlagInd{s} = (1:nx)';
        xTaVecInd{s} = (1:nxa(s))';
        xTuVecInd{s} = (1:nxu(s))';
        xTnVecInd{s} = (1:nxn(s))';
        xTalagVecInd{s} = (1:0)';
        xTulagVecInd{s} = (1:0)';
        xTnlagVecInd{s} = (1:nx)';
        aTiVecInd{s} = (1:nai(s))';
        aTnVecInd{s} = (1:nan(s))';
        uTiVecInd{s} = (1:nui(s))';
        uTnVecInd{s} = (1:nun(s))';
    else
        xTalagInd{s} = (1:nxa(s-1))';
        xTulagInd{s} = nxa(s-1)+(1:nxu(s-1))';
        xTnlagInd{s} = nxa(s-1)+nxu(s-1)+(1:nxn(s-1))';
        xTaVecInd{s} = sum(nxa(1:s-1))+(1:nxa(s))';
        xTuVecInd{s} = sum(nxu(1:s-1))+(1:nxu(s))';
        xTnVecInd{s} = sum(nxn(1:s-1))+(1:nxn(s))';
        aTiVecInd{s} = sum(nai(1:s-1))+(1:nai(s))';
        aTnVecInd{s} = sum(nan(1:s-1))+(1:nan(s))';
        uTiVecInd{s} = sum(nui(1:s-1))+(1:nui(s))';
        uTnVecInd{s} = sum(nun(1:s-1))+(1:nun(s))';
        if s == 2
            xTalagVecInd{s} = (1:nxa(s-1))';
            xTulagVecInd{s} = (1:nxu(s-1))';
            xTnlagVecInd{s} = nx+(1:nxn(s-1))';
        else
            xTalagVecInd{s} = sum(nxa(1:s-2))+(1:nxa(s-1))';
            xTulagVecInd{s} = sum(nxu(1:s-2))+(1:nxu(s-1))';
            xTnlagVecInd{s} = nx+sum(nxn(1:s-2))+(1:nxn(s-1))';
        end
    end
end

end

%% FUNCTION TO VECTORIZE FIXES AND SHOCKS NOT USED IN THE INVERSION
function [xaVec,xuVec,xalagVec,xulagVec,anVec,unVec,zbariVec] ...
    = stack_fixes_and_shocks(xa,xu,nxa,nxu,an,un,abari,ubari,S)
% This helper stacks the known shocks and variables across all periods.
% It stacks all the variables and shocks with known values into single
% vectors across all inversion periods, 1 to S. The known values are:
% anticipated model variable fixes, unanticipated model variable fixes,
% anticipated shocks not being used in the inversion, unanticipated shocks
% not being used in the inversion. It also creates lagged versions of the
% model variable fixes.
%
% INPUTS:
%   -> xa: cell array of anticipated fixes
%   -> xu: cell array of unanticipated fixes
%   -> nxa: number of anticipated model variable fixes in each period
%   -> nxu: number of unanticipated model variable fixes in each period
%   -> an: cell array of anticipated shocks not used in the fixes
%   -> un: cell array of unanticipated shocks not used in the fixes
%   -> abari: cell array of anticipated shocks being used in the fixes
%   -> ubari: cell array of unanticipated shocks being used in the fixes
%   -> S: horizon of the inversion
%
% OUTPUTS:
%   -> xaVec: time-stacked vector of anticipated fixes
%   -> xuVec: time-stacked vector of unanticipated fixes
%   -> xalagVec: time-stacked vector of lagged anticipated fixes
%   -> xulagVec: time-stacked vector of lagged unanticipated fixes
%   -> anVec: time-stacked vector of anticipated shocks held fixed
%   -> unVec: time-stacked vector of unanticipated shocks held fixed
%   -> zbariVec: time-stacked vector of shocks being used in the fix
%
% CALLS:
%   -> none

%% CREATE THE CONTEMPORANEOUSLY STACKED OUTPUT
% Use the MATLAB cell2mat command to convert the cell arrays of fixes and
% shocks being held fixed (ie those not being used to achieve the
% inversion) to numeric vectors.
xaVec = cell2mat(xa(:));
xuVec = cell2mat(xu(:));
anVec = cell2mat(an(:));
unVec = cell2mat(un(:));
abariVec = cell2mat(abari(:));
ubariVec = cell2mat(ubari(:));
zbariVec = [abariVec;ubariVec];

%% CREATE THE LAGGED STACKED OUTPUT
% Use the contemporaneously stacked output and the information about the
% number of each type of fixe across all periods of the inversion to create
% the lagged, stacked vector. Note that there is an assumption that all
% lagged model variables in the initial conditions are treated as non-fixed
% (so the lagged fixes in period 1 of the inversion are an empty 0*1
% vector).
xalagVec = xaVec(1:sum(nxa(1:S-1)));
xulagVec = xuVec(1:sum(nxu(1:S-1)));
xalagVec = xalagVec(:);
xulagVec = xulagVec(:);

end

%% CHECK FEASIBILITY OF THE REQUIRED INVERSIONS
function check_feasibility_of_inversion(...
    FTkTimesPHITa,PHITu,nxa,nxu,nai,nui,S)
% This helper checks whether or not both types of inversion are feasible.
% It checks whether or not both the anticipated and unanticipated
% inversions specified in the inputs are feasible in the model input. It
% also checks whether or not they are feasible given the model input and
% the shocks chosen to achieve each type of fix in the inversion. It will
% throw an exception if either of the inversion types input are infeasible.
%
% INPUTS:
%   -> FTkTimesPHITa: reordered anticipated shock loadings for each period
%   -> PHITu: reordered unanticipated shock loadings for each period
%   -> nxa: number of anticipated model variable fixes in each period
%   -> nxu: number of unanticipated model variable fixes in each period
%   -> nai: number of anticipated shocks in the inversion in each period
%   -> nui: number of unanticipated shocks in the inversion in each period
%   -> S: horizon of the inversion
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% SETUP MASTER EXCEPTION
% Setup a master exception for infeasibility to add any causes encountered
% below.
masterErrId = ['MAPS:',mfilename,':Infeasibility'];
InfeasibilityE = generate_MAPS_exception(masterErrId);

%% NUMERICALLY ADJUST THE LOADING MATRICES
% The solution to most LSS models is a numerical approximation. During the
% course of that approximation numerical values that should be zero
% analytically can instead be set to very small numbers (like 1e-15). The
% rank checks below rely on numbers that should be 0 being set precisely to
% 0, so adjust the loading matrices so that small numbers are set to 0.
for s = 1:S
    FTkTimesPHITa{s,s}(abs(FTkTimesPHITa{s,s})<1e-12) = 0;
    PHITu{s}(abs(PHITu{s})<1e-12) = 0;
end

%% CHECK FOR LINEAR INDEPENDENCE OF THE ANTICIPATED FIXES
% For each time period of the inversion check that the anticipated fixes
% are linearly independent in the model input.
badIndicesLogicals = false(1,S);
for s = 1:S
    if nxa(s)>0 && rank(FTkTimesPHITa{s,s}(1:nxa(s),:))<nxa(s)
        badIndicesLogicals(s) = true;
    end
end
badIndices = find(badIndicesLogicals);
if ~isempty(badIndices)
    errId = [masterErrId,':FixesNotLinearlyIndependent'];
    errArgs = {'anticipated' num2str(badIndices)};
    InfeasibilityE = generate_MAPS_exception_and_add_as_cause(...
        InfeasibilityE,errId,errArgs);
end

%% CHECK THAT THE ANTICIPATED FIX-SHOCK COMBINATION IS VALID
% For each time period of the inversion check that the anticipated shocks
% chosen in that period are valid instruments in the fix. Note that if a
% combination of fixes and shock usgages failed the above test, it would
% also fail this test so the time periods of failure from above are removed
% from the time periods of failures here to avoid repetition of (and 
% possible confusion due to) error messages.
badIndicesLogicals = false(1,S);
for s = 1:S
    if nxa(s)>0 && rank(FTkTimesPHITa{s,s}(1:nxa(s),1:nai(s)))<nxa(s)
        badIndicesLogicals(s) = true;
    end
end
badIndicesLogicals(badIndices) = false;
badIndices = find(badIndicesLogicals);
if ~isempty(badIndices)
    errId = [masterErrId,':InvalidShocksInstruments'];
    errArgs = {'anticipated' num2str(badIndices)};
    InfeasibilityE = generate_MAPS_exception_and_add_as_cause(...
        InfeasibilityE,errId,errArgs);
end

%% CHECK FOR LINEAR INDEPENDENCE OF THE UNANTICIPATED FIXES
% For each time period of the inversion check that the unanticipated fixes
% are linearly independent in the model input.
badIndicesLogicals = false(1,S);
for s = 1:S
    if nxu(s)>0 && rank(PHITu{s}(nxa(s)+1:nxa(s)+nxu(s),:))<nxu(s)
        badIndicesLogicals(s) = true;
    end
end
badIndices = find(badIndicesLogicals);
if ~isempty(badIndices)
    errId = [masterErrId,':FixesNotLinearlyIndependent'];
    errArgs = {'unanticipated' num2str(badIndices)};
    InfeasibilityE = generate_MAPS_exception_and_add_as_cause(...
        InfeasibilityE,errId,errArgs);
end

%% CHECK THAT THE UNANTICIPATED FIX-SHOCK COMBINATION IS VALID
% For each time period of the inversion check that the unanticipated shocks
% chosen in that period are valid instruments in the fix. Note that if a
% combination of fixes and shock usgages failed the above test, it would
% also fail this test so the time periods of failure from above are removed
% from the time periods of failures here to avoid repetition of (and 
% possible confusion due to) error messages.
badIndicesLogicals = false(1,S);
for s = 1:S
    if nxu(s)>0 && rank(PHITu{s}(nxa(s)+1:nxa(s)+nxu(s),1:nui(s)))<nxu(s)
        badIndicesLogicals(s) = true;
    end
end
badIndicesLogicals(badIndices) = false;
badIndices = find(badIndicesLogicals);
if ~isempty(badIndices)
    errId = [masterErrId,':InvalidShocksInstruments'];
    errArgs = {'unanticipated' num2str(badIndices)};
    InfeasibilityE = generate_MAPS_exception_and_add_as_cause(...
        InfeasibilityE,errId,errArgs);
end

%% THROW MASTER EXCEPTION
% If any of the above feasibility tests failed, then throw the master
% exception with the failures
if ~isempty(InfeasibilityE.cause)
    throw(InfeasibilityE)
end

end

%% FUNCTION TO COMPUTE REORDERED SYSTEM
function [BT,FTkTimesPHITa,PHITu] = reorder_model_solution(...
    B,F,PHI,Tx,Ta,Tu,S,nx)
% This helper reorders the solution to be consistent with the convention.
% It uses the reorder operators to reorder the model solution to be
% consistent with the inversion convention that model variables are
% reordered with anticipated fixes first, then unanticipated fixes, then
% variables not fixed and that both anticipated & unanticipated shocks are
% separately reordered so that the shocks being used to achieve the
% inversion are ordered first.
%
% INPUTS:
%   -> B: backward loadings in the model solution
%   -> F: forward loadings for anticipated shocks in the model solution
%   -> PHI: loadings on the shocks in the model solution
%   -> Tx: reorder operators for model variables in each period
%   -> Ta: reorder operators for anticipated shocks in each period
%   -> Tu: reorder operators for unanticipated shocks in each period
%   -> S: horizon of the inversion
%   -> nx: number of model variables in the model
%
% OUTPUTS:
%   -> BT: reordered solution backward loadings for each period
%   -> FTkTimesPHITa: reordered anticipated shock loadings for each period
%   -> PHITu: reordered unanticipated shock loadings for each period
%
% CALLS:
%   -> none

%% SETUP THE OUTPUT
% Setup the output as cell array row vectors with S elements.
BT = cell(S,1);
FTkTimesPHITa = cell(S,S);
PHITu = cell(S,1);

%% REORDER THE SOLUTION
% Reorder the backward, forward and shocks loadings in the model solution
% to be consistent with the reordered variables and shocks as given by the
% reorder operators. Note the assumption that the model variables are
% ordered in the same way as in the model in the initial conditions. Note
% also that the reordered loadings on the anticipated shocks form an upper
% triangular cell array.
for s = 1:S
    if s == 1
        BT{s} = Tx{s}*B*eye(nx,nx);
    else
        BT{s} = Tx{s}*B*(Tx{s-1}');
    end
    FT = Tx{s}*F*(Tx{s}');
    PHITu{s} = Tx{s}*PHI*(Tu{s}');
    for k = 0:S-s
        PHITa = Tx{s}*PHI*(Ta{s+k}');
        FTkTimesPHITa{s,s+k} = ((FT^k)*PHITa);
    end
end

end

%% FUNCTION TO COMPUTE REORDER OPERATORS
function [Tx,Ta,Tu] = compute_reorder_operators(nx,nxa,nxu,xInd,xaInd,...
    xuInd,nz,zInd,nai,aiInd,nan,anInd,nui,uiInd,nun,unInd,S)
% This helper computes the reorder operators to reorder the model.
% The reorder operators are consistent with the inversion convention that 
% model variables are reordered with anticipated fixes first, then 
% unanticipated fixes, then variables not fixed and that both anticipated & 
% unanticipated shocks are separately reordered so that the shocks being 
% used to achieve the inversion are ordered first, followed by the shocks
% not being used.
%
% INPUTS:
%   -> nx: number of model variables in the model
%   -> nxa: number of anticipated model variable fixes in each period
%   -> nxu: number of unanticipated model variable fixes in each period
%   -> xInd: all model variable indices (1 to nx)
%   -> xaInd: cell array of model indices of anticipated fixes
%   -> xuInd: cell array of model indices of unanticipated fixes
%   -> nz: number of shocks in the model
%   -> zInd: all model shock indices (1 to nz)
%   -> nai: number of anticipated shocks to be used in the inversion
%   -> aiInd: cell array of indices of anticipated shocks used
%   -> nan: number of anticipated shocks not used in the inversion
%   -> anInd: indices of anticipated shocks not used in the inversion
%   -> nui: number of unanticipated shocks to be used in the inversion
%   -> uiInd: cell array of indices of unanticipated shocks used
%   -> nun: number of unanticipated shocks not used in the inversion
%   -> unInd: indices of unanticipated shocks not used in the inversion
%   -> S: horizon of the inversion
%
% OUTPUTS:
%   -> Tx: reorder operators for model variables in each period
%   -> Ta: reorder operators for anticipated shocks in each period
%   -> Tu: reorder operators for unanticipated shocks in each period
%
% CALLS:
%   -> none

%% SETUP OUTPUT
% Setup the reorder operator outputs as row cell arrays across each of the
% periods in the inversion.
Tx = cell(1,S);
Ta = cell(1,S);
Tu = cell(1,S);

%% SETUP SOME MODEL ORDERED INDICES
% Setup some model ordered indices to be used below.
ixVec = xInd;
iaVec = zInd;
iuVec = zInd;

%% LOOP THROW INVERSION PERIODS
% In each iteration setup appropriately sized reordering oeprators. Setup a
% vector that will be used to create the reorder operators. For model
% variables, the convention used is that anticipated fixes are ordered
% first, then unanticipated fixes, then variables not fixed (ie those
% allowed to respond to the inversion). For each type of shock, the
% convention is that the shocks being used to achieve the inversion are
% ordered before those held fixed (not used) in the inversion. The vectors
% are set such that the indices of the vector point to the original model
% ordering of variables and shocks, while the values dictate the new
% ordering.
for s = 1:S
    Tx{s} = zeros(nx,nx);
    Ta{s} = zeros(nz,nz);
    Tu{s} = zeros(nz,nz);
    jxVec = zeros(nx,1);
    jaVec = zeros(nz,1);
    juVec = zeros(nz,1);
    for ixa = 1:nxa(s)
        jxVec(xaInd{s}(ixa)) = ixa;
    end
    for ixu = 1:nxu(s)
        jxVec(xuInd{s}(ixu)) = nxa(s)+ixu;
    end
    for iai = 1:nai(s)
        jaVec(aiInd{s}(iai)) = iai;
    end
    for ian = 1:nan(s)
        jaVec(anInd{s}(ian)) = nai(s)+ian;
    end
    for iui = 1:nui(s)
        juVec(uiInd{s}(iui)) = iui;
    end
    for iun = 1:nun(s)
        juVec(unInd{s}(iun)) = nui(s)+iun;
    end
    jxVec(jxVec==0) = (nxa(s)+nxu(s)+1:nx)';
    Tx{s}((ixVec-1)*nx+jxVec) = 1;
    Ta{s}((iaVec-1)*nz+jaVec) = 1;
    Tu{s}((iuVec-1)*nz+juVec) = 1;
end

end

%% FUNCTION TO CHECK THAT THE INVERSION IS VALID
function check_validity_of_inversion(...
    nxa,xaInd,nxu,xuInd,nai,aiInd,nui,uiInd,S)
% This helper checks that the inversion inputs meet the inversion rules.
% It checks that the number of variables to be fixed (in either type of
% fix) does not exceed the number of shocks being used. It also checks that
% there are no overlaps in the variables being fixed using anticipated and
% unanticipated shocks and in the shocks being used to achieve those fixes.
%
% INPUTS:
%   -> nxa: number of anticipated model variable fixes in each period
%   -> xaInd: cell array of model indices of anticipated fixes
%   -> nxu: number of unanticipated model variable fixes in each period
%   -> xuInd: cell array of model indices of unanticipated fixes
%   -> nai: number of anticipated shocks to be used in the inversion
%   -> aiInd: cell array of indices of anticipated shocks used
%   -> nui: number of unanticipated shocks to be used in the inversion
%   -> uiInd: cell array of indices of unanticipated shocks used
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% SETUP A MASTER EXCEPTION
% Setup a master exception to add cuases to as they are encountered below.
masterErrId = ['MAPS:',mfilename,':InvalidInversion'];
InversionValidityE = generate_MAPS_exception(masterErrId);

%% CHECK THAT THE NUMBER OF ANTICIPATED FIXES DOES NOT EXCEED SHOCKS USED
% Attempt to find periods where the number of anticipated fixes exceeds the
% number of shocks used to achieve those fixes. Throw an exception if any
% were found.
badAnticipatedFixPeriods = find(nai<nxa);
if ~isempty(badAnticipatedFixPeriods)
    errId = [masterErrId,':ShockFixMismatch'];
    errArgs = {'anticipated' num2str(badAnticipatedFixPeriods)};
    InversionValidityE = generate_MAPS_exception_and_add_as_cause(...
        InversionValidityE,errId,errArgs);
end

%% CHECK THAT THE NUMBER OF UNANTICIPATED FIXES DOES NOT EXCEED SHOCKS USED
% Attempt to find periods where the number of unanticipated fixes exceeds
% the number of shocks used to achieve those fixes. Throw an exception if
% any were found.
badUnanticipatedFixPeriods = find(nui<nxu);
if ~isempty(badUnanticipatedFixPeriods)
    errId = [masterErrId,':ShockFixMismatch'];
    errArgs = {'unanticipated' num2str(badUnanticipatedFixPeriods)};
    InversionValidityE = generate_MAPS_exception_and_add_as_cause(...
        InversionValidityE,errId,errArgs);
end

%% CHECK THAT ANTICIPATED FIXES AND UNANTICIPATED FIXES DO NOT OVERLAP
% Check that the variables being fixed using anticipated and unanticipated
% shocks do not overlap in any of the inversion periods. Throw an exception
% if they do.
for s = 1:S
    overlapInds = intersect(xaInd{s},xuInd{s});
    if ~isempty(overlapInds)
        errId = [masterErrId,':OverlappingInds'];
        errArgs = {'fixes' num2str(s) num2str(overlapInds')};
        InversionValidityE = generate_MAPS_exception_and_add_as_cause(...
            InversionValidityE,errId,errArgs);
    end
end

%% CHECK THAT ANTICIPATED SHOCKS AND UNANTICIPATED SHOCKS DO NOT OVERLAP
% Check that the shocks being used in the anticipated and unanticipated
% fixes do not overlap in any of the inversion periods. Throw an exception
% if they do. Note that relative to the cell above, the calculation here
% involves finding the set difference between the indices of all the shocks
% and the indices of the shocks not being used in the inversion (which is
% what is provided as input to this function).
for s = 1:S
    overlapInds = intersect(aiInd{s},uiInd{s});
    if ~isempty(overlapInds)
        errId = [masterErrId,':OverlappingInds'];
        errArgs = {'shocks' num2str(s) num2str(overlapInds')};
        InversionValidityE = generate_MAPS_exception_and_add_as_cause(...
            InversionValidityE,errId,errArgs);
    end
end

%% THROW ANY EXCEPTIONS ENCOUNTERED
% If the inversion was found to be invalid for any of the above reasons,
% throw the master exception.
if ~isempty(InversionValidityE.cause)
    throw(InversionValidityE)
end

end

%% FUNCTION TO UNPACK THE INVERSION INFO
function [xa,xaInd,xu,xuInd,an,anInd,un,unInd,abari,aiInd,ubari,uiInd] =...
    create_inversion_info(...
    Instructions,InstructionsInclude,Shocks,ShocksInclude,H,zInd)
% This helper creates all the insersion info needed for the algorithm.
% It unpacks the content of the inversion instructions and augments it to
% to include all data necessary for the inversion algorithm to operate
% efficiently. It also ensures that the algebra of the general case (which
% includes anticipated and unanticipated fixes that cover the whole
% forecast horizon) works by creating "empty" fix data for forecast
% quarters where an inversion was not requested, for forecast quarters
% where an inversion of a particular type (either anticipated or
% unanticipated) was not requested and zeros for shocks that were not
% included in the input.
%
% INPUTS:
%   -> Instructions: structure detailing the inversion required
%   -> InstructionsInclude: structure describing the content of the
%      inversion
%   -> Shocks: structure of existing shocks
%   -> ShocksInclude: structure describing the content of the shocks
%   -> H: numeric scalar representing the forecast horizon
%   -> zInd: vector of model shock indices from 1 to nz
%
% OUTPUTS:
%   -> xa: cell array of anticipated fix values
%   -> xaInd: cell array of anticipated fix indices
%   -> xu: cell array of unanticipated fix values
%   -> xuInd: cell array of unanticipated fix indices
%   -> an: cell array of values of anticipated shocks not used
%   -> anInd: cell array of indices of anticipated shocks not used
%   -> un: cell array of values of unanticipated shocks not used
%   -> unInd: cell array of indices of unanticipated shocks not used
%   -> abari: cell array of values of existing anticipated shocks used
%   -> aiInd: cell array of indices of anticipated shocks used
%   -> ubari: cell array of values of existing unanticipated shocks used
%   -> uiInd: cell array of indices of unanticipated shocks used
%
% CALLS:
%   -> augment_inversion_data (sub-function)
%   -> create_empty_inversion_data (sub-function)
%   -> reflect_shock_usage_indices_and_split_values_by_usage (sub-function)

%% AUGMENT OR CREATE INVERSION INFO DATA
% For each type of fix either augment the information using a helper
% function below if the fix type exists in the instructions input or create
% empty data if it does not. In the first case, the inversion info is
% extended with empty data to cover the whole forecast horizon (if
% necessary). In addition, the indices in the inversion instructions are
% validated against ...
fixTypes = {'Anticipated';'Unanticipated'};
nFixTypes = size(fixTypes,1);
for iType = 1:nFixTypes
    iFixType = [fixTypes{iType},'Fixes'];
    iFixTypeLower = [lower(fixTypes{iType}),'Fixes'];
    iShockTypeLower = lower(fixTypes{iType});
    if InstructionsInclude.(iFixTypeLower)
        Instructions.(iFixType) = augment_inversion_data(...
            Instructions.(iFixType),H);
    else
        Instructions.(iFixType) = create_empty_inversion_data(H);
    end
    if ~ShocksInclude.(iShockTypeLower)
        Shocks.(iShockTypeLower) = zeros(size(zInd,1),H);
    end
    Instructions.(iFixType) = ...
        reflect_shock_usage_indices_and_split_values_by_usage(...
        Instructions.(iFixType),Shocks.(iShockTypeLower),zInd,H);
end

%% UNPACK INVERSION INFO
% Unpack all the components of the inversion instructions augmented above
% and required for output to this function in the inversion.
xa = Instructions.AnticipatedFixes.modelVarValues;
xaInd = Instructions.AnticipatedFixes.modelVarIndices;
xu = Instructions.UnanticipatedFixes.modelVarValues;
xuInd = Instructions.UnanticipatedFixes.modelVarIndices;
an = Instructions.AnticipatedFixes.shockNonUsageValues;
anInd = Instructions.AnticipatedFixes.shockNonUsageIndices;
un = Instructions.UnanticipatedFixes.shockNonUsageValues;
unInd = Instructions.UnanticipatedFixes.shockNonUsageIndices;
abari = Instructions.AnticipatedFixes.shockUsagePreviousValues;
aiInd = Instructions.AnticipatedFixes.shockUsageIndices;
ubari = Instructions.UnanticipatedFixes.shockUsagePreviousValues;
uiInd = Instructions.UnanticipatedFixes.shockUsageIndices;

end

%% HELPER FUNCTION TO EXTEND THE INVERSION DATA TO END OF FORECAST HORIZON
function AugmentedInversionData = augment_inversion_data(InversionData,H)
% This helper augments the inversion data to cover the forecast horizon.
% If the input inversion data does not span the whole forecast horizon,
% this function augments it with empty (dummy) inversion data so that the
% information covers the whole forecast horizon.
%
% INPUTS:
%   -> InversionData: structure of inversion data
%   -> H: forecast horizon
%
% OUTPUTS:
%   -> AugmentedInversionData: structure of inversion data covering the
%      whole forecast horizon
%
% CALLS:
%   -> create_empty_inversion_data_cell (sub-function)

%% CREATE AN EMPTY CELL ARRAY
% Create a cell array of empty data to bridge the gap between the length of
% the inversion data and the forecast horizon. The logic relies on the
% fact that all components of the data structure must have the same length
% (which means that S could be computed from any of the 3 fields in the
% structure) and that the inversion data length cannot exceed the forecast
% horizon (these conditions were tested in the inversion instructions
% validation function).
S = size(InversionData.modelVarValues,2);
emptyDataCellToAugment = create_empty_inversion_data_cell(H-S);

%% AUGMENT THE DATA
% Augment each of the inversion data fields with the empty cell array
% computed above.
dataTypes = {'modelVarValues';'modelVarIndices';'shockUsageIndices'};
nDataTypes = size(dataTypes,1);
for iType = 1:nDataTypes
    iDataType = dataTypes{iType};
    AugmentedInversionData.(iDataType) = [InversionData.(iDataType) ...
        emptyDataCellToAugment];
end

end

%% HELPER FUNCTION TO CREATE EMPTY INVERSION DATA
function InversionData = create_empty_inversion_data(H)
% This helper creates inversion data for fix types not instructed.
% This function allows for both anticipated and unanticipated inversions
% but they are optional. If one of them is not included in the input, this
% function creates empty (dummy or non-) fix data. This can be safely used
% in the general inversion algebra that follows without requiring if
% statements.
%
% INPUTS:
%   -> H: forecast horizon
%
% OUTPUTS:
%   -> InversionData: structure of empty inversion data
%
% CALLS:
%   -> create_empty_inversion_data_cell (sub-function)

%% CREATE AN EMPTY CELL ARRAY
% Create a cell array of empty data covering the whole forecast horizon.
emptyDataCell = create_empty_inversion_data_cell(H);

%% CREATE THE DATA
% Create the expected inversion data fields with the empty data cell
% arrays.
dataTypes = {'modelVarValues';'modelVarIndices';'shockUsageIndices'};
nDataTypes = size(dataTypes,1);
for iType = 1:nDataTypes
    iDataType = dataTypes{iType};
    InversionData.(iDataType) = emptyDataCell;
end

end

%% HELPER FUNCTION TO CREATE EMPTY INVERSION CELL
function emptyDataCell = create_empty_inversion_data_cell(S)
% This helper creates an empty inversion data cell over a horizon S.
% This empty data can be used in the general inversion algebra for periods
% in which a particular fix type does not apply. The algebra follows
% through correctly, thereby avoiding the need for conditional execution
% using if statements.
%
% INPUTS:
%   -> S: horizon
%
% OUTPUTS:
%   -> emptyDataCell: 1*S cell array of empty zeros(0,1) data
%
% CALLS:
%   -> none

%% CREATE THE CELL ARRAY OF EMPTY DATA
% Create the cell array of empty data (zeros(0,1)) over the horizon S.
emptyDataCell = cell(1,S);
for s = 1:S
    emptyDataCell{s} = zeros(0,1);
end

end

%% HELPER FUNCTION TO REFLECT SHOCK INDICES & SPLIT VALUES BY USAGE
function InversionData = ...
    reflect_shock_usage_indices_and_split_values_by_usage(...
    InversionData,zf,zInd,H)
% This helper creates additional information required for the inversion.
% It creates the indices of the shocks not being used in the inversion
% (i.e. it reflects the indices of the shocks being used) and it extracts
% the values for the shocks not being used and the previous values for the
% shocks that are being used (i.e. the values for the shocks being used in
% the inversion before they get updated/computed in the inversion).
%
% INPUTS:
%   -> InversionData: structure of inversion data
%   -> zf: matrix of shocks over the forecast
%   -> zInd: vector of model shock indices from 1 to nz
%   -> H: forecast horizon
%
% OUTPUTS:
%   -> InversionData: Updated structure of inversion data
%
% CALLS:
%   -> reflect_shock_usage_indices (sub-function)
%   -> split_shock_values_by_usage (sub-function)

%% REFLECT SHOCK USAGE INDICES
% Unpack the cell array of index values for the shocks to be used in the
% particular inversion type and call a sub-function to compute the indices
% of the shocks not being used in that inversion type. Pack the result into
% the structure.
ziInd = InversionData.shockUsageIndices;
znInd = reflect_shock_usage_indices(ziInd,zInd,H);
InversionData.shockNonUsageIndices = znInd;

%% SPLIT SHOCK VALUES BY USAGE
% Split the forecast shocks matrix period-by-period into the values that
% correspond to the shocks not being used to implement the inversion and
% the values that are being used. Pack the result into the structure.
[zn,zbari] = split_shock_values_by_usage(zf,znInd,ziInd,H);
InversionData.shockNonUsageValues = zn;
InversionData.shockUsagePreviousValues = zbari;

end

%% FUNCTION TO REFLECT SHOCK SET
function znInd = reflect_shock_usage_indices(ziInd,zInd,S)
% This helper reflects the shock usage indices to shocks not being used.
% This is a helper function which reflects the set of shocks used
% to implement an inversion of a particular type (anticipated or
% unanticipated) to find the indices and the values of the shocks not
% being used to implement the inversion.
%
% INPUTS:
%   -> ziInd: cell array of shock indices organised by period
%   -> zInd: indices of the shocks in the model (ie 1:nz)
%   -> S: horizon of the inversion
%
% OUTPUTS:
%   -> znInd: cell array of indices for the shocks not being used
%
% CALLS:
%   -> none

%% SETUP OUTPUT
% The index output is a cell array of vectors of the same length as the
% number of inversion periods.
znInd = cell(1,S);

%% COMPUTE THE SHOCK INDICES PERIOD-BY PERIOD
% Loop through each period of the inversion. In each period remove the
% shock indices that are being used in the inversion in that period. Note
% that the final line in the loop is needed to ensures that any empty
% index vectors are set to dimensions 0*1 rather than 1*0.
for s = 1:S
    znInd{s} = zInd;
    znInd{s}(ziInd{s}) = [];
    znInd{s} = znInd{s}(:);
end

end

%% FUNCTION TO SPLIT SHOCKS BY USAGE
function [zn,zbari] = split_shock_values_by_usage(zbar,znInd,ziInd,S)
% This helper splits the existing shocks into those used and those not used
% It uses the indices computed above to split the existing shocks into two
% (possibly time-varying) groups. One group is the set of shocks that will
% not be affected by the inversion because they are not being used to
% implement. The other group is those that will have new values after the
% inversion because they are being used to implement it.
%
% INPUTS:
%   -> zbar: matrix of existing shocks
%   -> znInd: cell array of indices for shocks not being used
%   -> ziInd: cell array of indices for shocks being used
%   -> S: horizon of the inversion
%
% OUTPUTS:
%   -> zn: cell array of values for the shocks not being used
%   -> zbari: cell array of values for the shocks being used
%
% CALLS:
%   -> none

%% SETUP OUTPUT
% Setup the output as two 1 by S cell arrays (one for each period of the
% inversion).
zn = cell(1,S);
zbari = cell(1,S);

%% SPLIT THE SHOCKS BY USAGE
% Run through the inversion periods, splitting the shocks by usage in each
% period.
for s = 1:S
    zn{s} = zbar(znInd{s},s);
    zbari{s} = zbar(ziInd{s},s);
end

end

%% FUNCTION TO VALIDATE THE INVERSION INPUT INVERSION INSTRUCTIONS
function InstructionsInclude = validate_inversion_instructions(...
    Instructions,H,xInd,zInd,ShocksInclude,F)
% This helper function validates the inversion instructions input.
% It checks the validity of the input and throws an exception if any of
% its fields are not valid. If they are valid, it returns a structure of
% indicator variables describing the content of the inversion instructions.
%
% INPUTS:
%   -> Instructions: structure detailing the inversion required
%   -> H: numeric scalar representing the forecast horizon
%   -> xInd: vector of model variable indices
%   -> zInd: vector of model shock indices
%   -> ShocksInclude: structure describing the content of the shocks input
%   -> F: matrix of forward loadings on the anticipated shocks
%
% OUTPUTS:
%   -> InstructionsInclude: structure describing the content of the 
%      inversion input with true or false indicators for the following:
%       - anticipatedFixes (optional)
%       - unanticipatedFixes (optional)
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause
%   -> check_inversion_structure_input (sub-function)

%% CHECK INPUT
% Check that the inversion instructions input is a structure. If not, throw
% an exception.
if ~isstruct(Instructions)
    errId = ['MAPS:',mfilename,':BadInversionInstructionsShape'];
    generate_and_throw_MAPS_exception(errId);
end

%% SETUP A MASTER EXCEPTION
% Setup a master exception to add causes to as they are encountered below.
masterErrId = ['MAPS:',mfilename,':BadInversionInstructions'];
InstructionsE = generate_MAPS_exception(masterErrId);

%% CHECK INVERSION ALGORITHM FIELD
% Check that the input contains the expected inversion algorithm field and
% that it is of the right type.
algorithmField = 'algorithmIsMinimumDamage';
if ~isfield(Instructions,algorithmField)
    errId = [masterErrId,':MissingAlgorithmField'];
    InstructionsE = generate_MAPS_exception_and_add_as_cause(...
        InstructionsE,errId);
else
    if ~islogical(Instructions.(algorithmField)) || ...
            ~isscalar(Instructions.(algorithmField))
        errId = [masterErrId,':BadAlgorithmField'];
        InstructionsE = generate_MAPS_exception_and_add_as_cause(...
            InstructionsE,errId,{algorithmField});
    end
end

%% SETUP OUTPUTS
% Set the content of the output to false. These are then overwritten below
% if the inputs do include the specified information
InstructionsInclude.anticipatedFixes = false;
InstructionsInclude.unanticipatedFixes = false;

%% CHECK CONTENT OF FIXES
% Check both of the fix types, overwriting the output to true where they
% exist. Check that the content of the structure meets the requirements of
% the inversion algorithm and throw an error if not.
fixTypes = {'Anticipated';'Unanticipated'};
nInversionTypes = size(fixTypes,1);
for iItype = 1:nInversionTypes
    iInversionType = [fixTypes{iItype},'Fixes'];
    if isfield(Instructions,iInversionType)
        iFixType = fixTypes{iItype};
        InstructionsInclude.([lower(iFixType),'Fixes']) = true;
        try
            validate_fixes_structure_input(...
                Instructions.(iInversionType),H,xInd,zInd);
        catch StructInputEcause
            errId = [masterErrId,':BadFixesStruct'];
            StructInputE = generate_MAPS_exception_and_add_cause(...
                StructInputEcause,errId,{iInversionType});
            InstructionsE = addCause(InstructionsE,StructInputE);
        end
    end
end

%% CHECK THAT THERE ARE ANY FIXES
% Check that that the input contains either anticipated or unanticipated
% fixes.
if ~InstructionsInclude.anticipatedFixes && ...
        ~InstructionsInclude.unanticipatedFixes
    errId = [masterErrId,':MissingFixes'];
    InstructionsE = generate_MAPS_exception_and_add_as_cause(...
        InstructionsE,errId);
end

%% CHECK CONSISTENCY OF INPUTS
% Check consistency of the inversion instructions with the model and shocks 
% input. In particular, the inversion instructions should not contain 
% anticipated fixes if the model is backward-looking (the F matrix is all 
% zeros) and if the inversion instructions includes anticipated and/or 
% unanticipated fixes then the shocks structure input must include 
% anticipated and/or unanticipated shocks.
if InstructionsInclude.anticipatedFixes
    if ~any(any(F))
        errId = ['MAPS:',mfilename,':BadUseOfAnticipatedFixes'];
        InstructionsE = generate_MAPS_exception_and_add_as_cause(...
            InstructionsE,errId);
    end
    if ~ShocksInclude.anticipated
        errId = ['MAPS:',mfilename,':MissingShocks'];
        InstructionsE = generate_MAPS_exception_and_add_as_cause(...
            InstructionsE,errId,{'anticipated'});
    end
end
if InstructionsInclude.unanticipatedFixes
    if ~ShocksInclude.unanticipated
        errId = ['MAPS:',mfilename,':MissingShocks'];
        InstructionsE = generate_MAPS_exception_and_add_as_cause(...
            InstructionsE,errId,{'unanticipated'});
    end
end

%% THROW MASTER EXCEPTION IF ANY CAUSES WERE ADDED
% If any exceptions were added above throw the master exception.
if ~isempty(InstructionsE.cause)
    throw(InstructionsE);
end

end

%% HELPER FUNCTION TO CHECK THE SHAPE OF THE FIXES INPUT
function validate_fixes_structure_input(Fixes,H,xInd,zInd)
% This helper valdiates the inversion fixes info structure.
% It checks that the generic fixes structures in the inversion info input
% are formatted correctly.
%
% INPUTS:
%   -> Fixes: structure containing the fix information
%   -> H: numeric scalar representing the forecast horizon
%   -> xInd: vector of model variable indices
%   -> zInd: vector of model shock indices
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% CHECK THE SHAPE OF THE INPUT
% Check that that the input is a structure. If not throw an expcetion.
if ~isstruct(Fixes)
    errId = ['MAPS:',mfilename,':BadFixDataType'];
    generate_and_throw_MAPS_exception(errId);
end

%% SETUP A MASTER EXCEPTION
% Setup a master exception to add causes to as they are encountered below.
masterErrId = ['MAPS:',mfilename,':BadFixData'];
FixesE = generate_MAPS_exception(masterErrId);

%% CHECK THE FIXES FIELDS
% For each of the expected fields in the fix structure, check that the
% field is present, is a row cell array and contains finite numeric column
% vectors. If not add an exception to the master exception.
fixFields = {'modelVarIndices';'modelVarValues';'shockUsageIndices'};
nFixFields = size(fixFields,1);
for iField = 1:nFixFields
    iFixField = fixFields{iField};
    if ~isfield(Fixes,iFixField)
        errId = [masterErrId,':MissingField'];
        FixesE = generate_MAPS_exception_and_add_as_cause(...
            FixesE,errId,{iFixField});
    elseif ~iscell(Fixes.(iFixField)) || ndims(Fixes.(iFixField))~=2 || ...
            size(Fixes.(iFixField),1)~=1
        errId = [masterErrId,':BadField'];
        FixesE = generate_MAPS_exception_and_add_as_cause(...
            FixesE,errId,{iFixField});
    elseif size(Fixes.(iFixField),2)>H
        errId = [masterErrId,':HorizonInconsistentWithForecast'];
        FixesE = generate_MAPS_exception_and_add_as_cause(...
            FixesE,errId,{iFixField});
    elseif any(cell2mat(cellfun(@(x) (~isnumeric(x)|ndims(x)~=2|...
            size(x,2)~=1|~all(isfinite(x(:)))),...
            Fixes.(iFixField),'UniformOutput',false)))
        errId = [masterErrId,':BadFieldContent'];
        FixesE = generate_MAPS_exception_and_add_as_cause(...
            FixesE,errId,{iFixField});
    end
end

%% THROW THE EXCEPTION IF ANY CAUSES WERE ADDED
% Throw any exceptions already encountered.
if ~isempty(FixesE.cause)
    throw(FixesE)
end

%% CHECK CONSISTENCY OF FIXES FIELD DIMENSIONS
% If the input passed the above basic validation, then check that the
% number of columns in each of the expected fields are consistent with each
% other (which amounts to saying that the input is consistent with a
% particular inversion horizon).
fieldPairs = {'modelVarIndices' 'modelVarValues';
    'modelVarIndices' 'shockUsageIndices'};
nFieldPairs = size(fieldPairs,1);
for iPair = 1:nFieldPairs
    if isfield(Fixes,fieldPairs{iPair,1}) && ...
            isfield(Fixes,fieldPairs{iPair,2})
        if size(Fixes.(fieldPairs{iPair,1}),2) ~= ...
                size(Fixes.(fieldPairs{iPair,2}),2)
            errId = [masterErrId,':FixesFieldsDimsMismatch'];
            FixesE = generate_MAPS_exception_and_add_as_cause(...
                FixesE,errId,{fieldPairs{iPair,1} fieldPairs{iPair,2}});
        end
    end
end

%% CHECK CONSISTENCY OF FIX VALUES & INDICES
% Check that the fix values and indices have consistent dimensions in each
% period. If not, add an exception detailing the periods in which there is
% a mismatch.
fixPair = {'modelVarIndices' 'modelVarValues'};
if isfield(Fixes,fixPair{1}) && isfield(Fixes,fixPair{2})
    if any(cellfun('size',Fixes.(fixPair{1}),1)-...
            cellfun('size',Fixes.(fixPair{2}),1))
        mismatchPeriods = find(cellfun('size',Fixes.(fixPair{1}),1)-...
            cellfun('size',Fixes.(fixPair{2}),1));
        errId = [masterErrId,':FixValIndsMismatch'];
        FixesE = generate_MAPS_exception_and_add_as_cause(...
            FixesE,errId,{fixPair{1} fixPair{2} num2str(mismatchPeriods)});
    end
end

%% CHECK INDICES ARE VALID
% Check that the model variable and shock usage indices are valid. That is,
% that they are unique and that they are in the model indices. The check is
% made on a period-by-period basis with information about the offending
% time period and indices added to an exception.
indsToCheck = {'modelVarIndices';'shockUsageIndices'};
modelIndsToCheckAgainst = {xInd;zInd};
nIndsToCheck = size(indsToCheck,1);
for iCheck = 1:nIndsToCheck
    modelInds = modelIndsToCheckAgainst{iCheck};
    iIndsToCheck = Fixes.(indsToCheck{iCheck});
    S = size(iIndsToCheck,2);
    for s = 1:S
        if size(unique(iIndsToCheck{s}),1) ~= size(iIndsToCheck{s},1)
            if ~exist('NonUniqueIndsE','var')
                errId = [masterErrId,':NonUniqueInds'];
                NonUniqueIndsE = generate_MAPS_exception(...
                    errId,{iIndsToCheck});
                errCauseId = [masterErrId,':BadInds:Instance'];
            end
            [uniqueInds,indUniqueInds] = unique(iIndsToCheck{s});
            tempInds = zeros(size(iIndsToCheck{s}));
            tempInds(indUniqueInds) = uniqueInds;
            indReps = ~find(tempInds);
            errArgs = {num2str(iIndsToCheck{s}(indReps)') num2str(s)};
            NonUniqueIndsE = generate_MAPS_exception_and_add_as_cause(...
                NonUniqueIndsE,errCauseId,errArgs);
        end
        if any(~ismember(iIndsToCheck{s},modelInds))
            if ~exist('UnknownIndsE','var')
                errId = [masterErrId,':UnknownInds'];
                UnknownIndsE = generate_MAPS_exception(...
                    errId,{iIndsToCheck});
                errCauseId = [masterErrId,':BadInds:Instance'];
            end
            unknownIndsLogicals = ~ismember(iIndsToCheck{s},modelInd);
            errArgs = {num2str(iIndsToCheck{s}(unknownIndsLogicals)') ...
                num2str(s)};
            UnknownIndsE = generate_MAPS_exception_and_add_as_cause(...
                UnknownIndsE,errCauseId,errArgs);
        end
    end
    if exist('NonUniqueIndsE','var')
        FixesE = addCause(FixesE,NonUniqueIndsE);
    end
    if exist('UnknownIndsE','var')
        FixesE = addCause(FixesE,UnknownIndsE);
    end
end

%% THROW THE EXCEPTION IF ANY CAUSES WERE ADDED
% Throw any exceptions already encountered.
if ~isempty(FixesE.cause)
    throw(FixesE)
end

end