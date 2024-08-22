function cyclingDetected = check_OBC_guesses_for_cycling(OBChistory)
% This function checks to see if latest "guess" for occasionally binding
% constraint regimes matches any previously tried guesses. Since the
% previously tried guesses do not constitute an equilibrium, finding a
% match indicates that the solution algorithm is cycling through a
% determinstic sequence of guesses that do not constitute equilibria.
%
% INPUTS
%   -> OBChistory: N*iter numeric matrix indicating OBC regime
%      indicators (note that these are (row) vectorised for cases with
%      multiple OBCs).
%
% OUTPUTS


%% CHECK THERE IS ENOUGH HISTORY OF GUESSES TO TEST
if size(OBChistory,2)<2
    error('OBChistory matrix must have at least two columns (guesses).');
end

%% SPLIT SEQUENCE OF GUESSES INTO PREVIOUSLY TRIED GUESSES AND LATEST GUESS
previousTries = OBChistory(:,1:end-1);
latestGuess = OBChistory(:,end);

%% TEST WHETER LATEST GUESS HAS ALREADY BEEN TRIED
identicalGuesses = (previousTries==latestGuess);
cyclingTest = all(identicalGuesses);
[~,repeatGuessInds] = find(cyclingTest==true,1,'last');
if ~isempty(repeatGuessInds)
    cyclingDetected = true;
else
    cyclingDetected = false;
end
 
end


