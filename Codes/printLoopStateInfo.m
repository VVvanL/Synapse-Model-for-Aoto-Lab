function printLoopStateInfo

global timestep_internal
global timeSteps_internal

timestep_internal = timestep_internal + 1;
firstDec = floor(((timestep_internal/timeSteps_internal)*100-floor(timestep_internal/timeSteps_internal*100))*10);

fprintf('\b\b\b\b\b%3.0f.%1.0f', floor((timestep_internal/timeSteps_internal)*100), firstDec * (firstDec < 10));

