%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%copyright (c) 2012  Matthew Caudill
%
%this program is free software: you can redistribute it and/or modify
%it under the terms of the gnu general public license as published by
%the free software foundation, either version 3 of the license, or
%at your option) any later version.

%this program is distributed in the hope that it will be useful,
%but without any warranty; without even the implied warranty of
%merchantability or fitness for a particular purpose.  see the
%gnu general public license for more details.

%you should have received a copy of the gnu general public license
%along with this program.  if not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trials] = trialStruct_SquareMoveContinuous( stimType, table)
% TRIALSTRUCT creates a structure array based on the values from the Gui
% parameter table. The trial structure arrary contains fields that are
% similar to the table strings (white spaces are removed). The number of
% structures in the structure array matches the total number of trials and
% the parameters for each trial can be called out of the structure array
% using trials(i). This command returns a structure for trial number i
% For example lets say you run 12 gratings at 12 orientations and 0 repeats 
% the trial array structure will contain 144 structures. The variable that 
% changes fastest is the one listed last in the Gui table (if no 
% randomization selected). In this case orientation  would vary first then 
% the contrast.
% INPUTS
% STIMTYPE:         STRING SUPPLIED BY STIMGEN GUI
% TABLE:            CELL ARRAY SUPPLIED BY STIMGEN GUI TABLE
%
% OUTPUTS
% TRIALS:   A STRUCTURE ARRAY CONTAINING ALL STIMULUS TRIAL INFO

%%%% TESTING INPUT TABLE (MSC 4-3-12)
% if nargin<1
%   table = {'Spatial Frequency (cpd)', 0.04, .04, 0.04;...
%               'Temporal Frequency (cps)', 3, 1, 3;...
%               'Contrast (start,end,numsteps)', 1, 1, 1;...
%               'Orientation', 0, 30, 330;...
%               'Timing (delay,duration,wait) (s)', 1, 2, 1;...
%               'Blank', 1, [], [];... 
%               'Randomize', 1, [], [];...
%               'Interleave', 1, [], [];...
%               'Interleave Timing', 1, 2, 10;...
%               'Repeats', 0, [], [];...
%               'Initialization Screen (s)', 5, [], []};
%    stimType = 'Full-field Grating';
%   trials = trialStruct(stimType, table);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% CONSTRUCT FIELDNAMES FOR TRIALSTRUCT  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First locate all the special characters like (). We will break on these
% to construct the fieldnames of the structure. Fieldnames will follow the
% format 'Spatial_Frequency' etc since spaces are not allowed in fields

% Initialize a structure to hold parameters (see explanantion below)
params = struct();
% Initialize a structure to hold constants
constants = struct();

for i = 1:size(table,1)
    %Find the strings from the 1st column of the table breaking at the '('
    tableStrings{i} = strtok(table{i,1},'(');
    %remove trailing spaces
    tableStrings{i} = strtrim(tableStrings{i});
    %construct a fieldname by replacing the white space with '_'
    fieldname{i} = strrep(tableStrings{i},' ','_');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%CREATE A PARAMETERS AND CONSTANTS STRUCTURE %%%%%%%
    %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We will create 2 temporary structures. One will hold parameters to be
    % varied and the other will hold constants (we initialized them above
    % because we want to check whether they are empty later)
    
    % RULES FOR VALID PARAMETERS
    % If all columns after the first are numbers and if the numel of a row 
    % in the table >2 and the first and last values are not equal then we 
    % have a parameter to add to the params struct.
    % Note we meed to exclude the timing becasue it does not follow
    % start:step:end AND contrast becasue it will be varied logarithmically
    % not linearly and also does not follow start:step:end but rather
    % start,end,num_steps
    
    % Check if row contains strings, if so store to constants structure
    if ischar([table{i,2:end}])
        constants.(fieldname{i}) = {table{i,2:end}};
        
    % Check if row is the timing row if so add to constants structure    
    elseif strcmp(fieldname{i},'Timing')
        constants.(fieldname{i}) = horzcat(table{i,2:4});
    elseif strcmp(fieldname{i},'Timing_End')
        constants.(fieldname{i}) = horzcat(table{i,2:4});
        
    elseif strcmp(fieldname{i}, 'Interleave_Timing')
        constants.(fieldname{i}) = horzcat(table{i,2:4});
        interleaveTiming = horzcat(table{i,2:4});
    
    % Check if row is 'AV_V_A_ratio' for A-V conditioning test    
    elseif strcmp(fieldname{i},'omit_ratio')
        constants.(fieldname{i}) = horzcat(table{i,2:4});
    
%     % Get Square Position Vectors    
%     elseif strcmp(fieldname{i},'Square_PositionX_vector')
%         constants.(fieldname{i}) = table{i,2};
%         
%     elseif strcmp(fieldname{i},'Square_PositionY_vector')
%         constants.(fieldname{i}) = table{i,2};
        
    % Check number of elements in the row if less than 2 add to constants
    elseif numel(cell2mat(table(i,2:end))) <= 2 % 2 element rows
        constants.(fieldname{i}) = horzcat(table{i,2:4});
        
    % Check if row is Contrast row and the check whether exponent 1 equals
    % exponent 2
    % elseif strcmp(fieldname{i},'Contrast') && table{i,3} == table{i,4}
    elseif strcmp(fieldname{i},'Contrast') && table{i,2} == table{i,4}
        % The constant contrast will be the highest contrast or 1
        constants.(fieldname{i}) = min(table{i,2}^table{i,4},100)/100;
        
    % Else if row is contrast and exp1 ~= exp 2
    elseif strcmp(fieldname{i},'Contrast') && table{i,3} ~= table{i,4}
        % we will loop over the exponent range and construct the contrast
        for exponent = table{i,3}:table{i,4}
            contrast(exponent) = table{i,2}^exponent/100;
        end
        % we limit the contrast to values <= 1 and use one only once
        contrast = unique(min(contrast,1));
        params.Contrast = contrast;
                                
    % Check whether start and end vals are the same. If so add to constants
    elseif numel(cell2mat(table(i,2:end))) >2 && table{i,2}==table{i,4}
        constants.(fieldname{i}) = table{i,2};
    % Else if start and end not same add to the params structure
    elseif numel(cell2mat(table(i,2:end))) >2 && table{i,2}~=table{i,4}
        params.(fieldname{i}) = table{i,2}:table{i,3}:table{i,4};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCT TRIAL ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SqPosNum = length(constants.Square_PositionX_vector);

trial_arrays = [constants.Square_PositionX_vector, constants.Square_PositionY_vector, zeros(SqPosNum,1)];

% omitIdx = trial_arrays(:,1)==constants.Audio_TargetX & trial_arrays(:,2)==constants.Audio_TargetY;
% trial_arrays(audioIdx,3)=1;


%keyboard;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% RANDOMIZE/INTERLEAVE/REPEAT TRIAL ARRAYS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We will break the randomize and inerleave into cases. CASE 1:
    % Interleave and NO Randomize. CASE 2: Interleave and Randomize. CASE 
    % 3: Randomize and no interleave. Finally if asked for we repeat.
    %
    % special randomization for Receptive field mapping. 
    % Only for case 3: randomize and no interleave
    %
    
    %trial_arrays_cat = [];
    OmitFlag_cat = [];
    for i=1:(constants.Repeats+1)
        
        % RANDOMIZE AND NO INTERLEAVE

            
                new_trial_arrays = trial_arrays;

                fullnum = constants.omit_ratio(1);
                omitnum = constants.omit_ratio(2);
                allnum = fullnum + omitnum-1;
                if omitnum>0
                p=randperm(allnum)==allnum;
                else
                   p=false(allnum,1)'; 
                end
                OmitFlag = [p'; false];
                
        OmitFlag_cat = [OmitFlag_cat; OmitFlag];
        
    end
    
    trial_arrays = OmitFlag_cat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ADD ARRAYS TO TRIAL STRUCTURE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Make the trials array into a cell array
    trialCell=num2cell(trial_arrays);
%     % If we have interleaved then we add the new parameter Led_Condition to
%     % the params_fields and fieldname
%     %if constants.Interleave == 1
%         params_fields{end+1} = 'Audio_Cue';
%         fieldname{end+1} = 'Audio_Cue';
%     %end

    % params_fields = {'Square_PositionX','Square_PositionY','OmitFlag'};
    params_fields = {'OmitFlag'};
    % Call the funct cell2struct along the 2nd Dimension to make our struct
    trials = cell2struct(trialCell,params_fields,2);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% ADD CONSTANTS TO TRIAL STRUCTURE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get the fieldnames of the constants struct
constant_fields = fieldnames(constants);
%Now add the constants from the constants structure to the trial structure
for i=1:numel(constant_fields)
    [trials(:).(constant_fields{i})] =deal(constants.(constant_fields{i}));
    % [trials.(constant_fields{i})] =constants.(constant_fields{i});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% ADD STIMTYPE TO TRIAL STRUCTURE %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now add stimType to the structure
[trials(:).Stimulus_Type] = deal(stimType);

% For the blank trials we want to set the Stimulus_type to Blank (note
% blanks only exist if a parameter exist so we check for this
if ~isempty(fieldnames(params))
[trials(blankTrials).Stimulus_Type] = deal('Blank'); 
end

% keyboard;
% 
% % Since we added stimType last reorder the fields of the structure so
% % stimType appears first and all other fields appear in the order in which
% % they arrive in from the Gui table
% trials = orderfields(trials, ['Stimulus_Type',fieldname]);

% clean up the trials structure by removing unneccessary fields
% remove the interleave timing field since its values are stored in timing
if isfield(trials,'Interleave_Timing')
    trials = rmfield(trials, 'Interleave_Timing');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%% REPEAT CONSTANT TRIAL STRUCTURE %%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % If there are no parameters but the user has chosen to repeat the
% % stimulus we must increase the trial structure from a 1x1 structure to a
% % size of repeats x 1 structure of arrays. So check that there are no
% % parameters first then if so repeat the trial structure to the number of
% % times the user selected by appending the first trial.
% if isempty(fieldnames(params)) && trials.Repeats > 0
%     for i=1:trials.Repeats
%         trials=[trials;trials(1)];
%     end
% end


end


