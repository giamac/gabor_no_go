function data = nogo_staircase(subj_number)

trial_number = 40; %Number of trials in nogo_staircase

%Define Staircase parameters

initalMaskDuration = 0.128;
initalStimulusDuration = 0.1;
stepsize = 0.016;
blank = 0.032;

%Screen resolution
res = [1920 1080];

myHome = pwd;

%Create Subjectfolder and file

subjectsPath = ['subjects/' num2str(subj_number)]
if ~exist(subjectsPath, 'dir')
  mkdir(subjectsPath)
end

cd(subjectsPath);
file_name_mat = ['SubjNr_',num2str(subj_number),'_Staircase.mat'];
if exist(file_name_mat) == 2;
    display('logfile for current subject & run already exists');
    display('do you want to overwrite?');
    overwrite = input('enter = no / 5 = yes: ');
    if overwrite == 5;
        display(' '); display('will continue and overwrite...');
    else
        error('please check logfiles or specify a different run');
    end
end
cd(myHome);

Screen('Preference', 'SkipSyncTests', 1);
screenNum=0;
offx = 0; offy = 0;
maskDuration = initalMaskDuration;
stimulusDuration = initalStimulusDuration;
isi = 1;  % inter stimulus interval (0.5s default)
Trialtime = 2;
%%% offsets has the order x offset and y offset
offsets = [0 0];
maximum_value = 255;

%% Images
mask = imread('img/mask.png');
diamond_left = imread('img/left_diamond.png');
diamond_right = imread('img/right_diamond.png');


disp('Welcome to the No Go experiment');
%%% data.Subnum = input (['Enter subject number: ']);
data.Subnum = subj_number;
data.Date = date;
data.Data = [];
file_name_txt = ['SubjNr_',num2str(subj_number),'_Staircase.txt'];

cd(subjectsPath);
[fid message] = fopen(file_name_txt, 'w');
if fid == -1
    fprintf('Couldn''t open output file.\n%s\n', message);
end
fprintf(fid, 'Subject no: %d\n', data.Subnum);
fprintf(fid, 'trial\tresponse\tRT\tstart_trial\tcorrect\tCorr_answer\tPresentation_duration\r\n');
cd(myHome);

%% Add randomization of trials (left and right )

while KbCheck; end

[w,rect] = Screen('OpenWindow', screenNum, [255 255 255], [0 0 res(1) res(2)]);
% define window w and open screen
[xc,yc] = RectCenter(rect);

%polygon
polyWidth = 100;
polyHeight = 100;
xCoord = [xc, xc - (polyWidth/2), xc, xc + (polyWidth/2)]';
yCoord = [yc + (polyHeight /2), yc, yc - (polyHeight /2), yc]';
polyCoords = [xCoord yCoord];

%polygonOuter
polyWidthOuter = 160;
polyHeightOuter = 160;
xCoordOuter = [xc, xc - (polyWidthOuter/2), xc, xc + (polyWidthOuter/2)]';
yCoordOuter = [yc + (polyHeightOuter /2), yc, yc - (polyHeightOuter /2), yc]';
polyCoordsOuter = [xCoordOuter yCoordOuter];

%Right Rectangle
rightRect =  [(xc + (polyWidth/2) - 8) (yc-8) (xc + (polyWidth/2)) (yc+8)];

%left Rectangle
leftRect =  [(xc - (polyWidth/2)) (yc-8) (xc - (polyWidth/2) + 8) (yc+8)];
textSize = 30;
oldTextSize=Screen('TextSize', w, textSize); %%%,textSize
%%%keyboard

black = BlackIndex(w);  % Retrieves the CLUT color code for black.
white = WhiteIndex(w);  % Retrieves the CLUT color code for white.
grey  = GrayIndex(w);

%%% Prepare keys and trigger
KbName('UnifyKeyNames');
choice_key_names = {'b', 'z'}; % b = (blue)right, z(y) = (yellow)left
trigger_name = {'t'};
choice_keys = KbName(choice_key_names);
trigger_key = KbName(trigger_name);

%%% wait for scanner trigger (expects a "t" input!)
txt_color = black;
DrawFormattedText(w, 'Bitte warten Sie! Das Experiment startet gleich!', 'center', 'center', txt_color);
Screen(w, 'Flip');
WaitSecs(0.05);

waiting_for_scanner = 1;
%%% time_start_exp = GetSecs;
while waiting_for_scanner
    [key_is_down, ~, key_code] = KbCheck;
    if key_is_down && any(key_code(trigger_key))
        waiting_for_scanner = 0;
        time_start_exp = GetSecs;
    end
  end

%% Create trialList
%% 1 for left and 2 for right
totalTrials = repmat([1; 2],trial_number/2,1);
trials = zeros(trial_number);
trials = totalTrials(randperm(length(totalTrials)));

%% Staircase
correctInaRow = 0;

%Start the task here
for t = 1:trial_number
  Screen('ColorRange', w , white, 0, 1);
  start_time_trial(t) = GetSecs - time_start_exp;
  Screen('FillRect', w, white);
%  imageTexture = Screen('MakeTexture', w, diamond_left);
%  Screen('DrawTexture', w, imageTexture, [], [], 0);
%  Screen('Flip',w);
%  Screen('FramePoly', w, [0 0 0], polyCoords, 3);

%Display the stimulus
Screen('FillPoly', w, [0 0 0], polyCoords, 3);
if trials(t) == 1
  Screen('FillRect', w, white, leftRect );
elseif trials(t) == 2
  Screen('FillRect', w, white, rightRect );
end

Screen('Flip', w);
WaitSecs(stimulusDuration);

%% Blank time
Screen('FillRect', w, white);
Screen('Flip', w, white);
WaitSecs(blank);

%% Show mask
Screen('FillPoly', w, [0 0 0], polyCoordsOuter, 3);
%Screen('FillPoly', w, [255 255 255], polyCoords, 3);
Screen('Flip', w);
WaitSecs(maskDuration);
DrawFormattedText(w, '?', 'center', 'center', black); %draw question mark
Screen(w, 'Flip'); % display gabor by flipping
give_response = false;
      Starttime = GetSecs;
      while GetSecs < Starttime + Trialtime

          [keyIsDown, secs, keycode] = KbCheck;
          if keyIsDown && any(keycode(choice_keys)) %%% if one of the two keys are pushed
              down_key = find(keycode, 1);
              response(t) = down_key; %%% response is 66 for right and 90 for left
              resptime(t) = secs - Starttime;
              give_response = true;

              if down_key == 66 || down_key == 90 %%% blue = left
                  Screen(w, 'FillRect', white);
                  txt_color = [ 0 0.75*maximum_value 0]; %%% green
                  DrawFormattedText(w, '?', 'center', 'center', txt_color); %draw question mark
                  Screen(w, 'Flip'); % display gabor by flipping
                  WaitSecs(0.5);
              end
              break
          else
              response(t) = 1000;
              resptime(t) = 999;
          end
      end

      if  give_response == false

          txt_color = [0.75*maximum_value 0 0]; %%% red
          DrawFormattedText(w, '?', 'center', 'center', txt_color); %draw question mark
          Screen(w, 'Flip'); % display gabor by flipping
          WaitSecs(0.5);
      else
           txt_color = [0 0.75*maximum_value 0]; %%% red
           DrawFormattedText(w, '?', 'center', 'center', txt_color); %draw question mark
           Screen(w, 'Flip'); % display gabor by flipping
           WaitSecs(0.5);
      end
      presentation_duration(t) = stimulusDuration;

%Correct Answer for task
if trials(t) == 1
  Corr_answer(t) = 66;
else
  Corr_answer(t) = 90;
end
%%% dificulty level
if response(t) == Corr_answer(t)
  correct(t) = 1;
  correctInaRow = correctInaRow + 1;
  if correctInaRow == 3
    maskDuration = maskDuration + stepsize;
    stimulusDuration = stimulusDuration - stepsize;
    correctInaRow = 0;
  end
elseif response(t) ~= Corr_answer(t)
  correct(t) = 0;
  maskDuration = maskDuration - stepsize;
  stimulusDuration = stimulusDuration + stepsize;
  correctInaRow = 0;
end

%% Avoid presentation times below 0

if stimulusDuration < 0
  stimulusDuration = 0.016;
  maskDuration = 0.210;
end

%% Work with frame rates

%%%% record reaction time data into data
  data.Data = [data.Data resptime];

  %%% Write the trial information to the text file
  fprintf (fid, '%d\t%s\t%f\t%f\t%d\t%s\t%f\t\r\n', t, response(t), resptime(t), start_time_trial(t), correct(t), Corr_answer(t), presentation_duration(t));
  txt_color = [0.9*maximum_value];
  DrawFormattedText(w, '?', 'center', 'center', txt_color); %draw question mark
  WaitSecs(isi);
  end
  Screen('CloseAll');
  fclose(fid);
  cd(subjectsPath);
save(file_name_mat, 'response','resptime', 'start_time_trial', 'Corr_answer', 'noise', 'time_start_exp', 'correct');
cd(myHome);
end
