
function data = gabor_experiment_staircase(subj_number)
% call prep_gabor.m to prepare the gabor and display it similar to what we
% did in display_gabor.m
%%%difficulty_level_temp = [0 2.8:0.4:3.4];

difficulty_level_temp = [0 2.4:0.4:4.4];

%% Define Startpoint of Dicciulty Level
difficulty_startpoint = 2.4;

%% Define Stepsize for difficulty
difficulty_stepsize = 0.4;

%% Number of trials
trial_number = 10;

res = [800 600];
res=[1920 1080]; % screen resolution 800 x 600, adapt this resolution here (of the primary)

% % % difficulty_level_temp = [0 1.8 3.8];
% % % nr_of_trial_per_diffLevel = 2;
% % % res=[1280 768]; % screen resolution 800 x 600


%%% load timming trials
myHome = pwd;

%%% Change this, create directories instead of this

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
screenNum=2;
offx = 0; offy = 0;
Trialtime = 2;
isi = 1;  % inter stimulus interval (0.5s default)
%%% offsets has the order x offset and y offset
offsets = [0 0];
maximum_value = 255;

disp('Welcome to the gabor experiment');
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
fprintf(fid, 'trial\tresponse\tRT\tstart_trial\tcorrect\tCorr_answer\tDifficulty\r\n');
cd(myHome);

%%% rng('shuffle');

%%% define the gabor parameters: sd, sf, theta
stim_freq = 10; stim_ang = 3; stim_sd = 0.5;
trials_temp = [stim_sd stim_freq stim_ang; stim_sd stim_freq (360 - stim_ang)];

total_trails_nr = trial_number;

%%% create the gabor trials
%%% create total_trails_nr trials, altough at the end we will use only half of it. The other half are control trials.
trials = zeros(total_trails_nr,3);

%% Create shuffled order of angle

nr_angles = 2;

trials = repmat(trials_temp, (total_trails_nr/nr_angles), 1);
angle_perm = trials(:,3)(randperm(length(trials)));
trials(:,3) = angle_perm;

while KbCheck; end

[w,rect] = Screen('OpenWindow', screenNum);
% define window w and open screen
[xc,yc] = RectCenter(rect);

textSize = 30;
oldTextSize=Screen('TextSize', w, textSize); %%%,textSize
%%%keyboard

black = BlackIndex(w);  % Retrieves the CLUT color code for black.
white = WhiteIndex(w);  % Retrieves the CLUT color code for white.
grey  = GrayIndex(w);
%%% grey = (black + white) / 2;  % Computes the CLUT color code for gray.
inc = abs (white - grey);
%%% Prepare keys and trigger
KbName('UnifyKeyNames');
choice_key_names = {'b', 'z'}; % b = (blue)right, z(y) = (yellow)left
escape = KbName('ESCAPE');
trigger_name = {'t'};
left_name = {'b'};
right_name = {'z'};
choice_keys = KbName(choice_key_names);
trigger_key = KbName(trigger_name);
left_key = KbName(left_name);
right_key = KbName(right_name);

%% Training
txt_color = [200 200 200];
DrawFormattedText(w, 'In each trial you will see a circle with stripes that are either angled to the left or to the right. \n Your task will be to indicate the direction of the stripes.', 'center', 'center',  [0 0 0]);
Screen(w, 'Flip');
WaitSecs(0.05);
waiting_for_scanner = 1;
%%% time_start_exp = GetSecs;
while waiting_for_scanner
    [key_is_down, ~, key_code] = KbCheck;
    if key_is_down && any(key_code(trigger_key))
        waiting_for_scanner = 0;
    elseif key_code(escape)
            ShowCursor;
            sca;
            return
    end

end

% Example 1 - Left
DrawFormattedText(w, 'Example 1', 'center', 'center', [0 0 0]);
Screen(w, 'Flip');
WaitSecs(2);
maximum_value =256+600;
%%% [oldmaximumvalue, oldclampcolors, oldapplyToDoubleInputMakeTexture] = Screen('ColorRange', w , maximum_value, 0, 1);
Screen('ColorRange', w , maximum_value, 0, 1);
grey  = 256;
matsize = 400;
gabor = prep_gabor (0.5, 10, 357, matsize);
gabor = grey + inc * gabor*0.1 + 0.4;

[gw, gh] = size (gabor); % width and height of gabor
gabor_int_1 = max(max(gabor));
noiseimg=(50*randn(matsize+1, matsize+1) + 128);
noiseimg =  abs(noiseimg);
noiseimg = noiseimg/max(max(noiseimg));
noiseimg_scaled = imadjust(noiseimg,[0.1 1.0],[0.1 1]);
noiseimg_scaled =  noiseimg_scaled*1.4; %%% difficulty_level(t)

    %%% aplly gaussian filtering to the noise
    sigma = 200; imSize = matsize+1;  trim = 0.0001;
    s = sigma / imSize;                     % gaussian width as fraction of imageSize
    X = 1:imSize;                           % X is a vector from 1 to imageSize
    X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
    [Xm Ym] = meshgrid(X0, X0);
    gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
    gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)

    noiseimg_scaled = noiseimg_scaled .* gauss;


    gabor = gabor/max(max(gabor));
    gabor = (gabor + noiseimg_scaled);
    gabor = gabor/mean(mean(gabor))*gabor_int_1;

    %%% aplly gaussian filtering to the noise
    sigma = 70; imSize = matsize+1;  trim = 0.01;
    s = sigma / imSize;                     % gaussian width as fraction of imageSize
    X = 1:imSize;                           % X is a vector from 1 to imageSize
    X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
    [Xm Ym] = meshgrid(X0, X0);
    gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
    gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)
    gabor = gabor .* gauss;
    gabor = gabor + grey;

    gabortex = Screen('MakeTexture', w, gabor, [], [], []); % make gabor texture
    location = [xc - gw/2 + offsets(1,1), yc - gh/2 + offsets(1,2), xc + gw/2 + offsets(1,1), yc+ gh/2+ offsets(1,2)];

     %%% display gabor
    Screen(w, 'FillRect', grey);
    Screen(w, 'DrawTexture', gabortex, [], location); % draw gabor
    txt_color = [0.9*maximum_value 0.9*maximum_value 0.9*maximum_value];
    DrawFormattedText(w, '?', 'center', 'center', txt_color); %draw question mark
    Screen(w, 'Flip'); % display gabor by flipping

%Response
waitForResponse = 1
while waitForResponse
  [key_is_down, ~, key_code] = KbCheck;
  if key_is_down && any(key_code(left_key))
            Screen(w, 'FillRect', grey);
            Screen(w, 'DrawTexture', gabortex, [], location); % draw gabor
            txt_color = [ 0 0.75*maximum_value 0]; %%% green
            DrawFormattedText(w, 'Correct', 'center', res(2) - res(2)/3, txt_color); %draw question mark
            Screen(w, 'Flip'); % display gabor by flipping
            WaitSecs(0.5);
            waitForResponse=0;
        elseif key_code(escape)
            ShowCursor;
            sca;
            return
        end
    end

    % Example 2 - Right
    DrawFormattedText(w, 'Example 2', 'center', 'center', [0 0 0]);
    Screen(w, 'Flip')
    WaitSecs(2);
    maximum_value =256+600;
    %%% [oldmaximumvalue, oldclampcolors, oldapplyToDoubleInputMakeTexture] = Screen('ColorRange', w , maximum_value, 0, 1);
    grey  = 256;
    matsize = 400;
    gabor = prep_gabor (0.5, 10, 3, matsize);
    gabor = grey + inc * gabor*0.1 + 1;

    [gw, gh] = size (gabor); % width and height of gabor
    gabor_int_1 = max(max(gabor));
    noiseimg=(50*randn(matsize+1, matsize+1) + 128);
    noiseimg =  abs(noiseimg);
    noiseimg = noiseimg/max(max(noiseimg));
    noiseimg_scaled = imadjust(noiseimg,[0.1 1.0],[0.1 1]);
    noiseimg_scaled =  noiseimg_scaled*1.4; %%% difficulty_level(t)

        %%% aplly gaussian filtering to the noise
        sigma = 200; imSize = matsize+1;  trim = 0.0001;
        s = sigma / imSize;                     % gaussian width as fraction of imageSize
        X = 1:imSize;                           % X is a vector from 1 to imageSize
        X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
        [Xm Ym] = meshgrid(X0, X0);
        gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
        gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)

        noiseimg_scaled = noiseimg_scaled .* gauss;


        gabor = gabor/max(max(gabor));
        gabor = (gabor + noiseimg_scaled);
        gabor = gabor/mean(mean(gabor))*gabor_int_1;

        %%% aplly gaussian filtering to the noise
        sigma = 70; imSize = matsize+1;  trim = 0.01;
        s = sigma / imSize;                     % gaussian width as fraction of imageSize
        X = 1:imSize;                           % X is a vector from 1 to imageSize
        X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
        [Xm Ym] = meshgrid(X0, X0);
        gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
        gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)
        gabor = gabor .* gauss;
        gabor = gabor + grey;

        gabortex = Screen('MakeTexture', w, gabor, [], [], []); % make gabor texture
        location = [xc - gw/2 + offsets(1,1), yc - gh/2 + offsets(1,2), xc + gw/2 + offsets(1,1), yc+ gh/2+ offsets(1,2)];

         %%% display gabor
        Screen(w, 'FillRect', grey);
        Screen(w, 'DrawTexture', gabortex, [], location); % draw gabor
        txt_color = [0.9*maximum_value 0.9*maximum_value 0.9*maximum_value];
        DrawFormattedText(w, '?', 'center', 'center', txt_color); %draw question mark
        Screen(w, 'Flip'); % display gabor by flipping

    %Response
    waitForResponse = 1
    while waitForResponse
      [key_is_down, ~, key_code] = KbCheck;
      if key_is_down && any(key_code(right_key))
                Screen(w, 'FillRect', grey);
                Screen(w, 'DrawTexture', gabortex, [], location); % draw gabor
                txt_color = [ 0 0.75*maximum_value 0]; %%% green
                DrawFormattedText(w, 'Correct', 'center', res(2) - res(2)/3, [0 0 0]); %draw question mark
                Screen(w, 'Flip'); % display gabor by flipping
                WaitSecs(0.5);
                waitForResponse=0;
            elseif key_code(escape)
                ShowCursor;
                sca;
                return
            end
        end




%%% wait for scanner trigger (expects a "t" input!)
txt_color = [200 200 200];
DrawFormattedText(w, 'Bitte warten Sie! Das Experiment startet gleich!', 'center', 'center',  [0 0 0]);
Screen(w, 'Flip');
WaitSecs(0.05);


waiting_for_scanner = 1;
%%% time_start_exp = GetSecs;
while waiting_for_scanner
    [key_is_down, ~, key_code] = KbCheck;
    if key_is_down && any(key_code(trigger_key))
        waiting_for_scanner = 0;
        time_start_exp = GetSecs;

    elseif key_code(escape)
            ShowCursor;
            sca;
            return
    end

end
%% Parameters for noise
amount_noise = difficulty_startpoint;
%% Staircase
correctInaRow = 0;
%%%

for t = 1: total_trails_nr

    maximum_value =256+600;
    %%% [oldmaximumvalue, oldclampcolors, oldapplyToDoubleInputMakeTexture] = Screen('ColorRange', w , maximum_value, 0, 1);
    Screen('ColorRange', w , maximum_value, 0, 1);
    grey  = 256;

    start_time_trial(t) = GetSecs - time_start_exp;
    matsize = 400;
    gabor = prep_gabor (trials (t,1), trials (t,2), trials(t,3), matsize);
    gabor = grey + inc * gabor*0.1 + 1;

    [gw, gh] = size (gabor); % width and height of gabor
    gabor_int_1 = max(max(gabor));
    noiseimg=(50*randn(matsize+1, matsize+1) + 128);
    noiseimg =  abs(noiseimg);
    noiseimg = noiseimg/max(max(noiseimg));
    noiseimg_scaled = imadjust(noiseimg,[0.1 1.0],[0.1 1]);
    noiseimg_scaled =  noiseimg_scaled*amount_noise; %%% difficulty_level(t)
    disp(t);

        %%% aplly gaussian filtering to the noise
        sigma = 200; imSize = matsize+1;  trim = 0.0001;
        s = sigma / imSize;                     % gaussian width as fraction of imageSize
        X = 1:imSize;                           % X is a vector from 1 to imageSize
        X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
        [Xm Ym] = meshgrid(X0, X0);
        gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
        gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)

        noiseimg_scaled = noiseimg_scaled .* gauss;


        gabor = gabor/max(max(gabor));
        gabor = (gabor + noiseimg_scaled);
        gabor = gabor/mean(mean(gabor))*gabor_int_1;

        %%% aplly gaussian filtering to the noise
        sigma = 70; imSize = matsize+1;  trim = 0.01;
        s = sigma / imSize;                     % gaussian width as fraction of imageSize
        X = 1:imSize;                           % X is a vector from 1 to imageSize
        X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
        [Xm Ym] = meshgrid(X0, X0);
        gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
        gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)
        gabor = gabor .* gauss;
        gabor = gabor + grey;

        gabortex = Screen('MakeTexture', w, gabor, [], [], []); % make gabor texture
        location = [xc - gw/2 + offsets(1,1), yc - gh/2 + offsets(1,2), xc + gw/2 + offsets(1,1), yc+ gh/2+ offsets(1,2)];

         %%% display gabor
        Screen(w, 'FillRect', grey);
        Screen(w, 'DrawTexture', gabortex, [], location); % draw gabor
        txt_color = [0.9*maximum_value 0.9*maximum_value 0.9*maximum_value];
        DrawFormattedText(w, '?', 'center', 'center', txt_color); %draw question mark
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
                    Screen(w, 'FillRect', grey);
                    Screen(w, 'DrawTexture', gabortex, [], location); % draw gabor
                    txt_color = [ 0 0.75*maximum_value 0]; %%% green
                    DrawFormattedText(w, '?', 'center', 'center', txt_color); %draw question mark
                    Screen(w, 'Flip'); % display gabor by flipping
                    WaitSecs(0.5);
                end
                break
          elseif keycode(escape)
                ShowCursor;
                sca;
                return
            else
                response(t) = 1000;
                resptime(t) = 999;
            end
        end

        %%% display message if they do not give an answer
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
        noise(t) = amount_noise;
        %%% Save other stimulus info
        %%% corect answer for gabor
        if trials(t,3) == stim_ang
            Corr_answer(t) = 90;
        else
            Corr_answer(t) = 66;
        end
        %%% dificulty level
        if response(t) == Corr_answer(t)
          correct(t) = 1;
          correctInaRow = correctInaRow + 1;
          if correctInaRow == 3
            amount_noise = amount_noise + difficulty_stepsize;
            correctInaRow = 0;
          end
        elseif response(t) ~= Corr_answer(t)
          correct(t) = 0;
          amount_noise = amount_noise - difficulty_stepsize;
          correctInaRow = 0;
        end

    if amount_noise < 0
        amount_noise = 0;
    end




    %%%% record reaction time data into data
    data.Data = [data.Data resptime];

    %%% Write the trial information to the text file
    fprintf (fid, '%d\t%s\t%f\t%f\t%d\t%s\t%f\t\r\n', t, response(t), resptime(t), start_time_trial(t), correct(t), Corr_answer(t), noise(t));

    txt_color = [0.9*maximum_value];
    DrawFormattedText(w, '?', 'center', 'center', txt_color); %draw question mark
    WaitSecs(isi);
    end
    Screen('CloseAll');
    fclose(fid);




%%% save data in mat file
cd(subjectsPath);
save(file_name_mat, 'response','resptime', 'start_time_trial', 'Corr_answer', 'noise', 'time_start_exp', 'correct');
cd(myHome);

%% Some statistics

figure(1)
clf
stairs(noise);

noise_levels = unique(noise);

nCorrect = zeros(1,length(noise_levels));
ntrials = zeros(1,length(noise_levels));

for i=1:length(noise_levels)
  id = noise == noise_levels(i) & ~isnan(response);
  ntrials(i) = sum(id);
  nCorrect(i) = sum(correct(id));
end

pCorrect = nCorrect ./ ntrials;

figure(2)
clf
hold on
plot(noise_levels,pCorrect*100,'-','MarkerFaceColor','b');
%loop through each intensity so each data point can have it's own size.
for i=1:length(noise_levels);
   sz = ntrials(i)+2;
   plot(noise_levels(i),pCorrect(i)*100,'ko-','MarkerFaceColor','b','MarkerSize',sz);
end

set(gca,'XTick',noise_levels);
set(gca,'YLim',[0,100]);
xlabel('Noise');
ylabel('Percent Correct');
clear all
end
