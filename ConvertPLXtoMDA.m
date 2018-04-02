function ConvertPLXtoMDA
%% ConvertPLXtoMDA
%   File for converting PLX data files into MDA files in order to use them
%   with the MountainSort module in MountainLab.
%
%   NOTE: Make sure you're using an UNCUT file to go here.
%
%   - 01/22/2018 Created by GE; created using process described at:
%       https://github.com/mari-sosa/Mountainsort_for_snippets/blob/master/mountainsort_for_snippets.md
%
%% Identify PLX file
[fileName, filePath] = uigetfile('.plx','Identify .PLX File');
if fileName == 0
    disp('No file selected');
    return
end
plxFile = [filePath fileName];
origDir = cd;
cd(filePath);

%% Identify channels and waveforms
% Pull out relevant information from the recording file
%   sampRate = sample rate the recording was done at
%   nTrodes = number of channels per recording, tetrode = 4
%   preThresh = number of samples recorded prior to threshold crossing for
%       trace capture
%   ssnDur = duration of the recording in seconds
[~,~,sampRate,~,nTrodes,~,preThresh,~,~,~,~,ssnDur,~] = plx_information(plxFile);
% Create a blank vector where the spike waveforms can be filled in (below)
%   NOTE: Nans are used here because the blanks are filled in by connecting
%   the traces with a linear vector. This is done rather than using 0s
%   because the waveforms are sometimes DC shifted so using 0s makes it so
%   the edges have very sharp jumps back to 0 which MountainSort picks up
%   on as unit activity.
truncLim = round(sampRate*ssnDur);
mptRaw = nan(4,truncLim);
% Pull out the timestamp (spike) count in order to direct the subsequent
% for loop as well as identifying tetrode starting channels
[tsCountFl, ~, ~, ~] = plx_info(plxFile, 1);
tetChans = 1:nTrodes:(size(tsCountFl,2)-1);
tetNum = 0;
for chan = 1:size(tsCountFl,2)-1
    % Pull out waveform data from the recording. NOTE: The input of 0 means
    % pull the unsorted units so make sure you're only using an UNCUT file
    % in this analysis.
    %   numWFs = number of waveforms on that channel
    %   npw = number of points recorded per wave
    %   ts = timestamps associated with each waveform
    %   wave = waveform values in mV
    [numWFs, npw, ts, wave] = plx_waves_v(plxFile, chan, 0);
    if ismember(chan,tetChans)
        % Will only step here on the first channel per tetrode. This sets
        % up that tetrode's data
        % Identify the tetrode #
        tetNum = tetNum+1;
        % Copy mptRaw to create a new empty vector for filling in this
        % tetrode's data
        curTet = mptRaw;
        % Create a variable to extract the valley depth and index for the
        % plexon waveforms
        curValley = nan(length(ts),4);        
        curValleyInds = nan(length(ts),4);
        % Also create a variable in order to identify which channel on the
        % tetrode has the deepest valley
        curMinValley = nan(length(ts),1);
        curMinValleyNdx = nan(length(ts),1);
        tetChanNum = 1;
        disp(sprintf('----- Compiling Tetrode %02d -----',tetNum)); %#ok<DSPS>
    else
        % Every other channel in the tetrode just increases the tetrode
        % channel number 
        tetChanNum = tetChanNum + 1;
        % BUT, if the first channel had zero waveforms you need to make the
        % data structures that were otherwise made assuming there were no
        % waveforms
        if size(curValley,1)==1 && length(ts)>1
            curValley = nan(length(ts),4);
            curValleyInds = nan(length(ts),4);
            curMinValley = nan(length(ts),1);
            curMinValleyNdx = nan(length(ts),1);
        end
    end    
    % Now step through each waveform and place that waveform in it's
    % correct time spot
    for wf = 1:numWFs
        % Find the minimum value in the waveform and it's index
        [curValley(wf,tetChanNum), wfMin] = min(wave(wf,:));
        % Find the index for the timestamp associated with the waveform
        tsInd = round(ts(wf)*sampRate);
        % Fill in the waveform in it's correct position
        if (tsInd-preThresh)>1
            curTet(tetChanNum,(tsInd-preThresh):(tsInd+(npw-preThresh-1))) = wave(wf,:);
        else 
            % Only way it gets here is if the first spike was detected
            % prior to the # of samples specified by pre-threshold, 
            % so just consider the first waveform to be the start of the
            % recording.
            curTet(tetChanNum,1:npw) = wave(wf,:);
        end
        % Identify the waveform vally index
        curValleyInds(wf,tetChanNum) = (wfMin-preThresh) + round(ts(wf)*sampRate);
        
        if ismember(chan+1,tetChans)
            % Will only go here if it's the last channel in the tetrode.
            % This is where the deepest vally is identified.
            [curMinValley(wf),minTet] = min(curValley(wf,:));
            curMinValleyNdx(wf) = curValleyInds(wf,minTet);
        end
    end
    
    % Next, go through and connect all the waveforms and get rid of all the
    % nan values. 
    if ismember(chan+1,tetChans)
        % This is only done when you're on the final channel in a tetrode
        for wire = 1:size(curTet,1)
            % Identify when the transition from a waveform to a nan period
            % starts
            nanTrans = diff(isnan([0,curTet(wire,:),0]));
            nanBreakNdxs = [find(nanTrans==1)', find(nanTrans==-1)'];
            % Now step through and make sure our indexs are appropriate and
            % correct for the fact that our final value is identified as
            % outside our trial indices. Using a while here may not be the
            % best approach but it works.
            while nanBreakNdxs(end,2)>size(curTet,2)
                nanBreakNdxs(end,2) = nanBreakNdxs(end,2) - 1; 
            end
            % Now step through each nan period and connect the waveforms
            % with a linear interpolation between the end of one and the
            % start of the next.
            for nanBreak = 1:size(nanBreakNdxs,1)
                index1 = nanBreakNdxs(nanBreak,1);
                index2 = nanBreakNdxs(nanBreak,2);
                if nanBreak==1
                    startVal = 0;
                else
                    startVal = curTet(wire,index1-1);
                end
                if startVal == curTet(wire,index2)
                    curTet(wire,index1:index2) = ones(1,index2-index1+1)*startVal;
                elseif length(startVal:(curTet(wire,index2)-startVal)/(index2-index1):curTet(wire,index2))==(index2-index1)
                    curTet(wire,index1:index2) = [startVal, startVal:(curTet(wire,index2)-startVal)/(index2-index1):curTet(wire,index2)];
                else
                    curTet(wire,index1:index2) = startVal:(curTet(wire,index2)-startVal)/(index2-index1):curTet(wire,index2);
                end
            end
            finalNanChk = isnan(curTet(wire,:));
            curTet(wire,finalNanChk) = 0;
        end
        % Save it all!
        valleys = int32(curMinValleyNdx);
        writemda(valleys,sprintf('spike_times.tet%02d.mda',tetNum),'int32');
        writemda(curTet,sprintf('raw.tet%02d.mda',tetNum),'float32');
        disp(sprintf('::::: Tetrode %02d saved :::::',tetNum)); %#ok<DSPS>
        clear curTet curTS
    end
end

cd(origDir);