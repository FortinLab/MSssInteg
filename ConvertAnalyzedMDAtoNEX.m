function ConvertAnalyzedMDAtoNEX
%% ConvertAnalyzedMDAtoNEX
%   File for converting .mda data into a .nex file so that it can be opened
%   in Offline Sorter for spike verification.
%   Ugh... I really don't want to go through and comment this file more...
%   someone else do it, I'll buy you a candy bar. -GE
%   - 02/01/2018    Created by GE
%
%%
load('RosettaLC.mat');
origCD = cd;
%% Identify MDA and JSON files
% Load PARAMS file: the .json file created that contains the recordings
% parameters used to guide MountainSort
[jsonDataFile, jsonFilePath] = uigetfile('.json', 'Identify .json PARAMS file');
if jsonDataFile == 0
    disp('No PARAMS file selected');
    cd(origCD);
    return
else
    cd(jsonFilePath);
    text = fileread([jsonFilePath jsonDataFile]);
    params = jsondecode(text);
    fprintf('File %s Loaded\n', jsonDataFile);
end
% Identify Original .plx file used 
[origPLXfileName, origPLXfilePath] = uigetfile('.plx','Identify ORIGINAL .PLX File');
if origPLXfileName == 0
    disp('Original plx file not selected');
    return
else
    cd(origPLXfilePath);
    origPlxFile = [origPLXfilePath origPLXfileName];   
    [~,~,~,~,~,~,numSampsPreSpk,~,~,~,~,~,~] = plx_information(origPlxFile);
end
% Load the FIRINGS file: the MountainSort output of spike timestamps and
% cluster assignment
[tsDataFile, tsFilePath] = uigetfile('.mda', 'Identify .mda FIRINGS file');
if tsDataFile == 0
    disp('No FIRINGS file selected');
    cd(origCD);
    return
else
    mountainSpikes = readmda([tsFilePath tsDataFile]);
    fprintf('File %s Loaded\n', tsDataFile);   
    cd(tsFilePath);
end
% Load SPIKE TIMES file: the mda file created with the indexes of the
% original Plexon waveforms %%%Comment In If Used for Some Reason%%%
[tsPlxDataFile, tsPlxDataPath] = uigetfile('.mda', 'Identify .mda SPIKE TIMES (from PLX) file');
if tsPlxDataFile == 0
    disp('No SPIKE TIMES file selected');
    cd(origCD);
    return
else
    plxSpikes = readmda([tsPlxDataPath tsPlxDataFile]);
    fprintf('File %s Loaded\n', tsPlxDataFile);
    cd(tsPlxDataPath);
end
% Load RAW file: the original continuous file constructed from the Plexon
% waveforms
[rawDataFile, rawFilePath] = uigetfile('.mda', 'Identify .mda RAW file');
if rawDataFile == 0
    disp('No RAW file selected');
    cd(origCD);
    return
else
    raw = readmda([rawFilePath rawDataFile]);
    fprintf('File %s Loaded\n', rawDataFile);
    cd(rawFilePath);
end

rawDataFileNameSplit = strsplit(rawDataFile, '.');
tetName = rawDataFileNameSplit{2};
%%%%% Include file check here %%%%%
%% Create Empty Data Structures
% First identify the number of templates (units) that are present
templateNums = unique(mountainSpikes(3,:));
waveforms = cell(1,length(templateNums)+1); % +1 to include unsorted (see below)
spikeTimes = waveforms;

%% Identify waveforms excluded by MountainSort
msSpks = unique(mountainSpikes(2,:));
skpdSpks = plxSpikes;
for spk = 1:length(plxSpikes)
    if sum((msSpks>(plxSpikes(spk)-numSampsPreSpk)) & (msSpks<(plxSpikes(spk)+(params.clip_size-(numSampsPreSpk+1)))))>=1
        skpdSpks(spk) = nan;
    end
end
skpdSpks(isnan(skpdSpks)) = [];
skpdSpks(skpdSpks==0) = [];
if size(skpdSpks,2)>size(skpdSpks,1) % Dunno why it randomly transposes shit sometimes... but it does it here with some files.
    skpdSpks = skpdSpks';
end
skpdSpkWvfms = repmat({nan(1,params.clip_size)}, [length(skpdSpks), 4]);
for wvfm = 1:length(skpdSpks)
    curSkpdWvfmNdx = skpdSpks(wvfm);
    for chan = 1:4
        while length(raw)-(curSkpdWvfmNdx+(params.clip_size-(numSampsPreSpk+1)))<0
            curSkpdWvfmNdx = curSkpdWvfmNdx-1;
        end
        skpdSpkWvfms{wvfm,chan} = raw(chan,curSkpdWvfmNdx-numSampsPreSpk:curSkpdWvfmNdx+(params.clip_size-(numSampsPreSpk+1)));
    end
end
disp('Excluded waveforms ID''d');

if ~isempty(templateNums)
    %% Identify waveforms in duplicate clusters
    [counts,bins] = histcounts(mountainSpikes(2,:), msSpks);
    dupeWvfmsNdxs = bins(counts>1);
    dupeRemovalLog = zeros(1,size(mountainSpikes,2));
    dupeWaveforms = repmat({nan(1,params.clip_size)}, [length(dupeWvfmsNdxs), 4]);
    for dp = 1:length(dupeWvfmsNdxs)
        curDupeNdx = dupeWvfmsNdxs(dp);
        for chan = 1:4
            dupeWaveforms{dp,chan} = raw(chan,curDupeNdx-numSampsPreSpk:curDupeNdx+(params.clip_size-(numSampsPreSpk+1)));
        end
        tempDupeRemovalLog = mountainSpikes(2,:)==curDupeNdx;
        dupeRemovalLog = dupeRemovalLog + tempDupeRemovalLog;
    end
    disp('Duplicate waveforms ID''d');

    %% Combine Skipped waveforms and Duplicate waveforms into unsorted unit
    tempSpkTms = unique([skpdSpks/params.samplerate;(dupeWvfmsNdxs/params.samplerate)']);
    tempSpkTms = [tempSpkTms,(1:length(tempSpkTms))'];
    sortedSpkTms = sortrows(tempSpkTms);
    spikeTimes{1} = [sortedSpkTms(:,1)', nan];

    wvfmSortingKey = sortedSpkTms(:,2);
    tempWvfms = [skpdSpkWvfms; dupeWaveforms];
    waveforms{1} = [tempWvfms(wvfmSortingKey,:);...
        repmat({nan(1,params.clip_size)},[1,4]);];
    dupeRemovalLog = logical(dupeRemovalLog);
    mountainSpikes(:,dupeRemovalLog) = [];
        disp('Unsorted unit ID''d');

    %% Now identify waveforms by template identity
    for uni = 1:length(templateNums)
        curUni = templateNums(uni);
        curUniSpks = mountainSpikes(3,:)==curUni;
        curUniNdxs = mountainSpikes(2,curUniSpks);
        curUniSpkTimes = curUniNdxs/params.samplerate;
        tempUniWaveforms = repmat({nan(1,params.clip_size)}, [length(curUniSpkTimes), 4]);
        for chan = 1:4        
            for wvfmNum = 1:length(curUniNdxs)
                curNdx = curUniNdxs(wvfmNum);
                tempUniWaveforms{wvfmNum,chan} = raw(chan,curNdx-numSampsPreSpk:curNdx+(params.clip_size-(numSampsPreSpk+1)));
            end
        end
        waveforms{uni+1} = tempUniWaveforms;
        spikeTimes{uni+1} = curUniSpkTimes;
        fprintf('Unit %02d ID''d\n', uni);
    end
else
    spikeTimes{1} = [skpdSpks/params.samplerate; nan];
    waveforms{1} = [skpdSpkWvfms; repmat({nan(1,params.clip_size)},[1,4]);];
    disp('All waveforms unsorted... stupid MountainSort :P');
end
%% Identify the largest waveform
% By finding the largest waveform and duplicating it across channels
% ensures that Offline Sorter treats all the the channels as having the
% save voltage scaling/settings
if isnan(spikeTimes{1})
    disp('No spikes at all on this tet');
    return
end
uniMaxWvfmNdxs = cellfun(@(a)find(max(max(cellfun(@(b)max(b),a)))==cellfun(@(b)max(b),a),1,'first'),waveforms);
uniMinWvfmNdxs = cellfun(@(a)find(min(min(cellfun(@(b)min(b),a)))==cellfun(@(b)min(b),a),1,'first'),waveforms);
maxWvfms = nan(length(waveforms),params.clip_size);
minWvfms = nan(length(waveforms),params.clip_size);
for uni = 1:length(waveforms)
    maxWvfms(uni,:) = waveforms{uni}{uniMaxWvfmNdxs(uni)};
    minWvfms(uni,:) = waveforms{uni}{uniMinWvfmNdxs(uni)};
end
maxWvfmVal = max(maxWvfms(find(max(max(maxWvfms,[],2))==max(maxWvfms,[],2),1,'first'),:));
minWvfmVal = min(minWvfms(find(min(min(minWvfms,[],2))==min(minWvfms,[],2),1,'first'),:));
artMaxWvfm = minWvfmVal:(maxWvfmVal-minWvfmVal)/(params.clip_size-1):maxWvfmVal;

artMaxWvfmSpkTime = (max(cell2mat(spikeTimes))+1);

spikeTimes{1}(end) = artMaxWvfmSpkTime;
waveforms{1}(end,:) = repmat({artMaxWvfm}, [1,4]);
disp('Standard waveform added');

%% Create the .nex file
nexFile = nexCreateFileData(params.samplerate);
for uni = 1:length(waveforms)
    for chan = 1:4
        if uni==1            
            nexFile = nexAddNeuron(nexFile, spikeTimes{uni}', sprintf('%s_%d', tetName, chan),...
                chan, uni-1);
            nexFile = nexAddWaveform(nexFile, params.samplerate,...
                spikeTimes{uni}', cell2mat(waveforms{uni}(:,chan))',...
                sprintf('%s_wf_%d', tetName, chan), numSampsPreSpk/params.samplerate,...
                params.clip_size, chan, uni-1);
        else
            nexFile = nexAddNeuron(nexFile, spikeTimes{uni}', sprintf('%s_%s_%d', tetName, RosettaLC{uni-1}, chan),...
                chan, uni-1); %#ok<*USENS>
            nexFile = nexAddWaveform(nexFile, params.samplerate,...
                spikeTimes{uni}', cell2mat(waveforms{uni}(:,chan))',...
                sprintf('%s_%s_wf_%d', tetName, RosettaLC{uni-1}, chan), numSampsPreSpk/params.samplerate,...
                params.clip_size, chan, uni-1);
        end
    end
end
%%
writeNexFile(nexFile, [tetName '.nex']);
fprintf('%s.nex Saved\n',tetName);
cd(origCD);
