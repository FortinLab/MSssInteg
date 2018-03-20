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
%% Identify PLX file if not fed in
[fileName, filePath] = uigetfile('.plx','Identify .PLX File');
if fileName == 0
    disp('No file selected');
    return
end
plxFile = [filePath fileName];
origDir = cd;
cd(filePath);

%% Identify channels and waveforms
[~,~,sampRate,~,nTrodes,~,preThresh,~,~,~,~,ssnDur,~] = plx_information(plxFile);
truncLim = round(sampRate*ssnDur);
mptRaw = nan(4,truncLim);
[tsCountFl, ~, ~, ~] = plx_info(plxFile, 1);
tetChans = 1:nTrodes:(size(tsCountFl,2)-1);
tetNum = 0;
for chan = 1:size(tsCountFl,2)-1
    [numWFs, npw, ts, wave] = plx_waves_v(plxFile, chan, 0);
    if ismember(chan,tetChans)
        tetNum = tetNum+1;
        curTet = mptRaw;
        curValley = nan(length(ts),4);        
        curValleyInds = nan(length(ts),4);
        curMinValley = nan(length(ts),1);
        curMinValleyNdx = nan(length(ts),1);
        tetChanNum = 1;
        disp(sprintf('----- Compiling Tetrode %02d -----',tetNum)); %#ok<DSPS>
    else
        tetChanNum = tetChanNum + 1;
    end    
    for wf = 1:numWFs
        [curValley(wf,tetChanNum), wfMin] = min(wave(wf,:));
        tsInd = round(ts(wf)*sampRate);
        curTet(tetChanNum,(tsInd-8):(tsInd+(npw-preThresh-1))) = wave(wf,:);
        curValleyInds(wf,tetChanNum) = (wfMin-preThresh) + round(ts(wf)*sampRate);
        if ismember(chan+1,tetChans)
            [curMinValley(wf),minTet] = min(curValley(wf,:));
            curMinValleyNdx(wf) = curValleyInds(wf,minTet);
        end
    end
    
    if ismember(chan+1,tetChans)
        for wire = 1:size(curTet,1)
            nanTrans = diff(isnan([0,curTet(wire,:),0]));
            nanBreakNdxs = [find(nanTrans==1)', find(nanTrans==-1)'];
            while nanBreakNdxs(end,2)>size(curTet,2)
                nanBreakNdxs(end,2) = nanBreakNdxs(end,2) - 1; %because the above method gives us our final value as outside 
            end
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
        valleys = int32(curMinValleyNdx);
        writemda(valleys,sprintf('spike_times.tet%02d.mda',tetNum),'int32');
        writemda(curTet,sprintf('raw.tet%02d.mda',tetNum),'float32');
        disp(sprintf('::::: Tetrode %02d saved :::::',tetNum)); %#ok<DSPS>
        clear curTet curTS
    end
end

cd(origDir);