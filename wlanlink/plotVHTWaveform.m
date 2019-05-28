function plotVHTWaveform(x,cfgVHT,varargin)
%plotVHTWaveform Displays a VHT waveform with the fields highlighted

narginchk(2,3);
if nargin>2
    titlestr = varargin{1};
else
    titlestr = 'Waveform with Highlighted Fields';
end

sr = wlanSampleRate(cfgVHT);
txTime = (numel(x)/sr)*1e6; % Duration to plot in microseconds
tick = (1/sr)*1e6; % Microseconds per sample
hf = figure;
hp = plot(0:tick:txTime-tick,abs(x(1:round(txTime*sr*1e-6),:)));
cy = ylim(gca);
ylim(gca,[cy(1) cy(2)+1]);
ylabel('Magnitude (V)');
xlabel('Time (microseconds)');
title(titlestr);
vhtPlotFieldColored(cfgVHT,hf,hp,0); % 5 us offset for plot
end

% Overlays field names on a plot
function vhtPlotFieldColored(cfg,hf,hp,varargin)

narginchk(3,4);

offset = 0;
if nargin>3
    offset = varargin{1};
end

% Field durations in Microseconds
Fs = wlanSampleRate(cfg);
ind = wlanFieldIndices(cfg);
Tlstf = double((ind.LSTF(2)-ind.LSTF(1))+1)/Fs*1e6;
Tlltf = double((ind.LLTF(2)-ind.LLTF(1))+1)/Fs*1e6;
Tlsig = double((ind.LSIG(2)-ind.LSIG(1))+1)/Fs*1e6;
Tvhtsiga = double((ind.VHTSIGA(2)-ind.VHTSIGA(1))+1)/Fs*1e6;
Tvhtstf = double((ind.VHTSTF(2)-ind.VHTSTF(1))+1)/Fs*1e6;
Tvhtltf = double((ind.VHTLTF(2)-ind.VHTLTF(1))+1)/Fs*1e6;
Tvhtsigb = double((ind.VHTSIGB(2)-ind.VHTSIGB(1))+1)/Fs*1e6;
Tvhtdata = double((ind.VHTData(2)-ind.VHTData(1))+1)/Fs*1e6;
% Tagc = double((ind.DMGAGC(2)-ind.DMGAGC(1))+1)/Fs*1e6;
% Tagcsf  = double((ind.DMGAGCSubfields(1,2)-ind.DMGAGCSubfields(1,1))+1)/Fs*1e6;
% Ttrn = double((ind.DMGTRN(2)-ind.DMGTRN(1))+1)/Fs*1e6;
% Ttrnsf  = double((ind.DMGTRNSubfields(1,2)-ind.DMGTRNSubfields(1,1))+1)/Fs*1e6;
% N = cfg.TrainingLength/4;

figure(hf);
ax = gca;
hold(ax,'on');
legendTxt = [];
lh = [];

% Set first color to use for plotting
ListColors = colormap('colorcube');
colidx = 39; % Start within the color map

xlimits = xlim(ax);
Tcum = 0;
if offset<Tlstf
    % LSTF
    legendTxt = [legendTxt {'LSTF'}];   
    plotLine(Tlstf);
    if Tlstf>=xlimits(2)
        return
    end
end
Tcum = Tlstf;
if offset<(Tlstf+Tlltf)
    % LLTF
    legendTxt = [legendTxt {'LLTF'}];
    plotLine(Tlltf);
    if (Tcum+Tlltf)>=xlimits(2)
        return
    end
end
Tcum = double(ind.LLTF(2))/Fs*1e6;
if offset<(Tcum+Tlsig)
    % Header
    legendTxt = [legendTxt {'LSIG'}];
    plotLine(Tlsig);
    if (Tcum+Tlsig)>=xlimits(2)
        return
    end
end
Tcum = double(ind.LSIG(2))/Fs*1e6;
if offset<(Tcum+Tvhtsiga)
    % Data
    legendTxt = [legendTxt {'VHTSIGA'}];
    plotLine(Tvhtsiga);
    if (Tcum+Tvhtsiga)>=xlimits(2)
        return
    end
end
Tcum = double(ind.VHTSIGA(2))/Fs*1e6;
if offset<(Tcum+Tvhtstf)
    % Data
    legendTxt = [legendTxt {'VHTSTF'}];
    plotLine(Tvhtstf);
    if (Tcum+Tvhtstf)>=xlimits(2)
        return
    end
end
Tcum = double(ind.VHTSTF(2))/Fs*1e6;
if offset<(Tcum+Tvhtltf)
    % Data
    legendTxt = [legendTxt {'VHTLTF'}];
    plotLine(Tvhtltf);
    if (Tcum+Tvhtltf)>=xlimits(2)
        return
    end
end
Tcum = double(ind.VHTLTF(2))/Fs*1e6;
if offset<(Tcum+Tvhtsigb)
    % VHTSIGB
    legendTxt = [legendTxt {'VHTSIGB'}];
    plotLine(Tvhtsigb);
    if (Tcum+Tvhtsigb)>=xlimits(2)
        return
    end
end
Tcum = double(ind.VHTSIGB(2))/Fs*1e6;
if offset<(Tcum+Tvhtdata)
    % VHTData
    legendTxt = [legendTxt {'VHTData'}];
    plotLine(Tvhtdata);
    if (Tcum+Tvhtdata)>=xlimits(2)
        return
    end
end
Tcum = double(ind.VHTData(2))/Fs*1e6;
% if offset<(Tcum+Tagc)
%     % AGC
%     for i = 1:N*4
%         if (Tcum+i*Tagcsf)>=xlimits(2)
%             return
%         end
%         legendTxt = [legendTxt {['AGC-SF' num2str(i)]}]; %#ok<AGROW>
%         plotLine(Tagcsf);
%         Tcum = double(ind.DMGAGCSubfields(i,2))/Fs*1e6;
%     end
% end
% if offset<(Tcum+Ttrn)
%     % TRN
%     for i = 1:N*5
%         if mod(i-1,5)==0
%             % TRN-CE
%             if (Tcum+Tce)>=xlimits(2)
%                 return
%             end
%             legendTxt = [legendTxt {'TRN-CE'}]; %#ok<AGROW>
%             plotLine(Tce);
%             trnCEIdx = mod(i-1,5)+(i-1)/5+1;
%             Tcum = double(ind.DMGTRNCE(trnCEIdx,2))/Fs*1e6;
%         else
%             % TRN-SF
%             if (Tcum+Ttrnsf)>=xlimits(2)
%                 return
%             end
%             trnSFIdx = mod(i-1,5)+floor((i-1)/5)*4;
%             legendTxt = [legendTxt {['TRN-SF' num2str(trnSFIdx)]}]; %#ok<AGROW>
%             plotLine(Ttrnsf);
%             Tcum = double(ind.DMGTRNSubfields(trnSFIdx,2))/Fs*1e6;
%         end
%     end
% end
hold(ax,'off');
xlim([offset xlimits(2)]);
legend(lh,legendTxt,'location','best')
delete(hp) % Remove original waveform from plot

function plotLine(Tfield)
    % Get portion of waveform for current field and replot it with desired
    % color
    a = get(ax,'children');
    replotIdx = (a(end).XData>=Tcum)&(a(end).XData<(Tcum+Tfield));
    col = ListColors(mod(colidx,size(ListColors,1))+1,:); % Get color
    hl = plot(a(end).XData(replotIdx),a(end).YData(replotIdx),'Color',col);
    colidx = colidx+2; % Increment color index for next field
    lh = [lh hl]; % Store handle for legend
    
    % Plot field boundary line
    plot(ax,[Tcum+Tfield Tcum+Tfield],ylim(gca),'k:');
end

end