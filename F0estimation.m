% based on https://ccrma.stanford.edu/~jos/pasp/Fundamental_Frequency_Estimation.html

pause on;

[y, fs] = wavread('jobs.wav');

order = 10; % for LPC

window_size = 512;  
window = hann(window_size);
% hop_size = window_size - order;
hop_size = window_size / 2;


y_length = length(y);
y_start = 1;
y_end = length(y) - window_size;

pitchIndex = 1;
pitches = []; % SLOOOWW

for frame_start = y_start:hop_size:y_end
    
    t_frame = frame_start:frame_start+window_size-1; % times for this frame
    y_frame = y(t_frame.').*window; % signal for this frame

    Y_frame = fft(y_frame);
    Y_pos = Y_frame(1:ceil(end/2));
    Y_pos_magnitude = abs(Y_pos);
    Y_pos_magnitude_norm = Y_pos_magnitude./max(Y_pos_magnitude(:)); 
    
    % lpc (linear predictive coding)
    A = lpc(y_frame, order);
    B = [1];
    y_inv_frame = filter(A,B,y_frame);
    
    % y inverse filtered spectrum
    Y_inv_frame = fftshift(fft(y_inv_frame));
    Y_inv_frame_magnitude = abs(Y_inv_frame);
    Y_inv_pos_magnitude = Y_inv_frame_magnitude(1:ceil(end/2));
    Y_inv_pos_magnitude_norm = Y_inv_pos_magnitude./max(Y_inv_pos_magnitude(:));
    Y_inv_frame_magnitude_max = max(Y_inv_frame_magnitude(:));
    Y_inv_frame_magnitude_dB = 20*log10(Y_inv_frame_magnitude/Y_inv_frame_magnitude_max);
    Y_inv_frame_magnitude_dB_pos = Y_inv_frame_magnitude_dB(1:end/2+1);

    Y_inv_treated = Y_inv_pos_magnitude;
    ramp = ((length(Y_inv_treated):-1:1)*(1./length(Y_inv_treated))).';
    Y_inv_treated = Y_inv_treated .* ramp;
    Y_inv_treated = Y_inv_treated./max(Y_inv_treated(:));
    
    Y_pos_magnitude_norm_treated = Y_pos_magnitude_norm.*ramp;
    
    % peak finding
    K = 16; % max number of peaks to find
    
    findFrom = Y_pos_magnitude_norm;
    higherThanNeighbors = findFrom > circshift(findFrom, [1, 0])...
        & findFrom > circshift(findFrom, [-1, 0]);
    higherThanNeighbors = higherThanNeighbors.*findFrom;
    lowestFreq = 80;
    lowestBin = ceil(lowestFreq / (fs/window_size));
    higherThanNeighbors_allowed = higherThanNeighbors;
    higherThanNeighbors_allowed(1:lowestBin) = 0; % zero out lows
    
    [ascendingNeighbors,indicies] = sort(higherThanNeighbors_allowed);
    ascendingNeighbors = ascendingNeighbors(end-K:end);
    indicies = indicies(end-K:end);
 
    % parabolic fit peaks
    %  based on: https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html
    oldValues = findFrom;
    bins = 1:length(oldValues);
    
    alpha = circshift(oldValues, [1 0])';
    beta = oldValues';
    gamma = circshift(oldValues, [-1 0])';
    
    p = (1/2).*(alpha-gamma)./(alpha-2.*beta+gamma);
    k_beta = bins;
    k_star = k_beta+p;
    
    newPeaks = beta-(1/4).*(alpha-gamma).*p;
    
    newPeaks = newPeaks(indicies);
    newBins = k_star(indicies);
    
    % refinements on peak finding
    %   restrict frequency upper and lower bound
    %   parabolic fit peak interpolation
    
    % autocorrelation... not using this
%     lags = round(length(Y_frame_magnitude_pos)/8);
%     [Y_acf, Y_lags, Y_bounds] = autocorr(Y_frame_magnitude_pos, lags);
%     Y_acf_max = max(Y_acf(:));
%     exponent = 4;
%     Y_acf_pow = (Y_acf/Y_acf_max).^(1/exponent)*Y_acf_max;
    
    figure(1);
    subplot(2,1,1);
    viewExponent = 2; % just for visualizing
    plot(findFrom.^(1/viewExponent), 'Marker', '.');
    hold on;
    stem(higherThanNeighbors_allowed.^(1/viewExponent), ...
        'g', 'Marker', 'none');
    stem(indicies, ascendingNeighbors.^(1/viewExponent), ...
        'b');
    stem(newBins, newPeaks.^(1/viewExponent), ...
        'r', 'Marker', '*');
    hold off;
    xlabel('bins');
    ylabel('amplitude squared')
    xlim([1 window_size / 4]);
    ylim([0 1]);
    legend('spectrum','peaks','largest peaks', 'interpolated peaks');
    
    subplot(2,1,2);
    diff = sort(newBins') - circshift(sort(newBins'),[1 0]);
    theMean = mean(diff(2:end));
    theMedian = median(diff(2:end));
    X = 2:102;
    hist(diff(2:end), X);
    
    [histN, histX] = hist(diff(2:end), X); % get hist values
    [histNmax, maxIndex] = max(histN); % where is the max value
    maxBinCenter = histX(maxIndex);
    binWidth = 1;
    minSearch = maxBinCenter-binWidth;
    maxSearch = maxBinCenter+binWidth;
    
    indicesToAverage = (diff>minSearch)&(diff<maxSearch);
    theAverage = sum(diff.*indicesToAverage)/sum(indicesToAverage);
   
    hold on;
%     stem(theMean, 1, 'b', 'Marker', '*');
%     stem(theMedian, 1, 'g', 'Marker', '*');
    stem(theAverage, 2, 'r', 'Marker', '*','LineWidth', 4);
    hold off;
    xlim([0, 20]);
    xlabel('bins');
    ylabel('# of occurences');
    title('histogram - spacings between largest peaks');
    legend('bins', 'average of all items in max bin');
    
    pitches(pitchIndex) = mean(diff(2:end)) * fs/window_size;
    pitchIndex = pitchIndex+1;
    
    F0 = theAverage * fs/window_size;
%     F0 = theMedian * fs/window_size;
    display(num2str(F0));
%     display(num2str(sum(indicesToAverage)));
    
    pause(window_size / fs * 4);
    
    
end

figure(2)
semitones = log2(pitches./16.352)*12;

plot(semitones)