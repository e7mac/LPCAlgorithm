pause on;

[y, fs, nbits] = wavread('jobs.wav');
y = y(:,1); % mono, left channel only
y = y/max(abs(y(:)))*1.0; %normalize
% y = resample(y, 441, 80); % resample to 44100
% fs = fs*441/80; % new sampling rate

y_full = zeros(size(y)); % pre-allocate zeros
y_inv_full = zeros(size(y));
residue_full = zeros(size(y));
glottal_full = zeros(size(y));

lengthOfThroat = 5.5; %inches
speedOfSoundInches = 13397; % inches per second
samplesPerInch = fs / speedOfSoundInches;
order = round(samplesPerInch*lengthOfThroat);

order = 30; % lpc order
window_size = 256;
% hop_size = window_size-order;
hop_size = window_size/2;

y_length = length(y);
y_start = 1;
y_end = length(y) - window_size;

theWindow = hann(window_size).^0.5; % sum of squared windows overlappoing equals 1

for loops = 1:1
    for frame_start = y_start:hop_size:y_end
        
        t_frame = frame_start:frame_start+window_size-1; % times for this frame
        t_frame_sec = t_frame/fs; % times for this frame
        y_frame = y(t_frame).*theWindow; % signal for this frame
        
        A = lpc(y_frame, order);
        B = [1];
        y_inv_frame = filter(A,B,y_frame);
        residue = y_frame - y_inv_frame; % the residue
        glottal = filter(B, A, residue); % actual glottal pulse?
        
        % write into files
        y_full(t_frame) = y_full(t_frame) + y_frame.*theWindow;
        y_inv_full(t_frame) = y_inv_full(t_frame) + y_inv_frame.*theWindow;
        residue_full(t_frame) = residue_full(t_frame) + residue.*theWindow;
        glottal_full(t_frame) = glottal_full(t_frame) + glottal.*theWindow;
        
        bins = 0:length(y_frame)/2;
        
        % y spectrum
        Y_frame = fft(y_frame);
        Y_frame_magnitude = abs(Y_frame);
        Y_frame_magnitude_max = max(max(Y_frame_magnitude(:), 1));
        Y_frame_magnitude_dB = 20*log10(Y_frame_magnitude/Y_frame_magnitude_max);
        Y_frame_magnitude_dB_pos = Y_frame_magnitude_dB(1:end/2+1);
        
        % y inverse filtered spectrum
        Y_inv_frame = fft(y_inv_frame);
        Y_inv_frame_magnitude = abs(Y_inv_frame);
        Y_inv_frame_magnitude_max = max(max(Y_inv_frame_magnitude(:), 1));
        Y_inv_frame_magnitude_dB = 20*log10(Y_inv_frame_magnitude/Y_inv_frame_magnitude_max);
        Y_inv_frame_magnitude_dB_pos = Y_inv_frame_magnitude_dB(1:end/2+1);
        
        % residue spectrum
        RESIDUE = fft(residue);
        RESIDUE_magnitude = abs(RESIDUE);
        RESIDUE_magnitude_max = max(max(RESIDUE_magnitude(:), 1));
        RESIDUE_magnitude_dB = 20*log10(RESIDUE_magnitude/RESIDUE_magnitude_max);
        RESIDUE_magnitude_dB_pos = RESIDUE_magnitude_dB(1:end/2+1);
        
        % actual glottal pulse?
        GLOTTAL = fft(glottal);
        GLOTTAL_magnitude = abs(GLOTTAL);
        GLOTTAL_magnitude_max = max(max(GLOTTAL_magnitude(:), 1));
        GLOTTAL_magnitude_dB = 20*log10(GLOTTAL_magnitude/GLOTTAL_magnitude_max);
        GLOTTAL_magnitude_dB_pos = GLOTTAL_magnitude_dB(1:end/2+1);
        
        figure(1);
        subplot(4,4,1:4);
        t = (1:length(y))/fs; % time values
        plot(t, y);
        title('Signal');
        xlabel('time (seconds)');
        ylabel('amplitude');
        axis tight;
        ylim([-1, 1]);

        rectangle_bounds = [frame_start/fs, -1, (hop_size-1)/fs, 2];
        rectangle('Position', rectangle_bounds, 'EdgeColor', 'r');
        
        % current frame
        subplot(4,4,5);
        plot(t_frame_sec, y_frame);
        title('Current Frame');
        xlabel('time (seconds)');
        ylabel('amplitude');
        axis tight;
        ylim([-1, 1]);
        
        subplot(4,4,9);
        plot(bins, Y_frame_magnitude_dB_pos);
        title('Current Frame Spectrum');
        xlabel('bin');
        ylabel('amplitude(dB)');
        axis tight;
        ylim([-60 0]);
        
        % current frame, inverse filtered
        subplot(4,4,6);
        plot(t_frame_sec, y_inv_frame);
        title('w/o Formant');
        xlabel('time (seconds)');
        ylabel('amplitude');
        axis tight;
        ylim([-1, 1]);
        
        subplot(4,4,10);
        plot(bins, Y_inv_frame_magnitude_dB_pos);
        title('Spectrum, w/o Formant');
        xlabel('bin');
        ylabel('amplitude(dB)');
        axis tight;
        ylim([-60 0]);
        
        % current frame, residue
        subplot(4,4,7);
        plot(t_frame_sec, residue);
        title('Residue');
        xlabel('time (seconds)');
        ylabel('amplitude');
        axis tight;
        ylim([-1, 1]);
        
        subplot(4,4,11);
        plot(bins, RESIDUE_magnitude_dB_pos);
        title('Spectrum, Residue');
        xlabel('bin');
        ylabel('amplitude(dB)');
        axis tight;
        ylim([-60 0]);
        
        % current frame, actual glottal?
        subplot(4,4,8);
        plot(t_frame_sec, glottal);
        title('Glottal Pulse estimate?');
        xlabel('time (seconds)');
        ylabel('amplitude');
        axis tight;
        ylim([-1, 1]);
        
        subplot(4,4,12);
        plot(bins, GLOTTAL_magnitude_dB_pos);
        title('Spectrum, Glottal Pulse estimate?');
        xlabel('bin');
        ylabel('amplitude(dB)');
        axis tight;
        ylim([-60 0]);
        
        % coefficients
        subplot(4,4,13:16);
        A_plot = (A + 1) /2; % process A coefficients for plotting
        A_plot = (A_plot < 0.0)*0.0 + (A_plot >= 0.0).*(A_plot); %hacky?
        coeffNums = 1:length(A_plot);
        plot(coeffNums, A_plot, coeffNums, -A_plot);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        xlim([min(coeffNums)-1, max(coeffNums)+1]);
        ylim([-1, 1]);

%         pause(hop_size/fs); % wait
%         drawnow; % flush
%         sound(y_frame, fs); % let's hear

    end
end

% normalize
y_full       = 0.9 * y_full./max(y_full(:));
y_inv_full   = 0.9 * y_inv_full./max(y_inv_full(:));
residue_full = 0.9 * residue_full./max(residue_full(:));
glottal_full = 0.9 * glottal_full./max(glottal_full(:));

wavwrite(y_full, fs, 'jobs_reconstructed.wav');
wavwrite(y_inv_full, fs, 'jobs_inv.wav');
wavwrite(residue_full, fs, 'jobs_residue.wav');
wavwrite(glottal_full, fs, 'jobs_glottal.wav');

