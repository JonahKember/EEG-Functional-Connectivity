function [adjacency_tensor, xs, filt_params] = across_trial_pli(data,samp_rate,filt_band)
%%
% Inputs:
%   data             = [Channels x Samples x Trials] Time-locked EEG data.
%   samp_rate        = Sampling rate of data, in Hz.
%   filt_band        = [Min, Max] Band-pass filter range.  
%
% Output:
%   adjacency_tensor = [Channel x Channel x Sample], PLI adjacency matrices over time. 
%   xs               = Structure with time-series used in analyses (raw, filtered, phase, hilbert-transformed).
%   filt_params      = Structure specifying filter parameters (order, transition-width, coefficients, normalized frequencies).
%%

tic
fprintf('\nCalculating across-trial phase-lag index...\n')

% Define variables.
n_channels = size(data,1);
n_samples = size(data,2);
n_trials = size(data,3);

xs = [];
xs.raw = data;

% Estimate filter parameters.
filt_params = [];
filt_params.transition = mean(filt_band) * 0.2;                                                                   
filt_params.frequencies = [filt_band(1) - filt_params.transition, filt_band(1), filt_band(2), filt_band(2) + filt_params.transition];
filt_params.order = kaiserord(filt_params.frequencies, [0 1 0], [0.1 0.05 0.1], samp_rate);                                                                          
filt_params.coefficients = fir1(filt_params.order, filt_band*(2/samp_rate), 'bandpass');

% Z-score each channel time-series, then calculate the instantaneous phase time-series using the Hilbert transform.
for trial = 1:n_trials
    xs.norm(:,:,trial) = zscore(xs.raw(:,:,trial)')';
    for channel = 1:n_channels
        xs.filtered(channel,:,trial)  = filtfilt(filt_params.coefficients, 1, squeeze(xs.norm(channel,:,trial)));
    end
    xs.hilbert(:,:,trial) = hilbert(xs.filtered(:,:,trial));
    xs.phase(:,:,trial) = angle(xs.hilbert(:,:,trial));
end

% Define a function to calculate the phase-lag index. 
pli_fx = @(phase_set_i,phase_set_j) abs(mean(sign(phase_set_i - phase_set_j),2));

% Calculate channel-by-channel PLI values to create adjacency matrices. 
adjacency_tensor = zeros(n_channels,n_channels,n_samples);
channel_pairs =  nchoosek(1:n_channels,2);

for ij = 1:length(channel_pairs)
    pair = channel_pairs(ij,:);
    pli = pli_fx(squeeze(xs.phase(pair(1),:,:)),squeeze(xs.phase(pair(2),:,:)));
    adjacency_tensor(pair(1),pair(2),:) = pli;
end

% Fill in symmetric values.
for sample = 1:n_samples
    adjacency_tensor(:,:,sample) = adjacency_tensor(:,:,sample)' + adjacency_tensor(:,:,sample);
end
toc

end
