function [b_hat] = receive(r,plot_flag)
% [b_hat] = receive(r,plot_flag)
% Receiver program for part 1 of the project. The program should produce the
% received information bits given the samples of the received signal (sampled 
% at fs Hz.)
%
% Input:
%   r = samples of the received signal
%   plot_flag = flag for plotting [0: don't plot, 1: plot] 
%
% Output:
%   b_hat = vector containing the received information bits
%
% Rev. C (VT 2016)

%********** Begin program EDIT HERE

% Complete the code below:     

%1. filter with Rx filter (Matched filter)               

Ns = 2; % Length of the transmit pulse, defining how many samples represent each symbol
rolloffactor=0.5; % Roll-off factor for the pulse shaping filter
span=6; % Span of the filter in symbol durations
pulse = rcosdesign(rolloffactor,span,Ns(1));       % Designing a square root raised cosine filter
y = conv(r, pulse, 'same');        % Convolve the received signal with the matched filter

%2. Sample filter output

Rs=1000; % defining the rate at which symbols are transmitted
fs=Ns * Rs; % Calculate the sampling frequency as symbol rate times the number of samples per symbol
Ts = 1 / fs;
samplingPoints = 1:Ns:length(y);   % Points at which to sample the filtered signal
y_sampled = y(samplingPoints); %  Sample the filter output at these points


%3. Make decision on which symbol was transmitted

symbols = [2, 4]; % Define your symbol amplitudes
% % Calculate decision boundaries (midpoints between each pair of symbols)
boundaries = (symbols(1:end-1) + symbols(2:end)) / 2;
a_hat = zeros(size(y_sampled)); % Preallocate received symbols vector
for i = 1:length(y_sampled)
    % Determine which symbol the sampled value is closest to
    [~, index] = min(abs(symbols - y_sampled(i)));
    a_hat(i) = symbols(index);    % Assign the symbol to a_hat
end


%4. Convert symbols to bits

b_hat = zeros(size(a_hat));% Initialize the bit vector

% Define the decision boundary
decisionBoundary = 3; % Threshold between symbols to decide between bit 0 and bit 1


% Loop through received symbols and convert them to bits based on the decision boundary
for i = 1:length(a_hat)
    if a_hat(i) < decisionBoundary
        b_hat(i) = 0; % Symbol is closer to 2, so it's bit 0
    else
        b_hat(i) = 1; % Symbol is closer to 4, so it's bit 1
    end
end
%********** DON'T EDIT FROM HERE ON
% plot Rx signals
PlotSignals(plot_flag, 'Rx', r, y, y_sampled)
%********** End program
end