% Muhammad Hilmi

%% System parameters
K = 64; % number of subcarriers
N = 64; % number of samples equal to K.
L = floor(0.1*K); % cyclic prefix samples "approx 10%"
Ts = 0.001; % signaling time
Tobs = (N*Ts)/(N+L); % observation time
Tcp = Ts - Tobs; % length of cyclic prefix
fDelta = 1/Tobs; % separation of carriers
fsamp = N*fDelta; % Choosen
% fsampMin = K*fDelta; % the lowest rate we can sample according to theory. Not used. Here same as fsamp. 
sampleTimes = 0:Ts/(N+L):Ts - Ts/(N+L); % times where the N samples are found (including cyclic prefix).
N_analog = 1000*(N+L); % number of samples per Ts for the simulated analog signal.
sampleIndexes = 1:floor(N_analog/(L+N)):(L+N)*floor(N_analog/(L+N)); % sample times
% indexes for the analog signal where the N samples are.





fc = (100)*N_analog; % Carrier frequency

% Channel model: sum_k channelGains(k) \delta(t - channelTaus(k))

Nchannel = 4; % number of taps or time delays in the channel. 
% Precisely: length(channelTaus).

channelGains = [0.9, 0.3 0.03 0.02];
channelTaus = [0 0.02 0.3 0.33]*Tcp; % We have a very generous choice of Tcp,
% we see that Th is approximately (1/3)Tcp. But let's be safe and we chose
% Tcp approximately 10% of Tobs. 

channelSampleTimes = floor(channelTaus*(N_analog/Ts)); % the sample times of 
% the analog signal correspoinding to channelTaus.

SamplesDuringTcp = floor(Tcp*(N_analog/Ts)); % simply the amount of samples
% of the analog signal during a time interval of length Tcp.




%% We produce a bunch of subcarrierers (K such to be precise)
% Culd be simplified into one struct, since k and M is the the same for all
% carriers. 
subcarrier = cell(K,1);

for k = 1:K
    % 64-QAM for all of them
    carrier.k = 6;
    carrier.M = 2^carrier.k; %sqrt(M) must be integer!
    carrier.sqrtM = sqrt(carrier.M);
    
    
    % Constellation same for all, Dmin = 2.
    carrier.Ichoices = -carrier.sqrtM+1:2:carrier.sqrtM-1; % x-coordinates
    carrier.Qchoices = -carrier.sqrtM+1:2:carrier.sqrtM-1; % y-coordinates
    
    % Constellation matrix
    carrier.C = ones(carrier.sqrtM,1)*carrier.Ichoices + ...
        carrier.Qchoices'*ones(1,carrier.sqrtM)*1i;
    
    % Priors
    % carrier.p = rand(1,carrier.M); random prior
    % carrier.p = (1:carrier.M); % some specific choice of prior
    carrier.p = ones(1,carrier.M); % flat prior
    carrier.p = carrier.p/norm(carrier.p,1); % sum over the prior should be 1.
    carrier.P = reshape(carrier.p, [carrier.sqrtM, carrier.sqrtM]); % Prior in matrix form;
    
    
    % Cumulative distribution 
    carrier.pCum = carrier.p*tril(ones(carrier.M, carrier.M))';

    % Shifted cumulative distribution 
    carrier.pCumShift = [0, carrier.pCum(1:carrier.M-1)];

    % save to array of subcarriers
    subcarrier{k,1} = carrier;
end


%% We start the main loop

% number of loops (messages to be sent for each carrier)
Nmes = 1; % change to higher number later when one iteration of the loop works!

% These guys are used in the symbol detection later. 
NrSymbolErrorsML = 0;
NrSymbolErrorsMAP = 0;


for iterations = 1:Nmes

    %% Transmitter
    
    % choose messages and symbols
    messages = zeros(K,1);
    a = zeros(K,1);
    
    % create the random symbol to be sent for each carrier and put it in 
    % the vector a
    %for k = 1:K
    for k = 1:K
        % Draw a message from the prior
        M = subcarrier{k,1}.M;
        probVec = rand*ones(1,M);
        messageChoosen = find((ones(1,M) - sign((probVec - ...
            subcarrier{k,1}.pCum).*(probVec - subcarrier{k,1}.pCumShift))));
        messages(k) = messageChoosen;
        a(k) = subcarrier{k,1}.C(messageChoosen);
    end
    % the a-vector contains the symbols to be sent. 
    

    % b) Create the vector g of length K, where the k'th element is gk from the OFDM compendium.
    % ...
    if mod(K,2) == 1        %odd
        krc = (K-1)/2;
    elseif mod(K,2) == 0    %even
        krc = (K-2)/2; 
    end

    k = 0:K-1;
    g = k - krc;

    % a) Construct the matrix Qt and the matrix Qt'; Pages 20-21 OFDM
    % lecture notes.
    % Represent both of these as sparse matrices to save memory. 
    % ...
    % g_bin is length-K, values in 0..N-1

    g_bin = mod(g, N);

    g0   = g(1);      
    gK_1 = g(end);

    % range for every Xm
    start_m_1 = 0;
    stop_m_1 = gK_1;
    start_m_2 = gK_1 + 1; 
    stop_m_2 = g0 + N - 1;
    start_m_3 = g0 + N;
    stop_m_3 = N - 1;   

    m_first_range = start_m_1:stop_m_1;
    a_first_range = m_first_range-g0;

    m_third_range = start_m_3:stop_m_3;
    a_third_range = m_third_range-(g0+N);
    
    if start_m_2 <= stop_m_2
        len2 = stop_m_2 - start_m_2 + 1;
        a_second_range = -1 * ones(1, len2); 
        a_idx = [a_first_range, a_second_range, a_third_range];
        
    elseif start_m_2 > stop_m_2
        a_idx = [a_first_range, a_third_range];
    end

    a_idx = a_idx.';

    % construct Qt from a_idx
    valid = (a_idx >= 0) & (a_idx < K);
    
    row_idx = find(valid);      
    col_idx = a_idx(valid)+1;
    vals    = ones(length(row_idx), 1);
    
    Qt = sparse(row_idx, col_idx, vals, N, K);
    Qr = Qt.';

    % c) using g and fDelta, create f, which is the vector of subcarrier
    % frequencies. f(k) is equal to fk in the compendium. 
    % ...
    f = g_bin(:) * fDelta; 
    f_new = g(:) * fDelta;
    
    % to save computations, a) b) and c) could be done outside this loop.
    
    
    % d) sample in time domain by performing IDFT on N*Qt*a; (use ifft(...))
    % ...
    X = N * Qt * a;
    x = ifft(X);
   

    % e) Add the cyclic prefix
    % ...
    x_cp = [x(end-L+1:end); x];

    % f) Now we simulate an analogue signal with an upsampled discrete signal
    % that is interpolated using for example the sinc-function. See f) in 
    % the problem description.
    % ... 
    
    Tsamp_orig   = Ts / (L + N);          % Ts/(L+N)
    fsamp_orig   = 1 / Tsamp_orig;        % original sampling rate
    t_samples    = (0:L+N-1).' * Tsamp_orig;  % t_m = m * Ts/(L+N)
    
    fsampAnalog  = N_analog / Ts;         % dense "analog" sampling rate
    tAnalog      = (0:N_analog-1).' / fsampAnalog;   % 0 ... Ts
    
    xAnalog      = zeros(N_analog, 1);    % s_I(t) in Eq. (3.8)
    
    N_interpolation = 8;                  % truncate sinc to ±8 original samples
    
    for m = 1:(L+N)
        t_m = t_samples(m);               % time of sample m
    
        % only use tAnalog within ± N_interpolation * Tsamp_orig
        t_left  = t_m - N_interpolation * Tsamp_orig;
        t_right = t_m + N_interpolation * Tsamp_orig;
    
        idx = (tAnalog >= t_left) & (tAnalog <= t_right);
        if ~any(idx)
            continue;
        end
    
        tau = tAnalog(idx) - t_m;         % t - m Ts/(L+N)
    
        % g_i(t) = sinc interpolation filter
        xAnalog(idx) = xAnalog(idx) + x_cp(m) * sinc(fsamp_orig * tau);
    end


    % % %% SANITY CHECK. See if we get back what we should.
    % Are the samples of xAnalog correct at the N+L places or at least
    % very close? Does the signal look smooth and is only defined during
    % the interval [0,Ts] etc. You may plot things to do this
    % investigation. 
    % ...

    % Re-sample xAnalog at the original sampling instants using sampleIndexes
    x_rec = xAnalog(sampleIndexes);      % length N+L, should match x_cp
    
    % Interpolation error
    err = x_rec - x_cp;
    maxErr = max(abs(err));
    rmsErr = sqrt(mean(abs(err).^2));
    
    fprintf('Max interpolation error: %g\n', maxErr);
    fprintf('RMS interpolation error: %g\n', rmsErr);


    % g) Modulate to bandpass. Simply multiply xAnalog exp(1i*2*pi*fc*t)
    % (elementwise!) and then take the real part of that signal.
    %tAnalog = linspace(0, Ts, N_analog).'; 
    xBP_complex = xAnalog .*exp(1i*2*pi*fc.*tAnalog);
    xBP         = real(xBP_complex);

    

    % Power allocation... Skipped in this very lab.

    %% Channel

    % h) We produce z by convolving the bandpass signal with the channel.
    % See the h) in the problem description.
    % ...
    
    % dicreate-time channel
    % h_BP(t) = sum (ap * dirac(t - tp)
    % h_BP(n) = sum (ap * dirac(n - np)
    h_BP = zeros(N_analog,1);
    for i = 1:Nchannel
        np = channelSampleTimes(i) + 1;
        h_BP(np) = h_BP(np) + channelGains(i);
    end
    
    count_nonzero = nnz(h_BP);

    % convolution xBP and h_BP
    z = conv(xBP, h_BP, 'full');
    z = z(1:N_analog);

    % Time axis for z
    Nz = length(z);
    t_z = (0:Nz-1).' / fsampAnalog;

    %% Receiver
    % Now it's time for the receiver. 

    % i) We ought to frequency shift the stuff back to baseband.
    % This is simply done by elementwise multiplication of (the real)
    % signal z with exp(-1i*2*pi*fc*t) to create a signal zFS. 
    % Then lowpass-filter zFS to create zLP. See i) in the probplem
    % description. 
    % ... 

    % frequency shift to baseband
    Nz = length(z);
    zFS = z .* exp(-1i*2*pi*fc.*tAnalog(1:Nz));

    % Sinc LPF
    NLP = 200;
    B = (K/2)*fDelta;
    t_sinc = (-NLP:NLP)/fsampAnalog;
    sinctrunc = (2*B/fsampAnalog) * sinc(2*B*t_sinc);
    sinctrunc = sinctrunc / sum(sinctrunc);

    % convolution z with sinc LPF
    zLP = conv(zFS, sinctrunc, 'full');
    zLP = zLP(NLP+1:NLP+Nz);
    zLP = 2 * zLP;

    % 1) Spectrum of bandpass received signal z(t)
    plot_fft(z, fsampAnalog, 'Spectrum of z(t) - bandpass after channel');
    
    % 2) Spectrum after mixing down: zFS(t) = z * exp(-j2πf_c t)
    plot_fft(zFS, fsampAnalog, 'Spectrum of zFS(t) - baseband + image at 2f_c');
    
    % 3) Spectrum after lowpass: zLP(t)
    plot_fft(zLP, fsampAnalog, 'Spectrum of zLP(t) - lowpass baseband');

    % Checking matrixsize
    %length_xBP = length(xBP)
    %length_hbp = length(h_bp)
    %length_z   = length(z)
    %length_t   = length(tAnalog)
    %N_analog
    %[idx, vals] = find(h_bp);
    %[idx vals]

    % j) Sample zLP (and remove cyclic prefix).
    % Simply pick the N samples during the observation
    % interval. Those should be at the indexes sampleIndexes(L+1:L+N). 
    % ...
    z_disc_cp = zLP(sampleIndexes);   % length = N + L
    y = z_disc_cp(L+1 : L+N);

    % k)  Run DFT using fft( ) for the vector produced in j). 
    R = fft(y); 

    % l) Channel Hk's. We assume we know the channel inpulse response. (or
    % has been estimated with perfect accuracy). Use the definition of the
    % channel impulse response and the vector f produced in c) to compute
    % the vector H, where H(k) is the Fourier transform of the
    % channel impulse response at frequency f(k) + fc. 
    % ...

    % Channel H_k: Fourier transform of h(t) at f(k) + fc
    H = zeros(K,1);   % H(k) for each active subcarrier
    
    for ell = 1:Nchannel
        % α_ℓ * exp( -j 2π (f(k) + fc) T_ℓ ), applied for all k at once
        H = H + channelGains(ell) * exp( -1i * 2*pi * (f_new + fc) * channelTaus(ell) );
    end
    
    % Now H is a K×1 vector, same ordering as f and as Qr*R

  

    % % %% SANITY CHECK. See if we get back what we should.
    % Here we should be careful before we proceed. If we define 
    % aTilde = (1/(N)*Qr*R)./H, we may then compute norm(a - aTilde). If 
    % norm(a - aTilde) is smaller than 1 = D_min/2, we will at least 
    % have no symbol errors  (sufficient
    % condition) due to poor filtering and numeric integration. 
    % Try to run som itererations of the loop at this point and
    % make sure the scaling of the preivous filtering is correct. 
    % 
    % A good starting point is to make sure this works without the channel.
    % So you might skip the channel and go directly from g) to i) and make
    % sure in this case that if aTilde = (1/(N)*Qr*R), then norm(a -
    % aTilde) is small (less than 1 preferably, but otherwise close to that). 
    % If it works without the channel it should work well with the channel
    % (and also changing the channel). 

    aTilde = (1/(N)*Qr*R)./H;
    norm(a - aTilde)

    figure;
    plot(abs(a - aTilde))

    X_dig    = N * Qt * a;
    x_dig    = ifft(X_dig);
    xcp_dig  = [x_dig(end-L+1:end); x_dig];
    y_dig    = xcp_dig(L+1:L+N);
    R_dig    = fft(y_dig);
    
    aTilde_d = (1/N) * (Qr * R_dig);
    norm_dig = norm(a - aTilde_d);
    fprintf('Pure-digital ||a - aTilde_d|| = %.3g\n', norm_dig);

    
    % m) We now introduce additional noise (besides imperfect integration and filtering).  
    % This should capture all sorts of additional noise sources.  N0 = 2000;
    % Create Rnoisy from R by adding complex Gaussian noise with variance N0 to R.
    % ...

    

    % n)  Symbol detection will now follow.
    % 
    %%% * Do the ML and increase NrSymbolErrorsML 
    % 
    %%% * Do the MAP and increase NrSymbolErrorsMAP
    % ...


    % Congratulations, you are done with one iteration!

end


% visualization check
visualization(X,x, N, fDelta, L, Tobs, x_cp, ...
    tAnalog, xAnalog, sampleTimes, Ts, xBP, t_z, z, channelTaus)


function plot_fft(sig, fs, ttl)
%PLOT_FFT  Plot the magnitude spectrum of a signal in dB
%
%   plot_fft(sig, fs, ttl)
%
%   sig : input time-domain signal (vector)
%   fs  : sampling frequency in Hz
%   ttl : plot title (string)

    sig = sig(:);                 % ensure column
    N = length(sig);              % signal length

    SIG = fftshift(fft(sig));     % centered FFT
    faxis = linspace(-fs/2, fs/2, N);   % frequency axis

    figure;
    plot(faxis, 20*log10(abs(SIG) + 1e-12), 'LineWidth', 1.2);
    grid on;

    xlabel('Frequency [Hz]');
    ylabel('Magnitude [dB]');
    title(ttl);

end




