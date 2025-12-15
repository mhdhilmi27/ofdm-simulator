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
    

    % a) Construct the matrix Qt and the matrix Qt'; Pages 20-21 OFDM
    % lecture notes.
    % Represent both of these as sparse matrices to save memory. 
    % ...
    
    
    
    % b) Create the vector g of length K, where the k'th element is gk from the OFDM compendium.
    % ...
  
    
    
    % c) using g and fDelta, create f, which is the vector of subcarrier
    % frequencies. f(k) is equal to fk in the compendium. 
    % ...
    
    
    
    % to save computations, a) b) and c) could be done outside this loop.
   
    
    
    
    % d) sample in time domain by performing IDFT on N*Qt*a; (use ifft(...))
    % ...

   

    % e) Add the cyclic prefix
    % ...



    % f) Now we simulate an analogue signal with an upsampled discrete signal
    % that is interpolated using for example the sinc-function. See f) in 
    % the problem description.
    % ... 
    


    % % %% SANITY CHECK. See if we get back what we should.
    % Are the samples of xAnalog correct at the N+L places or at least
    % very close? Does the signal look smooth and is only defined during
    % the interval [0,Ts] etc. You may plot things to do this
    % investigation. 
    % ...



    % g) Modulate to bandpass. Simply multiply xAnalog exp(1i*2*pi*fc*t)
    % (elementwise!) and then take the real part of that signal.



    % Power allocation... Skipped in this very lab.

    %% Channel

    % h) We produce z by convolving the bandpass signal with the channel.
    % See the h) in the problem description.
    % ...



    %% Receiver
    % Now it's time for the receiver. 

    % i) We ought to frequency shift the stuff back to baseband.
    % This is simply done by elementwise multiplication of (the real)
    % signal z with exp(-1i*2*pi*fc*t) to create a signal zFS. 
    % Then lowpass-filter zFS to create zLP. See i) in the probplem
    % description. 
    % ... 
    


    % j) Sample zLP (and remove cyclic prefix).
    % Simply pick the N samples during the observation
    % interval. Those should be at the indexes sampleIndexes(L+1:L+N). 
    % ...



    % k)  Run DFT using fft( ) for the vector produced in j). 
    


    % l) Channel Hk's. We assume we know the channel inpulse response. (or
    % has been estimated with perfect accuracy). Use the definition of the
    % channel impulse response and the vector f produced in c) to compute
    % the vector H, where H(k) is the Fourier transform of the
    % channel impulse response at frequency f(k) + fc. 
    % ...
  

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




