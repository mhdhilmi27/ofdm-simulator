function visualization(X,x, N, fDelta, L, Tobs, x_cp, ...
    tAnalog, xAnalog, sampleTimes, Ts, xBP, t_z, z, channelTaus)

    %% X_m vs m
    m = 0:N-1;                 
    figure;
    stem(m, abs(X), 'filled');
    xlabel('FFT bin index m');
    ylabel('|X_m|');
    title('Frequency-domain OFDM symbol  |X_m| vs m');
    grid on;
    

    %% X_m vs bin frequency f_m
    f_bins = m.' * fDelta;    
    figure;
    stem(f_bins, abs(X), 'filled');
    xlabel('Bin frequency f_m [Hz]');
    ylabel('|X_m|');
    title('Frequency-domain OFDM symbol  |X_m| vs f_m');
    grid on;
    

    %% x_n
    % Time axes
    n_x   = 0 : N-1;                 
    n_cp  = -L : N-1;                
    
    Tsamp = Tobs / N;                
    t_x   = n_x.'  * Tsamp;          
    t_cp  = n_cp.' * Tsamp;          
    
    % x_n
    subplot(2,1,1);
    stem(t_x, real(x), 'filled', 'MarkerSize', 3); hold on;
    stem(t_x, imag(x), 'r.', 'MarkerSize', 8);
    xlabel('Time [s]');
    ylabel('x[n]');
    title('OFDM symbol WITHOUT CP');
    legend('Re\{x\}','Im\{x\}');
    grid on;
    xlim([min(t_cp) max(t_cp)]);   
    
    % x_cp
    subplot(2,1,2);
    stem(t_cp, real(x_cp), 'filled', 'MarkerSize', 3); hold on;
    stem(t_cp, imag(x_cp), 'r.', 'MarkerSize', 8);
    xlabel('Time [s]');
    ylabel('x_{cp}[n]');
    title('OFDM symbol WITH CP');
    legend('Re\{x_{cp}\}','Im\{x_{cp}\}');
    grid on;
    xlim([min(t_cp) max(t_cp)]);   

    %% analog signal
    % Plot real part: analog vs original discrete samples
    figure;
    subplot(2,1,1);
    plot(tAnalog, real(xAnalog), 'LineWidth', 1); hold on;
    stem(sampleTimes, real(x_cp), 'r', 'filled');
    xlabel('Time [s]');
    ylabel('Re\{x\}');
    title('Real part: xAnalog vs original x_{cp} samples');
    legend('xAnalog (oversampled)', 'x_{cp} samples');
    grid on;
    xlim([0 Ts]);
    
    % Plot im part: analog vs original discrete samples
    subplot(2,1,2);
    plot(tAnalog, imag(xAnalog), 'LineWidth', 1); hold on;
    stem(sampleTimes, imag(x_cp), 'r', 'filled');
    xlabel('Time [s]');
    ylabel('Im\{x\}');
    title('Im part: xAnalog vs original x_{cp} samples');
    legend('xAnalog (oversampled)', 'x_{cp} samples');
    grid on;
    xlim([0 Ts]);

    %% Bandpass signal
    figure;
    plot(tAnalog, xBP);
    xlabel('Time [s]');
    ylabel('x_{BP}(t)');
    title('Bandpass OFDM signal (real-valued)');
    xlim([0 Ts]);
    grid on;

    % z(t) plot
    figure;
    plot(t_z, z);
    title('Received bandpass signal z(t)');
    xlabel('Time [s]');
    ylabel('z(t)');
    grid on;
    xlim([0 Ts + max(channelTaus)]);
end