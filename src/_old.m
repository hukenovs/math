
%% ------------------------------------------------------------------------
%figure(2) % Quad Matrix !!!
%for n1 = 1:N1
%    subplot(2,1,1)
%    plot(1:N2, real(SRW(1:N2, n1)), '-', 'LineWidth', 1, 'Color', Col(mod(n1-1, 8)+1));
%    %SRW(k2, n1) 
%    grid on; hold on; axis tight; 
%    title(['Input signal']) 
%end
%%legend('show', 'Location', 'eastoutside', 'Orientation', 'vertical'); hold off;
%
%figure(2) % Plot loaded data in Time Domain
%for n1 = 1:N1
%    subplot(2,1,2)
%    plot(1:N2, imag(SRW(1:N2, n1)), '-', 'LineWidth', 1, 'Color', Col(mod(n1-1, 8)+1));
%    %SRW(k2, n1) 
%    grid on; hold on; axis tight;
%    title(['Input signal']) 
%end
%legend('show', 'Location', 'eastoutside', 'Orientation', 'vertical'); hold off;
%
%for i = 1:N1
%    for j = 1:N2
%        LRW((i-1)*N2+j, 1) = real(SRW(j, i)) / 1024; 
%        LRW((i-1)*N2+j, 2) = imag(SRW(j, i)) / 1024; 
%    end
%end



%figure(2) % Quad Matrix !!!
%    plot(1:NFFT, LRW(1:NFFT,2), '-', 'LineWidth', 1);
%    grid on; axis tight; 
%title(['WW']);


%RTL_WW = load ("rtl_rotate.dat");
%RTL_WW = RTL_WW / 2^14;

%for i = 1:NFFT
%    DFFT1(:,n1) = fft(DIN(:,n1));
%end

% CHECK ROTATOR
%figure(3) % Quad Matrix !!!
%    plot(1:NFFT, LRW - RTL_WW(1:NFFT,:), '-', 'LineWidth', 1);
%    grid on; axis tight; 
%title(['WW']);


%RWW(:,1) = real(Difft);
%RWW(:,2) = imag(Difft);




%% -----------------------------------------------------------------------
%figure(1) % Plot loaded data in Time Domain
%for i = 1:2
%    %subplot(1,1,1)
%    plot(t(1:NFFT), TX(1:NFFT,i), '-', 'LineWidth', 1, 'Color',[2-i 0  i-1])
%    grid on; hold on; axis tight;
%    title(['Input signal'])  
%end
%hold off;
%%% -----------------------------------------------------------------------
%figure(2) % Plot loaded data in Time Domain
%for i = 1:2
%    subplot(2,1,1)
%    plot(t(1:NFFT), DX(1:NFFT,i), '-', 'LineWidth', 1, 'Color',[2-i 0  i-1])
%    grid on; hold on; axis tight;
%    title(['DFT algorithm result'])  
%end
%hold off;
%%% -----------------------------------------------------------------------
%figure(2) % Plot loaded data in Time Domain
%for i = 1:2
%    subplot(2,1,2)
%    plot(t(1:NFFT), IX(1:NFFT,i), '-', 'LineWidth', 1, 'Color',[2-i 0  i-1])
%    grid on; hold on; axis tight;
%    title(['FFT algorithm result'])
%end
%hold off;
%% -----------------------------------------------------------------------


%function y = fft_decomposition(x, M)
%% y = fft_decomposition(x, M)
%% Computes FFT by decomposing into smaller FFTs.
%%
%% Inputs:
%% x is a 1D array of the input data.
%% M is the size of one of the FFTs to use.
%%
%% Outputs:
%% y is the FFT of x.  It has been computed using FFTs of size M and
%% length(x)/M.
%%
%% Note that this implementation doesn't explicitly use the 2D array U; it
%% works on samples of x in-place.
%
%q = 1; % Offset because MATLAB starts at one.  Set to 0 for C code.
%x_original = x;
%P = length(x);
%if mod(P,M)~=0, error('Invalid block size.'); end;
%N = P/M;
%
%% step 2: FFT-N on rows of U.
%for m = 0 : M-1
%    x(q+(m:M:P-1)) = fft(x(q+(m:M:P-1)));
%end;
%
%% step 3: Twiddle factors.
%for m = 0 : M-1
%    for n = 0 : N-1
%        x(m+n*M+q) = x(m+n*M+q) * exp(-2*pi*j*m*n/P);
%    end;
%end;
%
%% step 4:  FFT-M on columns of U.
%for n = 0 : N-1
%    x(q+n*M+(0:M-1)) = fft(x(q+n*M+(0:M-1)));
%end;
%
%% step 5:  Re-arrange samples for output.
%y = zeros(size(x));
%for m = 0 : M-1
%    for n = 0 : N-1
%        y(m*N+n+q) = x(m+n*M+q);
%    end;
%end;
%
%err = max(abs(y-fft(x_original)));
%fprintf( 1, 'The largest error amplitude is %g\n', err);
%return;
%% End of fft_decomposition().

%% ---- DFT Calculate (simple algorithm) ---- 
%for k = 1:NFFT
%    WX = 0;
%    for n = 1:NFFT
%        % Twiddle:
%        WW(n, k) = exp(-j*2*pi*(n-1)*(k-1)/NFFT); 
%        % ephi(n, k) = arg(WW(n, k));
%        % erho(n, k) = abs(WW(n, k));
%        WX = WX + DAT_IN(n) * WW(n, k);
%    end
%    DAT_IN(k) = WX;
%    WX = 0;
%end
%
%DX(:,1) = real(DAT_IN);
%DX(:,2) = imag(DAT_IN);

