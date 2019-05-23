%% ------------------------------------------------------------------------
%
% Title       : test_fir.m
% Author      : Alexander Kapitanov 
% E-mail      : sallador@bk.ru 
% Version     : 1.0  
%
% -------------------------------------------------------------------------
%
% Description : 
%    Create your FIR filters and use it in FPGA projects!
%
% -------------------------------------------------------------------------
%
% Version     : 1.0 
% Date        : 2017.03.20 
%
%% ------------------------------------------------------------------------ 
%
%   GNU GENERAL PUBLIC LICENSE
% Version 3, 29 June 2007
%
%   Copyright (c) 2018 Kapitanov Alexander
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
% APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT 
% HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY 
% OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, 
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM 
% IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF 
%  ALL NECESSARY SERVICING, REPAIR OR CORRECTION. 
%
%% ------------------------------------------------------------------------
  
set(0, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 11, 'DefaultTextFontName', 'Times New Roman'); 

close all;
clear all;
clf;

%% ------------------------------------------------------------------------

% if IS_XXX = 'y' - use section, else if IS_XXX = 'n' - don't use section
IS_PLOT1    = 'Y'; % Main info: FIR: Impulse/Freq responce, Quant. errors
IS_PLOT2    = 'Y'; % Signal data plot via FIR filter

IS_COE      = 'N'; % create *.COE Xilinx file
IS_TAPS     = 'N'; % create *.TAPS file for Graychip
IS_HDR      = 'N'; % create *.H file (header)

IS_HLS      = 'N'; % use example for Vivado HLS

f1          = 125000; % First freq (passband)
Fs          = 1000000; % Sampling freq

% f2 = 110000; % Second freq (stopband)

N           = 256; % Filter order
COE_WIDTH   = 24; % Real width for FIR coeffs (Implementation)

BETA        = 5; % Beta (Kaiser)
WIND        = kaiser(N, BETA); % KAISER WINDOW IS USED!
%WIND       = hamming(N, 'symmetric'); 

NFFT        = 2^12; % Number of FFT points (may leave this)

%% ------------------------------------------------------------------------
% ---- Set freqs for signal
F1          = 10; 
F2          = 100; 
F3          = 1400; 
F4          = 0; 

Asig        = (2^15)-1;
Fsig        = 4;
Fm          = NFFT/8;
Fnz         = 350;
Anz         = 1;

SNR         = 30; % Set noise params

%% ------------------------------------------------------------------------
% Calculate filter
%delta_f = f2-f1;
f   =  [f1]/(Fs/2);
t   = 1:N;

% Filter type: 'low', 'high', 'stop', 'pass';
Hc = fir1(N-1, f, 'low', WIND);

Hc_n = (Hc/max(Hc))*(2^(COE_WIDTH-1)-1);
Hc_r = round(Hc_n);
Hf = 20 * log10(abs(fftshift(fft(Hc_n, NFFT))));
Sp_err = 20*log10(abs(fftshift(fft(Hc_r, NFFT))));

% Hph = arg(fftshift(fft(Hc_n, NFFT)));
Hph = (fftshift(fft(Hc_n, NFFT)));
ff = -0.5:1/NFFT:0.5-1/NFFT;

%freqz(Hc_r);
%freqz(Hc_n);

% Plot FIR in T/F Domains
if (IS_PLOT1 == 'Y')
  figure(1)
    subplot(2,2,1);
    plot(t, Hc_r, '*-', 'LineWidth', 1, 'Color',[1 0 0]);
    axis tight;
    title(['Filter IR, Order = ',num2str(N)]);
    xlabel ('Time');
    ylabel ('Magnitude');
    grid on;

    %  subplot(2,2,2);
    %  plot(ff, Hf-max(Hf), '-', 'LineWidth', 1, 'Color',[0 0 1]);
    %  %axis tight;
    %  axis([0 ff(NFFT) -120 0])
    %  title('FIR filter Spectrum');
    %  xlabel ('Freq ( x rad / samples)');
    %  ylabel ('Magnitude (dB)');  
    %  grid on;
 
    subplot(2,2,2);
    plot(ff, Hf-max(Hf), '-', 'LineWidth', 1, 'Color',[1 0 0]);
    hold on
    
    plot(ff, Sp_err-max(Sp_err), '-', 'LineWidth', 1, 'Color',[0 0 1]);
    grid on;
    axis([0 ff(NFFT) -120 0]) 
    title('FIR Spectrum');
    xlabel ('Freq ( x rad / samples)');
    ylabel ('Magnitude (dB)');  
    legend ('Double dtype', ['CoeWidth: ',num2str(COE_WIDTH)], 'location', 'east'); 
end

%% ------------------------------------------------------------------------
% Find calc error
q1 = abs(Hc_n - Hc_r);
Hc_err = 20*log10(abs(fftshift(fft(q1, NFFT))));

if (IS_PLOT1 == 'Y')
  figure(1)
    subplot(2,2,3);
    plot(q1, '*-', 'LineWidth', 1, 'Color',[0 1 0]);
    axis tight;
    title('Quantization Error (Time)');
    xlabel ('Time');
    ylabel ('Magnitude');  
    grid on; 

    subplot(2,2,4);
    plot(ff(NFFT/2:NFFT), Hc_err(NFFT/2:NFFT), '-', 'LineWidth', 1, 'Color',[1 0 1]);
    title('Quantization Error (Freq)');
    % axis([0 0.5 -60 20])
    axis tight;
    xlabel ('Freq ( x rad / samples)');
    ylabel ('Magnitude (dB)');  
    grid on;
end

%% -----------------------------------------------------------------------
% Find spectrun error for implementation

% Sp_err = 20*log10(abs(fftshift(fft(Hc_r, NFFT))));

%figure(2)
%  plot(ff, Hf-max(Hf), '-', 'LineWidth', 1, 'Color',[1 0 0]);
%  hold on
%  
%  plot(ff, Sp_err-max(Sp_err), '-', 'LineWidth', 1, 'Color',[0 0 1]);
%  grid on;
%  axis([0 ff(NFFT) -120 0]) 
%  title('Frequency Spectrum for Implementation');
%  xlabel ('Freq ( x rad / samples)');
%  ylabel ('Magnitude (dB)');  
%  legend ('Red - Double', 'Blue - Quantized', 'location', 'northeast');
%% -----------------------------------------------------------------------

% Create signal for FIR filter !
SEED = 1;

% Signal #1
SIGNAL = sin(2*pi*[1:NFFT]*(F1/NFFT))+sin(2*pi*[1:NFFT]*(F2/NFFT)) + ...
  + sin(2*pi*[1:NFFT]*(F3/NFFT))+sin(2*pi*[1:NFFT]*(F4/NFFT));

% Signal #2

SIGNAL = zeros(1,NFFT);
for i = 1:NFFT
  SIGNAL(1,i) = round(Asig * cos((Fsig*i + (Fm / NFFT)*i*i/2) * 2*pi/NFFT) * sin(i* pi / NFFT)) + ...
  + round((Asig/Anz) * cos(Fnz * i * pi / NFFT));
% SIGNAL(1,i) = sin(600 * i * pi / NFFT);

end
  
Dt = awgn(SIGNAL, SNR, 0, SEED);     
SIGNAL = Dt;  

FFT_SIG = 20 * log10(abs(fftshift(fft(SIGNAL, NFFT))));
SIGFIL = filter2 (Hc_r, SIGNAL, "same");
% SIGFIL = filter(Hc_r, 1, SIGNAL);

NSig_n = SIGNAL / max(SIGNAL);
FSig_n = SIGFIL / max(SIGFIL);

%% -----------------------------------------------------------------------
figure(2)
  subplot(2,2,1);
  plot(NSig_n, '-', 'LineWidth', 1, 'Color',[1 0 0]);
  grid on;  
  %axis tight;
  axis([0 NFFT/1 -1 1]) 
  title('Signal w/ harmonics');
  xlabel('Time');
  ylabel('Normalized');
  
  subplot(2,2,2)
  plot(FSig_n, '-', 'LineWidth', 1, 'Color',[0 0 1]);
  grid on; 
  %axis tight;
  axis([0 NFFT/1 -1 1]) 
  title('Filtered Signal');
  xlabel('Time');
  ylabel('Normalized');

 
XFIL = 20 * log10(abs(fftshift(fft(Hc_r, NFFT))));
FFT_OUT = 20 * log10(abs(fftshift(fft(SIGFIL, NFFT))));
% FFT_OUT = 20 * log10( 10.^(FFT_SIG/20) .* 10.^(XFIL/20));
% FFT_OUT = FFT_SIG + XFIL;
figure(2)
  subplot(2,2,3);
  plot(ff, FFT_SIG-max(FFT_SIG), '-', 'LineWidth', 2, 'Color',[1 0 0]);
  hold on; 
  axis([0 ff(NFFT) -100 0]) 
  
  plot(ff, XFIL-max(XFIL), '-', 'LineWidth', 2, 'Color',[0 0 1]);
  hold off;
  %axis tight;
  axis([0 ff(NFFT) -100 0]) 
  title('Input spectrum');
  xlabel ('Freq ( x rad / samples)');
  ylabel ('Magnitude (dB)');  
  grid on; 
  
  subplot(2,2,4);
  plot(ff, FFT_OUT-max(FFT_OUT), '-', 'LineWidth', 2, 'Color',[1 0 0]);
  axis([0 ff(NFFT) -100 0]) 
  title('Output spectrum');
  xlabel ('Freq ( x rad / samples)');
  ylabel ('Magnitude (dB)');  
  grid on;

%[n, w, beta, ftype] = kaiserord ([1000, 1200], [1, 0], [0.05, 0.05], 11025);
% b = fir1 (n, w, kaiser (n+1, beta), ftype, "noscale");
% freqz (b, 1, [], 11025);  
  
%radix=10;
%coefdata=-200,1200,2047,1200,-200;

%% -----------------------------------------------------------------------
% Load coefs to files: .coe, .taps, .h

if (IS_COE == 'Y') 
  fid = fopen ("fir_filter.coe", "w");
  fprintf(fid, "Radix = 10;\n");
  fprintf(fid, "Coefficient_Width = %d;\n", COE_WIDTH);
  fprintf(fid, "Coefdata =\n");
  for i = 1:N
    if (i == N)
      fprintf(fid, "%d;\n", Hc_r(1,i));
    else
      fprintf(fid, "%d,\n", Hc_r(1,i));
    end
  end
  fclose(fid);
end

if (IS_TAPS == 'Y') 
  fid = fopen ("fir_filter.taps", "w");
  for i = 1:N
    fprintf(fid, "%d\n", Hc_r(1,i));
  end
  fclose(fid);
end

if (IS_HDR == 'Y') 
  fid = fopen ("fir_filter.h", "w");
  fprintf(fid, "const int BL = %d;\n",  N);
  fprintf(fid, "const int B[%d] = {\n", N);
  for i = 1:N
    if (i == N)
      fprintf(fid, "%d}\n", Hc_r(1,i));
    else
      fprintf(fid, "%d,\n", Hc_r(1,i));
    end
  end
  fclose(fid);
end
%% -----------------------------------------------------------------------
if (IS_HLS == 'Y')
  Asig = (2^14)-1;
  Fsig = 1;
  Fm = NFFT/2;

  % Signal #1
  Dat = zeros(NFFT,1);
  for i = 1:NFFT
    Dat(i,1) = round(Asig * cos((Fsig*i + (Fm / NFFT)*i*i/2) * 2*pi/NFFT) * sin(i* pi / NFFT));
  end

  % Signal #2
  MD = 128;
  for i = 1:NFFT
    if (mod(i,MD) == 1)
      Dat(i,1) = round(Asig);
    else
      Dat(i,1) = 0;
    end
  end

  % Adding noise
  SNR = -10;
  SEED = 1;

  DatX = awgn(Dat, SNR, 0, SEED);     
  DSV = round(DatX);
  %% ------------------------------------------------------------------------

  % Save data to file
  % fid = fopen ('din.dat', 'w');
  % for i = 1:NFFT
  %     fprintf(fid, '%d %d\n', i, DSV(i));
  % end
  % fclose(fid);

 
  % FIR Filter
  COEF = Hc_r / 2^18;
  Dout = filter(COEF, 1, DSV);
  DFIR = round(Dout/2^0);

 %% -----------------------------------------------------------------------
  figure(3) % Plot loaded data in Time Domain
  subplot(2,1,1)  
  plot(DSV(1:NFFT), '-', 'LineWidth', 1, 'Color',[1 0  0])

  hold on
  grid on
  axis tight     
  title(['Input Data'])
  %% ------------------------------------------------------------------------
  figure(3) % Plot loaded data in Time Domain
  subplot(2,1,2)
  % plot(Hc_r(1:N), '-', 'LineWidth', 1, 'Color',[0 0  1])
  plot(Dout(1:NFFT), '-', 'LineWidth', 1, 'Color',[0 0  1])
  hold on
  grid on
  axis tight     
  %title(['FIR IR'])
  title(['Output Data'])
  %% ------------------------------------------------------------------------


%  figure(3) % Plot loaded data in Time Domain
%  subplot(2,2,3)
%  plot(Dout(1:NFFT), '-', 'LineWidth', 1, 'Color',[1 0  0])
%  hold on
%  grid on
%  axis tight     
%  title(['Output Data'])
  %% ------------------------------------------------------------------------
  fid = fopen ('fir_gold.dat', 'w');
  for i = 1:NFFT
    fprintf(fid, '%d\n', DFIR(i));
  end
  fclose(fid);

  fid = fopen ('fir_data.dat', 'w');
  for i = 1:NFFT
    fprintf(fid, '%d\n', DSV(i));
  end
  fclose(fid);

  fid = fopen ('fir_coe.dat', 'w');
  for i = 1:N
    fprintf(fid, '%d\n', Hc_r(i));
  end
  fclose(fid);
end
