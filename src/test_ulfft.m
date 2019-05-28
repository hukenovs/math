%% ------------------------------------------------------------------------
%
% Title       : test_ulfft.m
% Author      : Alexander Kapitanov 
% E-mail      : sallador@bk.ru 
% Version     : 1.0
%
% -------------------------------------------------------------------------
%
% Description : 
%    Testing Ultra Long FFT model
%
% -------------------------------------------------------------------------
%
% Date        : 2018.04.15 
%
%% ------------------------------------------------------------------------ 
%
% GNU GENERAL PUBLIC LICENSE
% Version 3, 29 June 2007
%
% Copyright (c) 2018 Kapitanov Alexander
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

% Preparing to work
close all;
clear all;
set(0, 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 12, 'DefaultTextFontName', 'Times New Roman'); 


%% ---- Settings ----------------------------------------------------------
N1      = 16;           % Colomn
N2      = 8;            % Rows
NFFT    = N1 * N2;      % Total FFT length

Asig    = 1;
Fsig    = 4;
B       = 0.4;

%% ---- Input data for FFT calculate --------------------------------------
Dre = zeros(NFFT, 1);
Dim = zeros(NFFT, 1);

for i = 0:NFFT-1
  Dre(i+1,1) = Asig * cos(Fsig*i* 2*pi/NFFT);
  Dim(i+1,1) = Asig * sin(Fsig*i* 2*pi/NFFT);

  Dre(i+1,1) = Asig * cos((Fsig*i + B*i*i/2) * 2*pi/NFFT) * sin(i * pi / NFFT);
  Dim(i+1,1) = Asig * sin((Fsig*i + B*i*i/2) * 2*pi/NFFT) * sin(i * pi / NFFT);
  
  if (i == 4)
    Dre(i+1,1) = Asig;
    Dim(i+1,1) = 0;
  else
    Dre(i+1,1) = 0;
    Dim(i+1,1) = 0;
  end
end

% Adding noise to real signal
SNR = 90;

DatRe = awgn(Dre, SNR, 0, 1);
DatIm = awgn(Dim, SNR, 0, 1);

DAT_IN = DatRe + 1j*DatIm;

IR0(:,1) = DatRe;
IR0(:,2) = DatIm;

TX(:,1) = real(DAT_IN);
TX(:,2) = imag(DAT_IN);

%% ------------------------------------------------------------------------
%  ---- FFT CALCULATE (ETALON FFT ALGORITHM) ------------------------------
%  ------------------------------------------------------------------------
DFFT0 = fft(DAT_IN);

IX(:,1) = real(DFFT0);
IX(:,2) = imag(DFFT0);

%% ------------------------------------------------------------------------
%  ---- ULTRA LONG FFT ALGORITM FROM 2D FFT + ROTATION --------------------
%  ------------------------------------------------------------------------

% ---- STAGE 0: CREATE MATRIX FOR INPUT DATA ----
DIN = zeros(N1,N2);
for n1 = 1:N1
  for n2 = 1:N2
    DIN(n1, n2) = DAT_IN((n2-1)*N1+n1,1);
  end 
end

DIN_L = zeros(NFFT,1);
for n1 = 1:N1
  for n2 = 1:N2
    DIN_L((n1-1)*N2+n2, 1) = DIN(n1, n2);
  end 
end

IR1(:,1) = real(DIN_L);
IR1(:,2) = imag(DIN_L);

% ---- STAGE 1: CALCULATE FIRST FFT (A) ----
DFFT1 = zeros(N1,N2);
for n1 = 1:N1
  DFFT1(n1,:) = fft(DIN(n1,:));
end

DFFT1_L = zeros(NFFT,1);
for n1 = 1:N1
  for n2 = 1:N2
    DFFT1_L((n1-1)*N2+n2, 1) = DFFT1(n1, n2);
  end 
end

IR2(:,1) = real(DFFT1_L);
IR2(:,2) = imag(DFFT1_L);

% ---- STAGE 2: SHUFFLE FFT RESULT ----
SHF1_L = zeros(NFFT,1);
for n1 = 1:N1
  for n2 = 1:N2
    SHF1_L((n2-1)*N1+n1,1) = DFFT1(n1, n2);
  end 
end

IR3(:,1) = real(SHF1_L);
IR3(:,2) = imag(SHF1_L);

% SHF1 = zeros(N1,N2);
% for n1 = 1:N1
%   for n2 = 1:N2
%     SHF1(n1, n2) = SHF1_L((n1-1)*N2+n2,1);
%   end 
% end

% ---- STAGE 3: ROTATE AFTER SHUFFLER ----
SHW = zeros(N1,N2);
for n1 = 1:N1
  for n2 = 1:N2
    SHW(n1, n2) = exp(-1j*2*pi*(n1-1)*(n2-1)/(N1*N2)); 
  end
end

TWD_L = zeros(NFFT,1);
for n1 = 1:N1
  for n2 = 1:N2
    TWD_L((n2-1)*N1+n1, 1) = SHW(n1, n2);
  end 
end

IR4(:,1) = real(TWD_L);
IR4(:,2) = imag(TWD_L);

% ---- STAGE 4: MULTIPLY DATA ----
MULT = zeros(N1,N2);
for n1 = 1:N1
  for n2 = 1:N2
    MULT(n1, n2) = DFFT1(n1, n2) * SHW(n1, n2);
  end
end

MLT_L = zeros(NFFT,1);
for n1 = 1:N1
  for n2 = 1:N2
    MLT_L((n2-1)*N1+n1, 1) = MULT(n1, n2);
  end 
end

IR5(:,1) = real(MLT_L);
IR5(:,2) = imag(MLT_L);

% ---- STAGE 5: CALCULATE SECOND FFT (B) ----
DFFT2 = zeros(N1,N2);
for n2 = 1:N2
    DFFT2(:,n2) = fft(MULT(:,n2));
end

DFFT2_L = zeros(NFFT,1);
for n1 = 1:N1
  for n2 = 1:N2
    DFFT2_L((n2-1)*N1+n1, 1) = DFFT2(n1, n2);
  end 
end

IR6(:,1) = real(DFFT2_L);
IR6(:,2) = imag(DFFT2_L);

% ---- STAGE 4: CREATE MATRIX FOR OUTPUT DATA ----
DOUT = zeros(NFFT,1);
for n1 = 1:N1
  for n2 = 1:N2
    DOUT((n1-1)*N2+n2, 1) = DFFT2(n1, n2);
  end 
end

IR7(:,1) = real(DOUT);
IR7(:,2) = imag(DOUT);

%% ------------------------------------------------------------------------
%  ---- PLOT DATA: INPUT SIGNAL, FFT0, FFT1, TWIDDLES etc -----------------
%  ------------------------------------------------------------------------
figure('Name', 'Ultra-long FFT', 'NumberTitle','off');
for i = 1:2
  subplot(4,2,1)
  plot(IR0(:, i), '-', 'LineWidth', 1, 'Color', [2-i 0  i-1]);
  grid on; hold on; axis tight;
  title(['Input signal']) 

  subplot(4,2,2)
  plot(IR1(:, i), '-', 'LineWidth', 1, 'Color', [2-i 0  i-1]);
  grid on; hold on; axis tight;
  title(['Shuffle data']) 
  
  subplot(4,2,3)
  plot(IR2(:, i), '-', 'LineWidth', 1, 'Color', [2-i 0  i-1]);
  grid on; hold on; axis tight;
  title(['FFT0 (N1)'])
 
  subplot(4,2,4)
  plot(IR3(:, i), '-', 'LineWidth', 1, 'Color', [2-i 0  i-1]);
  grid on; hold on; axis tight;
  title(['FFT0 (Shuffle)'])
 
  subplot(4,2,5)
  plot(IR4(:, i), '-', 'LineWidth', 1, 'Color', [2-i 0  i-1]);
  grid on; hold on; axis tight;
  title(['Twiddle (Cordic)'])
 
  subplot(4,2,6)
  plot(IR5(:, i), '-', 'LineWidth', 1, 'Color', [2-i 0  i-1]);
  grid on; hold on; axis tight;
  title(['Rotate after Cordic'])

  subplot(4,2,7)
  plot(IR6(:, i), '-', 'LineWidth', 1, 'Color', [2-i 0  i-1]);
  grid on; hold on; axis tight;
  title(['FFT1 (N2)']) 
  
  subplot(4,2,8)
  plot(IR7(:, i), '-', 'LineWidth', 1, 'Color', [2-i 0  i-1]);
  grid on; hold on; axis tight;
  title(['Output signal'])
end

Col = ['k' 'r' 'b' 'g' 'm', 'c', 'y', 'k'];
Val = ['-' '--' ':' '-.'];
%% ------------------------------------------------------------------------