%% -----------------------------------------------------------------------
%
% Title       : test_fast_conv.m
% Author      : Alexander Kapitanov	
% Company     : AO "Insys"
% E-mail      : sallador@bk.ru 
% Version     : 1.0	 
%
% ------------------------------------------------------------------------
%
% Description :  
%  
%  Math model for calculate Fast Convolution by using FFT algorithm
%  Calculate linear convolution via circular convolution.
%
%  The linear convolution of an N-point vector, x, and an L-point vector, y, 
%  has length N + L - 1.
%
% ------------------------------------------------------------------------
%
% Version     : 1.0 
% Date        : 2017.05.14 
%
%% ----------------------------------------------------------------------- 
%
%	GNU GENERAL PUBLIC LICENSE
% Version 3, 29 June 2007
%
%	Copyright (c) 2018 Kapitanov Alexander
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
%% -----------------------------------------------------------------------	   

% Preparing to work
close all;
clear all;

set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 14, 'DefaultTextFontName', 'Times New Roman'); 

%% -------------------------------------------------------------------------- %%
% ---------------- CHANGE INPUT SIGNAL AND FFT PARAMETERS -------------------- %
%% -------------------------------------------------------------------------- %%

% Number of FFT length: (Filter length = NFFT / 2)
NFFT = 2^8;            % Number of FFT points
MT = 8;                % M multiplier for total length of data vector

% Signal parameters:
Asig = 2^15-1;         % Signal magnitude
Fsig = 0.3722;         % Signal frequency
B = 0.125;             % Base for chirp signal

SNR = -10;             % AWGN

% Filter parameters (FIR) and beta for Kaiser windowing:
Fcut = 190;           % First freq (passband)
Fs = 1000;            % Sampling freq
BETA = 7;             % Add window: Beta (Kaiser)

%% -------------------------------------------------------------------------- %%
% ---------------- 0: CREATE INPUT DATA FOR CPP/RTL -------------------------- %
%% -------------------------------------------------------------------------- %%
NT = MT/2;             % N multiplier (don't change)
NSIG = NFFT * MT;      % Total signal length
SEED = 1;

for i = 1:NSIG
  
  % Chirp:
  Dre(i,1) = round(Asig * cos((Fsig*i + B*i*i/2) * 1*pi/NFFT) * sin(2*i * pi / NSIG));
  Dim(i,1) = round(Asig * sin((Fsig*i + B*i*i/2) * 1*pi/NFFT) * sin(2*i * pi / NSIG));

  % Harmonic w/ modulation:
  % Dre(i,1) = Asig * cos(Fsig * i* 2*pi/NFFT) * sin(2*i * pi / NSIG);
  % Dim(i,1) = Asig * sin(Fsig * i* 2*pi/NFFT) * sin(2*i * pi / NSIG);

  % Harmonic
  % Dre(i,1) = Asig * cos(Fsig * i* 2*pi/NFFT);
  % Dim(i,1) = Asig * sin(Fsig * i* 2*pi/NFFT);
 
end

DatRe = awgn(Dre, SNR, 0, SEED);     
DatIm = awgn(Dim, SNR, 0, SEED);     

DSVRe = round(DatRe);
DSVIm = round(DatIm);

%% -------------------------------------------------------------------------- %%
% ---------------- 1:  INPUT DATA (TEST MATH SIGNAL) ------------------------- % 
%% -------------------------------------------------------------------------- %%
fid = fopen ("di_re.dat", "w");
for i = 1:NFFT*MT
    fprintf(fid, "%d \n", DSVRe(i));
end
fclose(fid);
%
fid = fopen ("di_im.dat", "w");
for i = 1:NFFT*MT
    fprintf(fid, "%d \n", DSVIm(i));
end
fclose(fid);

Din(:,1) = DSVRe;
Din(:,2) = DSVIm;
Dcm(:,1) = Din(:,1) + j*Din(:,2);

% Create two parts of input signal: normal data and shifted on NFFT/2
for i = 1:NFFT*NT
  Din1(i,1) = Din(i,1);
  Din1(i,2) = Din(i,2);
  Din2(i,1) = Din(i+NFFT/2,1);
  Din2(i,2) = Din(i+NFFT/2,2);
end

Dcm1(:,1) = Din1(:,1) + j*Din1(:,2);
Dcm2(:,1) = Din2(:,1) + j*Din2(:,2);

%% -------------------------------------------------------------------------- %%
% ---------------- 2:  SUPPORT FUNCTION FOR FC-FILTER ------------------------ % 
%% -------------------------------------------------------------------------- %%

N = NFFT/2;           % Number of taps for FIR filter = NFFT/2

WIND = kaiser(N, BETA);
f =  [Fcut]/(Fs/2);
Hc = fir1(N-1, f, 'low', WIND);
%Hc2 = [Hc zeros(1,N)];

DFc(:,1) = conv(Din(1:NT*NFFT,1), Hc);
DFc(:,2) = conv(Din(1:NT*NFFT,2), Hc);

SF = fft(Hc, NFFT);
SxFFT(:,1) = real(SF);
SxFFT(:,2) = imag(SF);
SFCm = SF';

SXFFT(:,1) = round (2^14 * SxFFT(:,1));
SXFFT(:,2) = round (2^14 * SxFFT(:,2));

fid = fopen ("sf_x64.dat", "w");
for i = 1:NFFT/2
    fprintf(fid, "%d %d %d %d\n",  SXFFT(i+NFFT/2,2), SXFFT(i+NFFT/2,1), SXFFT(i,2), SXFFT(i,1));
end
fclose(fid);

%% -------------------------------------------------------------------------- %%
% ---------------- 3:  CREATE TWO PARTS OF SPECTRUM -------------------------- % 
%% -------------------------------------------------------------------------- %%
for j = 1:NT
  Dfft1(:,j) = fft(Dcm1((j-1)*NFFT+1 : NFFT*j, 1), NFFT);
  Dfft2(:,j) = fft(Dcm2((j-1)*NFFT+1 : NFFT*j, 1), NFFT);
end

Dfft_l1 = reshape(Dfft1, NT*NFFT, 1);
Dfft_l2 = reshape(Dfft2, NT*NFFT, 1);

DatFFT1(:,1) = real(Dfft_l1);
DatFFT1(:,2) = imag(Dfft_l1);
DatFFT2(:,1) = real(Dfft_l2);
DatFFT2(:,2) = imag(Dfft_l2);

%% -------------------------------------------------------------------------- %%
% ---------------- 4:  CALCULATE FC FILTER (FFT-IFFT) ------------------------ % 
%% -------------------------------------------------------------------------- %%

for j = 1:NT
  DCM1(:,j) = Dfft1(:,j) .* SFCm;
  DCM2(:,j) = Dfft2(:,j) .* SFCm;
end

DCm_l1 = reshape(DCM1, NT*NFFT, 1);
DCm_l2 = reshape(DCM2, NT*NFFT, 1);    

DatCM1(:,1) = real(DCm_l1);
DatCM1(:,2) = imag(DCm_l1);
DatCM2(:,1) = real(DCm_l2);
DatCM2(:,2) = imag(DCm_l2);  

for j = 1:NT
  Difft1(:,j) = ifft(DCm_l1((j-1)*NFFT+1 : NFFT*j, 1), NFFT);
  Difft2(:,j) = ifft(DCm_l2((j-1)*NFFT+1 : NFFT*j, 1), NFFT);
end

Difft_l1 = reshape(Difft1, NT*NFFT, 1);
Difft_l2 = reshape(Difft2, NT*NFFT, 1);

DatIFFT1(:,1) = real(Difft_l1);
DatIFFT1(:,2) = imag(Difft_l1);
DatIFFT2(:,1) = real(Difft_l2);
DatIFFT2(:,2) = imag(Difft_l2);

%% -------------------------------------------------------------------------- %%
% ---------------- 5:  CALCULATE Output parts of data ------------------------ % 
%% -------------------------------------------------------------------------- %%

for i = 1:NFFT
  for j = 1:NT
    if (i <= NFFT/2)
      Dout1(i+(j-1)*NFFT,1) = 0;
      Dout1(i+(j-1)*NFFT,2) = 0;
    else
      Dout1(i+(j-1)*NFFT,1) = DatIFFT2(i-NFFT/2 + (j-1)*NFFT,1);
      Dout1(i+(j-1)*NFFT,2) = DatIFFT2(i-NFFT/2 + (j-1)*NFFT,2);  
    end
  end
end

for i = 1:NFFT
  for j = 1:NT
    if (i <= NFFT/2)
      Dout2(i+(j-1)*NFFT,1) = DatIFFT1(i + (j-1)*NFFT,1);
      Dout2(i+(j-1)*NFFT,2) = DatIFFT1(i + (j-1)*NFFT,2);
    else
      Dout2(i+(j-1)*NFFT,1) = 0;
      Dout2(i+(j-1)*NFFT,2) = 0;  
    end
  end
end

%% -------------------------------------------------------------------------- %%
% ---------------- 6:  PLOT: INPUT, SF SPEC, FFT/IFFT, OUTPUT ---------------- % 
%% -------------------------------------------------------------------------- %%

figure(1) % Plot loaded data in Time Domain
for i = 1:2
  subplot(5,2,i)
  plot(Din1(1:NT*NFFT,i), '--', 'LineWidth', 1, 'Color',[2-i 0  i-1])
  grid on
  hold on
  plot(Din2(1:NT*NFFT,i), '-.', 'LineWidth', 1, 'Color',[i-1 0  2-i])
  axis tight      
  title(['Input Data in Time Domain'])   
end

figure(1) % Plot FFT data in Frequency Domain
for i = 1:2
  subplot(5,2,i+2)
  plot(DatFFT1(1:NT*NFFT,i), '--', 'LineWidth', 1, 'Color',[2-i 0  i-1])
  grid on
  hold on 
  plot(DatFFT2(1:NT*NFFT,i), '-.', 'LineWidth', 1, 'Color',[i-1 0  2-i])
  axis tight       
  title(['Input Spectrum'])    
end

figure(1) % Plot FFT data in Frequency Domain
  subplot(5,2,5)
  plot(Hc, '-', 'LineWidth', 1, 'Color',[1 0  0])
  grid on
  axis tight       
  title(['Impulse responce'])    

figure(1) % Plot FFT data in Frequency Domain
for i = 1:2
  subplot(5,2,6)
  plot(SxFFT(1:NFFT,i), '-', 'LineWidth', 1, 'Color',[2-i 0  i-1])
  hold on
  grid on
  axis tight       
  title(['Freq responce'])    
end

figure(1) % Plot loaded data in Time Domain
for i = 1:2
  subplot(5,2,i+6)
  plot(DatCM1(1:NT*NFFT,i), '--', 'LineWidth', 1, 'Color',[2-i 0  i-1])
  grid on
  hold on
  plot(DatCM2(1:NT*NFFT,i), '-.', 'LineWidth', 1, 'Color',[i-1 0  2-i])
  axis tight      
  title(['Complex mult (w/ SF)'])   
end  

figure(1) % Plot loaded data in Time Domain
for i = 1:2
  subplot(5,2,i+8)
  plot(Dout1(1:NT*NFFT,i), '--', 'LineWidth', 1.5, 'Color',[2-i 0  i-1])
  grid on
  hold on
  plot(Dout2(1:NT*NFFT,i), '-.', 'LineWidth', 1.5, 'Color',[i-1 0  2-i])
  axis tight      
  title(['Output Data'])   
end

for i = 1:NFFT
  for j = 1:NT 
    if (i <= NFFT/2)
      Dout(i+(j-1)*NFFT, 1) = DatIFFT1(i+(j-1)*NFFT, 1);
      Dout(i+(j-1)*NFFT, 2) = DatIFFT1(i+(j-1)*NFFT, 2);
    else
      Dout(i+(j-1)*NFFT, 1) = DatIFFT2(i+(j-1)*NFFT - NFFT/2, 1);
      Dout(i+(j-1)*NFFT, 2) = DatIFFT2(i+(j-1)*NFFT - NFFT/2, 2);
    end
  end
end

figure(2)
for i = 1:2
  subplot(2,2,i)
  plot(Dout(1:NT*NFFT,i), '-', 'LineWidth', 1, 'Color',[2-i 0  i-1])
  grid on
  axis tight      
  title(['Test Data FFT-iFFT alg.'])   
end

figure(2)
for i = 1:2
  subplot(2,2,i+2)
  plot(DFc(1:NT*NFFT,i), '-', 'LineWidth', 1, 'Color',[2-i 0  i-1])
  grid on
  axis tight      
  title(['Fast Conv Data (conv func)'])   
end