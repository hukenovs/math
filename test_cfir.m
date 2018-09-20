%% -----------------------------------------------------------------------
%
% Title       : test_cfir.m
% Author      : Alexander Kapitanov	
% Company     : AO "Insys"
% E-mail      : sallador@bk.ru 
% Version     : 1.0	 
%
% ------------------------------------------------------------------------
%
% Description : 
%    Top level for testing Complex FIR filters: time/freq methods
%
% ------------------------------------------------------------------------
%
% Version     : 1.0 
% Date        : 2017.07.01 
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

% Settings
NFFT = 2^8;             % Sampling Frequency
tt = 1:NFFT;            % Time vector #2
VAL_SHIFT = 0;

Asig = (2^15)-1;
Fsig = 8;
F0 = 0;

Fd = 1;
Fm = NFFT/1;
B = Fm / NFFT;
Ffm = 1;
% For testing FORWARD and INVERSE FFT: FWT

STAGE = log2(NFFT);

%% -------------------------------------------------------------------------- %%
% ---------------- 0: CREATE INPUT DATA FOR CPP/RTL -------------------------- % 
%% -------------------------------------------------------------------------- %%
F = 1;

for i = 1:NFFT
  if (i < NFFT/Fd)
    Dre(i,1) = round(Asig * cos(F0 + (Fsig*i + B*i*i/2) * 2*Fd*pi/NFFT) * abs(sin(Fd * i * Ffm * pi / NFFT)));
    Dim(i,1) = round(Asig * sin(F0 + (Fsig*i + B*i*i/2) * 2*Fd*pi/NFFT) * abs(sin(Fd * i * Ffm * pi / NFFT)));    
  else
    Dre(i,1) = 0;
    Dim(i,1) = 0;
  end  
    Dre(i,1) = round(Asig * cos(F0 + (Fsig*i + B*i*i/2) * 2*pi/NFFT) * sin(i * Ffm * pi / NFFT));
    Dim(i,1) = round(Asig * sin(F0 + (Fsig*i + B*i*i/2) * 2*pi/NFFT) * sin(i * Ffm * pi / NFFT));  
end

for i = 1:NFFT
  if (i > VAL_SHIFT)
    Xre(i,1) = Dre(i-VAL_SHIFT, 1);
    Xim(i,1) = Dim(i-VAL_SHIFT, 1);
  else
    Xre(i,1) = Dre(NFFT-VAL_SHIFT+i, 1);
    Xim(i,1) = Dim(NFFT-VAL_SHIFT+i, 1);   
  end
end
Dre = Xre;
Dim = Xim;

% Adding noise to real signal 
SNR = -70;
SEED = 1;

DatRe = awgn(Dre, SNR, 0, SEED);     
DatIm = awgn(Dim, SNR, 0, SEED);     

DatRe = awgn(Dre, SNR, 0, SEED);     
DatIm = awgn(Dim, SNR, 0, SEED);     

Mre = max(abs(DatRe));
Mim = max(abs(DatIm));
Mdt = max(Mre, Mim);

DSVRe = round(((2^15 - 1)/Mdt)*DatRe);
DSVIm = round(((2^15 - 1)/Mdt)*DatIm);

Din(:,1) = DSVRe;
Din(:,2) = DSVIm;
Dcm(:,1) = Din(:,1) + j*Din(:,2);

if (1)
  figure(1) % Plot loaded data in Freq Domain
  for i = 1:2
    subplot(2,1,i)
    plot(tt(1:NFFT), Din(1:NFFT,i), '-', 'LineWidth', 1, 'Color',[2-i 0  i-1])
    grid on
    axis tight 
    title(['Input data for test module'])  
  end
end

Dfft = fft(Dcm) / NFFT; 
Dabs = (1/(NFFT)) * abs(Dfft).^2; 

DatFFT(:,1) = real(Dfft);
DatFFT(:,2) = imag(Dfft);

%% -------------------------------------------------------------------------- %%
% ---------------- 2:  SF DATA CMP (SUPPLY FUNCTION) ------------------------- % 
%% -------------------------------------------------------------------------- %%
Aof = (2^15)-1;
Fof = 1;
Fm = NFFT/1;
B = Fm / NFFT;

for i = 1:NFFT
    SF(i,1) = round(Aof * cos(F0 + (Fof*i + B*i*i/2) * 2*pi/NFFT) * sin(i * Ffm * pi / NFFT));
    SF(i,2) = round(Aof * sin(F0 + (Fof*i + B*i*i/2) * 2*pi/NFFT) * sin(i * Ffm * pi / NFFT));   
  if (i < NFFT/Fd)
    SF(i,1) = round(Aof * cos(F0 + (Fof*i + B*i*i/2) * 2*Fd*pi/NFFT) * abs(sin(Fd * i * Ffm * pi / NFFT)));
    SF(i,2) = round(Aof * sin(F0 + (Fof*i + B*i*i/2) * 2*Fd*pi/NFFT) * abs(sin(Fd * i * Ffm * pi / NFFT)));    
  else
    SF(i,1) = 0;
    SF(i,2) = 0;
  end      
end

for i = 1:NFFT
    SF(i,1) = round(Aof * cos(F0 + (Fof*i + B*i*i/2) * 2*pi/NFFT) * sin(i * Ffm * pi / NFFT));
    SF(i,2) = round(Aof * sin(F0 + (Fof*i + B*i*i/2) * 2*pi/NFFT) * sin(i * Ffm * pi / NFFT));   
end

% Combine new data to matrix
SFCm(:,1) = SF(:,1) + j*SF(:,2);

SF_FFT = fft(SFCm) / NFFT; 
Sabs = (1/(NFFT)) * abs(SF_FFT).^2; 

% Preparing data to integer:
SfFFT(:,1) = real(SF_FFT);
SfFFT(:,2) = -imag(SF_FFT);

Sfft = SfFFT(:,1) + j*SfFFT(:,2);

for i = 1:NFFT
    CMsig(i,1) = Dfft(i,1) * Sfft(i,1);
end

CMDat(:,1) = real(CMsig);
CMDat(:,2) = imag(CMsig);

Difft = ifft(CMsig) / NFFT; 
Dabs = (1/(NFFT)) * abs(Difft).^2; 

DatIFFT(:,1) = real(Difft);
DatIFFT(:,2) = imag(Difft);

figure(2) % Plot FFT data in Frequency Domain
for i = 1:2
    subplot(2,2,i)
    plot(tt(1:NFFT), DatIFFT(1:NFFT,i), '-', 'LineWidth', 1, 'Color',[2-i 0  i-1])
    grid on
    hold on
    axis tight    
    title(['Output data (Freq method)'])   
end

%% -------------------------------------------------------------------------- %%
% ---------------- 3:  FIR FILTER COMBINATION -------------------------------- % 
%% -------------------------------------------------------------------------- %%

Dat1 = DSVRe + DSVIm;
Dat2 = DSVIm - DSVRe;
Dat3 = DSVRe;

Cm1 = -SF(:,2);
Cm2 = SF(:,1);
Cm3 = SF(:,1) - SF(:,2);

Dx1 = filter2(Cm1, Dat1);
Dx2 = filter2(Cm2, Dat2);
Dx3 = filter2(Cm3, Dat3);

Ore = - Dx1 + Dx3;
Oim = Dx2 + Dx3;

Ore = Ore / NFFT / NFFT / NFFT;
Oim = Oim / NFFT / NFFT / NFFT;

for i = 1:NFFT
  if (i > NFFT/2)
    DTO(i,1) = Ore(i-NFFT/2, 1);
    DTO(i,2) = Oim(i-NFFT/2, 1);
  else
    DTO(i,1) = Ore(NFFT-NFFT/2+i, 1);
    DTO(i,2) = Oim(NFFT-NFFT/2+i, 1);   
  end
end


figure(2) % Plot FFT data in Frequency Domain
for i = 1:2
    subplot(2,2,i+2)
    plot(tt(1:NFFT), DTO(1:NFFT,i), '-', 'LineWidth', 1, 'Color',[2-i 0  i-1])
    grid on
    axis tight    
    title(['Output data (Time method)'])   
end

