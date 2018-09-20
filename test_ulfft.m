%% -----------------------------------------------------------------------
%
% Title       : test_ulfft.m
% Author      : Alexander Kapitanov	
% Company     : AO "Insys"
% E-mail      : sallador@bk.ru 
% Version     : 1.0	 
%
%-------------------------------------------------------------------------
%
% Description : 
%    Top level for testing Ultra Long FFT model
%
% ------------------------------------------------------------------------
%
% Version     : 1.0 
% Date        : 2018.04.15 
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


%% ---- Settings ----
N1 = 16;          % Colomn
N2 = 16;          % Rows
NFFT = N1 * N2;   % Total FFT length     
t = 0:NFFT-1;     % Time vector #1

Asig = (2^14)-1;
Fsig = 4;
F0 = 0;
Fm = NFFT/1;
B = Fm / NFFT;
Ffm = 1;

%% ---- Input data for FFT calculate ----
for i = 1:NFFT

  Dre(i,1) = i-1;
  Dim(i,1) = i-1;   
  
  Dre(i,1) = Asig * cos(Fsig*i* 2*pi/NFFT);
  Dim(i,1) = Asig * sin(Fsig*i* 2*pi/NFFT); 
  
%  Dre(i,1) = round(Asig * cos(F0 + (Fsig*i + B*i*i/2) * 2*pi/NFFT) * sin(i * Ffm * pi / NFFT));
%  Dim(i,1) = round(Asig * sin(F0 + (Fsig*i + B*i*i/2) * 2*pi/NFFT) * sin(i * Ffm * pi / NFFT)); 
  
  if (i == 254)
    Dre(i,1) = Asig;
    Dim(i,1) = 0;
  else
    Dre(i,1) = 0;
    Dim(i,1) = 0;   
  end  

end

% Adding noise to real signal 
SNR = 50;
SEED = 1;

DatRe = awgn(Dre, SNR, 0, SEED);     
DatIm = awgn(Dim, SNR, 0, SEED);     

DAT_IN = DatRe + j*DatIm; 


IR0(:,1) = DatRe;
IR0(:,2) = DatIm;


TX(:,1) = real(DAT_IN);
TX(:,2) = imag(DAT_IN);


%% -----------------------------------------------------------------------
%% ---- FFT CALCULATE (ETALON FFT ALGORITHM) -----------------------------
%% -----------------------------------------------------------------------
DFFT0 = fft(DAT_IN);

IX(:,1) = real(DFFT0);
IX(:,2) = imag(DFFT0);


%% -----------------------------------------------------------------------
%% ---- ULTRA LONG FFT ALGORITM FROM 2D FFT + ROTATION -------------------
%% -----------------------------------------------------------------------

% ---- STAGE 0: CREATE MATRIX FOR INPUT DATA ----
for n1 = 1:N1
    for n2 = 1:N2
        DIN(n2, n1) = DAT_IN((n2-1)*N1+n1,1);
    end 
end

for k1 = 1:N1
    for k2 = 1:N2
        DIN_L(1, (k2-1)*N1+k1) = DIN(k1, k2);
    end 
end
IR1(:,1) = real(DIN_L);
IR1(:,2) = imag(DIN_L);

% ---- STAGE 1: CALCULATE FIRST FFT (A) ----
for n1 = 1:N1
    DFFT1(:,n1) = fft(DIN(:,n1));
end

for k1 = 1:N1
    for k2 = 1:N2
        DFFT1_L(1, (k2-1)*N1+k1) = DFFT1(k1, k2);
    end 
end

IR2(:,1) = real(DFFT1_L);
IR2(:,2) = imag(DFFT1_L);

% ---- STAGE 2: ROTATE AFTER SHUFFLER ----
fid = fopen ("WW_ULFFT.dat", "w");
for n1 = 1:N1
    for k2 = 1:N2
        SHW(k2, n1) = exp(-j*2*pi*(n1-1)*(k2-1)/(N1*N2)); 
        SRW(k2, n1) = round((Asig * SHW(k2, n1)));
        fprintf(fid, "%d+i*%d\t", real(SRW(k2, n1)), imag(SRW(k2, n1)));
        % DSH(k2, n1) = DFFT1(k2, n1) * SHW(k2, n1);
        DVH(k2, n1) = DFFT1(k2, n1);
    end
    fprintf(fid, "\n");
end
fclose(fid);


for n1 = 1:N1
    for k2 = 1:N2
        DSH(k2, n1) = DFFT1(k2, n1) * SHW(k2, n1);
    end
end



for k1 = 1:N1
    for k2 = 1:N2
        DCORD(1, (k2-1)*N1+k1) = SHW(k1, k2);
    end 
end

IR4(:,1) = real(DCORD);
IR4(:,2) = imag(DCORD);

DXH = DVH';
for k1 = 1:N1
    for k2 = 1:N2
        DSH_L(1, (k2-1)*N1+k1) = DXH(k1, k2);
    end 
end

IR3(:,1) = real(DSH_L);
IR3(:,2) = imag(DSH_L);

DXH = DSH';
for k1 = 1:N1
    for k2 = 1:N2
        DSH_L(1, (k2-1)*N1+k1) = DXH(k1, k2);
    end 
end

IR5(:,1) = real(DSH_L);
IR5(:,2) = imag(DSH_L);


% ---- STAGE 3: CALCULATE SECOND FFT (B) ----
for n2 = 1:N2
    DFFT2(:,n2) = fft(DSH(n2,:));
end

for k1 = 1:N1
    for k2 = 1:N2
        DFFT2_L(1, (k2-1)*N1+k1) = DFFT2(k1, k2);
    end 
end

IR6(:,1) = real(DFFT2_L);
IR6(:,2) = imag(DFFT2_L);

% ---- STAGE 4: CREATE MATRIX FOR OUTPUT DATA ----
for k1 = 1:N1
    for k2 = 1:N2
        DOUT((k1-1)*N2+k2, 1) = DFFT2(k1, k2);
    end 
end

IR7(:,1) = real(DOUT);
IR7(:,2) = imag(DOUT);

%% -----------------------------------------------------------------------
%% ---- PLOT DATA: INPUT SIGNAL, FFT0, FFT1, TWIDDLES etc ----------------
%% -----------------------------------------------------------------------

figure(1) % Plot loaded data in Time Domain
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
 

%% -----------------------------------------------------------------------
Col = ['k' 'r' 'b' 'g' 'm', 'c', 'y', 'k'];
Val = ['-' '--' ':' '-.'];
%% -----------------------------------------------------------------------
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
%q = 1;   % Offset because MATLAB starts at one.  Set to 0 for C code.
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

