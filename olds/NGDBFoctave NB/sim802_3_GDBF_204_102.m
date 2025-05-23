%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Identifiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
simName  = 'NGDBF_802.3_example';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNRvals           = 4;  % SNR values to simulate
F                 = 1000;      % Number of frame errors to observe
Fchan             = 250;    % Number of channel frames to generate per batch
R                 = 0.5;  % code rate
T                 = 500;    % Max iterations for decoding
w                 = 0.3;   % Syndrome weight parameter
Ymax              = 2.95;    % Channel sample saturation magnitude
theta             = -1.5;   % Flipping threshold
eta               = 1.0;     % Perturbation noise scale parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pth4 = (fullfile(pwd, 'related_variables/alists'));
H_matrix_mat_fl_nm = '204.102.3.6.16';
load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
H=full(h);
H(H ~= 0) = 1;
% H=readalist('802_3_H.alist');
[M, N] = size(H);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GDBFerrors       = zeros(length(SNRvals), 1);
GDBFframeerrors  = zeros(length(SNRvals), 1);
GDBFframeundetected  = zeros(length(SNRvals), 1);
NFrames          = zeros(length(SNRvals), 1);

sdx              = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (sdx <= length(SNRvals))
  SNR       = SNRvals(sdx);
  N0        = (R*10^(SNR/10))^(-1);
  csigma    = sqrt(N0/2);
  nsigma    = eta*csigma;    

  NFrames(sdx)=0;
  keepGoing = 1;
  
  while (keepGoing)
		c        = ones(Fchan,N);                      % All-zero codeword
		noisevec = csigma*randn(Fchan,N);              % Matrix of AWGN noise
		y        = min(max(c + noisevec,-Ymax),Ymax);  % Clipped channel samples

		% Outputs from decoder
		%        d: decisions (bipolar format)
		%   failed: failed check nodes
		%        S: syndrome (bipolar format)
		%       Sb: syndrom (binary format)
		%        E: flip function values 
		[d, failed, S, Sb, E]  = decodeGDBFvec(y, H, N, Fchan, M, T, w, theta, nsigma);
		
		NFrames(sdx)         = NFrames(sdx)         + Fchan;
		GDBFerrors(sdx)      = GDBFerrors(sdx)      + sum(sum(d<0));    % BER (total)
		GDBFframeerrors(sdx) = GDBFframeerrors(sdx) + sum(sum(d'<0)>0);  % FER (total)
		GDBFframeundetected(sdx) = GDBFframeundetected(sdx) + sum(failed>0);  % FER (undetected)
		
		if (GDBFframeerrors(sdx) > F)
			keepGoing = 0;
		end

		fprintf(1,"SNR=%1.2f: %d Frame Errors in %d Frames, FER=%e, BER=%e, undetected FER=%e\n",...
						SNR, GDBFframeerrors(sdx), NFrames(sdx),...
						GDBFframeerrors(sdx)/NFrames(sdx),...
						GDBFerrors(sdx)/(NFrames(sdx)*N),GDBFframeundetected(sdx)/NFrames(sdx));
		

  end
  sdx = sdx+1;
end

