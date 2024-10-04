%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Identifiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simName  = 'NGDBF_Tanner_code_155_example';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



F                 = 1;      % Number of frame errors to observe
Fchan             = 1;    % Number of channel frames to generate per batch
R                 = 62/155;  % 62 or 64 info bits?
T                 = 100;     % Max iterations for decoding
w                 = 0.5;     % Syndrome weight parameter
Ymax              = 2.5;     % Channel sample saturation magnitude
theta             = -0.5;    % Flipping threshold
eta               = 0.9;     % Perturbation noise scale parameter

SNRvals           = 3:0.2:6;  % SNR values to simulate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=readalist('Tanner_155_64_uni-kl.alist');
[M, N] = size(H);
H = binGaussElim(H);
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
        GDBFframeerrors(sdx) = GDBFframeerrors(sdx) + sum(sum(d<0,2)>0);  % FER (total)
        GDBFframeundetected(sdx) = GDBFframeundetected(sdx) + sum(failed>0);  % FER (undetected)

        if (GDBFframeerrors(sdx) > F)
            keepGoing = 0;
        end


    end

    semilogy(SNRvals, GDBFframeerrors./NFrames, '-o')
    grid on
    xlim([SNRvals(1) SNRvals(end)])
    ylim([1e-5 1])
    xlabel("ebn0 (dB)")
    ylabel("FER")
    pause(0.1)
    fprintf(1,"SNR=%1.2f: %d Frame Errors in %d Frames, FER=%e, BER=%e, undetected FER=%e\n",...
        SNR, GDBFframeerrors(sdx), NFrames(sdx),...
        GDBFframeerrors(sdx)/NFrames(sdx),...
        GDBFerrors(sdx)/(NFrames(sdx)*N),GDBFframeundetected(sdx)/NFrames(sdx));



    sdx = sdx + 1;


end

