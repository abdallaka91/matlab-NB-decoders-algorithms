function constellation = getQAMConstellation(M)
    % Generate message symbols
    messageSymbols = (0:M-1).';
    
    % Modulate using QAM
    constellation = qammod(messageSymbols, M, 'UnitAveragePower', true);
end
