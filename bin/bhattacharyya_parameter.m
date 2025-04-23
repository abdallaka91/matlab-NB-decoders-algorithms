function Z = bhattacharyya_parameter(snr)
    % Here we use a simplified version of the Bhattacharyya parameter.
    % You can use a more detailed formula depending on the channel model.
    % For example, in an AWGN channel:
    
    % Assuming a simple AWGN channel model, Bhattacharyya parameter Z can be
    % approximated as Z = 0.5 * erfc(sqrt(snr) / 2), where erfc is the complementary error function.
    
    Z = 0.5 * erfc(sqrt(snr) / 2);
end