function [constellation_points, gray_labels, I1, Q1] = qam_constellation(M)

    n = sqrt(M);
    
    I = 0:n-1;
    Q = 0:n-1;
    
    I1 = 2*I - n + 1;
    Q1 = 2*Q - n + 1;
    
    [I, Q] = meshgrid(I1, Q1);
    constellation_points = I(:) + 1j*Q(:);
        [~, gray_labels] = comm.internal.utilities.bin2gray(0:(M-1), "QAM", M);
        gray_labels = gray_labels';
    

end

