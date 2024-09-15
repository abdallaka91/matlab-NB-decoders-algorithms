function graySymbols = bin2gray(decimalSymbols)

    graySymbols = bitxor(decimalSymbols, bitshift(decimalSymbols, -1));
    
end
