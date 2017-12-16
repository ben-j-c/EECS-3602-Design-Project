function H = t1_freq(h, w)
    N = length(h);
    a = zeros(1,(N - 1)/2 + 1);
    a(1) = h((N-1)/2+1);
    
    for i = 1:((N-1)/2)
        a(i + 1) = 2*h((N-1)/2 - i + 1);
    end
    
    n = 0:((N-1)/2);
    cosW = cos((n')*(w));
    
    H = exp(-1j.*(w)*(N-1)/2).*(a*cosW);
end

