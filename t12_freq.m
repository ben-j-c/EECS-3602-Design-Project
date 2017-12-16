function H = t12_freq(h, w)
    N = length(h);
    
    if(bitand(N,1))
        a = zeros(1,(N - 1)/2 + 1);
        a(1) = h((N-1)/2+1);
        
        for i = 1:((N-1)/2)
            a(i + 1) = 2*h((N-1)/2 - i + 1);
        end
        
        n = 0:((N-1)/2);
        cosW = cos((n')*(w));
        
        H = exp(-1j.*(w)*(N-1)/2).*(a*cosW);
    else
        b = zeros(1,N/2);
        for i = 1:(N/2)
            b(i) = 2*h(N/2 - i + 1);
        end
        
        n = 1:(N/2);
        cosW = cos((n' - 0.5)*(w));
        
        H = exp(-1j.*(w)*(N-1)/2).*(b*cosW);
    end
end

