function result = filter_length(rippleA, rippleB, beta)
    result = (-20.*log10(sqrt(rippleA.*rippleB)) - 13)./(14.6.*beta) + 1;
end