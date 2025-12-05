function x = CRT(remainders, moduli)
    % Check lengths
    assert(length(remainders) == length(moduli), 'Remainders and moduli must be same length');
    
    N = prod(moduli);         % product of all moduli
    x = 0;
    
    for i = 1:length(moduli)
        ni = moduli(i);
        Ri = N / ni;          % partial product
        % Modular inverse of Ri modulo ni
        [~, inv] = gcd(Ri, ni); 
        inv = mod(inv, ni);   % make sure inverse is positive
        x = x + remainders(i) * Ri * inv;
    end
    
    x = mod(x, N);  % final result modulo product of moduli
end