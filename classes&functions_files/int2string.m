function result = int2string(n)
  n_digits = floor(log10(n)+1);
  for k = 1:n_digits
    result(n_digits-k+1) = char('0'+ mod(n,10));
    n = floor(n/10);
  endfor
endfunction