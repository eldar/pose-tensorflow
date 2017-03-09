%
% function res = padZeros(str, npad)
%
function res = padZeros(str, npad)

  n = length(str);
  
  assert(n <= npad);

  res = [repmat('0', 1, npad - n) str];

