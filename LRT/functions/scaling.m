function sc=scaling(ksp)

tmp = dimnorm(ifft2c(ksp),3);

tmpnorm = dimnorm(tmp, 4);

tmpnorm2 = sort(tmpnorm(:), 'ascend');

% match convention used in BART

p100 = tmpnorm2(end);
p90 = tmpnorm2(round(.9 * length(tmpnorm2)));
p50 = tmpnorm2(round(.5 * length(tmpnorm2)));
if (p100 - p90) < 2 * (p90 - p50)
    sc = p90;
else
    sc = p100;
end
fprintf('\nScaling: %f\n\n', sc);
end