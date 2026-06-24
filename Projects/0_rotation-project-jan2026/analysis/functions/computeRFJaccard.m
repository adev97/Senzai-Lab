function J = computeRFJaccard(maskA, maskB)

% Ensure logical
maskA = logical(maskA);
maskB = logical(maskB);

intersection = sum(maskA(:) & maskB(:));
unionArea    = sum(maskA(:) | maskB(:));

if unionArea == 0
    J = NaN;   % no RF detected in either unit
else
    J = intersection / unionArea;
end
