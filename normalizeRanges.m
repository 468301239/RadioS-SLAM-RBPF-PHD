function normalizedRanges = normalizeRanges(ranges)
    normalizedRanges=ranges;
    normalizedRanges(normalizedRanges<0)=0;
end