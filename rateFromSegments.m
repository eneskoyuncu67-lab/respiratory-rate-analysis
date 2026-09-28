function [rateBpm, nBreaths] = rateFromSegments(inhaleOnsetTimes, segments)
% Count breaths whose onset falls inside any of the given [start end]
% segments, divide by total segment duration -> breaths per minute.
    nBreaths = 0;
    totalDuration = 0;
    for r = 1:size(segments, 1)
        startT = segments(r, 1);
        endT = segments(r, 2);
        inSeg = inhaleOnsetTimes >= startT & inhaleOnsetTimes <= endT;
        nBreaths = nBreaths + sum(inSeg);
        totalDuration = totalDuration + (endT - startT);
    end
    if totalDuration > 0
        rateBpm = nBreaths / (totalDuration / 60);
    else
        rateBpm = NaN;
    end
end
