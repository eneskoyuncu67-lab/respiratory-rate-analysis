function plotViolin(xCenter, values, faceColor)
% Draws one violin shape centered at xCenter, using a kernel density
% estimate (ksdensity) mirrored left/right, scaled so its max half-width
% is 0.35 (so it doesn't overlap the neighboring condition's violin).
    if numel(values) < 2
        return
    end
    [f, xi] = ksdensity(values);
    f = f / max(f) * 0.35;
    fill(xCenter + [f, -fliplr(f)], [xi, fliplr(xi)], faceColor, ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none');
end
