function B = linearize_bands(HIM)
    [xx,yy,ll] = size(HIM);
    num_pixels = xx * yy;
    B = reshape(HIM,num_pixels,ll);
end
