function edgecurves = extractedge(inpic,scale,threshold,shape)

if (nargin < 4)
    shape = 'same';
end

Lv_inpic = Lv(discgaussfft(inpic,scale),shape);
Lvv_inpic = Lvvtilde(discgaussfft(inpic,scale),shape);
Lvvv_inpic = Lvvtilde(discgaussfft(inpic,scale),shape);

%Need to threshold Lv and find points of Lvvv < 0.
Lv_thresh = (Lv_inpic > threshold) - 0.5; 

%Based on the sign of maskpic1 -> if 0 we need to make it -0.5
Lvvv_thresh = (Lvvv_inpic < 0) - 0.5; 

% Extracts level curves in a given image and rejects points based on the sign of the second input argument
% zerocrosscurves(zeropic, maskpic1)
edgecurves = zerocrosscurves(Lvv_inpic,Lvvv_thresh);

% Thresholds these curves with respect to the sign of another image.
% Returns a new set of curves containing only those points in the polygons for which the mask value is non-negative.
% thresholdcurves(curves, maskpic2)
edgecurves = thresholdcurves(edgecurves, Lv_thresh);

end