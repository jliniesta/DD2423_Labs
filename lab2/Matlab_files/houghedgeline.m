function [linepar, acc] = houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, verbose)

magnitude = Lv(discgaussfft(pic,scale),'same');

curves = extractedge(pic,scale,gradmagnthreshold,'same');

[linepar, acc] = houghline(curves, magnitude, nrho, ntheta, gradmagnthreshold, nlines, verbose);

figure('Name','Hough transform','NumberTitle','off');
overlaycurves(pic, linepar);
axis([1 size(pic, 2) 1 size(pic, 1)]);    
title('Hough transform')

figure('Name','Hough transform - Accumulator peaks','NumberTitle','off');
subplot(1, 2, 1);
    showgrey(acc)
    title('Accumulator')  
subplot(1, 2, 2);
    overlaycurves(pic, linepar);
    axis([1 size(pic, 2) 1 size(pic, 1)]);    
    title('Hough transform')



end