function SurfPlot(X, Y, Z, Xlabel, Ylabel, Zlabel, Title, filename)

    figure;
    surf(X, Y, Z);
    xlabel(Xlabel);
    ylabel(Ylabel);
    zlabel(Zlabel);
    title(Title);
    saveas(gcf,['plots/' filename '_surf.png'])

end
