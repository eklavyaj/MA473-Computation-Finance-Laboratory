function LinePlot(X, Y, Xlabel, Ylabel, Title, filename)

    figure;
    plot(X, Y);
    xlabel(Xlabel);
    ylabel(Ylabel);
    title(Title);
    saveas(gcf, ['plots/' filename '_line.png']);

end