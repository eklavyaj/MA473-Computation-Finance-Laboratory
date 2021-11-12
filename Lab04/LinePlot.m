function LinePlot(X, Y1, Y2, l1, l2, Xlabel, Ylabel, Title, filename)

    figure;
    plot(X, Y1); hold on;
    plot(X, Y2); hold off;
    legend(l1, l2);
    xlabel(Xlabel);
    ylabel(Ylabel);
    title(Title);
    saveas(gcf, strcat('plots/', filename, '_line.png'));

end