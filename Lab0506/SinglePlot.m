function SinglePlot(X, Y1, l1, Xlabel, Ylabel, Title, filename)

    figure;
    plot(X, Y1); hold on;
    legend(l1);
    xlabel(Xlabel);
    ylabel(Ylabel);
    title(Title);
    saveas(gcf, strcat('plots/', filename, '_line.png'));

end