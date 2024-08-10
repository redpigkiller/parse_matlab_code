function fig = snr_plot(id, snr, yval, ylgd, x_text, y_text, t_text, plot_settings)
%SNR_PLOT This function plot the figure for multiple data along snr
%   id:     figure id
%   snr:    snr array
%   yval:   value for each snr, different data should be located on the
%           first dimension
%   ylgd:   the display name of the different data
%   x_text: xlabel
%   y_text: ylabel
%   t_text: title
%   plot_settings:  struct for plot settings
%       is_ylog:    0/1
%       is_grid:    0/1
%       color:      color map
%       marker:     marker
%       lines:      line styles
%       lgd_loc:    location of legend (0 for no legend)
%

if nargin < 8
    is_ylog = 1;
    is_grid = 1;
    cmap = ['b', 'k', 'r', 'g', 'c', 'm', 'y'];
    makr = ['o', '+', '*', '.', 'x', '-'];
    lins = ["-", "--", ":", "-."];
    lgd_loc = {'Location', 'southwest'};
else
    is_ylog = plot_settings.is_ylog;
    is_grid = plot_settings.is_grid;
    cmap = plot_settings.color;
    makr = plot_settings.marker;
    lins = plot_settings.lines;
    lgd_loc = plot_settings.lgd_loc;
end
if nargin < 7
    t_text = sprintf('value of %s', inputname(3));
end
if nargin < 6
    y_text = inputname(3);
end
if nargin < 5
    x_text = inputname(2);
end

if size(yval, 1) ~= length(ylgd)
    error("different number of data y");
end

fig = figure(id);
for q = 1 : length(ylgd)
    if is_ylog == 1
	    semilogy(snr, yval(q, :), ...
		    'DisplayName', [inputname(4), ' = ', num2str(ylgd(q))], ...
		    'Color', cmap(q), ...
		    'Marker', makr(1), ...
		    'LineStyle', lins(1));
    else
        plot(snr, yval(q, :), ...
		    'DisplayName', [inputname(4), ' = ', num2str(ylgd(q))], ...
		    'Color', cmap(q), ...
		    'Marker', makr(1), ...
		    'LineStyle', lins(1));
    end
    hold on;
end
hold off;
xlabel(x_text);
ylabel(y_text);
title(t_text);

if is_grid == 1
    grid on;
end

if iscell(lgd_loc)
    legend(lgd_loc{:});
end

end

