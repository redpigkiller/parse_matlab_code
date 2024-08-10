function percent = parfor_progress(N)
%PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
%   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   PARFOR_PROGRESS(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   PARFOR_PROGRESS updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   PARFOR_PROGRESS(0) deletes parfor_progress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      parfor_progress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress;
%      end
%      parfor_progress(0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/

	narginchk(0, 1);
	
	if nargin < 1
    	N = -1;
	end
	
	percent = 0;
	% Width of progress bar
	w = 70;
	
	if N > 0
		if ~exist('_tmp', 'dir')
			mkdir('_tmp');
		end

    	f = fopen('_tmp/parfor_progress.txt', 'w');
		if f < 0
        	error('Do you have write permissions for %s?', pwd);
		end
		
		% Save N at the top of progress.txt
    	fprintf(f, '%d\n', N);
    	fclose(f);
    	
    	if nargout == 0
        	disp(['  0.00%[>', repmat(' ', 1, w), ']']);
    	end
	elseif N == 0
    	delete('_tmp/parfor_progress.txt');
    	percent = 100;
    	
    	if nargout == 0
        	disp([repmat(char(8), 1, (w+12)), newline, '100.00%[', repmat('=', 1, w+1), ']']);
    	end
	else
    	if ~exist('_tmp/parfor_progress.txt', 'file')
        	error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
    	end
    	
    	f = fopen('_tmp/parfor_progress.txt', 'a');
    	fprintf(f, '1\n');
    	fclose(f);
    	
    	f = fopen('_tmp/parfor_progress.txt', 'r');
    	progress = fscanf(f, '%d');
    	fclose(f);
    	percent = (length(progress)-1)/progress(1)*100;
    	
    	if nargout == 0
			% 7 characters wide, (including percentage mark)
        	perc = sprintf('%6.2f%%', percent);
			len = round(percent*w/100);
        	disp([repmat(char(8), 1, (w+12)), newline, perc, ...
				'[', repmat('=', 1, len), '>', ...
				repmat(' ', 1, w - len), ']']);
    	end
	end
end