function [answer, IsCancelled] = request_user_input_plot_estimate()
    % %%%%%%%%%%%%%%%%%%%%
    % get results
    % %%%%%%%%%%%%%%%%%%%%
    dir_res_now = get_subfolders_r('results', 2);

    if dir_res_now.count == 0
		answer      = struct;
		IsCancelled = 1;
		fprintf(1, "Can not find any simulation result.\n");
		return;
    end

    % %%%%%%%%%%%%%%%%%%%%
    % verify the results
    % %%%%%%%%%%%%%%%%%%%%
    dir_now = struct;

    count   = 0;
    name    = {};
    array   = [];
    for i_arr = 1 : dir_res_now.count
        exp_struct = dir_res_now.array(i_arr);
        if exp_struct.count ~= 0
            % ---------- choose a specific extension ----------
            res_exp_struct = struct;

            count1  = 0;
            name1   = {};
            array1  = [];

            for ii_arr = 1 : exp_struct.count
                subdir = fullfile('results', dir_res_now.name{i_arr}, exp_struct.name{ii_arr});
                subdir = fullfile(subdir, "*.est.mat");
                res_files = dir(subdir);
                if ~isempty(res_files)
                    count1          = count1 + 1;
                    name1{end+1}    = exp_struct.name{ii_arr}; %#ok<AGROW> 
                    array1          = [array1; exp_struct.array(ii_arr)]; %#ok<AGROW> 
                end
            end
            res_exp_struct.count    = count1;
            res_exp_struct.name     = name1;
            res_exp_struct.array    = array1;
            % ---------- end ----------
            if count1 ~= 0
                count       = count + 1;
                name{end+1} = dir_res_now.name{i_arr}; %#ok<AGROW> 
                array       = [array; res_exp_struct]; %#ok<AGROW> 
            end
        end
    end
    
    dir_now.count   = count;
    dir_now.name    = name;
    dir_now.array   = array;

    % %%%%%%%%%%%%%%%%%%%%
    % default answer
    % %%%%%%%%%%%%%%%%%%%%
	defAns = struct;
    sim_ind = 1;
	defAns.sim_name = dir_now.name(sim_ind);
    defAns.exp_name = dir_now.array(sim_ind).name(1);

    % %%%%%%%%%%%%%%%%%%%%
    % get MD5 of dir_now
    % %%%%%%%%%%%%%%%%%%%%
    md5         = GetMD5(dir_now, 'Array');
    defVal_path = fullfile('_tmp', 'defaultInput_plot_estimate.mat');
    if isfile(defVal_path)
		defVal      = load(defVal_path).defVal;
        md5_prev    = defVal.md5;

        if strcmp(md5, md5_prev)
            sim_ind = defVal.sim_ind;
			defAns  = defVal.defAns;
        else
            sim_ind = 1;
			defAns.sim_name = dir_now.name(sim_ind);
	        defAns.exp_name = dir_now.array(sim_ind).name(1);
        end
    end

	% ask for user input
	prompt  = {};
	title   = 'Plot setting';
	formats = {};
	options = struct;
	
	prompt(1,:)             = {'Enter simulation data name:', 'sim_name', []};
	formats(1,1).type       = 'list';
	formats(1,1).style      = 'popupmenu';
	formats(1,1).size       = -1;
	formats(1,1).format     = 'text';
	formats(1,1).items      = dir_now.name;
	formats(1,1).callback   = @(~,~,h,k)(set(h(k+1), ...
        'String', dir_now.array(get(h(k),'Value')).name, 'Value', 1));
	
	prompt(2,:)             = {'experiment', 'exp_name',[]};
	formats(2,1).type       = 'list';
	formats(2,1).style      = 'popupmenu';
	formats(2,1).size       = -1;
	formats(2,1).format     = 'text';
	formats(2,1).items      = dir_now.array(sim_ind).name;
	
	options.Resize          = 'off';
	options.AlignControls   = 'on';
	
	[answer, IsCancelled] = inputsdlg(prompt, title, formats, defAns, options);
	if IsCancelled == 0
		answer.sim_name = answer.sim_name{1};
		answer.exp_name = answer.exp_name{1};

        % record defAns
        defVal      = struct;
        defVal.md5  = md5;
        % compute sim_ind
        sim_ind     = find(ismember(dir_now.name, answer.sim_name));
        % save to defVal
        defVal.sim_ind  = sim_ind;
        defVal.defAns   = answer;
        if ~exist('_tmp', 'dir');	mkdir('_tmp'); end
        save(defVal_path, 'defVal');
	end
end

