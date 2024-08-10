function [answer, IsCancelled] = request_user_input_plot_crb()
    % %%%%%%%%%%%%%%%%%%%%
    % get results
    % %%%%%%%%%%%%%%%%%%%%
    dir_res_now_est = get_subfolders_r('results', 2);

    if dir_res_now_est.count == 0
		answer      = struct;
		IsCancelled = 1;
		fprintf(1, "Can not find any simulation result.\n");
		return;
    end

    dir_res_now_crb = get_subfolders_r('results_crb', 2);
    if dir_res_now_crb.count == 0
		answer      = struct;
		IsCancelled = 1;
		fprintf(1, "Can not find any simulation result.\n");
		return;
    end

    % %%%%%%%%%%%%%%%%%%%%
    % verify the results
    % %%%%%%%%%%%%%%%%%%%%
    dir_now_crb = struct;
    dir_now_est = struct;

    count_crb   = 0;
    name_crb    = {};
    array_crb   = [];

    count_est   = 0;
    name_est    = {};
    array_est   = [];

    for i_arr = 1 : dir_res_now_crb.count
        exp_struct = dir_res_now_crb.array(i_arr);
        if exp_struct.count ~= 0
            % ---------- choose a specific extension ----------
            res_exp_struct = struct;

            count1  = 0;
            name1   = {};
            array1  = [];

            for ii_arr = 1 : exp_struct.count
                subdir = fullfile('results_crb', dir_res_now_crb.name{i_arr}, exp_struct.name{ii_arr});
                subdir = fullfile(subdir, "*.crb.mat");
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
                count_crb       = count_crb + 1;
                name_crb{end+1} = dir_res_now_crb.name{i_arr}; %#ok<AGROW> 
                array_crb       = [array_crb; res_exp_struct]; %#ok<AGROW> 
            end
        end

        exp_struct = dir_res_now_est.array(i_arr);
        if exp_struct.count ~= 0
            % ---------- choose a specific extension ----------
            res_exp_struct = struct;

            count1  = 0;
            name1   = {};
            array1  = [];

            for ii_arr = 1 : exp_struct.count
                subdir = fullfile('results', dir_res_now_est.name{i_arr}, exp_struct.name{ii_arr});
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
                count_est       = count_est + 1;
                name_est{end+1} = dir_res_now_est.name{i_arr}; %#ok<AGROW> 
                array_est       = [array_est; res_exp_struct]; %#ok<AGROW> 
            end
        end
    end
    
    dir_now_crb.count   = count_crb;
    dir_now_crb.name    = name_crb;
    dir_now_crb.array   = array_crb;
    dir_now_est.count   = count_est;
    dir_now_est.name    = name_est;
    dir_now_est.array   = array_est;

    % %%%%%%%%%%%%%%%%%%%%
    % default answer
    % %%%%%%%%%%%%%%%%%%%%
	defAns = struct;
    sim_ind_crb = 1;
    sim_ind_est = 1;
	defAns.sim_name_crb     = dir_now_crb.name(sim_ind_crb);
    defAns.exp_name_crb     = dir_now_crb.array(sim_ind_crb).name(1);
	defAns.sim_name_est     = dir_now_est.name(sim_ind_est);
    defAns.exp_name_est     = dir_now_est.array(sim_ind_est).name(1);
    defAns.dataMetric   = 'DFO';

    % %%%%%%%%%%%%%%%%%%%%
    % get MD5 of dir_now
    % %%%%%%%%%%%%%%%%%%%%
    md5         = GetMD5(dir_now_est, 'Array');
    defVal_path = fullfile('_tmp', 'defaultInput_plot_crb.mat');
    if isfile(defVal_path)
		defVal      = load(defVal_path).defVal;
        md5_prev    = defVal.md5;

        if strcmp(md5, md5_prev)
            sim_ind_crb = defVal.sim_ind_crb;
            sim_ind_est = defVal.sim_ind_est;
			defAns  = defVal.defAns;
        else
            sim_ind_crb = 1;
            sim_ind_est = 1;
	        defAns.sim_name_crb     = dir_now_crb.name(sim_ind_crb);
            defAns.exp_name_crb     = dir_now_crb.array(sim_ind_crb).name(1);
	        defAns.sim_name_est     = dir_now_est.name(sim_ind_est);
            defAns.exp_name_est     = dir_now_est.array(sim_ind_est).name(1);
            defAns.dataMetric   = 'DFO';
        end
    end

	% ask for user input
	prompt  = {};
	title   = 'Plot setting';
	formats = {};
	options = struct;
	
	prompt(1,:)             = {'Enter CRB data name:', 'sim_name_crb', []};
	formats(1,1).type       = 'list';
	formats(1,1).style      = 'popupmenu';
	formats(1,1).size       = -1;
	formats(1,1).format     = 'text';
	formats(1,1).items      = dir_now_crb.name;
	formats(1,1).callback   = @(~,~,h,k)(set(h(k+1), ...
        'String', dir_now_crb.array(get(h(k),'Value')).name, 'Value', 1));
	
	prompt(2,:)             = {'CRB experiment', 'exp_name_crb',[]};
	formats(2,1).type       = 'list';
	formats(2,1).style      = 'popupmenu';
	formats(2,1).size       = -1;
	formats(2,1).format     = 'text';
	formats(2,1).items      = dir_now_crb.array(sim_ind_crb).name;

    prompt(3,:)             = {'Enter estimation data name:', 'sim_name_est', []};
	formats(3,1).type       = 'list';
	formats(3,1).style      = 'popupmenu';
	formats(3,1).size       = -1;
	formats(3,1).format     = 'text';
	formats(3,1).items      = dir_now_est.name;
	formats(3,1).callback   = @(~,~,h,k)(set(h(k+1), ...
        'String', dir_now_est.array(get(h(k),'Value')).name, 'Value', 1));
	
	prompt(4,:)             = {'estimation experiment', 'exp_name_est',[]};
	formats(4,1).type       = 'list';
	formats(4,1).style      = 'popupmenu';
	formats(4,1).size       = -1;
	formats(4,1).format     = 'text';
	formats(4,1).items      = dir_now_est.array(sim_ind_est).name;

    prompt(5,:) = {'Enter data type: ', 'dataMetric', []};
    formats(5,1).type = 'list';
    formats(5,1).style = 'radiobutton';
	formats(5,1).size = -1;
    formats(5,1).format = 'text';
    formats(5,1).items = {'DFO' 'AOA'};
	
	options.Resize          = 'off';
	options.AlignControls   = 'on';
	
	[answer, IsCancelled] = inputsdlg(prompt, title, formats, defAns, options);
    if IsCancelled == 0
		answer.sim_name_crb = answer.sim_name_crb{1};
		answer.sim_name_est = answer.sim_name_est{1};

        % record defAns
        defVal      = struct;
        defVal.md5  = md5;
        % compute sim_ind
        sim_ind_crb     = find(ismember(dir_now_crb.name, answer.sim_name_crb));
        sim_ind_est     = find(ismember(dir_now_est.name, answer.sim_name_est));
        % save to defVal
        defVal.sim_ind_crb  = sim_ind_crb;
        defVal.sim_ind_est  = sim_ind_est;
        defVal.defAns   = answer;
        if ~exist('_tmp', 'dir');	mkdir('_tmp'); end
        save(defVal_path, 'defVal');
    end
end

