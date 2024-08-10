function [answer, IsCancelled] = request_user_input_sim_estimate(default_Mc, default_exp_name)
    % %%%%%%%%%%%%%%%%%%%%
    % get gen_cases
    % %%%%%%%%%%%%%%%%%%%%
    dir_now = get_subfolders_r('_gendata', 2);

    if dir_now.count == 0
		answer      = struct;
		IsCancelled = 1;
		fprintf(1, "Can not find any simulation source.\n" + ...
            "Please execute gen_cases.m first!\n");
		return;
    end

    % %%%%%%%%%%%%%%%%%%%%
    % default answer
    % %%%%%%%%%%%%%%%%%%%%
	defAns                      = struct;
	defAns.sim_src_name         = dir_now.name(1);
	defAns.Mc                   = default_Mc;
	defAns.exp_name             = default_exp_name;
	defAns.flag_collect_mode    = false;
	defAns.flag_skip_simulated  = true;
	defAns.flag_use_parfor      = true;
	defAns.flag_use_gpu         = false;

    % %%%%%%%%%%%%%%%%%%%%
    % get MD5 of dir_now
    % %%%%%%%%%%%%%%%%%%%%
    md5 = GetMD5(dir_now, 'Array');

    % try to load last user input if nothing has changed
	defVal_path = fullfile('_tmp', 'defaultInput_sim_estimate_channel.mat');
    if isfile(defVal_path)
		defVal      = load(defVal_path).defVal;
        md5_prev    = defVal.md5;

        if strcmp(md5, md5_prev)
            defAns = defVal.defAns;
        else
            defAns.exp_name             = defVal.defAns.exp_name;
			defAns.flag_collect_mode    = defVal.defAns.flag_collect_mode;
			defAns.flag_skip_simulated  = defVal.defAns.flag_skip_simulated;
			defAns.flag_use_parfor      = defVal.defAns.flag_use_parfor;
			defAns.flag_use_gpu         = defVal.defAns.flag_use_gpu;
        end
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%% ask for user input %%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	prompt  = {};
	title   = 'Simulation setting';
	formats = {};
	options = struct;
	
	prompt(1,:)         = {'Enter simulation data source:', 'sim_src_name', []};
	formats(1,1).type   = 'list';
	formats(1,1).style  = 'popupmenu';
	formats(1,1).size   = -1;
	formats(1,1).format = 'text';
	formats(1,1).items  = dir_now.name;
	formats(1,1).span   = [1 4];

	prompt(2,:)         = {'number of simulation trials', 'Mc',[]};
	formats(2,1).type   = 'edit';
	formats(2,1).size   = -1;
	formats(2,1).format = 'integer';
	formats(2,1).span   = [1 4];
	
	prompt(3,:)         = {'experiment', 'exp_name',[]};
	formats(3,1).type   = 'edit';
	formats(3,1).size   = -1;
	formats(3,1).format = 'text';
	formats(3,1).span   = [1 4];

    prompt(4,:)         = {'collect mode', 'flag_collect_mode',[]};
	formats(4,1).type   = 'check';

	prompt(5,:)         = {'skip simulated', 'flag_skip_simulated',[]};
	formats(4,2).type   = 'check';

	prompt(6,:)         = {'use parfor', 'flag_use_parfor',[]};
	formats(4,3).type   = 'check';
	
	prompt(7,:)         = {'use GPU', 'flag_use_gpu',[]};
	formats(4,4).type   = 'check';

	options.Resize          = 'off';
	options.AlignControls   = 'on';
	
	[answer, IsCancelled] = inputsdlg(prompt, title, formats, defAns, options);
    if IsCancelled == 0
		answer.sim_src_name = answer.sim_src_name{1};

        % record defAns
        defVal          = struct;
        defVal.md5      = md5;
        defVal.defAns   = answer;

        if ~exist('_tmp', 'dir');	mkdir('_tmp'); end
        save(defVal_path, 'defVal');
    end
end

