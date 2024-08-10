function oupt_data = segment_save(inpt_data, seg_size)
%SEGMENT_SAVE Given input data, SEGMENT_SAVE will split the last dimension
% of input data into several part with each size = seg_size.

	Mc = size(inpt_data, ndims(inpt_data));
	otherdims = repmat({':'}, 1, ndims(inpt_data) - 1);
	
	% compute size of oupt_data
	if Mc < seg_size
		n_seg = 1;
	else
		n_seg = ceil(Mc / seg_size);
	end
	oupt_data = cell(n_seg, 1);
	
	% start segmentation
	for i_seg = 1 : n_seg
		if i_seg * seg_size <= Mc
			oupt_data{i_seg} = ...
				inpt_data(otherdims{:}, ...
				(i_seg - 1) * seg_size + 1 : i_seg * seg_size);
		else
			oupt_data{i_seg} = ...
				inpt_data(otherdims{:}, ...
				(i_seg - 1) * seg_size + 1 : Mc);
		end
	end
end