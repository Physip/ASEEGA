%    ___                         ______       __     __  ____              ______               __  
%   / _ | ___ ___ ___ ___ ____ _/ __/ /____ _/ /_   /  |/  (_)__________  / __/ /_______ ______/ /_ 
%  / __ |(_-</ -_) -_) _ `/ _ `/\ \/ __/ _ `/ __/  / /|_/ / / __/ __/ _ \_\ \/ __/ __/ // / __/ __/ 
% /_/ |_/___/\__/\__/\_, /\_,_/___/\__/\_,_/\__/__/_/  /_/_/\__/_/  \___/___/\__/_/  \_,_/\__/\__/  
%                   /___/                     /___/                                                 
% 
%
%      Name:	AseegaStat_MicroStruct.m
%   Purpose:	Analysis of sleep microstrutures
%   Comment: 
%
%        In:	aseega_analysis, structure array provided by Aseega EEG Analysis
%		flag_plot_fig,  1 = plot figures, 0 = no figure
%		flag_Excel_ouput, 1 = results in Excel file,  0 = no output file
%		flag_visual_scoring, 1 = use visual scoring, 0 = use Aseega scoring
%
%       Out:	Nb_events_per_epok, number of events detected, for each 30s-epoch
%
% Called by: 
%     Calls:	AseegaPlot_Hypno.m
%
%   Example:	load('GUEA-C4O2.AseegaAnalysis.A4R.mat')	% loads the Aseega structure result of the recording GUEA-C4O2.edf,
%								% this can also be achieved by a drag and drop of the matfile in the command tool
%		AseegaStat_MicroStruct(aseega_analysis,1,0,0);	% grab the results of interest, according to the Parameters section, and plot the results
%
%
%   History:	v1.0   24-02-2014 Creation
%		v2.0   25-04-2016 Add EEGWatch plot
%		v2.1   13-09-2016 CommandTool print debug    
%		v2.2   05-10-2016 User friendly syntax
%		v2.3   27-02-2017 Add Excel ouput file
%		v2.4   22-12-2017 can use visual scoring, if provided
%
%  Physip, Paris France
%  2002-2018




function [Nb_events_per_epok] = AseegaStat_MicroStruct(varargin)

me.version = '2.4';
me.date    = 'Physip, Dec 22nd, 2017';
fname = 'AseegaStat_MicroStruct';



%       ___                          __             
%      / _ \___ ________ ___ _  ___ / /____ _______ 
%     / ___/ _ `/ __/ _ `/  ' \/ -_) __/ -_) __(_-< 
%    /_/   \_,_/_/  \_,_/_/_/_/\__/\__/\__/_/ /___/ 
%                                                   
% 
% Comment/uncomment lines to change the parameters
% ===================================================

%% Choice of sleep microstructure:
type_microstructure = 'spindles';	%    - spindles
%type_microstructure = 'alpha_bursts';	%    - alpha_bursts

%% Sleep stages of interest (SOI), {'N2'}, {'N2';'N3'} ...
sleep_stages = {'N2'};			

%% Plot option
flag_SOI_restricted = 0;		% Plot for the whole revording
%flag_SOI_restricted = 1;		% Plot restricted to the events detected in SOI
	









					
%% - - - - - - - - - - - - DO NOT CROSS!  - - - - - - - - - - - - - - - - - -
%
% - - - - - - - - FOR EXPERIMENTED MATLAB USERS ONLY - - - - - - - - - - - -


if nargin < 4 || nargin > 4
	disp(['   Error using ',fname,'.m, incorrect number of input arguments.'])
	disp('   The correct syntax is:')
	disp('   >> AseegaStat_MicroStruct(aseega_analysis,flag_plot_fig,flag_Excel_ouput,flag_visual_scoring)')
	disp('   Thanks')
	disp(' ')
	return
else
	aseega_analysis		= varargin{1};
	flag_plot_fig		= varargin{2};
	flag_Excel_ouput	= varargin{3};
	flag_visual_scoring	= varargin{4};
end


%       ___             __         _    
%      / _ | ___  ___ _/ /_ _____ (_)__ 
%     / __ |/ _ \/ _ `/ / // (_-</ (_-< 
%    /_/ |_/_//_/\_,_/_/\_, /___/_/___/ 
%                      /___/            
%
% ==================================================
[~,name_subject,~] = fileparts(aseega_analysis.info.recording.filename);
name_subject_fig = strrep(name_subject,'_','\_');		% remplace  underscore by '\_' for display
event_results = aseega_analysis.micro.(type_microstructure);	% sub-structure containing the event analysis
event_name = strrep(type_microstructure,'_','\_');
FOI = event_results.adapted_frequency_band;			% frequency band in which the events are detected


% Index of epochs for which the subject is in the sleep stage of interest (SOI)
% =============================================================================
if flag_visual_scoring
	hypno_5stages = aseega_analysis.scoring.hypno_visual_1;	% Current five-state VISUAL hypnogram: W,R,N1,N2,N3
else
	hypno_5stages = aseega_analysis.scoring.hypno5;		% Current five-state ASEEGA hypnogram: W,R,N1,N2,N3
end
if isempty(hypno_5stages)
	if flag_visual_scoring
		disp('Sorry, visual scoring not found')
	else
		disp('Sorry, Aseega scoring not found')
	end
	return
end

[Nb_SOI, ~] = size(sleep_stages);			% Number of stages of interest

ptr_W = [];ptr_R = [];ptr_N1 = [];ptr_N2 = [];ptr_N3 = [];% Initialization
for st = 1: Nb_SOI
	if strcmp(sleep_stages(st),'W')	
		ptr_W = find(hypno_5stages == 10);		% W epochs pointer. W <-> 10 in the hypno5 case (cf. Aseega-Analysis-Format.pdf)
	elseif strcmp(sleep_stages(st),'R')
		ptr_R = find(hypno_5stages == 7);		% R epochs pointer. R <-> 7 in the hypno5 case (cf. Aseega-Analysis-Format.pdf)
	elseif strcmp(sleep_stages(st),'N1')
		ptr_N1 = find(hypno_5stages == 8);		% N1 epochs pointer. N1 <-> 8 in the hypno5 case (cf. Aseega-Analysis-Format.pdf)
	elseif strcmp(sleep_stages(st),'N2')
		ptr_N2 = find(hypno_5stages == 5);		% N2 epochs pointer. N2 <-> 5 in the hypno5 case (cf. Aseega-Analysis-Format.pdf)
	elseif strcmp(sleep_stages(st),'N3')
		ptr_N3 = find(hypno_5stages == 1);		% N3 epochs pointer. N3 <-> 1 in the hypno5 case  (cf. Aseega-Analysis-Format.pdf)
	end
end
ptr_SOI_epk = cat(2,ptr_W,ptr_R,ptr_N1,ptr_N2,ptr_N3);	% Index of epochs in SleepStage of Interest (SOI)
legend_SOI = '';
for i = 1 : Nb_SOI;
	legend_SOI = sprintf('%s %s',legend_SOI,sleep_stages{i});
end


% Identification of detected events in SOI
% ========================================
epok_lgth = aseega_analysis.info.analysis.nsecs_1epoch;
loc_events_epk = ceil(event_results.position/epok_lgth);		% localization of each detected event, in epoch
Ind_events_SOI = ismember(loc_events_epk,ptr_SOI_epk);			% localization among all detected events of the ones detected in SOI


% Event counting, whole file
% ==========================
Nb_events		= event_results.number;				% Nb of events detected in the recording
Nb_events_per_epok	= event_results.nbr_per_epoch;			% Nb of events per 30s-epoch
Nb_events_per_epok_of_SOI= event_results.nbr_per_epoch(ptr_SOI_epk);	% Nb of events per epoch spent in SOI
Nb_events_SOI 		= sum(Nb_events_per_epok_of_SOI);		% Nb of events detected in SOI


% Event counting, by sleep cycle
% Initialization to NaN to detect anomalies, there must remain no NaN at the end of the loop
% ==========================================================================================
aseega_sleep_cycle = aseega_analysis.sleepparams.sleep_onset_3S.sleep_cycles;	% sub-structure containing the sleep cycles analysis, using the "3S" sleep onset definition

Nb_sleep_cycles			= length(aseega_sleep_cycle.cycle);	% Nb of automatically detected sleep cycles
Nb_epochs_SOI_per_cycle  	= ones(Nb_sleep_cycles,1)*NaN;		% Nb of epochs of SOI in each cycle		
Nb_event_per_cycle		= ones(Nb_sleep_cycles,1)*NaN;		% Nb of events detected in each cycle	
Nb_event_per_cycle_of_SOI	= ones(Nb_sleep_cycles,1)*NaN;		
mat_Ind_events_cycle_SOI	= ones(Nb_sleep_cycles,Nb_events)*NaN;
Nb_events_cycle_per_epok_of_SOI	= ones(Nb_sleep_cycles,2)*NaN;		% Nb of events per cycle per epoch spent in SOI [µ,std]

for c = 1 : Nb_sleep_cycles
	
	ptr_epk_cur_cycle = (aseega_sleep_cycle.cycle{c}.start_e : aseega_sleep_cycle.cycle{c}.end_e);	% Index of epochs in cycle #'c'
	Ind_events_cur_cycle = ismember(loc_events_epk,ptr_epk_cur_cycle);				% localization among all events detected of the ones detected in this cycle
	
	ptr_epk_cur_cycle_SOI = intersect(ptr_epk_cur_cycle,ptr_SOI_epk);				% Index of epochs in cycle #'c' of sleep stage SOI
	mat_Ind_events_cycle_SOI(c,:) = ismember(loc_events_epk,ptr_epk_cur_cycle_SOI);			% localization among all events detected in this cycle of the ones detected in SOI

	Nb_epochs_SOI_per_cycle(c) = length(ptr_epk_cur_cycle_SOI);					% Nb of epochs of SOI in the current cycle
	Nb_event_per_cycle(c) = sum(Ind_events_cur_cycle);						% Nb of spindles detected in the current cycle, whatever the sleep stage
	Nb_event_per_cycle_of_SOI(c) = sum(mat_Ind_events_cycle_SOI(c,:));				% Nb of spindles detected in the current cycle, in SOI	

	Nb_events_cycle_per_epok_of_SOI(c,1) = mean(event_results.nbr_per_epoch(ptr_epk_cur_cycle_SOI));
	Nb_events_cycle_per_epok_of_SOI(c,2) = std(event_results.nbr_per_epoch(ptr_epk_cur_cycle_SOI));
end






%          _____                              __  ______          __  ___  _          __          
%         / ___/__  __ _  __ _  ___ ____  ___/ / /_  __/__  ___  / / / _ \(_)__ ___  / /__ ___ __ 
%        / /__/ _ \/  ' \/  ' \/ _ `/ _ \/ _  /   / / / _ \/ _ \/ / / // / (_-</ _ \/ / _ `/ // / 
%        \___/\___/_/_/_/_/_/_/\_,_/_//_/\_,_/   /_/  \___/\___/_/ /____/_/___/ .__/_/\_,_/\_, /  
%                                                                            /_/          /___/   

disp(' ');
fprintf('%s: version %s (%s)',fname,me.version,me.date)
disp(' ');disp(' ');
disp('________________________________________________________________________________________')
fprintf('Recording: %s, %s detection in frequency band [%0.2g-%0.2g] Hz\n',name_subject,event_name,FOI(1),FOI(2))
disp('________________________________________________________________________________________')

disp(' ')
disp('-----------------------------')
disp('* Analysis for the whole file')
disp('-----------------------------')
for iter = 1 : 2
	switch iter
		% For all the epochs
		case 1;	Nb_events_per_epok_iter = Nb_events_per_epok;
			Nb_events_iter = Nb_events;
			Ind_events_OK_iter = logical(ones(Nb_events,1));
			sleepstage_choice = '';
			fprintf('Considering all the %0.4g epochs of the recording, whatever the sleep stage:\n',length(hypno_5stages))
		% Only for the epochs of SOI
		case 2; Nb_events_per_epok_iter = Nb_events_per_epok_of_SOI;
			Nb_events_iter = Nb_events_SOI;
			Ind_events_OK_iter = Ind_events_SOI;
			sleepstage_choice = sprintf(' in stage%s',legend_SOI);
			fprintf('Considering only the %0.4g epochs labelled%s:\n',length(ptr_SOI_epk),legend_SOI)
	end
	fprintf('  - %0.5g %s detected%s \n',Nb_events_iter,event_name,sleepstage_choice)
	fprintf('  - number/epoch%s: mean = %0.2f/epoch, std = %0.2g \n',sleepstage_choice,mean(Nb_events_per_epok_iter),std(Nb_events_per_epok_iter))
	fprintf('  - duration%s: mean = %0.2f s, std = %0.2g \n',sleepstage_choice,mean(event_results.duration(Ind_events_OK_iter)),std(event_results.duration(Ind_events_OK_iter)))
	fprintf('  - maximum amplitude%s: mean = %0.1f µV, std = %0.1g \n',sleepstage_choice,mean(event_results.maximum_amplitude(Ind_events_OK_iter)),std(event_results.maximum_amplitude(Ind_events_OK_iter)))
	fprintf('  - power%s: mean = %0.1f µV^2, std = %0.2g \n',sleepstage_choice,mean(event_results.power(Ind_events_OK_iter)),std(event_results.power(Ind_events_OK_iter)))
	fprintf('  - frequency%s: mean = %0.1f Hz, std = %0.1g \n',sleepstage_choice,mean(event_results.frequency(Ind_events_OK_iter)),std(event_results.frequency(Ind_events_OK_iter)))
	fprintf('  - frequency instability%s: mean = %0.2f Hz^2, std = %0.2g \n',sleepstage_choice,mean(event_results.frequential_instability(Ind_events_OK_iter)),std(event_results.frequential_instability(Ind_events_OK_iter)))
	fprintf('  - frequential purity%s: mean = %0.0f%%, std = %0.0f \n',sleepstage_choice,mean(event_results.frequential_purity(Ind_events_OK_iter)),std(event_results.frequential_purity(Ind_events_OK_iter)))
	fprintf('  - temporal instability%s: mean = %0.1f µV/s, std = %0.2g \n',sleepstage_choice,mean(event_results.temporal_instability(Ind_events_OK_iter)),std(event_results.temporal_instability(Ind_events_OK_iter)))
	disp(' ')
end


disp(' ')
disp('--------------------------')
disp('* Analysis by sleep cycles')
disp('--------------------------')
for c = 1 : Nb_sleep_cycles
	
	Ind_events_OK_cycle = logical(mat_Ind_events_cycle_SOI(c,:));
	
	disp(['Sleep cycle #',num2str(c)])
	fprintf('  %0.5g %s detected, including %0.5g in stage %s \n',Nb_event_per_cycle(c),event_name,Nb_event_per_cycle_of_SOI(c),legend_SOI)
	fprintf('  Considering only the %0.4g epochs labelled%s in this sleep cycle: \n',Nb_epochs_SOI_per_cycle(c),legend_SOI)
	fprintf('  - number/epoch%s: mean = %0.2f/epoch, std = %0.2g \n',sleepstage_choice,Nb_events_cycle_per_epok_of_SOI(c,1),Nb_events_cycle_per_epok_of_SOI(c,2))
	fprintf('  - duration%s: mean = %0.2f s, std = %0.2g \n',sleepstage_choice,mean(event_results.duration(Ind_events_OK_cycle)),std(event_results.duration(Ind_events_OK_cycle)))
	fprintf('  - maximum amplitude%s: mean = %0.1f µV, std = %0.1g \n',sleepstage_choice,mean(event_results.maximum_amplitude(Ind_events_OK_cycle)),std(event_results.maximum_amplitude(Ind_events_OK_cycle)))
	fprintf('  - power%s: mean = %0.1f µV^2, std = %0.2g \n',sleepstage_choice,mean(event_results.power(Ind_events_OK_cycle)),std(event_results.power(Ind_events_OK_cycle)))
	fprintf('  - frequency%s: mean = %0.1f Hz, std = %0.1g \n',sleepstage_choice,mean(event_results.frequency(Ind_events_OK_cycle)),std(event_results.frequency(Ind_events_OK_cycle)))
	fprintf('  - frequency instability%s: mean = %0.2f Hz^2, std = %0.2g \n',sleepstage_choice,mean(event_results.frequential_instability(Ind_events_OK_cycle)),std(event_results.frequential_instability(Ind_events_OK_cycle)))
	fprintf('  - frequential purity%s: mean = %0.0f%%, std = %0.0f \n',sleepstage_choice,mean(event_results.frequential_purity(Ind_events_OK_cycle)),std(event_results.frequential_purity(Ind_events_OK_cycle)))
	fprintf('  - temporal instability%s: mean = %0.1f µV/s, std = %0.2g \n\n',sleepstage_choice,mean(event_results.temporal_instability(Ind_events_OK_cycle)),std(event_results.temporal_instability(Ind_events_OK_cycle)))
end
disp(' ')



%              ____             __            __            __  
%             / __/_ _________ / / ___  __ __/ /____  __ __/ /_ 
%            / _/ \ \ / __/ -_) / / _ \/ // / __/ _ \/ // / __/ 
%           /___//_\_\\__/\__/_/  \___/\_,_/\__/ .__/\_,_/\__/  
%                                             /_/               
if flag_Excel_ouput
%keyboard
	% Output file name
	xls_outpufile = sprintf('%s_Aseega_Microstruct_%s_Band_%0.5g-%0.5gHz_Stage%s.xls',name_subject,type_microstructure,FOI(1),FOI(2),legend_SOI(2:length(legend_SOI)));

	% Script version
	xlswrite(xls_outpufile, {sprintf('%s: version %s (%s)',fname,me.version,me.date)}, 1, 'B1')

	% Sheet title
	if flag_visual_scoring
		sheet_title = {sprintf('%s -  Sleep microstructure  - %s (%0.5g-%0.5gHz) detected in stages%s (visual scoring)',name_subject,type_microstructure,FOI(1),FOI(2),legend_SOI)};
	else
		sheet_title = {sprintf('%s -  Sleep microstructure  - %s (%0.5g-%0.5gHz) detected in stages%s (Aseega scoring)',name_subject,type_microstructure,FOI(1),FOI(2),legend_SOI)};
	end
	xlswrite(xls_outpufile, sheet_title, 1, 'B4')

	cur_line = 6;
	
	% Analysis for the whole file
	% ---------------------------
	xlswrite(xls_outpufile,{'**********************************'}, 1, sprintf('B%0.5g',cur_line));cur_line = cur_line + 1;
	xlswrite(xls_outpufile,{sprintf('All-night analysis')}, 1, sprintf('B%0.5g',cur_line));cur_line = cur_line + 1;
	
	for iter = 1 : 2
		switch iter
			% For all the epochs
			case 1;	Nb_events_per_epok_iter = Nb_events_per_epok;
				Nb_events_iter = Nb_events;
				Ind_events_OK_iter = logical(ones(Nb_events,1));
				sleepstage_choice = '';
				xlswrite(xls_outpufile,{sprintf('Considering all the %0.4g epochs of the recording, whatever the sleep stage:',length(hypno_5stages))}, 1, sprintf('C%0.5g',cur_line))
			% Only for the epochs of SOI
			case 2; Nb_events_per_epok_iter = Nb_events_per_epok_of_SOI;
				Nb_events_iter = Nb_events_SOI;
				Ind_events_OK_iter = Ind_events_SOI;
				sleepstage_choice = sprintf(' in stage%s',legend_SOI);
				xlswrite(xls_outpufile,{sprintf('Considering only the %0.4g epochs labelled%s:',length(ptr_SOI_epk),legend_SOI)}, 1, sprintf('C%0.5g',cur_line))
		end
		
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.5g',Nb_events_iter)}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' %s detected%s',event_name,sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2f',mean(Nb_events_per_epok_iter))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = number/epoch%s: mean',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2g',std(Nb_events_per_epok_iter))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = number/epoch%s: std',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2f',mean(event_results.duration(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = duration%s (s): mean',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2g',std(event_results.duration(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = duration%s (s): std',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1f',mean(event_results.maximum_amplitude(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = maximum amplitude%s (µV): mean',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1g',std(event_results.maximum_amplitude(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = maximum amplitude%s (µV): std',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))
		
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1f',mean(event_results.power(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = power%s (µV^2): mean',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2g',std(event_results.power(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = power%s (µV^2): std',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1f',mean(event_results.frequency(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency%s (Hz): mean',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1g',std(event_results.frequency(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency%s (Hz): std',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2f',mean(event_results.frequential_instability(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency instability%s (Hz^2): mean',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2g',std(event_results.frequential_instability(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency instability%s (Hz^2): std',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.0f',mean(event_results.frequential_purity(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequential purity%s (%%): mean',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.0f',std(event_results.frequential_purity(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequential purity%s (%%): std',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1f',mean(event_results.temporal_instability(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency instability%s (µV/s): mean',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2g',std(event_results.temporal_instability(Ind_events_OK_iter)))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency instability%s (µV/s): std',sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))

		cur_line = cur_line + 2;
	end
	cur_line = cur_line + 1;

	
	% Analysis by sleep cycles
	% ---------------------------
	xlswrite(xls_outpufile,{'**********************************'}, 1, sprintf('B%0.5g',cur_line));	cur_line = cur_line + 1;
	xlswrite(xls_outpufile,{sprintf('Analysis by sleep cycles')}, 1, sprintf('B%0.5g',cur_line));	cur_line = cur_line + 1;
	%sleepstage_choice = sprintf(' in stage%s',legend_SOI);
	
	for c = 1 : Nb_sleep_cycles

		Ind_events_OK_cycle = logical(mat_Ind_events_cycle_SOI(c,:));

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('Sleep cycle # %0.3g:',c)}, 1, sprintf('C%0.5g',cur_line));
		
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.5g',Nb_event_per_cycle(c))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' %s detected in the whole cycle',event_name)}, 1, sprintf('E%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.5g',Nb_event_per_cycle_of_SOI(c))}, 1, sprintf('D%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' %s detected%s',event_name,sleepstage_choice)}, 1, sprintf('E%0.5g',cur_line))

 		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('Considering only the %0.4g epochs labelled%s in this sleep cycle:',Nb_epochs_SOI_per_cycle(c),legend_SOI)}, 1, sprintf('E%0.5g',cur_line))

		cur_line = cur_line + 1;
 		xlswrite(xls_outpufile,{sprintf('%0.2f',Nb_events_cycle_per_epok_of_SOI(c,1))}, 1, sprintf('E%0.5g',cur_line))
 		xlswrite(xls_outpufile,{sprintf(' = number/epoch%s: mean',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))
 		cur_line = cur_line + 1;
 		xlswrite(xls_outpufile,{sprintf('%0.2g',Nb_events_cycle_per_epok_of_SOI(c,2))}, 1, sprintf('E%0.5g',cur_line))
 		xlswrite(xls_outpufile,{sprintf(' = number/epoch%s: std',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2f',mean(event_results.duration(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = duration%s (s): mean',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2g',std(event_results.duration(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = duration%s (s): std',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1f',mean(event_results.maximum_amplitude(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = maximum amplitude%s (µV): mean',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1g',std(event_results.maximum_amplitude(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = maximum amplitude%s (µV): std',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))
		
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1f',mean(event_results.power(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = power%s (µV^2): mean',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2g',std(event_results.power(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = power%s (µV^2): std',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1f',mean(event_results.frequency(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency%s (Hz): mean',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1g',std(event_results.frequency(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency%s (Hz): std',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2f',mean(event_results.frequential_instability(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency instability%s (Hz^2): mean',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2g',std(event_results.frequential_instability(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency instability%s (Hz^2): std',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.0f',mean(event_results.frequential_purity(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequential purity%s (%%): mean',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.0f',std(event_results.frequential_purity(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequential purity%s (%%): std',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))

		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.1f',mean(event_results.temporal_instability(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency instability%s (µV/s): mean',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))
		cur_line = cur_line + 1;
		xlswrite(xls_outpufile,{sprintf('%0.2g',std(event_results.temporal_instability(Ind_events_OK_cycle)))}, 1, sprintf('E%0.5g',cur_line))
		xlswrite(xls_outpufile,{sprintf(' = frequency instability%s (µV/s): std',sleepstage_choice)}, 1, sprintf('F%0.5g',cur_line))


% 		fprintf('  - duration%s: mean = %0.2f s, std = %0.2g \n',sleepstage_choice,mean(event_results.duration(Ind_events_OK_cycle)),std(event_results.duration(Ind_events_OK_cycle)))
% 		fprintf('  - maximum amplitude%s: mean = %0.1f µV, std = %0.1g \n',sleepstage_choice,mean(event_results.maximum_amplitude(Ind_events_OK_cycle)),std(event_results.maximum_amplitude(Ind_events_OK_cycle)))
% 		fprintf('  - power%s: mean = %0.1f µV^2, std = %0.2g \n',sleepstage_choice,mean(event_results.power(Ind_events_OK_cycle)),std(event_results.power(Ind_events_OK_cycle)))
% 		fprintf('  - frequency%s: mean = %0.1f Hz, std = %0.1g \n',sleepstage_choice,mean(event_results.frequency(Ind_events_OK_cycle)),std(event_results.frequency(Ind_events_OK_cycle)))
% 		fprintf('  - frequency instability%s: mean = %0.2f Hz^2, std = %0.2g \n',sleepstage_choice,mean(event_results.frequential_instability(Ind_events_OK_cycle)),std(event_results.frequential_instability(Ind_events_OK_cycle)))
% 		fprintf('  - frequential purity%s: mean = %0.0f%%, std = %0.0f \n',sleepstage_choice,mean(event_results.frequential_purity(Ind_events_OK_cycle)),std(event_results.frequential_purity(Ind_events_OK_cycle)))
% 		fprintf('  - temporal instability%s: mean = %0.1f µV/s, std = %0.2g \n\n',sleepstage_choice,mean(event_results.temporal_instability(Ind_events_OK_cycle)),std(event_results.temporal_instability(Ind_events_OK_cycle)))
	end

	

end




%          _____              __   _          __  ___  __     __  
%         / ___/______ ____  / /  (_)______ _/ / / _ \/ /__  / /_ 
%        / (_ / __/ _ `/ _ \/ _ \/ / __/ _ `/ / / ___/ / _ \/ __/ 
%        \___/_/  \_,_/ .__/_//_/_/\__/\_,_/_/ /_/  /_/\___/\__/  
%                    /_/                                          


if flag_plot_fig

	fig_size = [2 4 25  15];	% figure position and size
	Hfig_size = [2 4 25 20];	% figure position and size (big)
	size_police = 8;		
	Nb_epk = length(hypno_5stages);	% Number of epochs in recording
	
		  

	% Events per epochs
	% =================
	figure('units','centimeters','PaperUnits','centimeters','name','Event detection - Epoch scale','numbertitle','off','color','w','position',fig_size);

	subplot(211)
		plot(Nb_events_per_epok);
		set(gca, 'xlim',[1 Nb_epk]);grid; ylabel('Nb / epoch')
		if flag_visual_scoring
			title(sprintf('%s -  Detection of sleep %s (N = %0.5g) the whole recording (visual stage scoring)',name_subject_fig,event_name,Nb_events),'fontweight','bold')
		else
			title(sprintf('%s -  Detection of sleep %s (N = %0.5g) the whole recording (Aseega stage scoring)',name_subject_fig,event_name,Nb_events),'fontweight','bold')
		end
		plot_scale = axis;
	subplot(212)
		Plot_Nb_events_per_epok_SOI = zeros(size(Nb_events_per_epok));
		Plot_Nb_events_per_epok_SOI(ptr_SOI_epk) = Nb_events_per_epok(ptr_SOI_epk);
		plot(Plot_Nb_events_per_epok_SOI)
		if flag_visual_scoring
			set(gca, 'xlim',[1 Nb_epk],'ylim',[0 plot_scale(4)]);grid; ylabel(sprintf('Nb / epoch of visual stage %s',legend_SOI))
		else
			set(gca, 'xlim',[1 Nb_epk],'ylim',[0 plot_scale(4)]);grid; ylabel(sprintf('Nb / epoch of stage %s',legend_SOI))
		end
		xlabel('Time (30s-epochs)')


	% Event characteristics
	% =====================
	figure('units','centimeters','PaperUnits','centimeters','name','Event detection - Event scale','numbertitle','off','color','w','position',Hfig_size);
	Nb_subplots = 7;	% Seven characteristics for each spindle
	
	
	% General plot or restricted to the SOI
	% =====================================
	if flag_SOI_restricted
		ptr2plot = Ind_events_SOI;
		sleepstage_choice = sprintf('in stage%s',legend_SOI);
		cut_Nb_events = Nb_events_SOI;
	else
		ptr2plot = (1:length(event_results.duration));
		sleepstage_choice = 'in whole recording';
		cut_Nb_events = Nb_events;
	end
	
	
	% Plot
	% ====
	abscisse = ceil(event_results.position(ptr2plot)/epok_lgth);
	x_lim_epok = [0  abscisse(length(abscisse))];


	subplot(Nb_subplots,1,1)
		duration = event_results.duration(ptr2plot);
		plot(abscisse,duration,'.')
		title(sprintf('Aseega %s detections %s: Event analysis (N = %0.5g)  -  \\sigma =  [%0.3g - %0.3g] Hz',event_name,sleepstage_choice,cut_Nb_events,FOI(1),FOI(2)),'fontweight','bold');
		ylabel('Duration (s)');
		set(gca,'xlim',x_lim_epok,'xticklabel',[]);grid
		xlabel(['Duration:    \mu = ',num2str(0.01*round(100*mean(duration))),'      med = ',num2str(0.01*round(100*median(duration))),'      \sigma = ',num2str(0.01*round(100*std(duration)))],'fontsize',size_police);
	subplot(Nb_subplots,1,2)
		maximum_amplitude = event_results.maximum_amplitude(ptr2plot);
		plot(abscisse,maximum_amplitude,'.')	
		ylabel('Amp (µV)')
		set(gca,'xlim',x_lim_epok,'xticklabel',[]);grid
		xlabel(['Ampl. max:    \mu = ',num2str(0.01*round(100*mean(maximum_amplitude))),'      med = ',num2str(0.01*round(100*median(maximum_amplitude))),'      \sigma = ',num2str(0.01*round(100*std(maximum_amplitude)))],'fontsize',size_police);
	subplot(Nb_subplots,1,3)
		power = event_results.power(ptr2plot);
		plot(abscisse,power,'.')	
		ylabel('Pwr (µV ^2)');	set(gca,'xticklabel',[]);xlabel('')
		set(gca,'xlim',x_lim_epok,'xticklabel',[]);grid
		xlabel(['Power:    \mu = ',num2str(0.01*round(100*mean(power))),'      med = ',num2str(0.01*round(100*median(power))),'      \sigma = ',num2str(0.01*round(100*std(power)))],'fontsize',size_police);
	subplot(Nb_subplots,1,4)
		frequency = event_results.frequency(ptr2plot);
		plot(abscisse,frequency,'.')	
		ylabel('Freq. (Hz)');	
		set(gca,'xlim',x_lim_epok,'xticklabel',[]);grid
		xlabel(['Frequency:    \mu = ',num2str(0.01*round(100*mean(frequency))),'      med = ',num2str(0.01*round(100*median(frequency))),'      \sigma = ',num2str(0.01*round(100*std(frequency)))],'fontsize',size_police);
	subplot(Nb_subplots,1,5)
		frequential_instability = event_results.frequential_instability(ptr2plot);
		plot(abscisse,frequential_instability,'.')	
		ylabel('Inst freq (Hz ^2)');	
		set(gca,'xlim',x_lim_epok,'xticklabel',[]);grid
		xlabel(['Frequential instability:    \mu = ',num2str(0.01*round(100*mean(frequential_instability))),'      med = ',num2str(0.01*round(100*median(frequential_instability))),'      \sigma = ',num2str(0.01*round(100*std(frequential_instability)))],'fontsize',size_police);
	subplot(Nb_subplots,1,6)
		frequential_purity = event_results.frequential_purity(ptr2plot);
		plot(abscisse,frequential_purity,'.')	
		ylabel('Purity (%)');	
		set(gca,'xlim',x_lim_epok,'xticklabel',[]);grid
		xlabel(['Frequency purity:    \mu = ',num2str(0.01*round(100*mean(frequential_purity))),'      med = ',num2str(0.01*round(100*median(frequential_purity))),'      \sigma = ',num2str(0.01*round(100*std(frequential_purity)))],'fontsize',size_police);
	subplot(Nb_subplots,1,7)
		temporal_instability = event_results.temporal_instability(ptr2plot);
		plot(abscisse,temporal_instability,'.')	
		ylabel('Inst temp. (µ/s)');	
		set(gca,'xlim',x_lim_epok);grid
		xlabel(['Temporal instability:    \mu = ',num2str(0.01*round(100*mean(temporal_instability))),'      med = ',num2str(0.01*round(100*median(temporal_instability))),'      \sigma = ',num2str(0.01*round(100*std(temporal_instability)))],'fontsize',size_police);
	
	
	% Hypnogram
	% =========
	AseegaPlot_Hypno(aseega_analysis)
	
end	

disp('End of Aseega result extraction.')
