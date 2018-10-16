%
%     ___                         ___  __     __     __ __                   
%    / _ | ___ ___ ___ ___ ____ _/ _ \/ /__  / /_   / // /_ _____  ___  ___  
%   / __ |(_-</ -_) -_) _ `/ _ `/ ___/ / _ \/ __/  / _  / // / _ \/ _ \/ _ \ 
%  /_/ |_/___/\__/\__/\_, /\_,_/_/  /_/\___/\__/__/_//_/\_, / .__/_//_/\___/ 
%                    /___/                    /___/    /___/_/               
%
%      Name:	AseegaPlot_Hypno.m
%   Purpose:	Plot the EEGWatch image and Aseega scoring, or visual scoring if available
%   Comment: 
%
%        In:	aseega_analysis, structure array provided by Aseega EEG Analysis (aseega_analysis)
%       Out:	
% Called by: 
%     Calls: 
%
%   Syntax:	AseegaPlot_Hypno(aseega_analysis)
%            
%   History:	v1.0   25-04-2016 Creation
%		v1.1   05-10-2016 User friendly syntax
%		v1.3   22-12-2017 enable plot visual scoring, if included in the Aseega matlab result file 
%		v1.3.1 17-06-2018 bugfix 
%
%  Physip, Paris France
%  2002-2018




function [] = AseegaPlot_Hypno(varargin)

me.version = '1.3.1';
me.date    = 'June 17th, 2017';
fname = 'AseegaPlot_Hypno';

%%
%       ___                          __             
%      / _ \___ ________ ___ _  ___ / /____ _______ 
%     / ___/ _ `/ __/ _ `/  ' \/ -_) __/ -_) __(_-< 
%    /_/   \_,_/_/  \_,_/_/_/_/\__/\__/\__/_/ /___/ 
%                                                   
% 
% Comment/uncomment lines to change the parameters
% ===================================================
    
% Choice of Sleep Onset latency definition, for the per-cycle analysis (see Aseega-Analysis-Format.pdf for more definitions)
SOL_def = 'sleep_onset_3S';			% R&K definition: the first 3 epochs of stage S1 or the 1st epoch of any other sleep 
%SOL_def = 'sleep_onset_1';			% AASM definition: first sleep epoch 
%SOL_def = 'sleep_onset_5m';  			% Carskadon & R definition: first 5 consecutive minutes of sleep 
%SOL_def = 'sleep_onset_10m';			% Carskadon & R definition: first 10 consecutive minutes of sleep 










					
%% - - - - - - - - - - - - DO NOT CROSS!  - - - - - - - - - - - - - - - - - -
%
%% - - - - - - - - FOR EXPERIMENTED MATLAB USERS ONLY - - - - - - - - - - - -

if nargin == 0
	disp(['   Error using ',fname,'.m, incorrect number of input arguments.'])
	disp('   The correct syntax is:')
	disp('   >> AseegaPlot_Hypno(aseega_analysis)')
	disp('   Thanks')
	disp(' ')
	return
else
	aseega_analysis = varargin{1};
end

% disp(' ');
% fprintf('%s: version %s (%s)',fname,me.version,me.date)
% disp(' ');disp(' ');

name_subject = aseega_analysis.info.recording.patient_id;
name_subject_fig = strrep(name_subject,'_','\_');	% remplace  underscore by '\_' for display

flag_hyp_Aseega = 0;
flag_hyp_visuel = 0;
Nb_scorers = 1;

% Aseega hypnogram, if available
if isfield(aseega_analysis.scoring,'hypno5')
	hypno_5stages = aseega_analysis.scoring.hypno5;		% Current five-state hypnogram: W,R,N1,N2,N3
	if ~isempty(hypno_5stages)
		Nb_epk = length(hypno_5stages);			% Number of epochs in recording
		flag_hyp_Aseega = 1;
	end
end

if ~flag_hyp_Aseega
	disp('Aseega hypnogram not found ...')
	hypno_5stages = NaN*ones(1,2000);
end


% Visual hypnogram, if available
if isfield(aseega_analysis.scoring,'hypno_visual_1')
	hypno_visual = aseega_analysis.scoring.hypno_visual_1;
	if ~isempty(hypno_visual)
		Nb_epk = length(hypno_visual);			% Number of epochs in recording
		Nb_scorers = Nb_scorers + 1;
		flag_hyp_visuel = 1;
		if ~flag_hyp_Aseega
			hypno_5stages = hypno_5stages(1:Nb_epk);
		end
		hypno_5stages = [hypno_5stages;hypno_visual];	% Hypno auto (eventually NaN) + hypno visual
	end
end

% Exit if no hypnogram found
if ~flag_hyp_Aseega && ~flag_hyp_visuel
	disp('No hypnogram found in this file')
	return
end


% Sleep cycles
if flag_hyp_Aseega
	cycle_analysis = aseega_analysis.sleepparams.(SOL_def).sleep_cycles.cycle;
	cycles_Nb = length(cycle_analysis);
	cycle_boundaries = zeros(cycles_Nb,2);
	for icycle = 1 : cycles_Nb
		cycle_boundaries(icycle,1) = cycle_analysis{icycle}.start_e;
		cycle_boundaries(icycle,2) = cycle_analysis{icycle}.end_e;
	end
end

% Graphics
fig_size = [2 6 25  6+6*Nb_scorers];		% figure position and size
colormap_MM =[	0.4938    0         0.6000	% static colormap
		0         0         0.6250
		0         0.2500    1.0000
		0         0.6250    1.0000
		0         0.8125    1.0000
		0         1.0000    0.8125
		0         0.9000         0
		0.6000    1.0000         0
		1.0000    1.0000         0
		1.0000    0.6250         0
		1.0000    0.3750         0
		0.9000         0         0
		0.6250         0         0 ];





figure('units','centimeters','PaperUnits','centimeters','name','Aseega scoring & EEGWatch','numbertitle','off','color','w','position',fig_size);

%%
%     ______________      __     __      __  
%    / __/ __/ ___/ | /| / /__ _/ /_____/ /  
%   / _// _// (_ /| |/ |/ / _ `/ __/ __/ _ \ 
%  /___/___/\___/ |__/|__/\_,_/\__/\__/_//_/ 
%                                            
subplot(1+Nb_scorers,1,1)

	% data check
	if (field_nonexistent_or_empty('aseega_analysis.macro.eegwatch'))	% if analysis file does not contains EEGWatch data
		not_avail_text = '[ ANALYSIS NOT AVAILABLE ]';
		size_fft_DC = 128 * 2;
	else
		not_avail_text = '';
		size_fft_DC = size(aseega_analysis.macro.eegwatch, 1) * 2;
	end

	% analysis parameters
	fech		= 100;		% static sampling frequency
	imax		= size_fft_DC / 2;
	ylims		= [0  imax];
	yticklabels	= char(num2str(round([0;10;20;30;40] * fech/100)), 'Hz');
	yticks		= imax * [0  0.2  0.4  0.6  0.8  1];
	
	% x-axis managment
	datenum_start  = datenum(aseega_analysis.info.analysis.start_time);		% start time of the 1st drawn epoch in datenum format
	nepochs = size(aseega_analysis.macro.eegwatch, 2);
	hms_1epoch     = sprintf('00:00:%02d', aseega_analysis.info.analysis.nsecs_1epoch);% = '00:00:30'
	datenum_1epoch = datenum(hms_1epoch) - datenum('00:00:00');	% duration of one epoch in datenum format
	abscisses = datenum_start + (0 : nepochs - 1) * datenum_1epoch;
	ordonnees = (1 : imax);
	datenum_end = datenum_start + nepochs * datenum_1epoch;	% time of the right end of the graph in datenum format

	% plot EEGWatch
	if (isempty(not_avail_text))						% 1) Plot EEGWatch ...
		matRGB = egw2matRGB(aseega_analysis.macro.eegwatch, colormap_MM);
		cla
		image(abscisses, ordonnees, matRGB);		% | IMAGE(matRGB) |
		set(gca, 'YDir','normal','Layer', 'top')
		grid('on');
	else									% 2) ... or text "NOT AVAILABLE"
		xt = (abscisses(1) + abscisses(nepochs)) / 2;
		yt = (ordonnees(1) + ordonnees(imax))    / 2;
		ht = text(xt, yt, not_avail_text, 'Parent', gca);		% | XXX_NOT_AVAILABLE |
		set(ht,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
	end
	
	% x-axis
	set(gca,'FontSize',10);
	xlim = [datenum_start,datenum_end];	% graph x-axis must spread from epoch 1 to epoch N
	set(gca, 'XLim', xlim);
	datetick(gca, 'x', 15, 'keeplimits');	% automatically create 'HH:MM' labels on the X axis
	xlabel('Time (wall clock)')

	% y-axis
	set(gca, 'Ylim',	ylims,'YTick',yticks,'YTickLabel',yticklabels,'Ydir','normal','Unit','Normalized');
	ylabel('Frequency')
	
	% title
	title([name_subject_fig,' -  EEGWatch'],'fontweight','bold','FontName','Calibri','FontSize',12)
	
	
%%	
%    ___                            __    __         _               __  _                     _           
%   / _ | ___ ___ ___ ___ ____ _  _/_/ __/ /_  _  __(_)__ __ _____ _/ / | |  ___ _______  ____(_)__  ___ _ 
%  / __ |(_-</ -_) -_) _ `/ _ `/ / /  /_  __/ | |/ / (_-</ // / _ `/ /  / / (_-</ __/ _ \/ __/ / _ \/ _ `/ 
% /_/ |_/___/\__/\__/\_, /\_,_/ / /    /_/    |___/_/___/\_,_/\_,_/_/ _/_/ /___/\__/\___/_/ /_/_//_/\_, /  
%                   /___/       |_|                                  /_/                           /___/   


for i = 1 : Nb_scorers
	
	switch i
		case 1; title_part = 'Aseega scoring';
		case 2; title_part = 'Visual scoring';
	end
	cur_hypno = hypno_5stages(i,:);
	
	subplot(1+Nb_scorers,1,1+i)

		% Sleep stage indexes
		hypno_5stages_without_Art = cur_hypno;
		hypno_5stages_without_Art(cur_hypno == 5) = 2;
		hypno_5stages_without_Art(cur_hypno == 7) = 3;
		hypno_5stages_without_Art(cur_hypno == 8) = 4;
		hypno_5stages_without_Art(cur_hypno == 10)= 5;
		hypno_5stages_without_Art(cur_hypno == 11) = NaN;	

		% Plot
		plot(hypno_5stages_without_Art);hold on
		ptr_R = find(cur_hypno == 7);
		plot(ptr_R,3*ones(size(ptr_R)),'.r')
		ptr_Art =  find(cur_hypno == 11);
		plot(ptr_Art,6*ones(size(ptr_Art)),'*','color',0.7*[1 1 1])
		ptr_W = find(cur_hypno == 10);
		plot(ptr_W,5*ones(size(ptr_W)),'.','linewidth',2)

		title([name_subject_fig,' -  ',title_part],'fontweight','bold','FontName','Calibri','FontSize',12)
		ylabel ('Sleep stages')
		xlabel('Time (30s-epoch)')
		ymin = 0.5;
		ymax = 6.5;
		set(gca, 'xlim',[1 Nb_epk],'ytick',(1:6),'yticklabel',[' N3 ';' N2 ';'  R ';' N1 ';'  W ';'Art '],'ylim',[ymin ymax])

		% Mark the Aseega sleep cycles for visual check
		if i == 1 
			if flag_hyp_Aseega
				for icycle = 1 : cycles_Nb
					plot([cycle_boundaries(icycle,1) cycle_boundaries(icycle,1)],[ymin ymax],'r--','linewidth',2);
					plot([cycle_boundaries(icycle,2) cycle_boundaries(icycle,2)],[ymin ymax ],'r--','linewidth',2);
					text(mean([cycle_boundaries(icycle,1) cycle_boundaries(icycle,2)]), 5.5,sprintf('cycle #%0.5g',icycle),'color','r','horizontalalignment','center')
				end
			else
				text(Nb_epk/2,3.5,'Aseega hypnogram not found','fontname','Calibri','fontsize',12,'horizontalalignment','center','color',.75*[1 1 1])
			end
		end
end


