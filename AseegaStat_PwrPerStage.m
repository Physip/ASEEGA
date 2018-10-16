%    ___                         ______       __     ___             ___          ______                
%   / _ | ___ ___ ___ ___ ____ _/ __/ /____ _/ /_   / _ \_    ______/ _ \___ ____/ __/ /____ ____ ____  
%  / __ |(_-</ -_) -_) _ `/ _ `/\ \/ __/ _ `/ __/  / ___/ |/|/ / __/ ___/ -_) __/\ \/ __/ _ `/ _ `/ -_) 
% /_/ |_/___/\__/\__/\_, /\_,_/___/\__/\_,_/\__/__/_/   |__,__/_/ /_/   \__/_/ /___/\__/\_,_/\_, /\__/  
%                   /___/                     /___/                                         /___/       
% 
%
%      Name:	AseegaStat_PwrPerStage.m
%   Purpose:	Compute spectral power in a given frequency band of interest
%   Comment:	This script will only run for the five-state hypnogram: W,R,N1,N2,N3, since
%		sleep stage labeling in not dynamic. For more information on result format, please
%		check the Aseega-Analysis-Format pdf file. 
%           
%        In:	aseega_analysis, structure array provided by Aseega EEG Analysis
%		flag_cycleanalysis,  1 = per-cycle analysis, 0 = all-night analysis 
%		flag_plot_fig,  1 = plot figures, 0 = no figure
%		flag_Excel_ouput, 1 = results in Excel file,  = no output file
%		flag_visual_scoring, 1 = use visual scoring, 0 = use Aseega scoring
%
%       Out:	Pwr_FOI_SOI, Spectral power in frequency band of interest (FOI) for sleep stages of interest (SOI)
% Called by: 
%     Calls:	AseegaPlot_Hypno.m
%
%   History:	v1.0   24-02-2014 Creation
%		v2.0   25-04-2016 Add Excel ouput file, per-cycle analysis and EEGWatch plot
%		v2.1   05-10-2016 User friendly syntax
%		v2.2   01-03-2017 Excel Layout
%		v2.3   22-12-2017 can use visual scoring, if provided
%
%   Example:	load('GUEA-C4O2.AseegaAnalysis.A4R.mat')	% loads the Aseega structure result of the recording GUEA-C4O2.edf
%		AseegaStat_PwrPerStage(aseega_analysis,1,1,1,0);% grab the results of interest according to the Parameters section,...
%								% ... perform per-cycle analysis, plot the results and write an Excel ouput file
%            
%
%  Physip, Paris France
%  2002-2018




function [Pwr_FOI] = AseegaStat_PwrPerStage(varargin)

me.version = '2.3';
me.date    = 'Physip, Dec 22nd, 2017';
fname = 'AseegaStat_PwrPerStage';





%       ___                          __             
%      / _ \___ ________ ___ _  ___ / /____ _______ 
%     / ___/ _ `/ __/ _ `/  ' \/ -_) __/ -_) __(_-< 
%    /_/   \_,_/_/  \_,_/_/_/_/\__/\__/\__/_/ /___/ 
%                                                   
% 
% Comment/uncomment lines to change the parameters
% ===================================================

%% Choice of frequency bands:
type_frequency_bands = 'usual';	% usual bands: [0.1  4  8  12   16   50] Hz
%type_frequency_bands = 'adapted';	% auto-adjust: [0.7  4  x  x+4  x+8  50] Hz ((cf. Aseega-Analysis-Format.pdf)

%% Choice of EEG data:
%type_EEG_data = 'raw';			% results on raw EEG data	
type_EEG_data = 'clean';		% results on auto-cleaned EEG data	

%% Choice of units:
%type_unit = 'physical';			% µV^2	
type_unit = 'normalized';		% normalized: for a given epoch, sum(power in each band) = 1	

%% Choice of SOI
%sleep_stages = {'N2';'N3'};		% Sleep stages of interest (SOI)
sleep_stages = {'N2'};			% Sleep stages of interest (SOI)
%sleep_stages = {'N3'};		% Sleep stages of interest (SOI)

%% Choice of FOI
%freq_band = 'delta';			% Frequency band of interest (FOI) among {'delta','theta','alpha','sigma','beta'}
freq_band = 'alpha';			% See above, type_frequency_bands, for frequency band boudaries.
    
%% Choice of Sleep Onset latency definition, for the per-cycle analysis (see Aseega-Analysis-Format.pdf for more definitions)
SOL_def = 'sleep_onset_3S';			% R&K definition: the first 3 epochs of stage S1 or the 1st epoch of any other sleep 
%SOL_def = 'sleep_onset_1';			% AASM definition: first sleep epoch 
%SOL_def = 'sleep_onset_5m';  			% Carskadon & R definition: first 5 consecutive minutes of sleep 
%SOL_def = 'sleep_onset_10m';			% Carskadon & R definition: first 10 consecutive minutes of sleep 

					





%% - - - - - - - - - - - - DO NOT CROSS!  - - - - - - - - - - - - - - - - - -
%
% - - - - - - - - FOR EXPERIMENTED MATLAB USERS ONLY - - - - - - - - - - - -



if nargin < 5 || nargin > 5
	disp(['   Error using ',fname,'.m, incorrect number of input arguments.'])
	disp('   The correct syntax is:')
	disp('   >> AseegaStat_PwrPerStage(aseega_analysis,flag_cycleanalysis,flag_plot_fig,flag_Excel_ouput,flag_visual_scoring)')
	disp('   Thanks')
	disp(' ')
	return
else
	aseega_analysis		= varargin{1};
	flag_cycleanalysis	= varargin{2};
	flag_plot_fig		= varargin{3};
	flag_Excel_ouput	= varargin{4};
	flag_visual_scoring	= varargin{5};
end





%       ___             __         _    
%      / _ | ___  ___ _/ /_ _____ (_)__ 
%     / __ |/ _ \/ _ `/ / // (_-</ (_-< 
%    /_/ |_/_//_/\_,_/_/\_, /___/_/___/ 
%                      /___/            
%
% Will only work for AASL 5-state hypnogramm.
% To be customized for other number of sleep states
% ==================================================
[~,name_subject,~] = fileparts(aseega_analysis.info.recording.filename);
[Nb_SOI, ~] = size(sleep_stages);			% Number of stages of interest

if strcmp(type_frequency_bands,'usual')
    bd_Hz = [0.1  4  8  12  16  50];
else
    bd_Hz = aseega_analysis.macro.adapted_band_freqs;
end

legend_SOI = '';filename_SOI = '';
for i = 1 : Nb_SOI;
    legend_SOI = sprintf('%s %s',legend_SOI,sleep_stages{i});
    filename_SOI = sprintf('%s%s',filename_SOI,sleep_stages{i});
end



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
ptr_SOI = cat(2,ptr_W,ptr_R,ptr_N1,ptr_N2,ptr_N3);		% SleepStage of interest epochs pointer



% Index of EEG rythm of interest
% ==============================
if strcmp(freq_band,'delta')
	band = 1;
elseif strcmp(freq_band,'theta')
	band = 2;	
elseif strcmp(freq_band,'alpha')
	band = 3;
elseif strcmp(freq_band,'sigma')
	band = 4;
elseif strcmp(freq_band,'beta')
	band = 5;
end




%              ___                               __                     __         _    
%             / _ \___ ____  ____  ______ ______/ /__   ___ ____  ___ _/ /_ _____ (_)__ 
%            / ___/ -_) __/ /___/ / __/ // / __/ / -_) / _ `/ _ \/ _ `/ / // (_-</ (_-< 
%           /_/   \__/_/          \__/\_, /\__/_/\__/  \_,_/_//_/\_,_/_/\_, /___/_/___/ 
%                                    /___/                             /___/            


% Power in frequency bands of interest for these SOI, according to the chosen parameters
% ======================================================================================
Pwr_FOI = aseega_analysis.macro.(sprintf('power_in_%s_bands',type_frequency_bands)).(sprintf('%s_eeg',type_EEG_data)).(type_unit)(:,band);	% Pwr for all epochs
allnight_Pwr_FOI_SOI = Pwr_FOI(ptr_SOI);

if ~flag_cycleanalysis			% if no cycle analysis requested, the whole recording is considered as one unique cycle
	cycles_Nb = 1;
	cycle_boundaries = [1 length(hypno_5stages)];
else
	cycle_analysis = aseega_analysis.sleepparams.(SOL_def).sleep_cycles.cycle;
	cycles_Nb = length(cycle_analysis);
	
	cycle_boundaries = zeros(cycles_Nb,2);
	for icycle = 1 : cycles_Nb
		cycle_boundaries(icycle,1) = cycle_analysis{icycle}.start_e;
		cycle_boundaries(icycle,2) = cycle_analysis{icycle}.end_e;
	end
end

disp(' ');
fprintf('%s: version %s (%s)',fname,me.version,me.date)
disp(' ');
disp('______________________________________________________________________')

for icycle = 1 : cycles_Nb
	
	cur_ptr_SOI		= ptr_SOI(ptr_SOI >= cycle_boundaries(icycle,1) & ptr_SOI <= cycle_boundaries(icycle,2));	% Stages of interest during the ith cycle
	cur_Pwr_FOI_SOI		= Pwr_FOI(cur_ptr_SOI);									% Pwer in band of interest during stages of interest during the ith cycle
	cur_lgth_Pwr_FOI_SOI	= length(cur_Pwr_FOI_SOI);
	
	cur_mean_Pwr_FOI_SOI	= mean(cur_Pwr_FOI_SOI);
	cur_std_Pwr_FOI_SOI	= std(cur_Pwr_FOI_SOI);
	cur_median_Pwr_FOI_SOI	= median(cur_Pwr_FOI_SOI);
	

	%              ___      _      __              __            __  
	%             / _ \____(_)__  / /_  ___  __ __/ /____  __ __/ /_ 
	%            / ___/ __/ / _ \/ __/ / _ \/ // / __/ _ \/ // / __/ 
	%           /_/  /_/ /_/_//_/\__/  \___/\_,_/\__/ .__/\_,_/\__/  
	%                                              /_/               

	disp(' ');
	if ~flag_cycleanalysis
		disp(['Subject ',name_subject, ' power'])
	else
		disp(['Subject ',name_subject, ' power - cycle #',num2str(icycle),'/',num2str(cycles_Nb)])
	end
	disp(['   . mean   =  ',num2str(cur_mean_Pwr_FOI_SOI)])
	disp(['   . std    =  ',num2str(cur_std_Pwr_FOI_SOI)])
	disp(['   . median =  ',num2str(cur_median_Pwr_FOI_SOI)])



	%              ____             __            __            __  
	%             / __/_ _________ / / ___  __ __/ /____  __ __/ /_ 
	%            / _/ \ \ / __/ -_) / / _ \/ // / __/ _ \/ // / __/ 
	%           /___//_\_\\__/\__/_/  \___/\_,_/\__/ .__/\_,_/\__/  
	%                                             /_/               
	if flag_Excel_ouput

		% Output file name
		xls_outpufile = sprintf('%s_%sPwr_%sstage_%s_data_%s_units.xls',aseega_analysis.info.recording.filename,freq_band,filename_SOI,type_EEG_data,type_unit);

		% Script version
		xlswrite(xls_outpufile, {sprintf('%s: version %s (%s)',fname,me.version,me.date)}, 1, 'B1')

		% Sheet title
		if flag_visual_scoring
			sheet_title = {sprintf('%s -  Spectral Power in %s band [%0.5g - %0.5g] Hz for epochs in stages%s (visual scoring, %s data, %s units)',name_subject,freq_band,bd_Hz(band),bd_Hz(band+1),legend_SOI,type_EEG_data,type_unit)};
		else
			sheet_title = {sprintf('%s -  Spectral Power in %s band [%0.5g - %0.5g] Hz for epochs in stages%s (Aseega scoring, %s data, %s units)',name_subject,freq_band,bd_Hz(band),bd_Hz(band+1),legend_SOI,type_EEG_data,type_unit)};
		end
		xlswrite(xls_outpufile, sheet_title, 1, 'B3')

		% Cycle result column
		cur_column1 = char(3*icycle + 'A'-1);
		cur_column2 = char(3*icycle + 'A');
		
		if flag_cycleanalysis
			xlswrite(xls_outpufile,{sprintf('Cycle #%0.5g',icycle)}, 1, sprintf('%s5',cur_column1))
			xlswrite(xls_outpufile,{sprintf('%0.5g epochs',cur_lgth_Pwr_FOI_SOI)}, 1, sprintf('%s5',cur_column2))
		else
			xlswrite(xls_outpufile,{legend_SOI}, 1, sprintf('%s5',cur_column1))
			xlswrite(xls_outpufile,{sprintf('%0.5g epochs',cur_lgth_Pwr_FOI_SOI)}, 1, sprintf('%s5',cur_column2))
		end
			
		% Statistics writing
		xlswrite(xls_outpufile, {'Mean'   cur_mean_Pwr_FOI_SOI }, 1, sprintf('%s7',cur_column1))
		xlswrite(xls_outpufile, {'Std'    cur_std_Pwr_FOI_SOI },  1, sprintf('%s8',cur_column1))
		xlswrite(xls_outpufile, {'Median' cur_median_Pwr_FOI_SOI},1, sprintf('%s9',cur_column1))

		% Data writing
		xlswrite(xls_outpufile,{sprintf('Power')}, 1, sprintf('%s11',cur_column1))
		if ~isempty(cur_Pwr_FOI_SOI)
			xlswrite(xls_outpufile,cur_Pwr_FOI_SOI, 1, sprintf('%s11',cur_column2))
		else
			xlswrite(xls_outpufile,{sprintf('Empty')}, 1, sprintf('%s11',cur_column2))
		end

	end
end



%       ___  __     __
%      / _ \/ /__  / /_ 
%     / ___/ / _ \/ __/ 
%    /_/  /_/\___/\__/  
%                       


if flag_plot_fig

	fig_size = [2 3 25  15];				% figure position and size
	Nb_epk = length(hypno_5stages);				% Number of epochs in recording
	name_subject_fig = strrep(name_subject,'_','\_');	% remplace  underscore by '\_' for display


	% Spectral power in FOI
	% =====================
	figure('units','centimeters','PaperUnits','centimeters','name','Spectral power','numbertitle','off','color','w','position',fig_size);

	subplot(411)
		plot(aseega_analysis.macro.power_in_usual_bands.raw_eeg.physical(:,band));hold on
		plot(ptr_SOI,aseega_analysis.macro.power_in_usual_bands.raw_eeg.physical(ptr_SOI,band),'r.')
		set(gca, 'xlim',[1 Nb_epk]);grid; ylabel('Raw EEG (µV ^2)')
		title(sprintf('%s -  Spectral Power in %s frequency band [%0.5g - %0.5g] Hz (in red: epochs of interest)',name_subject_fig,freq_band,bd_Hz(band),bd_Hz(band+1)),'fontweight','bold')
	subplot(412)
		plot(aseega_analysis.macro.power_in_usual_bands.clean_eeg.physical(:,band));hold on
		plot(ptr_SOI,aseega_analysis.macro.power_in_usual_bands.clean_eeg.physical(ptr_SOI,band),'r.')
		set(gca, 'xlim',[1 Nb_epk]);grid; ylabel('Cleaned EEG (µV ^2)')
	subplot(413)
		plot(aseega_analysis.macro.power_in_usual_bands.raw_eeg.normalized(:,band));hold on
		plot(ptr_SOI,aseega_analysis.macro.power_in_usual_bands.raw_eeg.normalized(ptr_SOI,band),'r.')
		set(gca, 'xlim',[1 Nb_epk]);grid; ylabel('Raw EEG (%)')
	subplot(414)
		plot(aseega_analysis.macro.power_in_usual_bands.clean_eeg.normalized(:,band));hold on
		plot(ptr_SOI,aseega_analysis.macro.power_in_usual_bands.clean_eeg.normalized(ptr_SOI,band),'r.')
		set(gca, 'xlim',[1 Nb_epk]);grid; ylabel('Cleaned EEG (%)')
	xlabel(' Time in 30s-epochs (all-night recording)')

	% Spectral power values for the SOI
	% =================================
	figure('units','centimeters','PaperUnits','centimeters','name','Spectral power for specific sleep stages','numbertitle','off','color','w','position',fig_size);
	
	subplot(211)
		hist(allnight_Pwr_FOI_SOI,100)
		if flag_visual_scoring
			title(sprintf('%s -  Spectral Power in %s band [%0.5g - %0.5g] Hz for epochs in stages%s (visual scoring, %s data, %s units)',name_subject_fig,freq_band,bd_Hz(band),bd_Hz(band+1),legend_SOI,type_EEG_data,type_unit),'fontweight','bold')
		else
			title(sprintf('%s -  Spectral Power in %s band [%0.5g - %0.5g] Hz for epochs in stages%s (Aseega scoring, %s data, %s units)',name_subject_fig,freq_band,bd_Hz(band),bd_Hz(band+1),legend_SOI,type_EEG_data,type_unit),'fontweight','bold')
		end
		ylabel('Pwr histogram')
		xlabel(sprintf('Mean power = %0.5g  - Std = %0.5g  -  Median = %0.5g',mean(allnight_Pwr_FOI_SOI),std(allnight_Pwr_FOI_SOI),median(allnight_Pwr_FOI_SOI)),'color','r')
	subplot(212)
		plot(allnight_Pwr_FOI_SOI);grid
		set(gca,'xlim',[1 length(allnight_Pwr_FOI_SOI)])
        if strcmp(type_unit,'physical')
    		ylabel('Pwr values (µV ^2)')
        elseif strcmp(type_unit,'normalized')
    		ylabel('Pwr values (%)')
        end
	xlabel(sprintf('Epochs of %s (all-night recording)',legend_SOI))
	
	
	% Hypnogram
	% =========
	AseegaPlot_Hypno(aseega_analysis)
			
end	


disp('End of Aseega result extraction.')
