%    ____     __   __                          _     __           __  
%   / _(_)__ / /__/ / ___  ___  ___  _____ __ (_)__ / /____ ___  / /_ 
%  / _/ / -_) / _  / / _ \/ _ \/ _ \/ -_) \ // (_-</ __/ -_) _ \/ __/ 
% /_//_/\__/_/\_,_/ /_//_/\___/_//_/\__/_\_\/_/___/\__/\__/_//_/\__/  
%                             __       
%  ___  ____  ___ __ _  ___  / /___ __ 
% / _ \/ __/ / -_)  ' \/ _ \/ __/ // / 
% \___/_/    \__/_/_/_/ .__/\__/\_, /  
%                    /_/       /___/   
%
%
% Usage: ret = field_nonexistent_or_empty(field_string)
%
% Regarde si un certain champ d'une struct existe (récursivement) et s'il est vide.
% Passer dans 'field_string' le "chemin du champ", par exemple 'options.spectro.colormap'
%
% Retourne:
%   0 si ce champ existe et est non-vide
%   1 si ce champ existe mais est vide
%   2 si ce champ (ou un des sous-champs intermédiaires) n'existe pas dans cette variable
%   3 si la variable n'existe même pas.
%
% Exemple:   ret = field_nonexistent_or_empty('aseega_analysis.info.recording')
%
% Physip, Paris France
% 2002-2017

% 16.09.13  v1.0  création

function [ret] = field_nonexistent_or_empty(field_string)
	fname = mfilename;

	ret = 0;

	try
		value = evalin('caller', field_string);			% essaie d'accéder au champ

	catch								% err1) champ inexistant

		if (isempty(findstr(lasterr, 'Reference to non-existent field')) ...
		 || isempty(findstr(lasterr, 'Attempt to reference field of non-structure array')))
			ret = 2;

		elseif (findstr(lasterr, 'Undefined variable'))		% err2) variable inexistante
			ret = 3;
		else							% autres erreurs ...	
			error(sprintf('%s: %s', fname, lasterr));	% erreur !
		end

	end
					% si j'ai pu évaluer field_string
	if (ret == 0)			% <=> si ce champ existe :
		if (isempty(value))	% regarde maintenant s'il est vide,
			ret = 1;	% si oui retourne 1
		end
	end

	return

