%                   ___             __  ___  ________  
%  ___ ___ __    __|_  |__ _  ___ _/ /_/ _ \/ ___/ _ ) 
% / -_) _ `/ |/|/ / __//  ' \/ _ `/ __/ , _/ (_ / _  | 
% \__/\_, /|__,__/____/_/_/_/\_,_/\__/_/|_|\___/____/  
%    /___/                                             
%
% egw --> matRGB : transforme en couleurs les indices contenus dans egw
%
%
% Physip, Paris France
% 2002-2017


function [matRGB] = egw2matRGB(egw, colormap_DC)

  if (isempty(egw)), matRGB = ''; return; end

  Nfft2 = size(egw, 1);		% = Nfft_DC / 2
  nepq  = size(egw, 2);

  matRGB = zeros(Nfft2, nepq, 3);

  for epq = (1 : nepq)
    for ifreq = (1 : Nfft2)
      matRGB(ifreq,               epq, :) = colormap_DC(egw(ifreq, epq), :);
    end
  end

  return

