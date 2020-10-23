% fooof_check_settings() - Check a struct of settings for the FOOOF model.
%
% Usage:
%  >> settings = fooof_check_settings(settings)
%
% Inputs:
%   settings        = struct, can optionally include:
%       settings.peak_width_limts
%       settings.max_n_peaks
%       settings.min_peak_amplitude
%       settings.peak_threshold
%       settings.background_mode
%       settings.verbose
%
% Outputs:
%   settings        = struct, with all settings defined:
%       settings.peak_width_limts
%       settings.max_n_peaks
%       settings.min_peak_amplitude
%       settings.peak_threshold
%       settings.background_mode
%       settings.verbose
%
% Notes:
%   This is a helper function, probably not called directly by the user.
%   Any settings not specified are set to default values

function settings = fooof_check_settings(settings)

    % Set defaults for all settings
    defaults = struct(...
        'peak_width_limits', [0.5, 12], ...
        'max_n_peaks', Inf, ...
        'min_peak_amplitude', 0.0, ...
        'peak_threshold', 2.0, ...
        'background_mode', 'fixed', ...
        'verbose', true);

    % Overwrite any non-existent or nan settings with defaults
    for field = fieldnames(defaults)'
        if ~isfield(settings, field) || all(isnan(settings.(field{1})))
            settings.(field{1}) = defaults.(field{1});
        end
    end

end