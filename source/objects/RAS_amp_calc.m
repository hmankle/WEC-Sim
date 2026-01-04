function amplitudes = RAS_amp_array(Hs, spectrum, dOmega, debug)
% RAS_amp_array - Generate an array of random amplitudes (RAM)
%                 with truncated Rayleigh distribution for each frequency bin
%                 and plot all PDFs and sampled amplitudes if debug=true
%
% Inputs:
%   Hs       - Significant wave height [m]
%   spectrum - Spectral density array S(ω)
%   dOmega   - Frequency spacing [rad/s]
%   debug    - optional true/false; if true, plots PDFs and amplitudes
%
% Output:
%   amplitudes - array of random amplitudes, same length as spectrum

    if nargin < 5
        debug = false;
    end

    % Ensure column vectors
    spectrum = spectrum(:);
    
    N = length(spectrum);
    amplitudes = zeros(N,1);

    % Maximum amplitude (truncation)
    Amax = 2.5*(Hs/2);

    % Precompute PDFs for plotting
    if debug
        x_plot = linspace(0, Amax, 500);
        pdf_matrix = zeros(N, length(x_plot));
    end

    % Loop over frequency bins
    for i = 1:N
        var_i = 2 * spectrum(i) * dOmega(i);
        if var_i > 0 & isfinite(var_i)
            B = sqrt(var_i / 2);
            
            % Truncated Rayleigh via inverse CDF
            u = rand();
            amplitudes(i) = B .* sqrt(-2 .* log(1 - u .* (1 - exp(-Amax.^2 ./ (2 .* B.^2)))));
            
            % Store PDF for debug plotting
            if debug
                pdf_i = (x_plot ./ B.^2) .* exp(-x_plot.^2 ./ (2 .* B.^2));
                pdf_i = pdf_i ./ trapz(x_plot, pdf_i); % normalize for truncated range
                pdf_matrix(i,:) = pdf_i;
            end
        else
            amplitudes(i) = 0;
            if debug
                pdf_matrix(i,:) = zeros(size(x_plot));
            end
        end
    end

    % Debug plot
    if debug
        figure;
        hold on;
        % Plot all individual PDFs in light gray
        for i = 1:N
            plot(x_plot, pdf_matrix(i,:), 'Color', [0.7 0.7 0.7]);
        end
        % Plot summed PDF in black
        sum_pdf = sum(pdf_matrix,1);
        sum_pdf = sum_pdf ./ trapz(x_plot, sum_pdf); % normalize
        plot(x_plot, sum_pdf, 'k', 'LineWidth', 2);

        % Overlay sampled amplitudes
        stem(amplitudes, interp1(x_plot, sum_pdf, amplitudes,'linear','extrap'), ...
             'r', 'LineWidth', 1.5);

        % Vertical line showing truncation limit
        xline(Amax, '--b', '2.5*Hs', 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center');

        xlabel('Amplitude [m]');
        ylabel('PDF');
        title('Truncated Rayleigh PDFs per frequency bin and sampled amplitudes');
        legend('Individual PDFs','Sum of PDFs','Sampled amplitudes','Truncation','Location','northeast');
        grid on;
    end
disp('SUM OF AMPS !!!!!!\n')
disp(sum(amplitudes(:)))

end



% 
% function amplitude = RAS_amp_calc(Hs, spectrum, omega, dOmega, debug)
% % RAS_amp_calc - Generate a single random amplitude using the Random Amplitude Method
% %                with a truncated Rayleigh PDF (max amplitude = 2.5*Hs)
% %
% % Inputs:
% %   Hs       - Significant wave height [m]
% %   spectrum - Spectral density array S(ω)
% %   omega    - Frequency array [rad/s]
% %   dOmega   - Frequency spacing [rad/s]
% %   debug    - (optional) true/false; if true, plots PDF and sampled amplitude
% %
% % Output:
% %   amplitude - Randomly sampled wave amplitude [m]
% 
%     if nargin < 5
%         debug = false;
%     end
% 
%     % Total variance
%     variance = 2 * sum(spectrum .* dOmega);
%     if variance <= 0 || ~isfinite(variance)
%         error('Total variance is non-positive or invalid. Check spectrum and dOmega.');
%     end
% 
%     % Rayleigh scale parameter
%     B = sqrt(variance / 2);
% 
%     % Maximum amplitude (truncation)
%     Amax = 2.5 * Hs;
% 
%     % Truncated Rayleigh using inverse transform method
%     % CDF for Rayleigh: F(a) = 1 - exp(-a^2 / (2*B^2))
%     % Truncate at Amax: normalized CDF = F(a)/F(Amax)
%     u = rand(); % uniform [0,1]
% 
%     % Solve inverse CDF for truncated distribution
%     amplitude = B * sqrt(-2*log(1 - u * (1 - exp(-Amax^2/(2*B^2)))));
% 
% 
%     % Debug plot if requested
%     if debug
%         x = linspace(0, Amax, 500);
%         pdf = (x ./ B^2) .* exp(-x.^2/(2*B^2));
%         % Normalize PDF over truncated range
%         pdf = pdf / trapz(x, pdf);
% 
%         figure;
%         plot(x, pdf, 'LineWidth', 2); hold on;
%         stem(amplitude, interp1(x, pdf, amplitude), 'r', 'LineWidth', 2);
%         xlabel('Amplitude [m]');
%         ylabel('PDF');
%         title('Truncated Rayleigh PDF and sampled amplitude');
%         legend('PDF','Sampled amplitude','Location','northeast');
%         grid on;
%     end
% 
% end





% function [amplitude] = RAS_amp_calc(Hs,spectrum,omega,dOmega)
% %UNTITLED3 Summary of this function goes here
% %   Detailed explanation goes here
% 
% 
% variance = 2 .* spectrum .* dOmega;  % Variance for each frequency
% 
% % Ensure variance is non-negative (if numerical issues arise)
% variance(variance < 0) = 0;
% 
% % Define the range of random variables to evaluate the Rayleigh distribution
% x = linspace(0, 2.5*Hs, 2*length(omega));  % Range of random variable for Rayleigh distribution
% 
% % Preallocate array for the Rayleigh distribution
% rayleigh_pdf = zeros(size(x));
% 
% % Compute the Rayleigh distribution for each variance
% for i = 1:length(variance)
%     if variance(i) > 0
%         B = sqrt(variance(i) / 2);     % Correct Rayleigh scale parameter
%         rayleigh_pdf = rayleigh_pdf + raylpdf(x, B);
%     end
% end
% 
% disp([min(rayleigh_pdf), max(rayleigh_pdf)])
% 
% % Normalize the result
% rayleigh_pdf = rayleigh_pdf / length(variance);
% 
% % Step 1: Compute the cumulative distribution function (CDF)
% rayleigh_cdf = cumtrapz(x, rayleigh_pdf);  % Numerically integrate the PDF to get the CDF
% 
% disp(max(rayleigh_cdf))
% 
% % Normalize the CDF to ensure it goes from 0 to 1
% rayleigh_cdf = rayleigh_cdf / max(rayleigh_cdf);
% 
% % Step 2: Generate a uniformly distributed random number between 0 and 1
% u = rand();
% 
% % Step 3: Use inverse transform sampling to find the random value
% % Find the index of the closest value in the CDF to the random number u
% [~, idx] = min(abs(rayleigh_cdf - u));
% 
% % The corresponding value in x is the random value from the Rayleigh distribution
% amplitude = x(idx);
% 
% 
% end