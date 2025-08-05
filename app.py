import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from flask import Flask, render_template, request, jsonify, send_file
import io
import base64
import csv
from scipy import interpolate
import gc  # Garbage collector
import tracemalloc  # For memory profiling (optional)
import os  # Added to check environment variable
import logging  # Added for better error logging


app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 5 * 1024 * 1024  # 5MB max file size

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Check if we're in development mode
is_development = os.environ.get('APP_ENV') == 'development'

# Physical Constants
c = 299.792458  # Speed of light in nm/fs

def gaussian(x, mu, sigma):
    """Defines a Gaussian function."""
    if sigma < 1e-15:
        return np.where(np.abs(x - mu) < 1e-9, 1.0, 0.0)
    return np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def sech2(x, mu, sigma):
    """Defines a Sech^2 function with a FWHM corresponding to a Gaussian with std dev sigma."""
    if sigma < 1e-15:
        return np.where(np.abs(x - mu) < 1e-9, 1.0, 0.0)
    
    # The 'sigma' parameter is for a Gaussian profile. First, find the FWHM it represents.
    fwhm_target = sigma * (2 * np.sqrt(2 * np.log(2)))
    
    # The FWHM of a sech^2 pulse is related to its width parameter 'gamma' by:
    # FWHM = 2 * arccosh(sqrt(2)) * gamma
    # We use the exact constant instead of the approximation '1.76'.
    fwhm_to_gamma_sech = 2 * np.arccosh(np.sqrt(2)) # This is approx 1.7627
    
    # Calculate the correct gamma for the sech^2 shape to match the target FWHM.
    gamma = fwhm_target / fwhm_to_gamma_sech
    
    arg = (x - mu) / gamma
    return 1.0 / np.cosh(arg)**2

def get_fwhm(x, y):
    """Calculates the Full-Width at Half-Maximum (FWHM) of a signal."""
    if np.max(y) < 1e-9:
        return np.nan
    
    y_normalized = y / np.max(y)
    half_max = 0.5
    
    above_half_max = np.where(y_normalized > half_max)[0]
    
    if len(above_half_max) < 2:
        return np.nan
        
    first_idx = above_half_max[0]
    last_idx = above_half_max[-1]
    if first_idx > 0:
        x1, y1 = x[first_idx-1], y_normalized[first_idx-1]
        x2, y2 = x[first_idx], y_normalized[first_idx]
        left_cross = np.interp(half_max, [y1, y2], [x1, x2])
    else:
        left_cross = x[first_idx]
    if last_idx < len(x) - 1:
        x1, y1 = x[last_idx], y_normalized[last_idx]
        x2, y2 = x[last_idx+1], y_normalized[last_idx+1]
        right_cross = np.interp(half_max, [y2, y1], [x2, x1])
    else:
        right_cross = x[last_idx]
        
    return right_cross - left_cross

def get_pulse_center(x, y):
    """Find the center position of a pulse (peak position)."""
    if np.max(y) < 1e-9:
        return 0
    peak_idx = np.argmax(y)
    return x[peak_idx]

def parse_csv_spectrum(file_content, delimiter, skip_rows, x_multiplier, x_exponent):
    """Parse CSV file content and extract wavelength and intensity data."""
    try:
        # Decode file content if it's bytes
        if isinstance(file_content, bytes):
            file_content = file_content.decode('utf-8')
        
        # Parse CSV
        lines = file_content.strip().split('\n')
        
        # Skip header rows
        data_lines = lines[skip_rows:]
        
        wavelengths = []
        intensities = []
        
        for line in data_lines:
            if line.strip():  # Skip empty lines
                parts = line.split(delimiter)
                if len(parts) >= 2:
                    try:
                        # Apply multiplier and exponent to wavelength
                        wavelength = float(parts[0]) * x_multiplier
                        if x_exponent != 1.0:
                            wavelength = wavelength ** x_exponent
                        
                        intensity = float(parts[1])
                        
                        wavelengths.append(wavelength)
                        intensities.append(intensity)
                    except ValueError:
                        continue  # Skip lines that can't be parsed as numbers
        
        if len(wavelengths) < 2:
            raise ValueError("Not enough valid data points in CSV file")
        
        return np.array(wavelengths), np.array(intensities)
    
    except Exception as e:
        raise ValueError(f"Error parsing CSV file: {str(e)}")

def process_imported_spectrum(wavelengths, intensities, cropping_window):
    """Process imported spectrum: normalize, find peak, calculate FWHM, and crop."""
    # Normalize intensity
    intensities = intensities / np.max(intensities)
    
    # Find peak wavelength
    peak_idx = np.argmax(intensities)
    peak_wavelength = wavelengths[peak_idx]
    
    # Calculate FWHM
    fwhm = get_fwhm(wavelengths, intensities)
    
    if np.isnan(fwhm):
        fwhm = 10.0  # Default fallback
    
    # Crop data around peak
    crop_width = cropping_window * fwhm
    mask = np.abs(wavelengths - peak_wavelength) <= crop_width / 2
    
    cropped_wavelengths = wavelengths[mask]
    cropped_intensities = intensities[mask]
    
    # Clean up arrays
    del wavelengths
    del intensities
    
    return cropped_wavelengths, cropped_intensities, peak_wavelength, fwhm

def get_plot_metadata(fig):
    """Extract metadata from matplotlib figure for coordinate display."""
    metadata = {
        'subplots': [],
        'image_size_px': [fig.get_figwidth() * fig.dpi, fig.get_figheight() * fig.dpi]
    }
    
    # Get all axes in the figure
    axes = fig.get_axes()
    
    for ax in axes:
        # Get axis limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        
        # Get axis position in pixels
        bbox = ax.get_window_extent().bounds  # (x, y, width, height)
        
        # Convert y-coordinate to be from top of figure
        y0_from_top = metadata['image_size_px'][1] - (bbox[1] + bbox[3])
        
        subplot_info = {
            'xlim': list(xlim),
            'ylim': list(ylim),
            'pixel_bbox': [bbox[0], y0_from_top, bbox[2], bbox[3]]
        }
        metadata['subplots'].append(subplot_info)
    
    return metadata

def calculate_multi_pulse_properties(pulses, grid_exponent, axis_settings):
    """Calculate and plot properties for multiple laser pulses."""
    
    try:
        # Grid Setup
        num_points = 2**grid_exponent
        
        # Warn if grid size is very large
        if grid_exponent > 20:
            logger.warning(f"Large grid size ({num_points:,} points) may use significant memory")
        
        # Convert all pulse parameters to frequency domain
        pulse_params = []
        omega_centers = []
        omega_widths = []
        imported_spectra = []
        
        for pulse in pulses:
            if pulse.get('shape') == 'imported' and pulse.get('importedData'):
                # Process imported spectrum
                import_data = pulse['importedData']
                wavelengths, intensities = parse_csv_spectrum(
                    import_data['content'],
                    import_data['delimiter'],
                    import_data['skipRows'],
                    import_data['xMultiplier'],
                    import_data['xExponent']
                )
                
                cropped_wavelengths, cropped_intensities, peak_wavelength, fwhm = process_imported_spectrum(
                    wavelengths, intensities, import_data['croppingWindow']
                )
                
                # Convert to frequency domain
                omega0 = 2 * np.pi * c / peak_wavelength
                delta_omega_fwhm = (2 * np.pi * c / (peak_wavelength**2)) * fwhm
                
                imported_spectra.append({
                    'wavelengths': cropped_wavelengths,
                    'intensities': cropped_intensities,
                    'omega0': omega0,
                    'peak_wavelength': peak_wavelength,
                    'fwhm': fwhm
                })
                
                pulse_params.append({
                    'omega0': omega0,
                    'amplitude': pulse['amplitude'],
                    'phi0': pulse['phi0'],
                    'phi1': pulse['phi1'],
                    'phi2': pulse['phi2'],
                    'phi3': pulse['phi3'],
                    'phi4': pulse['phi4'],
                    'shape': 'imported',
                    'imported_index': len(imported_spectra) - 1
                })
                
                omega_centers.append(omega0)
                omega_widths.append(delta_omega_fwhm)
            else:
                # Standard gaussian or sech2 pulse
                lambda0 = pulse['lambda0']
                fwhm_lambda = pulse['fwhm']
                
                omega0 = 2 * np.pi * c / lambda0
                delta_omega_fwhm = (2 * np.pi * c / (lambda0**2)) * fwhm_lambda
                
                # FIX: Calculate sigma for the intensity profile directly from FWHM.
                # This sigma corresponds to a Gaussian intensity profile with the given FWHM.
                # It's used as a basis for both Gaussian and Sech^2 shapes.
                sigma_omega = delta_omega_fwhm / (2 * np.sqrt(2 * np.log(2)))
                
                # FIX: Removed incorrect sigma_E_omega calculation and pass the correct sigma.
                pulse_params.append({
                    'omega0': omega0,
                    'sigma_omega': sigma_omega,
                    'amplitude': pulse['amplitude'],
                    'phi0': pulse['phi0'],
                    'phi1': pulse['phi1'],
                    'phi2': pulse['phi2'],
                    'phi3': pulse['phi3'],
                    'phi4': pulse['phi4'],
                    'shape': pulse['shape']
                })
                
                omega_centers.append(omega0)
                omega_widths.append(delta_omega_fwhm)
        
        # Create frequency grid
        omega_center = np.mean(omega_centers)
        omega_span = 10 * max(omega_widths + [np.ptp(omega_centers), 0.1])
        omega_axis = np.linspace(omega_center - omega_span / 2, 
                                omega_center + omega_span / 2, num_points)
        
        # Calculate individual pulses and total field
        E_omega_total = np.zeros_like(omega_axis, dtype=complex)
        individual_pulses = []
        pulse_centers = []
        pulse_fwhms = []
        
        for i, params in enumerate(pulse_params):
            if params['shape'] == 'imported':
                # Handle imported spectrum
                imported_data = imported_spectra[params['imported_index']]
                
                # Convert wavelengths to omega for interpolation
                imported_omega = 2 * np.pi * c / imported_data['wavelengths']
                
                # Sort for interpolation (omega is inversely related to wavelength)
                sort_idx = np.argsort(imported_omega)
                imported_omega_sorted = imported_omega[sort_idx]
                imported_intensity_sorted = imported_data['intensities'][sort_idx]
                
                # Interpolate to our frequency grid
                interp_func = interpolate.interp1d(
                    imported_omega_sorted, 
                    imported_intensity_sorted,
                    bounds_error=False,
                    fill_value=0.0,
                    kind='cubic'
                )
                
                spectral_shape = interp_func(omega_axis)
                spectral_shape[spectral_shape < 0] = 0  # Remove any negative values from interpolation
                
                # Normalize to preserve energy
                if np.max(spectral_shape) > 0:
                    spectral_shape = spectral_shape / np.max(spectral_shape)
            else:
                # Select spectral shape function
                # FIX: Use the corrected 'sigma_omega' parameter instead of 'sigma_E_omega'
                if params['shape'] == 'sech2':
                    spectral_shape = sech2(omega_axis, params['omega0'], params['sigma_omega'])
                else:  # gaussian
                    spectral_shape = gaussian(omega_axis, params['omega0'], params['sigma_omega'])
            
            # Calculate phase
            phase = (params['phi0'] +
                    params['phi1'] * (omega_axis - params['omega0']) +
                    (params['phi2'] / 2) * (omega_axis - params['omega0'])**2 +
                    (params['phi3'] / 6) * (omega_axis - params['omega0'])**3 +
                    (params['phi4'] / 24) * (omega_axis - params['omega0'])**4)
            
            # The spectral_shape variable now correctly represents the INTENSITY profile.
            # E_omega is calculated from its square root.
            E_omega = np.sqrt(params['amplitude']) * np.sqrt(spectral_shape) * np.exp(-1j * phase)
            E_omega_total += E_omega
            
            # Calculate temporal intensity for individual pulse
            E_t = np.fft.ifft(np.fft.ifftshift(E_omega))
            I_t = np.fft.fftshift(np.abs(E_t)**2)
            
            # Rescale to match input amplitude
            if np.max(I_t) > 1e-9:
                I_t = I_t / np.max(I_t) * params['amplitude']
            
            individual_pulses.append(I_t)
            
            # Clean up intermediate arrays
            del E_omega
            del E_t
            del spectral_shape
        
        # Calculate spectral intensity
        spectral_intensity = np.abs(E_omega_total)**2
        
        # Define time and wavelength axes
        d_omega = omega_axis[1] - omega_axis[0]
        time_step = 2 * np.pi / (num_points * d_omega)
        time_axis = (np.arange(num_points) - num_points / 2) * time_step
        lambda_axis_plot = 2 * np.pi * c / omega_axis
        
        # Calculate total temporal intensity and autocorrelation
        E_t_total = np.fft.ifft(np.fft.ifftshift(E_omega_total))
        I_t = np.abs(E_t_total)**2
        del E_t_total  # Free memory
        
        I_t_fft = np.fft.fft(I_t)
        autocorr_fft = np.fft.ifft(np.abs(I_t_fft)**2)
        del I_t_fft  # Free memory
        
        autocorr = np.fft.fftshift(autocorr_fft).real
        del autocorr_fft  # Free memory
        
        I_t_total_shifted = np.fft.fftshift(I_t)
        del I_t  # Free memory
        
        # Normalize the total temporal intensity
        if np.max(I_t_total_shifted) > 1e-9:
            I_t_total_normalized = I_t_total_shifted / np.max(I_t_total_shifted)
        else:
            I_t_total_normalized = I_t_total_shifted
        
        # Normalize the spectral intensity
        if np.max(spectral_intensity) > 1e-9:
            spectral_intensity_normalized = spectral_intensity / np.max(spectral_intensity)
        else:
            spectral_intensity_normalized = spectral_intensity
        
        # Normalize the autocorrelation
        if np.max(autocorr) > 1e-9:
            autocorr_normalized = autocorr / np.max(autocorr)
        else:
            autocorr_normalized = autocorr
        
        # Calculate FWHMs and centers for zoom
        for i, I_t in enumerate(individual_pulses):
            fwhm = get_fwhm(time_axis, I_t)
            center = get_pulse_center(time_axis, I_t)
            pulse_fwhms.append(fwhm)
            pulse_centers.append(center)
        
        # Calculate the default shared x-axis limits for temporal plots
        default_temporal_xlim = None
        if any(not np.isnan(f) for f in pulse_fwhms):
            valid_centers = [c for c, f in zip(pulse_centers, pulse_fwhms) if not np.isnan(f)]
            valid_fwhms = [f for f in pulse_fwhms if not np.isnan(f)]
            
            if valid_centers:
                min_center = min(valid_centers)
                max_center = max(valid_centers)
                max_fwhm = max(valid_fwhms)
                
                zoom_margin = 3 * max_fwhm  # 6x FWHM total (3x on each side)
                x_min = min_center - zoom_margin
                x_max = max_center + zoom_margin
                
                default_temporal_xlim = (x_min, x_max)
        
        # Generate plots with zoom
        plt.rcParams['font.size'] = 10
        fig, axes = plt.subplots(4, 1, figsize=(10, 16))
        
        # Plot 1: Individual pulse temporal intensities (NOT NORMALIZED)
        ax_temp = axes[0]
        for i, I_t in enumerate(individual_pulses):
            fwhm = pulse_fwhms[i]
            label = f'Pulse {i+1}'
            if not np.isnan(fwhm):
                label += f', FWHM: {fwhm:.2f} fs'
            ax_temp.plot(time_axis, I_t, linewidth=2, label=label)
        
        # Set x-axis limits for temporal plots
        if axis_settings['temporal']['auto'] and default_temporal_xlim:
            ax_temp.set_xlim(default_temporal_xlim)
        elif not axis_settings['temporal']['auto']:
            ax_temp.set_xlim(axis_settings['temporal']['min'], axis_settings['temporal']['max'])
        
        ax_temp.set_title('Temporal Intensity of Individual Pulses', fontsize=14)
        ax_temp.set_xlabel('Time (fs)')
        ax_temp.set_ylabel('Intensity (a.u.)')
        ax_temp.grid(True, linestyle='--', alpha=0.6)
        ax_temp.legend()
        
        # Plot 2: Total temporal intensity (NORMALIZED)
        ax_interf = axes[1]
        fwhm_total = get_fwhm(time_axis, I_t_total_normalized)
        
        # Create title with FWHM
        title_total = 'Total Temporal Intensity (Interference)'
        if not np.isnan(fwhm_total):
            title_total += f', FWHM: {fwhm_total:.2f} fs'
        
        ax_interf.plot(time_axis, I_t_total_normalized, color='green', linewidth=2)
        
        # Use the same x-axis limits as individual pulses
        if axis_settings['temporal']['auto'] and default_temporal_xlim:
            ax_interf.set_xlim(default_temporal_xlim)
        elif not axis_settings['temporal']['auto']:
            ax_interf.set_xlim(axis_settings['temporal']['min'], axis_settings['temporal']['max'])
        
        ax_interf.set_title(title_total, fontsize=14)
        ax_interf.set_xlabel('Time (fs)')
        ax_interf.set_ylabel('Normalized Intensity')
        ax_interf.grid(True, linestyle='--', alpha=0.6)
        
        # Plot 3: Spectral intensity (NORMALIZED)
        ax_spec = axes[2]
        sort_indices = np.argsort(lambda_axis_plot)
        sorted_lambda = lambda_axis_plot[sort_indices]
        sorted_intensity_normalized = spectral_intensity_normalized[sort_indices]
        
        # Calculate spectral FWHM on normalized data
        spec_fwhm = get_fwhm(sorted_lambda, sorted_intensity_normalized)
        
        # Create title with FWHM
        title_spectral = 'Spectral Intensity'
        if not np.isnan(spec_fwhm):
            title_spectral += f', FWHM: {spec_fwhm:.2f} nm'
        
        ax_spec.plot(sorted_lambda, sorted_intensity_normalized, color='royalblue', linewidth=2)
        
        # Set x-axis limits for spectral plot
        if axis_settings['spectral']['auto']:
            if not np.isnan(spec_fwhm) and spec_fwhm > 0:
                peak_idx = np.argmax(sorted_intensity_normalized)
                peak_lambda = sorted_lambda[peak_idx]
                zoom_margin = 3 * spec_fwhm  # 6x FWHM total
                ax_spec.set_xlim(peak_lambda - zoom_margin, peak_lambda + zoom_margin)
        else:
            ax_spec.set_xlim(axis_settings['spectral']['min'], axis_settings['spectral']['max'])
        
        ax_spec.set_title(title_spectral, fontsize=14)
        ax_spec.set_xlabel('Wavelength (nm)')
        ax_spec.set_ylabel('Normalized Intensity')
        ax_spec.grid(True, linestyle='--', alpha=0.6)
        
        # Plot 4: Autocorrelation (NORMALIZED)
        ax_auto = axes[3]
        autoco_fwhm = get_fwhm(time_axis, autocorr_normalized)
        title_text = 'Intensity Autocorrelation Trace'
        if not np.isnan(autoco_fwhm):
            title_text += f', FWHM: {autoco_fwhm:.2f} fs'
        
        ax_auto.plot(time_axis, autocorr_normalized, color='crimson', linewidth=2)
        
        # Set x-axis limits for autocorrelation
        if axis_settings['autocorr']['auto']:
            if default_temporal_xlim:
                # Use the maximum absolute value from the shared limits
                max_abs = max(abs(default_temporal_xlim[0]), abs(default_temporal_xlim[1]))
                ax_auto.set_xlim(-max_abs, max_abs)
        else:
            ax_auto.set_xlim(axis_settings['autocorr']['min'], axis_settings['autocorr']['max'])
        
        ax_auto.set_title(title_text, fontsize=14)
        ax_auto.set_xlabel('Time Delay (fs)')
        ax_auto.set_ylabel('Normalized Intensity')
        ax_auto.grid(True, linestyle='--', alpha=0.6)
        
        plt.tight_layout()
        
        # Extract plot metadata before saving
        plot_metadata = get_plot_metadata(fig)
        
        # Save plot to memory buffer
        img_buffer = io.BytesIO()
        plt.savefig(img_buffer, format='png', dpi=100, bbox_inches='tight')
        plt.close(fig)  # Explicitly close the figure
        
        # Clear matplotlib's internal references
        plt.clf()
        plt.cla()
        
        img_buffer.seek(0)
        
        # Convert to base64 for embedding in HTML
        img_base64 = base64.b64encode(img_buffer.getvalue()).decode()
        
        # Clean up the buffer
        img_buffer.close()
        
        # Clean up large arrays
        del omega_axis
        del E_omega_total
        del spectral_intensity
        del spectral_intensity_normalized
        del time_axis
        del lambda_axis_plot
        del I_t_total_shifted
        del I_t_total_normalized
        del autocorr
        del autocorr_normalized
        del individual_pulses
        
        # Force garbage collection
        gc.collect()
        
        return img_base64, plot_metadata
        
    except Exception as e:
        # Clean up on error
        plt.close('all')
        gc.collect()
        logger.error(f"Error in calculate_multi_pulse_properties: {str(e)}", exc_info=True)
        raise e

# Intro page route
@app.route('/')
def intro_page():
    """Render the intro page."""
    return render_template('intro.html')

# Tool page route
@app.route('/tool')
def tool_page():
    """Render the calculator tool page."""
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    try:
        # Check if request has JSON data
        if not request.is_json:
            logger.error("Request is not JSON")
            return jsonify({'error': 'Request must be JSON'}), 400
        
        data = request.get_json()
        if data is None:
            logger.error("No JSON data in request")
            return jsonify({'error': 'No data provided'}), 400
        
        # Validate grid exponent
        try:
            grid_exponent = int(data.get('gridExponent', 16))
        except (ValueError, TypeError):
            logger.error("Invalid grid exponent value")
            return jsonify({'error': 'Invalid grid exponent value'}), 400
            
        # Only apply grid size limit in production mode
        if not is_development and grid_exponent > 22:
            logger.error(f"Grid exponent {grid_exponent} exceeds maximum allowed value of 22")
            return jsonify({'error': 'Grid exponent cannot exceed 22'}), 400
        
        # Parse pulses
        pulses = []
        pulses_data = data.get('pulses', [])
        if not pulses_data:
            logger.error("No pulses defined in request")
            return jsonify({'error': 'No pulses defined'}), 400
            
        for pulse_data in pulses_data:
            try:
                pulse = {
                    'lambda0': float(pulse_data.get('lambda0', 800)),
                    'fwhm': float(pulse_data.get('fwhm', 10)),
                    'amplitude': float(pulse_data.get('amplitude', 1)),
                    'phi0': float(pulse_data.get('phi0', 0)),
                    'phi1': float(pulse_data.get('phi1', 0)),
                    'phi2': float(pulse_data.get('phi2', 0)),
                    'phi3': float(pulse_data.get('phi3', 0)),
                    'phi4': float(pulse_data.get('phi4', 0)),
                    'shape': pulse_data.get('shape', 'gaussian')
                }
                
                # Add imported data if shape is 'imported'
                if pulse['shape'] == 'imported' and pulse_data.get('importedData'):
                    pulse['importedData'] = pulse_data['importedData']
                
                pulses.append(pulse)
            except (ValueError, TypeError) as e:
                logger.error(f"Error parsing pulse data: {str(e)}")
                return jsonify({'error': f'Invalid pulse data: {str(e)}'}), 400
        
        # Parse axis settings
        try:
            axis_settings = {
                'temporal': {
                    'auto': data.get('axisSettings', {}).get('temporal', {}).get('auto', True),
                    'min': float(data.get('axisSettings', {}).get('temporal', {}).get('min', -1000)),
                    'max': float(data.get('axisSettings', {}).get('temporal', {}).get('max', 1000))
                },
                'spectral': {
                    'auto': data.get('axisSettings', {}).get('spectral', {}).get('auto', True),
                    'min': float(data.get('axisSettings', {}).get('spectral', {}).get('min', 700)),
                    'max': float(data.get('axisSettings', {}).get('spectral', {}).get('max', 900))
                },
                'autocorr': {
                    'auto': data.get('axisSettings', {}).get('autocorr', {}).get('auto', True),
                    'min': float(data.get('axisSettings', {}).get('autocorr', {}).get('min', -1000)),
                    'max': float(data.get('axisSettings', {}).get('autocorr', {}).get('max', 1000))
                }
            }
        except (ValueError, TypeError) as e:
            logger.error(f"Error parsing axis settings: {str(e)}")
            return jsonify({'error': f'Invalid axis settings: {str(e)}'}), 400
        
        # Calculate and generate plots
        try:
            plot_base64, plot_metadata = calculate_multi_pulse_properties(pulses, grid_exponent, axis_settings)
        except Exception as e:
            logger.error(f"Error in calculation: {str(e)}", exc_info=True)
            return jsonify({'error': f'Calculation error: {str(e)}'}), 500
        
        # Force garbage collection after calculation
        gc.collect()
        
        return jsonify({
            'success': True,
            'plot_data': plot_base64,
            'plot_metadata': plot_metadata
        })
        
    except Exception as e:
        logger.error(f"Unexpected error in calculate route: {str(e)}", exc_info=True)
        gc.collect()  # Clean up on error
        return jsonify({'error': f'Unexpected error: {str(e)}'}), 500

@app.errorhandler(413)
def request_entity_too_large(error):
    return jsonify({'error': 'File size exceeds 5 MB limit'}), 413

# Optional: Add a route to check memory usage
@app.route('/memory-stats')
def memory_stats():
    import psutil
    import os
    
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    
    # Get current memory snapshot if tracemalloc is enabled
    stats = {}
    if tracemalloc.is_tracing():
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')[:10]
        stats['top_allocations'] = [
            f"{stat.filename}:{stat.lineno}: {stat.size / 1024 / 1024:.1f} MB"
            for stat in top_stats
        ]
    
    return jsonify({
        'rss_mb': memory_info.rss / 1024 / 1024,
        'vms_mb': memory_info.vms / 1024 / 1024,
        'available_mb': psutil.virtual_memory().available / 1024 / 1024,
        'percent': psutil.virtual_memory().percent,
        **stats
    })

if __name__ == '__main__':
    app.run(debug=True, port=5004)