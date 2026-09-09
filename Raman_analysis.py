# Raman_analysis.py - LIG Raman spectrum analysis (D, G and 2D peaks).
#
# Algorithm credits:
#
# Baseline correction - `baseline_als` implements Asymmetric Least Squares
#   (AsLS) smoothing from:
#     P. H. C. Eilers and H. F. M. Boelens, "Baseline Correction with
#     Asymmetric Least Squares Smoothing", Leiden University Medical Centre
#     report, 2005.
#   The sparse-matrix Python formulation used here follows the widely
#   circulated community implementation of that paper.
#
# Spike removal - `modified_z_score` / `fixer` implement the modified
#   z-score despiking algorithm from:
#     D. A. Whitaker and K. Hayes, "A simple algorithm for despiking Raman
#     spectra", Chemometrics and Intelligent Laboratory Systems, 179 (2018),
#     82-84.
#   The 0.6745 factor scales the median absolute deviation to an estimate of
#   the standard deviation; the threshold of 7 and the local-mean repair of
#   flagged points are as described in that paper.

import numpy as np
from scipy.signal import savgol_filter, find_peaks, peak_widths
from scipy import sparse
from scipy.sparse.linalg import spsolve
import plotly.graph_objects as go

def modified_z_score(ys):
    ysb = np.diff(ys)  # Differentiated intensity values
    median_y = np.median(ysb)  # Median of the intensity values
    median_absolute_deviation_y = np.median([np.abs(y - median_y) for y in ysb])  # median_absolute_deviation of the differentiated intensity values
    if not np.isfinite(median_absolute_deviation_y) or median_absolute_deviation_y == 0:
        # A flat or near-flat stretch has zero MAD. Dividing by it gives NaN,
        # abs(NaN) > threshold is False, and despiking silently switches off.
        # A constant signal has no spikes to find, so say so explicitly.
        return np.zeros_like(ysb, dtype=float)
    modified_z_scores = [0.6745 * (y - median_y) / median_absolute_deviation_y for y in ysb]  # median_absolute_deviation modified z scores
    return modified_z_scores

def fixer(y, ma):
    threshold = 7  # binarization threshold
    spikes = abs(np.array(modified_z_score(y))) > threshold
    y_out = y.copy()
    for i in np.arange(len(spikes)):
        if spikes[i] != 0:
            w = np.arange(max(0, i-ma), min(len(spikes), i+1+ma))  # Clamp to len(spikes) == len(y)-1 (np.diff) so spikes[w] stays in bounds
            we = w[spikes[w] == 0] if len(w) > 0 else []
            if len(we) > 0:
                y_out[i] = np.mean(y[we])
    return y_out

def baseline_als(y, lam, p, niter=100):
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L-2), dtype=float)
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z.tocsc(), w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

def full_analysis2(intensity, x_range, l=10000000, p=0.05, w=30, pol=2, distance=40, nb_of_max=3, ma=5):
    mix_spectrum = intensity.copy()

    # Create mask for filtering
    mask = x_range > 1000
    mix_spectrum = mix_spectrum[mask]
    x_range_filtered = x_range[mask]

    # Despike before baselining: cosmic-ray spikes would otherwise drag the
    # AsLS baseline and can be picked up as peaks (Whitaker & Hayes 2018)
    if len(mix_spectrum) > 2:
        mix_spectrum = fixer(mix_spectrum, ma)

    estimated_baselined = baseline_als(mix_spectrum, l, p)
    baselined_spectrum = mix_spectrum - estimated_baselined
    smoothed_spectrum = savgol_filter(baselined_spectrum, w, polyorder=pol, deriv=0)

    # Create search area mask
    search_mask = ((x_range_filtered > 1000) & (x_range_filtered < 1800)) | ((x_range_filtered > 2500) & (x_range_filtered < 2800))
    search_area = smoothed_spectrum[search_mask]
    x_area = x_range_filtered[search_mask]

    # Check if search area has data
    if len(search_area) == 0:
        print("Warning: No data found in search area")
        return x_range_filtered, x_area, search_area, smoothed_spectrum, np.array([]), np.array([]), np.array([]), estimated_baselined, baselined_spectrum

    peak_indices, _ = find_peaks(search_area, height=0, distance=distance, rel_height=0.5)

    if len(peak_indices) == 0:
        print("Warning: No peaks found")
        return x_range_filtered, x_area, search_area, smoothed_spectrum, np.array([]), np.array([]), np.array([]), estimated_baselined, baselined_spectrum

    peak_heights = search_area[peak_indices]

    # Select top peaks
    nb_of_max = min(nb_of_max, len(peak_heights))  # Ensure we don't ask for more peaks than available
    idn = np.argpartition(peak_heights, -nb_of_max)[-nb_of_max:]
    peaks = peak_indices[idn]

    results = peak_widths(search_area, peaks, rel_height=0.5)
    peak_positions = x_area[peaks]
    peak_widths_values = results[0] * (x_area[1] - x_area[0]) if len(x_area) > 1 else np.array([])  # Convert indices to actual width values

    return x_range_filtered, x_area, search_area, smoothed_spectrum, peaks, peak_positions, peak_widths_values, estimated_baselined, baselined_spectrum

def parse_two_column_txt(filepath):
    """
    Parses a text file with two numeric columns into two numpy arrays.

    The delimiter is detected automatically: tabs, commas, semicolons and runs
    of whitespace are all accepted. Columns beyond the first two are ignored.
    """
    column1_data = []
    column2_data = []

    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().replace('\t', ' ').replace(',', ' ').replace(';', ' ').split()
                if len(parts) >= 2:
                    try:
                        column1_data.append(float(parts[0]))
                        column2_data.append(float(parts[1]))
                    except ValueError:
                        print(f"Warning: Could not convert to float: {line.strip()}")
                else:
                    print(f"Warning: Skipping malformed line: {line.strip()}")
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return np.array([]), np.array([])
    except Exception as e:
        print(f"An error occurred: {e}")
        return np.array([]), np.array([])

    return np.array(column1_data), np.array(column2_data)

# Main plotting code
def plot_spectrum_analysis(filepath):
    # Parse the data
    x, y = parse_two_column_txt(filepath)

    if len(x) == 0 or len(y) == 0:
        print("Error: No data loaded from file")
        return

    # Call the full_analysis function
    results = full_analysis2(y, x)
    x_range, x_area, search_area, smoothed_spectrum, peaks, peak_positions, peak_widths_values, estimated_baseline, baselined_spectrum = results

    # Create the plot
    fig = go.Figure()

    # Add the original spectrum
    fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name='Original Spectrum', line=dict(color='blue', width=1)))

    # Add the baseline
    if len(estimated_baseline) > 0:
        fig.add_trace(go.Scatter(x=x_range, y=estimated_baseline, mode='lines', name='Estimated Baseline', line=dict(color='green', width=2)))

    # Add the baseline-corrected smoothed spectrum
    if len(smoothed_spectrum) > 0:
        fig.add_trace(go.Scatter(x=x_range, y=smoothed_spectrum, mode='lines', name='Smoothed Spectrum (Baseline Corrected)', line=dict(color='red', width=2)))

    print(f"Number of peaks found: {len(peaks)}")
    print(f"Peak positions: {peak_positions}")
    print(f"Peak widths: {peak_widths_values}")

    # Identify peaks by their characteristic positions for LIG
    peak_assignments = {}
    peak_intensities = {}

    if len(peaks) > 0:
        for i, peak_pos in enumerate(peak_positions):
            peak_height = search_area[peaks[i]]
            peak_intensities[i] = peak_height

            # Assign peaks based on typical Raman shifts for graphene
            if 1250 <= peak_pos <= 1450:  # D peak
                peak_assignments[i] = 'D'
                print(f"D peak identified at {peak_pos:.0f} cm⁻¹ with intensity {peak_height:.3f}")
            elif 1500 <= peak_pos <= 1650:  # G peak
                peak_assignments[i] = 'G'
                print(f"G peak identified at {peak_pos:.0f} cm⁻¹ with intensity {peak_height:.3f}")
            elif 2550 <= peak_pos <= 2850:  # 2D peak
                peak_assignments[i] = '2D'
                print(f"2D peak identified at {peak_pos:.0f} cm⁻¹ with intensity {peak_height:.3f}")
            else:
                peak_assignments[i] = f'Peak@{peak_pos:.0f}'

    # Calculate peak ratios for LIG characterization
    I_D = None
    I_G = None
    I_2D = None

    for i, assignment in peak_assignments.items():
        if assignment == 'D':
            I_D = peak_intensities[i]
        elif assignment == 'G':
            I_G = peak_intensities[i]
        elif assignment == '2D':
            I_2D = peak_intensities[i]

    # Print peak ratios
    print("\n--- Peak Ratio Analysis for LIG ---")
    if I_D is not None and I_G is not None:
        Id_Ig_ratio = I_D / I_G
        print(f"I_D/I_G ratio: {Id_Ig_ratio:.3f}")
    else:
        print("I_D/I_G ratio: Cannot calculate (D or G peak not identified)")

    if I_2D is not None and I_G is not None:
        I2d_Ig_ratio = I_2D / I_G
        print(f"I_2D/I_G ratio: {I2d_Ig_ratio:.3f}")
    else:
        print("I_2D/I_G ratio: Cannot calculate (2D or G peak not identified)")

    # Add FWHM measurements and peak labels for peaks
    if len(peaks) > 0 and len(peak_widths_values) > 0:
        for i in range(len(peaks)):
            peak_x = x_area[peaks[i]]
            peak_height = search_area[peaks[i]]
            peak_width = peak_widths_values[i]
            left = peak_x - peak_width / 2
            right = peak_x + peak_width / 2

            # Get peak assignment
            peak_label = peak_assignments.get(i, f'Peak{i+1}')

            # Add FWHM line
            fig.add_shape(
                type='line',
                x0=left, x1=right,
                y0=peak_height/2, y1=peak_height/2,
                line=dict(color='orange', width=3)
            )

            # Add FWHM annotation
            fig.add_annotation(
                x=(left + right) / 2,
                y=peak_height/2 + max(smoothed_spectrum) * 0.05,
                text=f"FWHM: {peak_width:.1f} cm⁻¹",
                showarrow=False,
                font=dict(family="Arial", size=9, color="orange")
            )

            # Add peak label annotation
            fig.add_annotation(
                x=peak_x,
                y=peak_height + max(smoothed_spectrum) * 0.08,
                text=f"{peak_label}<br>({peak_x:.0f} cm⁻¹)",
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=2,
                arrowcolor="red",
                font=dict(family="Arial", size=11, color="red"),
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="red",
                borderwidth=1
            )

            # Add vertical markers at FWHM boundaries
            fig.add_shape(
                type='line',
                x0=left, x1=left,
                y0=peak_height/2 - max(smoothed_spectrum) * 0.02,
                y1=peak_height/2 + max(smoothed_spectrum) * 0.02,
                line=dict(color='orange', width=3)
            )
            fig.add_shape(
                type='line',
                x0=right, x1=right,
                y0=peak_height/2 - max(smoothed_spectrum) * 0.02,
                y1=peak_height/2 + max(smoothed_spectrum) * 0.02,
                line=dict(color='orange', width=3)
            )

            # Mark peak position
            fig.add_trace(go.Scatter(
                x=[peak_x], y=[peak_height],
                mode='markers',
                marker=dict(color='red', size=8, symbol='diamond'),
                name=f'{peak_label} ({peak_x:.0f} cm⁻¹)',
                showlegend=True
            ))

    # Add ratio information to the plot title
    ratio_text = ""
    if I_D is not None and I_G is not None:
        ratio_text += f"I_D/I_G = {I_D/I_G:.3f}"
    if I_2D is not None and I_G is not None:
        if ratio_text:
            ratio_text += " | "
        ratio_text += f"I_2D/I_G = {I_2D/I_G:.3f}"

    title_text = 'LIG Raman Spectrum Analysis with Peak Ratios'
    if ratio_text:
        title_text += f'<br><sub>{ratio_text}</sub>'

    # Configure the layout
    fig.update_layout(
        xaxis=dict(title=dict(text='Wavenumber (cm⁻¹)', font=dict(size=14))),
        yaxis=dict(title=dict(text='Intensity (a.u.)', font=dict(size=14))),
        title=dict(text=title_text, font=dict(size=16)),
        legend=dict(x=0.02, y=0.98),
        width=1000,
        height=650
    )

    # Display the figure
    fig.show()

    return fig

# Example usage:
# Uncomment the line below to run the analysis
# fig = plot_spectrum_analysis("/content/control_01.txt")
