from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from numpy import prod, squeeze
from scipy.io import loadmat, savemat
from fooof import FOOOF

def unpack(var):
    while prod(var.shape) == 1:
        var = var[0]
    return squeeze(var)

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('inputfilename')
parser.add_argument('outputbasename')
parser.add_argument('--varname', default='timefreq', help='Name of the variable in the MAT-file.')
parser.add_argument('--freq_range', type=float, nargs=2, default=[1, 40], metavar=('lower_bound','upper_bound'), help='Frequency range')
parser.add_argument('--peak_width_limits', type=tuple, nargs=2, default=(1, 12), metavar=('lower_bound','upper_bound'), help='Limits on possible peak width, in Hz.')
parser.add_argument('--max_n_peaks', default=5, help='Maximum number of peaks to fit.')
parser.add_argument('--min_peak_height', default=0.1, help='Absolute threshold for detecting peaks (in log power).')
parser.add_argument('--peak_threshold', default=2, help='Relative threshold for detecting peaks (in SD).')
parser.add_argument('--aperiodic_mode', default='fixed', help='Approach to take for fitting the aperiodic component')
args = parser.parse_args()

print(args.freq_range)

data = loadmat(args.inputfilename)
channels = [dat[0] for dat in unpack(data['timefreq']['label'])]
freqs = unpack(data[args.varname]['freq']).astype('float')
psds = unpack(data[args.varname]['powspctrm']).astype('float')

fm = FOOOF(peak_width_limits=args.peak_width_limits, max_n_peaks=args.max_n_peaks, min_peak_height=args.min_peak_height, peak_threshold=args.peak_threshold, aperiodic_mode=args.aperiodic_mode, verbose=False)
fm.print_settings()

aperiodic = {}
peaks = {}
for ch in range(len(channels)):
    fm.fit(freqs,psds[ch,],args.freq_range)
    # fm.save_report(file_name=filename.stem.replace('freq','fooof_report'), file_path=filename.parent)
    aperiodic[channels[ch]] = fm.get_results().aperiodic_params
    peaks[channels[ch]] = fm.get_results().peak_params

savemat(args.outputbasename + '_aperiodic.mat',aperiodic)
savemat(args.outputbasename + '_peaks.mat',peaks)

