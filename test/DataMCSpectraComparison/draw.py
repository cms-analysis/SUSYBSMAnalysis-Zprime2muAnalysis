#!/usr/bin/env python

import sys, os, glob
from collections import defaultdict
from pprint import pprint
from optparse import OptionParser

# We have to optparse before ROOT does, or else it will eat our
# options (at least -h/--help gets eaten). So don't move this!
parser = OptionParser()
parser.add_option('-d', '--histo-dir', dest='histo_dir', default='data/Run2012MuonsOnly',
                  help='Directory containing the input files for the data. Default is %default. The files expected to be in this directory are ana_datamc_data.root, the ROOT file containing the input histograms, and ana_datamc_data.lumi, the log file from the output of LumiCalc. Optionally the directory can contain a link to a directory for MC histogram ROOT files; the link/directory must be named "mc".')
parser.add_option('--no-print-table', action='store_false', dest='print_table', default=True,
                  help='Do not print out the ASCII table of event counts in specified mass ranges.')
parser.add_option('--no-save-plots', action='store_false', dest='save_plots', default=True,
                  help='Do not save the plots drawn.')
parser.add_option('--luminosity', dest='override_int_lumi', type='float',
                  help='Set the integrated luminosity manually (in 1/pb) rather than attempting to get it from the LumiCalc log file.')
parser.add_option('--no-lumi-rescale', action='store_false', dest='rescale_lumi', default=True,
                  help='Do not rescale the luminosity.')
parser.add_option('--for-rescale-factors', action='store_true', dest='for_rescale_factors', default=False,
                  help='Just print the tables for the Z peak counts to determine the luminosity rescaling factors (implies --no-lumi-rescale and --no-save-plots).')
parser.add_option('--lumi_syst_frac', dest='lumi_syst_frac', type='float', default=0.06,
                  help='Set the systematic uncertainty for the luminosity (as a relative value). Default is %default.')
parser.add_option('--no-draw-zprime', action='store_false', dest='draw_zprime', default=True,
                  help='Do not draw the Z\' curve.')
parser.add_option('--no-poisson-intervals', action='store_false', dest='use_poisson_intervals', default=True,
                  help='Do not use Poisson 68% CL intervals for the error bars on the data points.')
parser.add_option('--no-overflow-in-last-bin', action='store_false', dest='put_overflow_in_last_bin', default=True,
                  help='Do not add the overflow amount to the last bin.')
parser.add_option('--no-joins', action='store_false', dest='do_joins', default=True,
                  help='Do not lump together the MC contributions from different samples.')
parser.add_option('--no-join-ttbar', action='store_false', dest='join_ttbar_and_other_prompt_leptons', default=True,
                  help='Do not lump ttbar and other prompt leptons contributions together.')
parser.add_option('--exclude-sample', action='append', dest='exclude_samples',
                  help='Specify a sample not to use (by name, e.g. wjets). To exclude more than one sample, give this option more than once.')
parser.add_option('--include-quantity', action='append', dest='include_quantities',
                  help='If specified, will override the default list of quantities to compare in favor of the specified one. Specify this option more than once to use multiple quantities to compare.')
parser.add_option('--include-dilepton', action='append', dest='include_dileptons',
                  help='If specified, will override the default list of dileptons in favor of the specified one. Specify this option more than once to use multiple dileptons.')
parser.add_option('--include-cutset', action='append', dest='include_cutsets',
                  help='If specified, will override the default list of cutsets in favor of the specified one. Specify this option more than once to use multiple cutsets.')
parser.add_option('--include-massrange', action='append', dest='include_mass_ranges_for_table',
                  help='If specified, will override the default list of mass ranges for the ASCII table in favor of the specified one. Specify this option more than once to use multiple mass ranges.')
parser.add_option('--plot-dir-tag', action='store',
                  help='Adds argument to plot_dir path, useful for tagging the current version.')
parser.add_option('--plot-size', default='900,600',
                  help='The canvas size for drawing the plots.')
parser.add_option('--no-guess-yrange', action='store_false', dest='guess_yrange', default=True,
                  help='Do not try to guess out the y-axis range for plots, using instead the fixed values in the script.')
options, args = parser.parse_args()
#pprint(options) ; raise 1

# Done with option parsing, now can import things that import ROOT.

import SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples as MCSamples
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import cumulative_histogram, get_integral, move_below_into_bin, move_above_into_bin, plot_saver, poisson_intervalize, real_hist_max, real_hist_min, set_zp2mu_style, ROOT

class Drawer:
    # Terminology:
    #
    # "cut set": the label for the set of cuts applied on the input
    # data. These correspond to substrings of the folder names in the
    # input file. E.g. "OurNewMuonsPlusMuonsMinus" folder -> "OurNew"
    # cut set.
    #
    # "join group": a set of MC samples that will be lumped together
    # in the plotted histograms (by color and in the legend) and in
    # the ASCII table.
    #
    # "nice-name": the string used for e.g. legend captions, usually
    # TLatex-ified. E.g. t#bar{t}.
    #
    # "sample name": the short name for the MC sample. For the sample
    # objects in MCSamples imported above, this is just
    # sample.name. E.g. zmumu.
    #
    # JMTBAD rest of it

    def __init__(self, options):
        self.histo_dir = options.histo_dir
        self.print_table = options.print_table
        self.save_plots = options.save_plots
        self.rescale_lumi = options.rescale_lumi
        self.lumi_syst_frac = options.lumi_syst_frac
        self.draw_zprime = options.draw_zprime
        self.use_poisson_intervals = options.use_poisson_intervals
        self.put_overflow_in_last_bin = options.put_overflow_in_last_bin
        self.do_joins = options.do_joins
        self.join_ttbar_and_other_prompt_leptons = options.join_ttbar_and_other_prompt_leptons
        self.guess_yrange = options.guess_yrange

        if not os.path.isdir(self.histo_dir):
            raise ValueError('histo_dir %s is not a directory' % self.histo_dir)

        self.data_fn = os.path.join(self.histo_dir, 'ana_datamc_data.root')
        if not os.path.isfile(self.data_fn):
            raise ValueError('data_fn %s is not a file' % self.data_fn)
        self.data_f = ROOT.TFile(self.data_fn)
        if not self.data_f.IsOpen() or self.data_f.IsZombie():
            raise ValueError('data_fn %s is not a ROOT file' % self.data_fn)
        
        # We look for a separate mc dir in our histo_dir first, for
        # the case where we want to use a specific set of MC files for
        # a given histo_dir. If none there, check if there is one
        # specified in options.  Else, use the general mc dir from the
        # current directory.
        self.mc_dir = os.path.join(self.histo_dir, 'mc')
        if not os.path.isdir(self.mc_dir):
            if hasattr(options, 'mc_dir'):
                self.mc_dir = options.mc_dir
            else:
                self.mc_dir = 'mc'
        if not os.path.isdir(self.mc_dir):
            raise ValueError('mc_dir %s is not a directory' % self.mc_dir)

        # If the options override the int_lumi, take that; otherwise
        # try to get the int_lumi from the file corresponding to
        # data_fn.
        if type(options.override_int_lumi) == float:
            self.int_lumi = options.override_int_lumi
        else:
            lumi_fn = self.data_fn.replace('.root', '.lumi')
            self.int_lumi = self.parse_lumi_from_log(lumi_fn)
            if self.int_lumi is None:
                raise ValueError('int_lumi could not be parsed from file %s' % lumi_fn)
        if self.int_lumi < 0:
            raise ValueError('int_lumi %f makes no sense' % self.int_lumi)

        # Use all samples specified in MCSamples that are not
        # requested to be dropped in the options. The sample objects
        # already-set values (like color, partial_weight) we won't
        # change, but the objects themselves may be modified. In
        # particular, we will use them as storage containers for
        # things like ROOT histograms and calculated things like scale
        # factors and event counts.
        self.samples = [s for s in MCSamples.samples if options.exclude_samples is None or s.name not in options.exclude_samples]
        self.hdata = None

        # Defaults for the dileptons, cutsets, and mass ranges for the
        # ASCII table.
        self.quantities_to_compare = ['DileptonMass', 'DimuonMassVertexConstrained']
        self.dileptons = ['MuonsPlusMuonsMinus', 'MuonsSameSign', 'MuonsAllSigns', 'MuonsElectronsOppSign', 'MuonsElectronsSameSign', 'MuonsElectronsAllSigns']
        self.cutsets = ['VBTF', 'OurNew', 'OurOld', 'Simple', 'EmuVeto', 'OurNoIso', 'OurMuPrescaled', 'VBTFMuPrescaled']
        self.mass_ranges_for_table = [(60,120), (120,200), (200,400), (400,600), (120,), (200,), (400,), (600,)]

        if options.include_quantities is not None:
            self.quantities_to_compare = options.include_quantities
        if options.include_dileptons is not None:
            self.dileptons = options.include_dileptons
        if options.include_cutsets is not None:
            self.cutsets = options.include_cutsets
        if options.include_mass_ranges_for_table is not None:
            self.mass_ranges_for_table = options.include_mass_ranges_for_table
        
        if options.for_rescale_factors:
            self.save_plots = False
            self.rescale_lumi = False
            self.mass_ranges_for_table = [(60,120)]
            self.dileptons = ['MuonsPlusMuonsMinus']

        self.setup_root()
        base = 'plots/datamc'
        if options.plot_dir_tag is not None:
            base += '_' + options.plot_dir_tag
        self.plot_dir_base = os.path.join(base, os.path.basename(self.histo_dir))
        os.system('mkdir -p %s' % self.plot_dir_base)
        width,height = options.plot_size.split(',')
        self.ps = plot_saver(self.plot_dir_base, size=(int(width), int(height)), pdf_log=True, pdf=True)

        self.table_sections = []
        self.table_rows = []
        
    def setup_root(self):
        set_zp2mu_style()
        ROOT.gStyle.SetPadLeftMargin(0.13)
        ROOT.gStyle.SetPadRightMargin(0.07)
        ROOT.TH1.AddDirectory(False)
        
    def get_join(self, sample_name):
        # If we're to merge the sample given with other histograms,
        # return a new nice-name and color for the sum. Otherwise,
        # return None.

        if self.do_joins:
            if 'qcd' in sample_name or sample_name in ('inclmu15', 'wmunu', 'wjets'):
                return 'jets', 4
            elif self.join_ttbar_and_other_prompt_leptons and sample_name in ('ttbar', 'tW', 'tbarW', 'ztautau', 'ww', 'wz', 'zz'):
                return 't#bar{t} + other prompt leptons', 2
            elif not self.join_ttbar_and_other_prompt_leptons and sample_name in ('tW', 'tbarW', 'ztautau', 'ww', 'wz', 'zz'):
                return 'other prompt leptons', 2
            elif 'dy' in sample_name or sample_name == 'zmumu':
                return '#gamma/Z #rightarrow #mu^{+}#mu^{-}', 7

        return None, None

    def is_join(self, sample_name):
        return self.get_join(sample_name)[0] is not None
    
    def get_join_name(self, sample_name):
        return self.get_join(sample_name)[0]

    def get_join_color(self, sample_name):
        return self.get_join(sample_name)[1]

    def get_color(self, sample):
        color = self.get_join_color(sample.name)
        return sample.color if color is None else color

    def get_rebin_factor(self, dilepton, quantity_to_compare):
        # For the combination of the arguments, figure out by which
        # factor to rebin the input histograms. E.g. for DileptonMass
        # the input is currently 1-GeV bins; here we change this to
        # 10-GeV bins.
        if dilepton == 'MuonsSameSign':
            if quantity_to_compare == 'DileptonMass':
                return 20
        if quantity_to_compare in ['DileptonMass', 'DimuonMassVertexConstrained', 'DileptonPt', 'LeptonPt']:
            return 10
#        if quantity_to_compare in ['RelCombIso', 'RelIsoSumPt']:
#            return 5
        if quantity_to_compare in ['DileptonPhi', 'DileptonRap', 'LeptonPhi', 'LeptonEta']:
            return 5
        return 1
        
    def rebin_histogram(self, h, cutset, dilepton, quantity_to_compare):
        # JMTBAD Make this more flexible to do arbitrary binning, e.g. by
        # mass resolution.
        ndim = h.GetDimension()
        if ndim == 1:
            h.Rebin(self.get_rebin_factor(dilepton, quantity_to_compare))
        elif ndim == 2:
            h.Rebin2D(self.get_rebin_factor(dilepton, quantity_to_compare))

    def get_x_axis_range(self, cutset, dilepton, quantity_to_compare):
        # For the given combination of the arguments, return the
        # desired restriction on the viewable x-axis range, if
        # any. E.g. for DileptonMass, only show from 60-2000 GeV on
        # the displayed plot.
        if dilepton == 'MuonsSameSign':
            if quantity_to_compare == 'DileptonMass':
                return 0,450
            if quantity_to_compare == 'DileptonPt':
                return 0,350
            if quantity_to_compare == 'DileptonRap':
                return -3,3
            if quantity_to_compare == 'LeptonPt':
                return 0,300
        if quantity_to_compare == 'RelCombIso':
            return 0,0.55
        if quantity_to_compare == 'RelIsoSumPt':
            return 0,0.1
        if 'Electron' in dilepton:
            return 100, 1000
#        if 'MuonsSameSign' in dilepton:
#            return 50, 700
        if quantity_to_compare in ['DileptonMass', 'DimuonMassVertexConstrained']:
            return 60,2000
        elif quantity_to_compare in ['DileptonPt', 'LeptonPt']:
            return 0, 700
        elif quantity_to_compare == 'LeptonEta':
            return -3,3
        return None

    def get_log_x(self, cutset, dilepton, quantity_to_compare):
        return quantity_to_compare in ['DimuonMassVtxConstrainedLog']
            
    def parse_lumi_from_log(self, log_fn):
        # JMTBAD magic, fragile parsing
        lumi_scale = 1 # assume /pb
        this = False
        for line in open(log_fn):
            if this:
                x = float(line.split()[-2])
                return x*lumi_scale
            if line == '---------------------------------------------------------------\n':
                this = True
            # lumi returned is expected to be in /pb; try to determine units from log file
            if 'Recorded' in line and 'Run' not in line:
                if '/fb' in line:
                    lumi_scale = 1000
                elif '/pb' not in line:
                    raise ValueError('cannot determine units from lumi log file: neither /fb nor /pb strings found')

    def get_lumi_rescale_factor(self, cutset, dilepton):
        # Get the cut set dependent factor by which we rescale the
        # luminosity.

        # If we're instructed not to rescale at all, don't.
        if not self.rescale_lumi:
            return 1.

        # Don't rescale e-mu plots. 
        if 'Electron' in dilepton:
            return 1.

        # These numbers are specified manually here, transcribing from
        # the 60-120 GeV tables (whether from the prescaled path) that
        # are produced when running this script with rescaling turned
        # off, rather than trying to be smart and getting them from
        # the histogram files.
        # Factors below were calculated for the MuonPhys JSON for 2012
        # released on November 2 and correspond to 15.640/fb.
        if cutset == 'VBTF':
            return 350613./337938.2
        elif cutset == 'OurNew':
            return 367659./357970.1
        elif cutset == 'OurOld':
            return 416016./396433.0
        elif cutset == 'OurNoIso':
            return 373475./362446.1
        elif cutset == 'OurMuPrescaled':
            return 16644./16460.2
        elif cutset == 'VBTFMuPrescaled':
            return 15522./15001.
        # If the cutset is not one of the above, don't rescale.
        return 1.

    def advertise_lines(self):
        s = []
        s.append('total lumi from data: %.f/pb' % self.int_lumi)
        s.append('rescaling mumu MC histograms by these cut-dependent factors:')
        for cutset in self.cutsets:
            s.append('%15s:%.10f' % (cutset, self.get_lumi_rescale_factor(cutset, '')))
        s.append('"joins" are:')
        for sample in self.samples:
            s.append('%10s -> %s' % (sample.name, self.get_join_name(sample.name)))
        return s

    def subtitleize(self, dilepton):
        return {
            'MuonsPlusMuonsMinus': '#mu^{+}#mu^{-}',
            'MuonsPlusMuonsPlus':  '#mu^{+}#mu^{+}',
            'MuonsMinusMuonsMinus': '#mu^{-}#mu^{-}',
            'MuonsSameSign': '#mu^{#pm}#mu^{#pm}',
            'MuonsAllSigns': '#mu#mu',
            'ElectronsPlusElectronsMinus': 'e^{+}e^{-}',
            'ElectronsPlusElectronsPlus': 'e^{+}e^{+}',
            'ElectronsMinusElectronsMinus': 'e^{-}e^{-}',
            'ElectronsSameSign': 'e^{#pm}e^{#pm}',
            'ElectronsAllSigns': 'ee',
            'MuonsPlusElectronsMinus': '#mu^{+}e^{-}',
            'MuonsMinusElectronsPlus': '#mu^{-}e^{+}',
            'MuonsPlusElectronsPlus': '#mu^{+}e^{+}',
            'MuonsMinusElectronsMinus': '#mu^{-}e^{-}',
            'MuonsElectronsOppSign': '#mu^{+}e^{-}/#mu^{-}e^{+}',
            'MuonsElectronsSameSign': '#mu^{#pm}e^{#pm}',
            'MuonsElectronsAllSigns': 'e#mu',
            }[dilepton]

    def titleize(self, quantity_to_compare):
        return {
            'DileptonMass': 'm(%s)%s',
            'DimuonMassVertexConstrained': 'm(%s)%s',
            'DimuonMassVtxConstrainedLog': 'm(%s)%s',
            'DileptonPt': '%s p_{T}%s',
            'DileptonRap': '%s rapidity%s',
            'LeptonPt': "%s leptons' p_{T}%s",
            'LeptonEta': "%s leptons' #eta%s",
            'RelIsoSumPt': "%s leptons' relative tk. iso.%s",
            'RelCombIso': "%s leptons' relative comb. iso.%s",
            }.get(quantity_to_compare, quantity_to_compare + ', %s, %s')

    def unitize(self, quantity_to_compare):
        return {
            'DileptonMass': ' [GeV]',
            'DimuonMassVertexConstrained': ' [GeV]',
            'DimuonMassVtxConstrainedLog': ' [GeV]',
            'DileptonPt': ' [GeV]',
            'LeptonPt': ' [GeV]',
            'LeptonEta': '',
            'LeptonPhi': ' [rad]',
            'DileptonPhi': ' [rad]',
            'DileptonRap': '',
            'RelCombIso': '',
            'RelIsoSumPt': '',
            }.get(quantity_to_compare, ' [XXX]')

    def get_dir_name(self, cutset, dilepton):
        return cutset + dilepton + 'Histos'

    def get_y_axis_range(self, dilepton, cumulative):
        return None

    def handle_overflows(self, h, range):
        if h.GetDimension() != 1:
            return
        if not self.put_overflow_in_last_bin:
            return
        if range is None:
            range = h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX()+1)
        move_above_into_bin(h, range[1])

    def prepare_mc_histograms(self, cutset, dilepton, quantity_to_compare, cumulative):
        # For each sample in the list of MC samples, we get the
        # specified histogram from the appropriate input file, clone
        # it, rebin/scale/otherwise manipulate it as necessary, and
        # store the result as The Histogram in the sample object.
        for sample in self.samples:
            mc_fn = os.path.join(self.mc_dir, 'ana_datamc_%s.root' % sample.name)
            f = ROOT.TFile(mc_fn)
            sample.histogram = getattr(f, self.get_dir_name(cutset, dilepton)).Get(quantity_to_compare).Clone()

            self.rebin_histogram(sample.histogram, cutset, dilepton, quantity_to_compare)

            sample.scaled_by = self.int_lumi * self.get_lumi_rescale_factor(cutset, dilepton) * sample.partial_weight
            sample.histogram.Scale(sample.scaled_by)

            if cumulative:
                sample.histogram = cumulative_histogram(sample.histogram)
            else:
                self.handle_overflows(sample.histogram, self.get_x_axis_range(cutset, dilepton, quantity_to_compare))

            sample.histogram.SetMarkerStyle(0)
            color = self.get_color(sample)
            sample.histogram.SetLineColor(color)
            if not sample.is_zprime:
                sample.histogram.SetFillColor(color)

    def prepare_data_histogram(self, cutset, dilepton, quantity_to_compare, cumulative):
        self.hdata = getattr(self.data_f, self.get_dir_name(cutset, dilepton)).Get(quantity_to_compare).Clone()
        self.rebin_histogram(self.hdata, cutset, dilepton, quantity_to_compare)
        range = self.get_x_axis_range(cutset, dilepton, quantity_to_compare)

        if cumulative:
            self.hdata = cumulative_histogram(self.hdata)
        else:
            self.handle_overflows(self.hdata, range)

        # if not self.put_overflow_in_last_bin:
            # Not so important if the MC histograms have entries past
            # the view range that isn't shown, but not so for
            # data. Check that there's no data off screen.
            # overflow_integral = get_integral(self.hdata, range[1], integral_only=True)
            # if overflow_integral > 0:
            #    raise ValueError('WARNING: in %s, data histogram has points in overflow (mass bins above %.f GeV)! integral = %f' % (cutset + dilepton, range[1], overflow_integral))

    def make_table(self, cutset, dilepton):
        # Make a nicely formatted ASCII table of event counts in the
        # specified mass ranges, along with uncertainties.
        for mass_range in self.mass_ranges_for_table:
            self.table_sections.append((cutset, dilepton, mass_range))
            
            self.table_rows.append('*'*(50+20*9) + '\n\n')
            self.table_rows.append('ANCHORMEcuts: %s  dilepton: %s  mass range: %s\n' % (cutset, dilepton, mass_range))

            # For all the MC samples, calculate the integrals over the
            # current mass range and store it in the sample object.
            for sample in self.samples:
                sample.integral = get_integral(sample.histogram, *mass_range, integral_only=True, include_last_bin=False)
                sample.raw_integral = sample.integral / sample.scaled_by

            # Header. (I hope you have a widescreen monitor, a small font, and good eyes.)
            self.table_rows.append('%50s%20s%20s%20s%20s%20s%20s%20s%20s%20s\n' % ('sample', 'weight', 'raw integral', 'integral', 'stat error', 'limit if int=0', 'syst error', 'syst(+)stat', 'lumi error', 'total'))

            # Print the row for the event count from the data (only
            # the integral and statistical uncertainty columns will be
            # filled).
            data_integral = get_integral(self.hdata, *mass_range, integral_only=True, include_last_bin=False)
            self.table_rows.append('%50s%20s%20i%20.6f%20.6f\n' % ('data', '-', int(data_integral), data_integral, data_integral**0.5))

            # As we loop over the MC samples, keep some running sums
            # of integrals and variances. Do one such set including
            # the whole MC, and a set for each join group.
            sum_mc = 0.
            var_sum_mc = 0.
            syst_var_sum_mc = 0.
            sums = defaultdict(float) # will initialize to 0. on first lookup
            var_sums = defaultdict(float)
            syst_var_sums = defaultdict(float)

            # Loop over the MC samples in order of decreasing
            # integral.
            for sample in sorted(self.samples, key=lambda x: x.integral, reverse=True):
                # Calculate the contributions to the variances for the
                # statistical and systematic uncertainties for this
                # sample.
                w = sample.scaled_by
                var = w * sample.integral # not w**2 * sample.integral because sample.integral is already I*w
                syst_var = (sample.syst_frac * sample.integral)**2

                # Add the integral and the variances to the whole-MC
                # values, if it is to be included (i.e. don't include
                # Z' samples).
                if not sample.is_zprime:
                    sum_mc += sample.integral
                    var_sum_mc += var
                    syst_var_sum_mc += syst_var

                # If this sample belongs to a join group, add the
                # integral and the variances to the values for the
                # join group.
                join_name = self.get_join_name(sample.name)
                if join_name is not None:
                    sums[join_name] += sample.integral
                    var_sums[join_name] += var
                    syst_var_sums[join_name] += syst_var

                # For integrals that turn out to be zero due to not
                # enough statistics, we will give also the 95% CL
                # upper limit.
                limit = '%.6f' % (3*w) if sample.integral == 0 else '-'

                # For this sample alone, determine the combined
                # statistical+systematic uncertainty, the uncertainty
                # due to luminosity, and the total of all three.
                syst_plus_stat = (var + syst_var)**0.5
                lumi_err = self.lumi_syst_frac * sample.integral
                tot_err = (var + syst_var + lumi_err**2)**0.5

                # Print this row of the table.
                self.table_rows.append('%50s%20.6f%20f%20.6f%20.6f%20s%20.6f%20.6f%20.6f%20.6f\n' % (sample.nice_name, w, sample.raw_integral, sample.integral, var**0.5, limit, syst_var**0.5, syst_plus_stat, lumi_err, tot_err))

            self.table_rows.append('-'*(50+20*9) + '\n')
            
            # Determine the uncertainties and print the rows for each
            # of the join groups. Sort this section by decreasing integral.
            join_integrals = sums.items()
            join_integrals.sort(key=lambda x: x[1])
            for join_name, join_integral in join_integrals:
                lumi_err = self.lumi_syst_frac * sums[join_name]
                syst_plus_stat = (var_sums[join_name] + syst_var_sums[join_name])**0.5
                tot_err = (var_sums[join_name] + syst_var_sums[join_name] + lumi_err**2)**0.5
                self.table_rows.append('%50s%20s%20s%20.6f%20.6f%20s%20.6f%20.6f%20.6f%20.6f\n' % (join_name, '-', '-', sums[join_name], var_sums[join_name]**0.5, '-', syst_var_sums[join_name]**0.5, syst_plus_stat, lumi_err, tot_err))

            self.table_rows.append('-'*(50+20*9) + '\n')

            # For the sum of all MC, determine the combined
            # statistical+systematic uncertainty, the uncertainty due
            # to luminosity, and the total of all three. Then print
            # the whole-MC row.
            syst_plus_stat = (var_sum_mc + syst_var_sum_mc)**0.5
            lumi_err = self.lumi_syst_frac * sum_mc
            tot_err = (var_sum_mc + syst_var_sum_mc + lumi_err**2)**0.5
            self.table_rows.append('%50s%20s%20s%20.6f%20.6f%20s%20.6f%20.6f%20.6f%20.6f\n' % ('sum MC (not including Z\')', '-', '-', sum_mc, var_sum_mc**0.5, '-', syst_var_sum_mc**0.5, syst_plus_stat, lumi_err, tot_err))

            self.table_rows.append('\n')
        self.table_rows.append('\n')

    def should_draw_zprime(self, dilepton):
        return self.draw_zprime and dilepton == 'MuonsPlusMuonsMinus'

    def get_zprime_histogram(self):
        # JMTBAD Extend to the rest of the Z' samples when there are any.
        try:
            from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import zssm1000
            return zssm1000.histogram        # Loaded/scaled already in prepare_mc_histograms since the loop over the samples list modified the original object that zssm1000 points to.
        except ImportError:
            pass

    def draw_legend(self, dilepton, cumulative, log_x):
        # Legend placement coordinates and sizes depend on factors set
        # elsewhere, too, so this is fragile.
        if dilepton == 'MuonsPlusMuonsMinus' and cumulative:
            legend = ROOT.TLegend(0.60, 0.69, 0.86, 0.88)
        elif log_x:
            legend = ROOT.TLegend(0.60, 0.55, 0.86, 0.88)
        else:
            legend = ROOT.TLegend(0.60, 0.69, 0.86, 0.88)

        legend.SetFillColor(0)
        legend.SetBorderSize(0)

        # Add an entry for the data points.
        entry = legend.AddEntry('data_marker', 'Data', 'EP')
        entry.SetMarkerStyle(20)
        entry.SetMarkerSize(0.8)
        entry.SetMarkerColor(ROOT.kBlack)

        # Add entries for the MC samples to the legend, respecting
        # join groups (i.e. don't add the same nice-name twice).
        legend_already = set()
        for sample in reversed(self.samples):
            if sample.is_zprime and not self.should_draw_zprime(dilepton):
                continue

            nice_name = sample.nice_name

            join_name = self.get_join_name(sample.name)
            if join_name is not None:
                if join_name in legend_already:
                    continue
                else:
                    legend_already.add(join_name)
                    nice_name = join_name
            legend.AddEntry(sample.histogram, nice_name, 'F') # Gets the color and fill style from the histogram.

        legend.SetTextSize(0.03)
        legend.Draw('same')
        ## "EP" in TLegend::AddEntry doesn't seem to work, so draw the error bar by hand
        ll = ROOT.TLine()
        if log_x:
            ll.DrawLineNDC(0.632, 0.845, 0.632, 0.875)
        else:
            ll.DrawLineNDC(0.632, 0.835, 0.632, 0.875)

        return legend

    def draw_data_on_mc(self, cutset, dilepton, quantity_to_compare, cumulative):
        # Make a Stack for the MC histograms. We draw it first so the
        # data points will be drawn on top later.  We assume that in
        # the list of samples, the join groups are contiguous already
        # so that like-colored histograms will be next to each other.
        
        s = ROOT.THStack('s', '')
        for sample in self.samples:
            # Don't stack the Z' samples.
            if sample.is_zprime:
                continue
            s.Add(sample.histogram) # Assumes they've already been prepared.

        # Figure out its titles.
        xtitle = self.titleize(quantity_to_compare) % (self.subtitleize(dilepton), self.unitize(quantity_to_compare))
        if cumulative:
            # E.g. Events >= m(mu+mu-).
            ytitle = 'Events #geq %s' % (self.titleize(quantity_to_compare) % (self.subtitleize(dilepton), ''))
        else:
            # E.g. Events/5 GeV.
            ytitle = 'Events / %i%s' % (self.get_rebin_factor(dilepton, quantity_to_compare), self.unitize(quantity_to_compare).translate(None, '()[]')) # Events/5 GeV. JMTBAD assumes original histogram had 1 GeV bins and was rebinned simply -- ignores rebin_histogram ability to have arb. bins
            if quantity_to_compare == 'DileptonRap':
                ytitle = 'Events / 0.5'
            elif quantity_to_compare == 'RelIsoSumPt':
                ytitle = 'Events / 0.02'
            elif quantity_to_compare == 'RelCombIso':
                ytitle = 'Events / 0.05'
            elif quantity_to_compare == 'DimuonMassVtxConstrainedLog':
                ytitle = 'Events / bin'
        s.SetTitle(';%s;%s' % (xtitle, ytitle))

        s.Draw('hist')

        # Must call Draw first or the THStack doesn't have a histogram/axes.
        s.GetXaxis().SetTitleOffset(0.9)
        s.GetXaxis().SetTitleSize(0.047)
        s.GetYaxis().SetTitleOffset(1.2)
        s.GetYaxis().SetTitleSize(0.047)

        # Set the x-axis range as specified. Then determine what
        # the real extrema of the data histogram and the MC stack are
        # over the xrange specified.
        xrange = self.get_x_axis_range(cutset, dilepton, quantity_to_compare)
        mymin, mymax = None, None
        if xrange is not None:
            s.GetXaxis().SetRangeUser(*xrange)
            self.hdata.GetXaxis().SetRangeUser(*xrange)
            mymin = real_hist_min(s.GetStack().Last(), user_range=xrange) * 0.7
            mymax = real_hist_max(s.GetStack().Last(), user_range=xrange, use_error_bars=False) * 1.05
            if self.hdata.GetEntries() > 0:
                rhm = real_hist_max(self.hdata, user_range=xrange)
                mymax = max(mymax, rhm)

        if self.guess_yrange:
            mymin = 0
            mymax = real_hist_max(self.hdata, user_range=xrange)

        # Can override the above fussing.
        yrange = self.get_y_axis_range(dilepton, cumulative)
        if yrange is not None:
            if yrange[0] is not None:
                mymin = yrange[0]
            if yrange[1] is not None:
                mymax = yrange[1]

        mymin = 0.05
    
        if mymin is not None: s.SetMinimum(mymin)
        if mymax is not None: s.SetMaximum(mymax)

        # Calculate (data-bckg)/bckg.  Do it before TH1 gets converted
        # to TGraphAsymmErrors by poisson_intervalize.
        if not cumulative:
            ifois = 0
            for sample in self.samples:
                # Don't add the Z' samples.
                if sample.is_zprime:
                    continue
                if ifois == 0:
                    mc_sum = sample.histogram.Clone()
                    ifois = 1
                else:
                    mc_sum.Add(sample.histogram, 1.)

            data_mc_diff = self.hdata.Clone()
            data_mc_diff.Divide(mc_sum)
            nbins = data_mc_diff.GetNbinsX()
            for ibin in range(1, nbins):
                f_bin = data_mc_diff.GetBinContent(ibin)
                data_mc_diff.SetBinContent(ibin, f_bin-1.)

        # Now draw the data on top of the stack.
        self.hdata.SetStats(0)
        data_draw_cmd = 'same p e'
        if self.use_poisson_intervals:
            self.hdata = poisson_intervalize(self.hdata, True)
            data_draw_cmd += ' z'
        if mymin is not None: self.hdata.SetMinimum(mymin)
        if mymax is not None: self.hdata.SetMaximum(mymax)
        self.hdata.SetMarkerStyle(20)
        self.hdata.SetMarkerSize(0.8)
        self.hdata.Draw(data_draw_cmd)

        # Draw the Z' curve separately. It is overlaid, not stacked
        # with the rest of the MC expectation.
        zp = self.get_zprime_histogram()
        if self.should_draw_zprime(dilepton) and zp:
            # JMTBAD Extend to the rest of the Z' samples when there are any.
            zp.SetTitle('')
            zp.SetLineWidth(2)
            if xrange is not None:
                zp.GetXaxis().SetRangeUser(*xrange)
                zp.SetMinimum(mymin)
                zp.SetMaximum(mymax)
            zp.SetStats(0)
            zp.Draw('hist same')

        # Use log(x) whenever needed
        log_x = self.get_log_x(cutset, dilepton, quantity_to_compare)
        if log_x:
            self.ps.c.SetLogx()

        # Adorn the plot with legend and labels.
        l = self.draw_legend(dilepton, cumulative, log_x)
        t = ROOT.TPaveLabel(0.20, 0.89, 0.86, 0.99, 'CMS 2012 preliminary   #sqrt{s} = 8 TeV    #int L dt = %.f pb^{-1}' % round(self.int_lumi), 'brNDC')
        t.SetTextSize(0.35)
        t.SetBorderSize(0)
        t.SetFillColor(0)
        t.SetFillStyle(0)
        t.Draw()

        # Done; save it!
        plot_fn = dilepton
        if cumulative:
            plot_fn += '_cumulative'
        self.ps.save(plot_fn)

        if log_x:
            self.ps.c.SetLogx(0)

        if not cumulative:
            data_mc_diff.SetMinimum(-1.)
            data_mc_diff.SetMaximum(1.)
            data_mc_diff.SetMarkerStyle(20)
            data_mc_diff.SetMarkerSize(0.8)
            data_mc_diff.SetTitle("(data-bckg)/bckg")
            data_mc_diff.SetStats(0)
            data_mc_diff.Draw("p e")
            if xrange is not None:
                l1 = ROOT.TLine(xrange[0], 0., xrange[1],  0.)
            else:
                l1 = ROOT.TLine(data_mc_diff.GetXaxis().GetXmin(), 0., data_mc_diff.GetXaxis().GetXmax(), 0.)
            l1.Draw()
            plot_fn += '_diff'
            self.ps.save(plot_fn, log=False, pdf_log=False)
            
    def finalize_table(self, dir_base):
        table_fn = os.path.join(dir_base, 'mass_counts.html')
        table_f = open(table_fn, 'wt')
        table_f.write('<html><body><pre>\n')
        table_f.write('\n'.join(self.advertise_lines()) + '\n')
        last_cutset = None
        last_dilepton = None
        anchors = []
        for cutset, dilepton, mass_range in self.table_sections:
            anchor = cutset+dilepton+str(mass_range[0])
            if len(mass_range) > 1:
                anchor += str(mass_range[1])
            anchors.append(anchor)
                
            cutset = '%12s' % cutset
            dilepton = '%25s' % dilepton
            mass_range = '%15s' % repr(mass_range)

            text = ''
            if cutset != last_cutset:
                text += cutset
                last_cutset = cutset
            else:
                text += ' '*12
            if dilepton != last_dilepton:
                text += dilepton
                last_dilepton = dilepton
            else:
                text += ' '*20
            text += mass_range
            
            table_f.write('<a href="#%s">%s</a>\n'% (anchor, text))

        for row in self.table_rows:
            if 'ANCHORME' in row:
                row = '<h4 id="%s">%s</h4>' % (anchors.pop(0), row.replace('ANCHORME', ''))
            table_f.write(row)

        table_f.write('</pre></body></html>\n')
        table_f.close()
        self.table_sections = []
        self.table_rows = []
        
    def go(self):
        for quantity_to_compare in self.quantities_to_compare:
            print quantity_to_compare
            
            if quantity_to_compare != 'DileptonMass':
                cutsets = ['OurNew']
            else:
                cutsets = self.cutsets
            
            for cutset in cutsets:
                print cutset
                
                # If the cut set doesn't exist in the input file, silently skip it.
                if not hasattr(self.data_f, self.get_dir_name(cutset, 'MuonsPlusMuonsMinus')):
                    continue

                # Directory structure example:
                # plots/datamc/lumi_mask_name/quantity_to_compare/cut_set/.
                plot_dir = self.plot_dir_base + '/%s/%s' % (quantity_to_compare, cutset)
                self.ps.set_plot_dir(plot_dir)

                # Depending on the quantity to compare and cut set, skip certain dileptons.
                if cutset == 'EmuVeto': # Only care about e-mu dileptons here.
                    dileptons = [x for x in self.dileptons if 'Electron' in x]
                elif 'MuPrescaled' in cutset: # Don't care about e-mu dileptons here.
                    dileptons = [x for x in self.dileptons if 'Electron' not in x]
                else:
                    dileptons = [x for x in self.dileptons if 'MuonsAllSigns' not in x]
#                    dileptons = self.dileptons

                # Also depending on the quantity to be compared, skip certain
                # dileptons.
                if 'Dimuon' in quantity_to_compare: # Only defined for mumu.
                    dileptons = [x for x in dileptons if 'Electron' not in x]

                self.ps.save_dir('mass_counts.html')

                for dilepton in dileptons:
                    print dilepton
                    
                    for cumulative in (False, True):
                        # Prepare the histograms. The MC histograms are stored in
                        # their respective sample objects, and the data histogram is
                        # kept in self.hdata.
                        self.prepare_mc_histograms(cutset, dilepton, quantity_to_compare, cumulative)
                        self.prepare_data_histogram(cutset, dilepton, quantity_to_compare, cumulative)

                        # Print the entries for the ASCII table for the current
                        # cutset+dilepton. Could extend this to support counts for
                        # ranges that aren't mass.
                        if not cumulative and quantity_to_compare in ['DileptonMass', 'DimuonMassVertexConstrained'] and self.print_table:
                            self.make_table(cutset, dilepton)

                        if self.save_plots:
                            self.draw_data_on_mc(cutset, dilepton, quantity_to_compare, cumulative)

                self.finalize_table(plot_dir)


d = Drawer(options)
print '\n'.join(d.advertise_lines())
d.go()
