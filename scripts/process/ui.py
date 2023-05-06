import collections
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg as tkagg
import nmrglue as ng
import numpy as np
import tkinter as tk
from tkinter import ttk


SCDIR = os.path.dirname(__file__)   # location of script

def euclid_dist(a, b):
    """Compute euclidean distance between two points."""
    return np.linalg.norm(a - b)

def empty_function(arg):
    pass

class Isotope:
    def __init__(self, label, ppm_range):
        self.label = label
        self.ppm_min = float(min(ppm_range))
        self.ppm_max = float(max(ppm_range))

class Trace:
    """Real or simulated 1D NMR spectrum, or list of peaks."""
    def __init__(self, line):
        self.line = line
        self.line.set(picker=True) # allow interaction

    def y(self, x):
        """Get amplitude of nearest point."""
        i = np.argmin([abs(a - x) for a in self.line.get_xdata()])
        return self.line.get_ydata()[i]

    def get_line(self):
        return self.line

    def get(self, *args):
        """Get properties from line."""
        return {arg: plt.getp(self.line, arg) for arg in args}

    def set(self, **kwargs):
        """Set line properties."""
        self.line.set(**kwargs)

    def handle_click(self, event):
        """Handle click and return nearest point."""
        event_on_line, indices = self.line.contains(event)
        indices = indices['ind']
        if event_on_line and indices:
            event_xy = np.array([event.x, event.y])
            indices = sorted(indices, key=lambda xy: euclid_dist(xy, event_xy))
            return indices[0]
        else:
            return None

class Peak:
    """Picked NMR peak"""
    def __init__(self, isotope, shift, height, l_fwhm=None,
                 g_fwhm=None, cpd=None):
        self.isotope = isotope
        self.x = shift    # Voigt curve ppm shift
        self.y = height   # Voigt curve height
        self.lw = l_fwhm  # Voigt curve Lorentzian full-width half maximum
        if g_fwhm is not None:
            self.cv = 'v'
            self.gw = g_fwhm  # Voigt curve Gaussian full-width half maximum
        else:
            self.cv = 'l'
            self.gw = 0 
        self.cpd = cpd
        self.trace = self.render_trace(resolution=100000)

    def get_x(self): return self.x

    def get_y(self): return self.y

    def get_lw(self): return self.lw

    def set_x(self, x): self.update(x=float(x))

    def set_y(self, y): self.update(y=float(y))

    def set_lw(self, lw): self.update(lw=float(lw))

    def set_gw(self, gw): self.update(gw=float(gw))

    def set_cpd(self, cpd): self.cpd = str(cpd)

    def get_curve_param_dict(self):
        pars = {
            'shift': self.x,
            'height': self.y,
            'l_fwhm': self.lw,
            'g_fwhm': self.gw,
        }

    def ng_params(self):
        if self.cv == 'v':
            pars = [self.x, self.lw, self.gw]
        else:
            pars = [self.x, self.lw]
        return pars

    def ng_params_i(self, res):
        x0 = self.isotope.ppm_min # may need to switch these two
        x1 = self.isotope.ppm_max
        def tx(x): return res*(x - x0)/(x1 - x0)
        return [[[tx(x) for x in self.ng_params()]]]

    def ng_amp(self):
        return [[self.y]]

    def render_trace(self, resolution, update=True):
        params = self.ng_params_i(resolution)
        amps = self.ng_amp()
        trace = ng.analysis.linesh.sim_NDregion(
            (resolution,), [self.cv], params, amps,
        )
        if update:
            self.trace = trace
        return trace

    def align_peak_to_curve(self):
        if isinstance(self.x, float) and isinstance(self.y, float):
            self.set_x(self.x)
            self.set_y(self.y)
        else:
            typestring = f"{type(self.x)} and {type(self.y)}"
            raise ValueError("invalid types for coordinates: " + typestring)

    def update(self, **kwargs):
        for k, v in kwargs.items():
            self.__setattr__(k, v)
        self.render_trace(len(self.trace))

    def to_record(self):
        record = {
            'cpd': self.cpd,    # compound
            'cv': self.cv,      # curveshape
            'ppm': self.x,      # ppm shift
            'amp': self.y,      # amplitude
            'fwhm': self.lw,    # Lorentzian FWHM
        }
        if self.gw is not None:
            record['fwhm2'] = self.gw    # Gaussian FWHM
        return record

class FloatText:
    def __init__(self, master):
        self.frame = tk.Frame(master)
        vcmd = (master.register(self.validate))
        self.entry = tk.Entry(self.frame, validate='all', validatecommand=(vcmd, '%P', '%V'))
        self.entry.pack()

    def validate(self, value, reason):
        if reason == 'key':
            if all(c in '0123456789.+-' for c in value) or value == '':
                return True
            else:
                return False
        elif reason == 'focusout':
            try:
                float(value)
                return True
            except:
                # print("Invalid float!")
                return False
        return True

class Window():
    def __init__(self, isotope, peaks, spec): # spec, pack_peaks_fn, unpack_peaks_fn, curve_fit_fn):
        # Set up plot
        self.isotope = isotope
        self.peaks = peaks
        self.res = 100000
        self.assignments = collections.defaultdict(list)
        self.fid = spec
        for pk in peaks:
            self.assignments[pk.cpd].append(pk)
        if len(peaks) > 0:
            self.active_peak = peaks[0]
        else:
            self.active_peak = None
        self.isotope = isotope
        
        # self.fig, self.plot_ax = plt.subplots(constrained_layout=True, figsize=(10, 10))
        self.fig = mpl.figure.Figure()
        self.plot_ax = self.fig.add_subplot(111)
        self.plot_ax.set_xlim(left=isotope.ppm_max, right=isotope.ppm_min)
        
        simspace = np.linspace(isotope.ppm_min, isotope.ppm_max, num=self.res)
        simspec = sum([peak.render_trace(self.res) for peak in peaks])
        l1 = self.plot_ax.plot(  # plot NMR spectrum
            spec.ppm, 
            spec.data, 
            c='black', 
            lw=1,
            zorder=1,
        )
        l2 = self.plot_ax.plot(  # plot simulated spectrum
            simspace, 
            simspec, 
            c='royalblue', 
            lw=1,
            zorder=2,
        )
        self.cpd_simspecs = {}
        for cpd, peaklist in self.assignments.items():
            l = self.plot_ax.plot( # plot simulated compound spectrum
                simspace, 
                sum([peak.trace for peak in peaklist]), 
                lw=1, 
                zorder=3
            )
            self.cpd_simspecs[cpd] = Trace(l[0])

        l3 = self.plot_ax.plot(   # plot active peak
            np.linspace(isotope.ppm_min, isotope.ppm_max, num=self.res), 
            self.active_peak.render_trace(self.res),
            c='deeppink',
            lw=1.5,
            zorder=4,
        )
        mask = (spec.ppm < self.isotope.ppm_max) & (spec.ppm > self.isotope.ppm_min)
        l4 = self.plot_ax.plot(  # plot residuals
            spec.ppm[mask], 
            spec.data[mask] - np.flip(sum([peak.render_trace(np.sum(mask), update=False) for peak in peaks])), 
            c='gray', 
            lw=1,
            zorder=0,
        )
        self.spec = Trace(l1[0]) # should not be edited
        self.simspec = Trace(l2[0])
        self.overlay = Trace(l3[0])
        self.residuals = Trace(l4[0])

        # Set up window
        self.root = tk.Tk()
        self.root.geometry("1280x800")
        self.root.title("Manual Peak Fitting Subroutine")
        self.root.config(bg="skyblue")
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=3)
        self.root.rowconfigure(0, weight=1)

        # Create frame for plot
        self.plot_frame = tk.Frame(self.root)
        self.plot_frame.grid(row=0, column=1, padx=5, pady=5, sticky=tk.E+tk.W+tk.N+tk.S)
        self.canvas = tkagg.FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar = tkagg.NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Create frames for contents
        self.options_frame = tk.Frame(self.root)
        self.options_frame.grid(row=0, column=0, padx=5, pady=5, sticky=tk.E+tk.W+tk.N+tk.S)
        self.options_frame.columnconfigure(0, weight=1)
        self.options_frame.rowconfigure(0, weight=1)
        self.options_frame.rowconfigure(1, weight=3)

        # Initialize peak selector
        self.selector_frame = tk.Frame(self.options_frame)
        self.selector_frame.grid(row=0, column=0, padx=5, pady=5)

        self.sel_label = tk.Label(self.selector_frame, text="Active Peak")
        self.sel_label.grid(row=0, column=0)
        self.selector = ttk.Combobox(
            self.selector_frame,
        )
        self.update_selector()
        self.selector.bind("<<ComboboxSelected>>", self.select_callback)
        # selector.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.selector.grid(row=0, column=1)

        self.fit_button = tk.Button(self.selector_frame, text="Curve Fit", command=self.refit_peaks)
        self.fit_button.grid(row=2, column=0)
        self.fit_progress = tk.Label(self.selector_frame, text="Fitting peaks...")
        self.fit_progress.grid(row=2, column=1)
        self.fit_progress.grid_forget()

        self.selector_frame.columnconfigure(0, weight=1)
        self.selector_frame.columnconfigure(1, weight=1)
        self.selector_frame.rowconfigure(0, weight=1)
        self.selector_frame.rowconfigure(1, weight=1)
        self.selector_frame.rowconfigure(2, weight=1)

        # sel_scroll = tk.Scrollbar(selector_frame, orient=tk.VERTICAL)
        # sel_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        # sel_scroll['command'] = selector.yview

        # Initialize dialogues
        self.dialogue_frame = tk.Frame(self.options_frame)
        self.dialogue_frame.grid(row=1, column=0, padx=5, pady=5, sticky=tk.E+tk.W+tk.N+tk.S)
        self.dialogue_nb = ttk.Notebook(self.dialogue_frame)
        self.edit_frame = tk.Frame(self.dialogue_nb)
        self.add_frame = tk.Frame(self.dialogue_nb)
        self.rem_frame = tk.Frame(self.dialogue_nb)
        self.dialogue_nb.add(self.edit_frame, text="Edit Peak")
        self.dialogue_nb.add(self.add_frame, text="Add Peak")
        self.dialogue_nb.add(self.rem_frame, text="Remove Peak")
        # dialogue_nb.pack()

        # Edit dialogue
        self.edit_title = tk.Label(self.edit_frame, text="Compound, Shift")
        self.edit_title.grid(row=0, column=0, columnspan=2)
        self.edit_var = [tk.StringVar() for i in range(4)]
        self.edit_elements = {
            'Compound': tk.Entry(self.edit_frame),
            'Shift (ppm)': FloatText(self.edit_frame),
            'Amplitude': FloatText(self.edit_frame),
            'Lorentzian FWHM (ppm)': FloatText(self.edit_frame),
        }
        for i, (k, v) in enumerate(self.edit_elements.items()):
            tk.Label(self.edit_frame, text=k).grid(row=i+1, column=0, sticky=tk.E)
            if i == 0:
                v.grid(row=1, column=1, sticky=tk.W)
                v.config(textvariable=self.edit_var[0])
            else:
                v.frame.grid(row=i+1, column=1, sticky=tk.W)
                v.entry.config(textvariable=self.edit_var[i])
        self.edit_elements['Compound'].bind("<FocusOut>", self.edit_cpd_callback)
        self.edit_elements['Shift (ppm)'].entry.bind("<FocusOut>", self.edit_shift_callback)
        self.edit_elements['Amplitude'].entry.bind("<FocusOut>", self.edit_amp_callback)
        self.edit_elements['Lorentzian FWHM (ppm)'].entry.bind("<FocusOut>", self.edit_lw_callback)
        self.select_callback(None)
        self.edit_frame.columnconfigure(0, weight=1)
        self.edit_frame.columnconfigure(1, weight=1)

        # Add dialogue
        self.add_var = [tk.StringVar(), tk.StringVar()]
        self.add_title = tk.Label(self.add_frame, text="Add Peak")
        self.add_title.grid(row=0, column=0, columnspan=2)
        tk.Label(self.add_frame, text="Compound").grid(row=1, column=0, sticky=tk.E)
        self.add_cpd = tk.Entry(self.add_frame, textvariable=self.add_var[0])
        self.add_cpd.grid(row=1, column=1, sticky=tk.W)
        tk.Label(self.add_frame, text="Shift (ppm)").grid(row=2, column=0, sticky=tk.E)
        self.add_shift = FloatText(self.add_frame)
        self.add_shift.frame.grid(row=2, column=1, sticky=tk.W)
        self.add_shift.entry.config(textvariable=self.add_var[1])
        self.add_submit = tk.Button(self.add_frame, text="Submit", command=self.add_button_callback)
        self.add_submit.grid(row=3, column=0, columnspan=2)
        self.add_frame.columnconfigure(0, weight=1)
        self.add_frame.columnconfigure(1, weight=1)

        # Remove dialogue
        self.remove_title = tk.Label(self.rem_frame, text="Remove Active Peak")
        self.remove_title.grid(row=0, column=0)
        self.remove_submit = tk.Button(self.rem_frame, text="Submit", command=self.remove_button_callback)
        self.remove_submit.grid(row=1, column=0)
        self.dialogue_nb.pack(fill=tk.BOTH, expand=True)
        self.rem_frame.columnconfigure(0, weight=1)

    def display(self):
        self.root.mainloop()

    def get_peaks(self):
        return self.peaks

    def update_cpd_simspec(self, cpd):
        if cpd not in self.cpd_simspecs:
            l = self.plot_ax.plot(
                np.linspace(self.isotope.ppm_min, self.isotope.ppm_max, num=self.res),
                sum([peak.trace for peak in self.peaks]),
                lw=1, 
                zorder=2
            )
            self.cpd_simspecs[cpd] = Trace(l[0])
        elif cpd in self.assignments and len(self.assignments.get(cpd, [])) > 0:
            self.cpd_simspecs[cpd].set(ydata=sum([pk.trace for pk in self.assignments[cpd]]))
        else:
            trace = self.cpd_simspecs.pop(cpd)
            trace.line.remove()

    def update_plot(self):
        if self.active_peak is not None:
            self.overlay.set(ydata=self.active_peak.trace)
            self.simspec.set(ydata=sum([peak.render_trace(self.res) for peak in self.peaks]))
            ppm = self.spec.get_line().get_xdata()
            mask = (ppm < self.isotope.ppm_max) & (ppm > self.isotope.ppm_min)
            trace = self.spec.get_line().get_ydata()[mask]
            sim2 = sum([peak.render_trace(len(trace), update=False) for peak in self.peaks])
            self.residuals.set(ydata=(trace - np.flip(sim2)))
            self.update_cpd_simspec(self.active_peak.cpd)
            self.canvas.draw()
            self.canvas.flush_events()

    def update_selector(self, initial_selected=None):
        """Update selector fields and optionally select initial peak."""
        self.selector.config(values=[f"{peak.cpd} {peak.x}" for peak in self.peaks])
        if initial_selected is not None:
            self.selector.current(newindex=initial_selected)
        else:
            self.selector.current(newindex=0)
        self.update_plot()

    def select_callback(self, event):
        """Callback for peak selector; gets the current peak.
        TODO: populate edit panel.
        """
        selected_idx = self.selector.current()
        self.active_peak = self.peaks[selected_idx]
        self.edit_title.config(text=f"Active Peak: {self.active_peak.cpd} {self.active_peak.x}")
        self.edit_var[0].set(self.active_peak.cpd)
        self.edit_var[1].set(str(self.active_peak.get_x()))
        self.edit_var[2].set(str(self.active_peak.get_y()))
        self.edit_var[3].set(str(self.active_peak.get_lw()))
        self.update_plot()

    def refit_peaks(self):
        self.fit_button["state"] = "disabled"
        idx = self.peaks.index(self.active_peak)
        self.fid.unpack_peaks(self.peaks)
        self.fit_progress.grid(row=2, column=1) # Notify progress during curve fitting
        self.fid.curve_fit()
        self.peaks = self.fid.pack_peaks()
        self.assignments = collections.defaultdict(list)
        for pk in self.peaks:
            self.assignments[pk.cpd].append(pk)
        if len(self.peaks) > 0:
            self.active_peak = self.peaks[0]
        else:
            self.active_peak = None
        for cpd in self.assignments:
            self.update_cpd_simspec(cpd)
        self.update_plot()
        self.update_selector(0)
        self.select_callback(None)
        self.fit_progress.grid_forget()
        self.fit_button["state"] = "normal"


    def edit_cpd_callback(self, event):
        """TODO: callback for edit peak compound entry:
        - update peak compound variable
        - update edit box title
        - update dropdown entry
        """
        if self.active_peak is not None:
            oldcpd = self.active_peak.cpd
            val = self.edit_elements["Compound"].get()
            self.active_peak.set_cpd(val)
            self.assignments[oldcpd].remove(self.active_peak)
            self.assignments[val].append(self.active_peak)
            self.update_cpd_simspec(oldcpd)
            self.update_selector(self.peaks.index(self.active_peak))
            self.edit_title.config(text=f"Active Peak: {self.active_peak.cpd} {self.active_peak.x}")
            self.update_plot()

    def edit_shift_callback(self, event):
        """TODO: callback for edit peak shift entry:
        - check for float
        - update peak shift variable
        - update edit box title
        - update dropdown entry
        """
        if self.active_peak is not None:
            try:
                val = float(self.edit_var[1].get())
                apk = self.active_peak
                apk.set_x(val)
                self.update_selector(self.peaks.index(apk))
                self.edit_title.config(text=f"Active Peak: {apk.cpd} {apk.x}")
                self.update_plot()
            except ValueError:
                pass

    def edit_amp_callback(self, event):
        """Change peak Amplitude."""
        if self.active_peak is not None:
            try:
                val = float(self.edit_var[2].get())
                self.active_peak.set_y(val)
                self.update_plot()
            except ValueError:
                pass

    def edit_lw_callback(self, event):
        """TODO: callback for edit peak L-FWHM entry:
        - check for float
        - update peak lw variable
        """
        if self.active_peak is not None:
            try:
                val = float(self.edit_var[3].get())
                self.active_peak.set_lw(val)
                self.update_plot()
            except ValueError:
                pass

    def add_button_callback(self):
        """TODO: callback for add peak button:
        - validate compound and shift not empty
        - check shift for float
        - store cpd and shift, clear entries
        - make new peak with default values
        - add peak to peaklist
        - update assignments and residual
        - update dropdown
        - select new peak in dropdown
        """
        cpd = self.add_var[0].get()
        try:
            shift = float(self.add_var[1].get())
        except ValueError:
            return
        self.add_var[0].set('')
        self.add_var[1].set('')
        if cpd != '':
            height = 100 #self.spec.y(shift)
            newpeak = Peak(self.isotope, shift, height,
                           l_fwhm=0.5, cpd=cpd)
            self.peaks.append(newpeak)
            self.assignments[newpeak.cpd].append(newpeak)
            self.update_selector(len(self.peaks) - 1)
            self.select_callback(None)
            self.update_plot()

    def remove_button_callback(self):
        """TODO: callback for remove peak button:
        - remove peak from peaklist
        - update assignments and residual
        - update dropdown
        - select next peak using formula
        """
        if len(self.peaks) > 1:
            rpk = self.active_peak
            peak_idx = self.peaks.index(rpk)
            self.peaks.remove(rpk)
            self.assignments[rpk.cpd].remove(rpk)
            self.update_selector(min(peak_idx, len(self.peaks) - 1))
            self.select_callback(None)
            self.update_plot()

if __name__ == "__main__":

    iso = Isotope('13C', (0., 200.))

    long_test_peaks = []
    with open(f'{SCDIR}/../../data/test/20210519_13CGlc/cfg_13C.txt', 'r') as rf:
        for line in rf:
            ls = line.strip().split('\t')
            sh = float(ls[0])
            nm = ls[1]
            ht = np.random.randint(100, 200)
            pk = Peak(iso, sh, ht, l_fwhm=0.5, g_fwhm=0.5, cpd=nm)
            long_test_peaks.append(pk)

    test_peaks = [
        Peak(iso, 25.9762, 200., l_fwhm=0.5, cpd="Acetate"), # g_fwhm=0.5, 
        Peak(iso, 18.78, 200., l_fwhm=0.6, cpd="Alanine"), # g_fwhm=0.6, 
        Peak(iso, 53.425, 200., l_fwhm=0.4, cpd="Alanine"), # g_fwhm=0.4, 
    ]

    t = Window(iso, test_peaks, None)
    # a = PlotWindow(iso, long_test_peaks, None)
    # plt.show()