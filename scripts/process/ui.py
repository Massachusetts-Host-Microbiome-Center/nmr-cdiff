import os

import matplotlib.pyplot as plt
import matplotlib.widgets as mpw
import matplotlib.colors as mpc
import nmrglue as ng
import numpy as np

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
        self.gw = g_fwhm  # Voigt curve Gaussian full-width half maximum
        self.cpd = cpd
        self.trace = self.render_trace(resolution=100000)

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def set_x(self, x):
        self.x = float(x)

    def set_y(self, y):
        self.y = float(y)

    def get_curve_param_dict(self):
        pars = {
            'shift': x,
            'height': y,
            'l_fwhm': lw,
            'g_fwhm': gw,
        }

    def ng_params(self):
        return [self.x, self.lw, self.gw]

    def ng_params_i(self, res):
        x0 = self.isotope.ppm_min # may need to switch these two
        x1 = self.isotope.ppm_max
        def tx(x): return res*(x - x0)/(x1 - x0)
        return [[[tx(x) for x in self.ng_params()]]]

    def ng_amp(self):
        return [[self.y]]

    def render_trace(self, resolution):
        params = self.ng_params_i(resolution)
        amps = self.ng_amp()
        self.trace = ng.analysis.linesh.sim_NDregion(
            (resolution,), ['v'], params, amps,
        )
        return self.trace

    def align_peak_to_curve(self):
        if isinstance(self.x, float) and isinstance(self.y, float):
            self.set_x(x)
            self.set_y(y)
        else:
            typestring = f"{type(self.x)} and {type(self.y)}"
            raise ValueError("invalid types for coordinates: " + typestring)

    def update(self, **kwargs):
        for k, v in kwargs.items():
            self.__setattr__(k, v)
        self.render_trace(len(self.trace))

class TextBoxObject():
    def __init__(self, parent, position, label, **kwargs):
        self.parent = parent
        self.ax = plt.axes(position)
        if 'initial' in kwargs:
            initial = kwargs.pop('initial')
            if isinstance(initial, float):
                initial = f"{initial:.4f}"
        else:
            initial = ''
        self.obj = mpw.TextBox(self.ax, label, initial=initial, **kwargs)
        self.obj.on_submit(self.on_changed)

    def get_val(self):
        return self.obj.text

    def set_text(self, val):
        if isinstance(val, float):
            val = f"{val:.4f}"
        self.obj.set_val(val)

    def on_changed(self, val):
        self.parent.update(self, val)

    def set_active(self, bool):
        self.obj.set_active(bool)

    def set_visible(self, bool):
        self.ax.set_visible(bool)

class RadioButtonsObject():
    def __init__(self, parent, position, labels):
        self.parent = parent
        self.pos_params = position
        self.ax = plt.axes(self.calc_position(len(labels)))
        self.map = {l: i for i, l in enumerate(labels)}
        self.obj = mpw.RadioButtons(self.ax, labels=labels)
        self.obj.on_clicked(self.on_changed)
        self.labels = labels
        self.idx_selected = 0

    def select(self, idx):
        self.obj.set_active(idx)
        self.idx_selected = idx

    def calc_position(self, numel):
        ht_total = min((numel + 1)*self.pos_params[3], 0.6)
        bottom = self.pos_params[1] - ht_total
        return [self.pos_params[0], bottom, self.pos_params[2], ht_total]

    def refresh(self, labels):
        self.obj._active = False
        self.ax.set_visible(False)
        self.ax = plt.axes(self.calc_position(len(labels)))
        self.map = {l: i for i, l in enumerate(labels)}
        self.obj = mpw.RadioButtons(self.ax, labels=labels)
        self.obj.on_clicked(self.on_changed)
        self.labels = labels
        self.idx_selected = 0

    def get_selected(self):
        return self.idx_selected

    def get_label(self):
        return self.obj.value_selected

    def set_label(self, idx, text):
        self.obj.labels[idx].set_text(text)
        self.obj.value_selected = text

    def on_changed(self, val):
        # stupid way to get selected index
        for i, p in enumerate(self.obj.circles):
            if mpc.to_rgba(p.get_facecolor(), 0) == mpc.to_rgba(self.obj.activecolor, 0):
                self.idx_selected = i
        self.parent.update(self, self.get_selected())

class ButtonObject():
    def __init__(self, parent, position, label, on_press, active=True, activators=None):
        self.parent = parent
        self.ax = plt.axes(position)
        self.obj = mpw.Button(self.ax, label=label)
        self.obj.on_clicked(on_press)
        self.set_active(active)
        self.activators = activators

    def set_active(self, bool):
        self.obj.set_active(bool)
        if bool:
            self.ax.set_facecolor([i/255. for i in (175, 219, 245)])
            self.obj.color = [i/255. for i in (175, 219, 245)]
            self.obj.hovercolor = [i/255. for i in (240, 255, 255)]
        else:
            self.ax.set_facecolor([0.85, 0.85, 0.85])
            self.obj.color = [0.85, 0.85, 0.85]
            self.obj.hovercolor = [0.95, 0.95, 0.95]

    def check_active(self):
        if self.activators is not None:
            activated = True
            for activator in self.activators:
                if activator.get_val() == '':
                    activated = False
            self.set_active(activated)

    def set_visible(self, bool):
        self.ax.set_visible(bool)

class SliderObject():
    def __init__(self, parent, position, label, vallim, valinit=0.5):
        self.parent = parent
        self.ax = plt.axes(position)
        self.obj = mpw.Slider(self.ax, label, *vallim, valinit=valinit)

class SliderTextObject():
    def __init__(self, parent, slider_position, text_position, label, vallim,
                 valinit=0.5):
        self.parent = parent
        self.val = valinit
        self.s_ax = plt.axes(slider_position)
        self.slider = mpw.Slider(self.s_ax, label, *vallim, valinit=valinit)
        self.t_ax = plt.axes(text_position)
        self.text = mpw.TextBox(self.t_ax, '', initial=valinit)
        self.slider.on_changed(self.on_slider_changed)
        self.text.on_submit(self.on_text_changed)

    def set_val(self, val):
        self.slider.set_val(val)
        self.text.set_val(val)

    def on_slider_changed(self, val):
        if self.val != val:
            self.val = val
            self.parent.update(self, self.val)
            self.text.set_val(f"{val:.4f}")

    def on_text_changed(self, val):
        fval = float(val)
        if f"{self.val:8.2f}" != f"{fval:.4f}":
            self.val = fval
            self.parent.update(self, self.val)
            self.slider.set_val(fval)

class PlotWindow():
    def __init__(self, isotope, peaks, spec):
        """Initialize plot window.
        
        Parameters:
        -- isotope : Isotope object matching labeled isotope of the NMR spectrum
        -- peaks : list of Peak objects detected in the spectrum
        -- spec : 2D array-like: (ppm, trace), currently ignored
        """
        self.res = 100000
        self.peaks = peaks
        if len(peaks) > 0:
            self.active_peak = peaks[0]
        else:
            self.active_peak = None
        self.isotope = isotope
        self.fig = plt.figure(figsize=(10, 10))
        self.text_canvas = plt.axes([0, 0, 1, 1])
        self.plot_ax = plt.axes([0.05, 0.4, 0.65, 0.55])
        spec = (
            np.linspace(isotope.ppm_min, isotope.ppm_max, num=20000),
            sum([peak.render_trace(20000) for peak in peaks])
        )
        simspec = sum([peak.render_trace(self.res) for peak in peaks])
        l1 = self.plot_ax.plot(
            spec[0], 
            spec[1], 
            c='slategray', 
            lw=2,
        ) # plot NMR spectrum
        l2 = self.plot_ax.plot(
            np.linspace(isotope.ppm_min, isotope.ppm_max, num=self.res), 
            simspec, 
            c='royalblue', 
            lw=1,
        )
        l3 = self.plot_ax.plot(
            np.linspace(isotope.ppm_min, isotope.ppm_max, num=self.res), 
            self.active_peak.render_trace(self.res),
            c='darkorange',
            lw=1,
        )
        self.spec = Trace(l1[0]) # should not be edited
        self.simspec = Trace(l2[0])
        self.overlay = Trace(l3[0])

        self.peak_selector = RadioButtonsObject(self, [0.71, 0.9, 0.25, 0.02], [f"{pk.cpd} {pk.x}" for pk in self.peaks])
        self.peak_title = self.text_canvas.text(0.375, 0.276, f"Active Peak: {self.peak_selector.get_label()}", ha='center', va='bottom', size='large')
        self.peak_cpd = TextBoxObject(self, [0.2, 0.225, 0.3, 0.04], "Peak compound ", initial=self.active_peak.cpd)
        self.peak_x = TextBoxObject(self, [0.2, 0.175, 0.3, 0.04], "Peak ppm shift ", initial=self.active_peak.x)
        self.peak_y = TextBoxObject(self, [0.2, 0.125, 0.3, 0.04], "Peak amplitude ", initial=self.active_peak.y)
        self.peak_lw = SliderTextObject(self, [0.2, 0.075, 0.3, 0.04], [0.5, 0.075, 0.08, 0.04], "Peak Lorentzian FWHM ", [0., 2.], valinit=self.active_peak.lw)
        self.peak_gw = SliderTextObject(self, [0.2, 0.025, 0.3, 0.04], [0.5, 0.025, 0.08, 0.04], "Peak Gaussian FWHM ", [0., 2.], valinit=self.active_peak.gw)
        self.add_peak = ButtonObject(self, [0.91, 0.9, 0.05, 0.04], "+", self.open_add_dialogue)
        self.remove_peak = ButtonObject(self, [0.85, 0.9, 0.05, 0.04], "–", self.open_remove_dialogue)
        cdialogue = TextBoxObject(self, [0.8, 0.125, 0.16, 0.04], "Peak compound ")
        sdialogue = TextBoxObject(self, [0.8, 0.075, 0.16, 0.04], "Chemical shift ")
        self.add_dialogue_elements = {
            'title': self.text_canvas.text(0.835, 0.175, 'Add Peak', ha='center', va='bottom', size='large'),
            'cpd': cdialogue,
            'shift': sdialogue,
            'cancel': ButtonObject(self, [0.71, 0.025, 0.07, 0.04], "Cancel", self.clear_add_dialogue),
            'submit': ButtonObject(self, [0.89, 0.025, 0.07, 0.04], "Submit", self.submit_add_dialogue, active=False, activators=[cdialogue, sdialogue]),
        }
        self.remove_dialogue_elements = {
            'title': self.text_canvas.text(0.835, 0.1, 'Remove Selected Peak?', ha='center', va='bottom', size='large'),
            'cancel': ButtonObject(self, [0.71, 0.025, 0.07, 0.04], "Cancel", self.clear_remove_dialogue),
            'submit': ButtonObject(self, [0.89, 0.025, 0.07, 0.04], "Submit", self.submit_remove_dialogue),
        }
        self.clear_add_dialogue(None)
        self.clear_remove_dialogue(None)

        self.plot_ax.set_xlim(left=isotope.ppm_max, right=isotope.ppm_min)
        ylim = self.plot_ax.get_ylim()

    def populate_panel(self, pk):
        self.peak_title.set_text(f"Active Peak: {self.peak_selector.get_label()}")
        self.peak_cpd.set_text(pk.cpd)
        self.peak_x.set_text(pk.x)
        self.peak_y.set_text(pk.y)
        self.peak_lw.set_val(pk.lw)
        self.peak_gw.set_val(pk.gw)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def refresh_labels(self):
        pk = self.active_peak
        i = self.peak_selector.get_selected()
        self.peak_selector.set_label(i, f"{pk.cpd} {pk.x}")
        self.peak_title.set_text(f"Active Peak: {self.peak_selector.get_label()}")

    def update(self, widget, val):
        if widget is self.peak_selector:
            pk = self.peaks[val]
            self.active_peak = None # prevents the following setters from updating the PlotWindow
            self.populate_panel(pk)
            self.active_peak = pk
            self.overlay.set(ydata=self.active_peak.trace)
        elif widget is self.add_dialogue_elements['cpd'] or widget is self.add_dialogue_elements['shift']:
            self.add_dialogue_elements['submit'].check_active()
        elif self.active_peak is not None:
            if widget is self.peak_cpd:
                self.active_peak.update(cpd=val)
                self.refresh_labels()
            else:
                if widget is self.peak_x:
                    self.active_peak.update(x=float(val))
                    self.refresh_labels()
                elif widget is self.peak_y:
                    self.active_peak.update(y=float(val))
                elif widget is self.peak_lw:
                    self.active_peak.update(lw=val)
                elif widget is self.peak_gw:
                    self.active_peak.update(gw=val)
                self.overlay.set(ydata=self.active_peak.trace)
                self.simspec.set(ydata=sum([peak.render_trace(self.res) for peak in self.peaks]))
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def open_add_dialogue(self, arg):
        self.clear_remove_dialogue(None)
        for k, ele in self.add_dialogue_elements.items():
            if k in ['cpd', 'shift']:
                ele.set_text('')
            if k == 'submit':
                ele.check_active()
            elif k != 'title':
                ele.set_active(True)
            ele.set_visible(True)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def clear_add_dialogue(self, arg):
        for k, ele in self.add_dialogue_elements.items():
            if k in ['cpd', 'shift']:
                ele.set_text('')
            if k != 'title':
                ele.set_active(False)
            ele.set_visible(False)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def submit_add_dialogue(self, arg):
        shiftstring = self.add_dialogue_elements['shift'].get_val()
        cpd = self.add_dialogue_elements['cpd'].get_val()
        if shiftstring != '' and cpd != '':
            shift = float(shiftstring)
            height = self.spec.y(shift)
            newpeak = Peak(self.isotope, shift, height,
                           l_fwhm=0.5, g_fwhm=0.5, cpd=cpd)
            self.clear_add_dialogue(arg)
            self.peaks.append(newpeak)
            self.peak_selector.refresh([f"{pk.cpd} {pk.x}" for pk in self.peaks])
            self.peak_selector.select(len(self.peaks) - 1)

    def open_remove_dialogue(self, arg):
        self.clear_add_dialogue(None)
        for k, ele in self.remove_dialogue_elements.items():
            if k != 'title':
                ele.set_active(True)
            ele.set_visible(True)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def clear_remove_dialogue(self, arg):
        for k, ele in self.remove_dialogue_elements.items():
            if k != 'title':
                ele.set_active(False)
            ele.set_visible(False)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def submit_remove_dialogue(self, arg):
        self.clear_remove_dialogue(arg)
        selected = self.peak_selector.get_selected()
        self.peaks.pop(selected) # bye bye!
        self.peak_selector.refresh([f"{pk.cpd} {pk.x}" for pk in self.peaks])
        self.peak_selector.select(min(selected, len(self.peaks) - 1))

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
        Peak(iso, 25.9762, 200., l_fwhm=0.5, g_fwhm=0.5, cpd="Acetate"),
        Peak(iso, 18.78, 200., l_fwhm=0.6, g_fwhm=0.6, cpd="Alanine"),
        Peak(iso, 53.425, 200., l_fwhm=0.4, g_fwhm=0.4, cpd="Alanine"),
    ]

    a = PlotWindow(iso, long_test_peaks, None)
    plt.show()