import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.backends.backend_tkagg as tkagg
import numpy as np
import tkinter as tk

class IntListText:
    def __init__(self, master):
        self.frame = tk.Frame(master)
        vcmd = (master.register(self.validate))
        self.entry = tk.Entry(self.frame, validate='all', validatecommand=(vcmd, '%P', '%V'))
        self.entry.pack()

    def validate(self, value, reason):
        if reason == 'key':
            if all(c in '0123456789 ' for c in value) or value == '':
                return True
            else:
                return False
        return True

class FormField():
    def __init__(self, master, row, cpd):
        self.master = master
        self.row = row
        self.cpd = cpd
        self.label = tk.Label(master, text=cpd)
        self.label.grid(row=row, column=0, sticky=tk.E)
        self.field = IntListText(master)
        self.var = tk.StringVar()
        self.field.frame.grid(row=row, column=1, sticky=tk.W)
        self.field.entry.config(textvariable=self.var)
        self.warning = tk.Label(master, text="Invalid assignments!", fg="red")
        self.warning.grid(row=row, column=2, sticky=tk.W)
        self.warning.grid_forget()
        self.master.rowconfigure(row, weight=1)

    def flag_invalid(self):
        self.warning.grid(row=self.row, column=2, sticky=tk.W)

    def flag_valid(self):
        self.warning.grid_forget()

    def get_contents(self):
        return self.var.get().strip().split()

    def get_values(self):
        return [int(i) for i in self.get_contents()]

    def validate_field(self, maxindex):
        try:
            indices = self.get_values()
            if any(i > maxindex for i in indices):
                self.flag_invalid()
                return False
            else:
                self.flag_valid()
                return True
        except:
            self.flag_invalid()
            return False


class ConflictWindow():
    def __init__(self, ppm, data, ppm_bounds, cid, subpeaks, colliding_refs, assignments_var):
        self.ppm = ppm
        self.data = data
        self.ppm_bounds = ppm_bounds
        self.cid = cid
        self.subpeaks = subpeaks
        self.colliding_refs = colliding_refs
        self.assignments_var = assignments_var

        # Set up plot
        self.fig = mpl.figure.Figure()
        self.plot_ax = self.fig.add_subplot(111)
        self.plot_ax.set_xlim(left=ppm_bounds[0], right=ppm_bounds[1])
        self.plot_ax.plot(ppm, data, lw=0.5)
        self.plot_ax.plot(ppm[subpeaks], data[subpeaks], ls='', marker='*')
        for ci, xi, yi in zip(range(len(subpeaks)), ppm[subpeaks], data[subpeaks]):
            self.plot_ax.annotate(ci, (xi,yi), textcoords='offset points', 
                                  xytext=(0,10), ha='center')

        # Set up window
        self.root = tk.Tk()
        self.root.geometry("1280x800")
        self.root.title(f"Resolve Peaks: Cluster {cid}")
        self.root.config(bg="skyblue")
        self.root.columnconfigure(0, weight=3)
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)

        # Create frame for plot
        self.plot_frame = tk.Frame(self.root)
        self.plot_frame.grid(row=0, column=0, padx=5, pady=5, sticky=tk.E+tk.W+tk.N+tk.S)
        self.canvas = tkagg.FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar = tkagg.NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Create frames for contents
        self.options_frame = tk.Frame(self.root)
        self.options_frame.grid(row=0, column=1, padx=5, pady=5, sticky=tk.E+tk.W+tk.N+tk.S)
        self.options_frame.columnconfigure(0, weight=1)
        self.options_frame.rowconfigure(0, weight=1)
        self.options_frame.rowconfigure(1, weight=2)

        # Initialize collision list
        self.list_frame = tk.Frame(self.options_frame)
        self.list_frame.grid(row=0, column=0, padx=5, pady=5)

        self.list_title = tk.Label(self.list_frame, text="Colliding Peaks")
        self.list_title.grid(row=0, column=0, columnspan=2)
        self.list_frame.rowconfigure(0, weight=1)

        self.list_items = []
        compounds = set([])
        for i, (idx, entry) in enumerate(colliding_refs.iterrows()):
            labels = (
                tk.Label(self.list_frame, text=entry["Shift"]), 
                tk.Label(self.list_frame, text=entry["Compound"]),
            )
            labels[0].grid(row=i+1, column=0, sticky=tk.E)
            labels[1].grid(row=i+1, column=1, sticky=tk.W)
            self.list_items.append(labels)
            self.list_frame.rowconfigure(i+1, weight=1)
            compounds.add(entry["Compound"])
        self.compounds = list(sorted(compounds))
        self.list_frame.columnconfigure(0, weight=1)
        self.list_frame.columnconfigure(1, weight=1)

        # Initialize assignments form
        self.form_frame = tk.Frame(self.options_frame)
        self.form_frame.grid(row=1, column=0, padx=5, pady=5, sticky=tk.E+tk.W+tk.N+tk.S)
        self.form_title = tk.Label(self.form_frame, text="Assign peaks to compounds by index:")
        self.form_title.grid(row=0, column=0, columnspan=3)
        self.form_frame.rowconfigure(0, weight=1)
        self.forms = []
        for i, cpd in enumerate(self.compounds):
            form = FormField(self.form_frame, i+1, cpd)
            self.forms.append(form)
            form.field.entry.bind("<FocusOut>", self.validate_fields_callback)
        self.submit = tk.Button(self.form_frame, text="Submit", command=self.submit_button_callback)
        self.submit.grid(row=len(self.compounds)+2, column=0, columnspan=3)
        self.form_frame.rowconfigure(len(self.compounds)+2, weight=1)
        self.form_frame.columnconfigure(0, weight=1)
        self.form_frame.columnconfigure(1, weight=1)
        self.form_frame.columnconfigure(2, weight=1)

        # Display
        self.root.mainloop()

    def validate_fields_callback(self, event):
        valid = [form.validate_field(len(self.subpeaks) - 1) for form in self.forms]
        if all(valid):
            self.submit["state"] = "normal"
            return True
        else:
            self.submit["state"] = "disabled"
            return False

    def submit_button_callback(self):
        all_entries_valid = self.validate_fields_callback(None)
        if all_entries_valid:
            for form in self.forms:
                indices = np.array(form.get_values())
                if len(indices) > 0:
                    self.assignments_var[form.cpd] = self.subpeaks[indices]
            self.root.destroy()


class FloatField():
    def __init__(self, master, row, label):
        self.master = master
        self.row = row
        self.label = tk.Label(master, text=label)
        self.label.grid(row=row, column=0, sticky=tk.E)
        self.field = tk.Entry(master)
        self.var = tk.DoubleVar()
        self.field.grid(row=row, column=1, sticky=tk.W)
        self.field.config(textvariable=self.var)
        self.warning = tk.Label(master, text="Invalid value!", fg="red")
        self.warning.grid(row=row, column=2, sticky=tk.W)
        self.warning.grid_forget()
        self.master.rowconfigure(row, weight=1)

    def flag_invalid(self):
        self.warning.grid(row=self.row, column=2, sticky=tk.W)

    def flag_valid(self):
        self.warning.grid_forget()

    def get_value(self):
        return self.var.get()

    def validate_field(self):
        try:
            self.var.get()
            self.flag_valid()
            return True
        except:
            self.flag_invalid()
            return False


class CalibrateWindow():
    def __init__(self, ppm, data, stack_var):
        self.ppm = ppm
        self.data = data
        self.refpeaks = [t for t in stack_var.refpeaks.itertuples(index=False, name=None)]
        self.refshifts = [float(pk[0]) for pk in self.refpeaks]
        self.ppm_bounds = stack_var.ppm_bounds
        self.stack_var = stack_var

        # Set up plot
        self.fig = mpl.figure.Figure()
        self.plot_ax = self.fig.add_subplot(111)
        self.plot_ax.set_xlim(left=self.ppm_bounds[0], right=self.ppm_bounds[1])
        self.plot_ax.plot(ppm, data, lw=0.5)

        # Set up window
        self.root = tk.Tk()
        self.root.geometry("1280x800")
        self.root.title(f"Find reference shift.")
        self.root.config(bg="skyblue")
        self.root.columnconfigure(0, weight=3)
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)

        # Create frame for plot
        self.plot_frame = tk.Frame(self.root)
        self.plot_frame.grid(row=0, column=0, padx=5, pady=5, sticky=tk.E+tk.W+tk.N+tk.S)
        self.canvas = tkagg.FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar = tkagg.NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Create frames for contents
        self.options_frame = tk.Frame(self.root)
        self.options_frame.grid(row=0, column=1, padx=5, pady=5, sticky=tk.E+tk.W+tk.N+tk.S)
        self.options_frame.columnconfigure(0, weight=1)
        self.options_frame.rowconfigure(0, weight=1)
        # self.options_frame.rowconfigure(1, weight=2)

        # Initialize assignments form
        self.form_frame = tk.Frame(self.options_frame)
        self.form_frame.grid(row=0, column=0, padx=5, pady=5, sticky=tk.E+tk.W+tk.N+tk.S)
        self.form_title = tk.Label(self.form_frame, text="PPM Shifts:")
        self.form_title.grid(row=0, column=0, columnspan=3)
        self.form_frame.rowconfigure(0, weight=1)
        self.reference_form = FloatField(self.form_frame, 1, "Reference Shift")
        self.reference_form.field.bind("<FocusOut>", self.validate_fields_callback)
        self.experimental_form = FloatField(self.form_frame, 2, "Experimental Shift")
        self.experimental_form.field.bind("<FocusOut>", self.validate_fields_callback)
        self.submit = tk.Button(self.form_frame, text="Submit", command=self.submit_button_callback)
        self.submit.grid(row=3, column=0, columnspan=3)
        self.form_frame.rowconfigure(3, weight=1)
        self.form_frame.columnconfigure(0, weight=1)
        self.form_frame.columnconfigure(1, weight=1)
        self.form_frame.columnconfigure(2, weight=1)

        # Display
        self.root.mainloop()

    def validate_fields_callback(self, event):
        if self.reference_form.validate_field() and self.experimental_form.validate_field():
            self.submit["state"] = "normal"
            return True
        else:
            self.submit["state"] = "disabled"
            return False

    def submit_button_callback(self):
        if self.validate_fields_callback(None):
            ref_shift = self.reference_form.get_value()
            exp_shift = self.experimental_form.get_value()
            self.stack_var.cf = ref_shift - exp_shift
            self.root.destroy()


if __name__ == "__main__":
    w = ConflictWindow()