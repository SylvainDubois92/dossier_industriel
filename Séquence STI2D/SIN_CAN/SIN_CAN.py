import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class CANConverterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Caractéristique convertisseur CAN")

        self.frame = ttk.Frame(self.root)
        self.frame.pack(fill=tk.BOTH, expand=True)

        self.bits_label = ttk.Label(self.frame, text="Nombre de bits:")
        self.bits_label.grid(row=0, column=0, padx=(5, 1), pady=5)

        self.bits_value_label = ttk.Label(self.frame, text="3")
        self.bits_value_label.grid(row=0, column=1, padx=(1, 1), pady=5)

        self.bits_slider = ttk.Scale(self.frame, from_=1, to=8, orient=tk.HORIZONTAL, command=self.update_plot)
        self.bits_slider.set(3)
        self.bits_slider.grid(row=0, column=2, padx=(1, 5), pady=5)

        self.ref_label = ttk.Label(self.frame, text="Tension de référence (V):")
        self.ref_label.grid(row=1, column=0, padx=(5, 1), pady=5)

        self.ref_entry = ttk.Entry(self.frame)
        self.ref_entry.insert(0, "8")
        self.ref_entry.grid(row=1, column=1, columnspan=2, padx=(1, 5), pady=5)
        self.ref_entry.bind('<Return>', self.update_plot)

        self.quantum_label = ttk.Label(self.frame, text="résolution : q = 8 / 2^3 = 1.0 V")
        self.quantum_label.grid(row=2, column=0, columnspan=3, padx=5, pady=5)

        self.fig, self.ax = plt.subplots(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.get_tk_widget().grid(row=3, column=0, columnspan=3, padx=5, pady=5)

        self.save_button = ttk.Button(self.frame, text="Enregistrer en SVG", command=self.save_figure)
        self.save_button.grid(row=4, column=0, columnspan=3, pady=10)

        self.update_plot()

    def update_plot(self, event=None):
        n = int(self.bits_slider.get())
        try:
            Uref = float(self.ref_entry.get())
        except ValueError:
            self.ref_entry.delete(0, tk.END)
            self.ref_entry.insert(0, "8")
            Uref = 8.0

        self.bits_value_label.config(text=str(n))

        levels = 2 ** n
        step = Uref / levels
        self.quantum_label.config(text=f"résolution : q = {Uref} / 2^{n} = {step:.2f} V")

        self.ax.clear()

        binary_values = [format(i, f'0{n}b') for i in range(levels)]
        analog_values = np.linspace(0, Uref, levels + 1)
        codes = list(range(levels)) + [levels - 1]

        self.ax.step(analog_values, codes, where='post', color='red')
        self.ax.axvline(Uref, color='gray', linestyle='--')
        self.ax.axvline(Uref - step, color='gray', linestyle='--')

        if n <= 4:
            for i, bin_val in enumerate(binary_values):
                self.ax.text(analog_values[i], codes[i] + 0.1, f'{bin_val}', ha='center')

        # Adding full scale label with double arrow, slightly adjusted to the right
        self.ax.annotate(
            '',
            xy=(Uref - 0.5, 0),
            xytext=(Uref - 0.5, levels - 1),
            arrowprops=dict(arrowstyle='<->', lw=1.5),
            va='center',
            ha='right'
        )
        self.ax.text(Uref - 0.5, (levels - 1) / 2, 'Pleine échelle', va='center', ha='right', rotation=90)

        self.ax.set_xlim(0, Uref)
        self.ax.set_ylim(-0.5, levels - 0.5)
        self.ax.set_xticks(np.arange(0, Uref + step, step))
        self.ax.set_yticks(range(levels))

        self.ax.set_xlabel('Tension Analogique (V)')
        self.ax.set_ylabel('Code Numérique')
        self.ax.set_title('Caractéristique d\'un convertisseur CAN')
        self.ax.grid(True)

        # Adding the legend with number of bits and reference voltage
        self.ax.legend([f'Nombre de bits : {n}', f'Tension de référence : {Uref} V'], loc='upper left')

        self.canvas.draw()

    def save_figure(self):
        self.fig.savefig("convertisseur_CAN.svg", format="svg")
        print("Figure sauvegardée en tant que convertisseur_CAN.svg")

if __name__ == "__main__":
    root = tk.Tk()
    app = CANConverterApp(root)
    root.mainloop()
