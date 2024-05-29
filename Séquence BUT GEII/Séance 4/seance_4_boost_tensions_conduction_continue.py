import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk

def generate_square_wave(alpha, period=1, duration=1.5):
    t = np.arange(0, duration * period, 0.01)
    duty_cycle = alpha * period
    signal = np.zeros_like(t)
    signal[t % period < duty_cycle] = 1
    return t, signal

def generate_voltage_signal(command_signal, Vi, Vo, alpha, period):
    voltage_signal = np.zeros_like(command_signal)
    for i in range(len(command_signal)):
        if command_signal[i] == 1:
            voltage_signal[i] = Vi
        else:
            voltage_signal[i] = Vi - Vo
    return voltage_signal

def update_plot():
    alpha = alpha_slider.get()
    Vi = Vi_slider.get()
    Vo = Vi / (1 - alpha)  # Calcul de Vo pour un hacheur boost idéal
    t, command_signal = generate_square_wave(alpha)
    voltage_signal = generate_voltage_signal(command_signal, Vi, Vo, alpha, period)
    alpha_Vi_signal = alpha * Vi * np.ones_like(t)
    ax1.clear()
    ax1.plot(t, command_signal, 'b-', label='Signal de commande')
    ax1.set_ylabel("Signal commande interrupteur ")
    ax1.set_title('Commande avec rapport cyclique de {:.2f}'.format(alpha))
    ax1.set_xlabel('Temps')
    ax2.clear()
    ax2.plot(t, voltage_signal, 'r-', label='Tension $V_L$')
    ax2.set_ylabel('$V_L$ (V)')
    ax2.set_title('Tension $V_L$ aux bornes de la bobine')
    ax2.set_xlabel('Temps')
    ax2.text(0.5, 0.9, 'Max $V_L$ = {:.2f} V\nMin $V_L$ = {:.2f} V'.format(np.max(voltage_signal), np.min(voltage_signal)), 
             transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))
    ax3.clear()
    ax3.plot(t, Vo * np.ones_like(t), 'g-', label='Tension Vo')
    ax3.set_xlabel('Temps')
    ax3.set_ylabel('Tension Vo en sortie de filtre (V)')
    ax3.set_title('Tension Vo en sortie du hacheur')
    ax3.text(0.5, 0.9, 'Vo = {:.2f} V'.format(Vo), transform=ax3.transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))
    ax1.legend(loc='upper right')
    ax2.legend(loc='upper right')
    ax3.legend(loc='upper right')
    plt.tight_layout(pad=2.0, h_pad=2.5)  # Adjust spacing between plots
    canvas.draw()

# Paramètres initiaux
alpha_initial = 0.5  # Rapport cyclique initial
Vi_initial = 20      # Valeur initiale de Vi en volts
period = 1           # Période du signal
duration = 1.5       # Durée d'affichage

# Création de la figure et des axes
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 8))

# Génération du signal initial
t, command_signal = generate_square_wave(alpha_initial, period, duration)
Vo_initial = Vi_initial / (1 - alpha_initial)
voltage_signal = generate_voltage_signal(command_signal, Vi_initial, Vo_initial, alpha_initial, period)
alpha_Vi_signal = alpha_initial * Vi_initial * np.ones_like(t)

# Tracé du signal initial
ax1.plot(t, command_signal, 'b-', label='Signal de commande')
ax1.set_ylabel('Signal')
ax1.set_title('Signal de commande carré avec rapport cyclique de {:.2f}'.format(alpha_initial))
ax1.set_xlabel('Temps')
ax2.plot(t, voltage_signal, 'r-', label='Tension $V_L$')
ax2.set_ylabel('$V_L$ (V)')
ax2.set_title('Tension $V_L$ en fonction du signal de commande')
ax2.set_xlabel('Temps')
ax2.text(0.5, 0.9, 'Max $V_L$ = {:.2f} V\nMin $V_L$ = {:.2f} V'.format(np.max(voltage_signal), np.min(voltage_signal)), 
         transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))
ax3.plot(t, Vo_initial * np.ones_like(t), 'g-', label='Tension Vo')
ax3.set_xlabel('Temps')
ax3.set_ylabel('Tension Vo en sortie de filtre (V)')
ax3.set_title('Tension Vo en sortie du hacheur')
ax3.text(0.5, 0.9, 'Vo = {:.2f} V'.format(Vo_initial), transform=ax3.transAxes,
         fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))
ax1.legend(loc='upper right')
ax2.legend(loc='upper right')
ax3.legend(loc='upper right')

plt.tight_layout(pad=2.0, h_pad=2.5)  # Adjust spacing between plots

# Création de la fenêtre Tkinter
root = tk.Tk()
root.title("Simulation hacheur boost, filtre idéal")

# Ajout du canevas matplotlib dans la fenêtre Tkinter
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

# Ajout du frame pour les sliders
slider_frame = tk.Frame(root)
slider_frame.pack(side=tk.RIGHT, fill=tk.Y)

# Ajout du slider pour régler le rapport cyclique
alpha_slider = tk.Scale(slider_frame, from_=0, to=0.99, resolution=0.01, orient=tk.HORIZONTAL, label="Rapport cyclique",
                         command=lambda event: update_plot())
alpha_slider.set(alpha_initial)
alpha_slider.pack(side=tk.TOP, fill=tk.X)

# Ajout du slider pour régler Vi
Vi_slider = tk.Scale(slider_frame, from_=0, to=100, resolution=1, orient=tk.HORIZONTAL, label="Vi (V)",
                     command=lambda event: update_plot())
Vi_slider.set(Vi_initial)
Vi_slider.pack(side=tk.TOP, fill=tk.X)

tk.mainloop()
