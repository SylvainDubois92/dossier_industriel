import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk

def generate_square_wave(alpha, period=1, duration=1.5):
    t = np.arange(0, duration * period, 0.001)
    duty_cycle = alpha * period
    signal = np.zeros_like(t)
    signal[t % period < duty_cycle] = 1
    return t, signal

def generate_voltage_signal(command_signal, Vi):
    voltage_signal = command_signal * Vi
    return voltage_signal

def generate_inductor_current(t, alpha, Vi, L, period, R=1):
    Vo = alpha * Vi  # Tension de sortie moyenne
    delta_i = alpha * (Vi - Vo) * period / L  # Variation maximale du courant
    i_L = np.zeros_like(t)
    Vo_R = Vo / R  # Vo / R
    for i in range(len(t)):
        t_mod = t[i] % period
        if t_mod < alpha * period:
            i_L[i] = max(0, -delta_i / 2 + ((Vi - Vo) / L) * t_mod + Vo_R)
        else:
            i_L[i] = max(0, delta_i / 2 - (Vo / L) * (t_mod - alpha * period) + Vo_R)
    return i_L, delta_i, Vo_R

def update_plot():
    alpha = alpha_slider.get()
    Vi = Vi_slider.get()
    L = L_slider.get()
    frequency = frequency_slider.get()
    period = 1 / frequency
    
    t, command_signal = generate_square_wave(alpha, period)
    voltage_signal = generate_voltage_signal(command_signal, Vi)
    alpha_Vi_signal = alpha * Vi * np.ones_like(t)
    inductor_current, delta_i, Vo_R = generate_inductor_current(t, alpha, Vi, L, period)
    
    ax1.clear()
    ax1.plot(t, command_signal, 'b-', label='Signal de commande')
    ax1.set_ylabel("Signal commande interrupteur ")
    ax1.set_title('Commande avec rapport cyclique de {:.2f}'.format(alpha))
    ax1.set_xlabel("Temps")
    ax1.set_ylim(-0.5, 1.5)
    
    ax2.clear()
    ax2.plot(t, voltage_signal, 'r-', label='Tension Vd')
    ax2.set_ylabel('Vd (V)')
    ax2.set_title('Tension Vd aux bornes de la diode')
    ax2.set_xlabel("Temps")
    ax2.set_ylim(-0.5 * Vi, 1.5 * Vi)
    
    ax3.clear()
    ax3.plot(t, alpha_Vi_signal, 'g-', label='Tension Vo')
    ax3.set_xlabel('Temps')
    ax3.set_ylabel('Tension Vo en sortie de filtre (V)')
    ax3.set_title('Tension Vo en sortie du hacheur')
    ax3.set_ylim(-0.5 * Vi, 1.5 * Vi)
    
    ax4.clear()
    ax4.plot(t, inductor_current, 'm-', label='Courant dans la bobine')
    ax4.set_ylabel('Courant (A)')
    ax4.set_title('Courant dans la bobine')
    ax4.set_xlabel("Temps")
    ax4.set_ylim(0, 40)  # Ajustez cette valeur en fonction de la variation attendue du courant
    ax4.text(0.5, 0.5, 'ΔI = {:.2f} A\nVo/R = {:.2f} A'.format(delta_i, Vo_R), transform=ax4.transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))
    
    ax1.legend(loc='upper right')
    ax2.legend(loc='upper right')
    ax3.legend(loc='upper right')
    ax4.legend(loc='upper right')
    
    plt.tight_layout(pad=2.0, h_pad=2.5)  # Adjust spacing between plots
    canvas.draw()

# Paramètres initiaux
alpha_initial = 0.5  # Rapport cyclique initial
Vi_initial = 20      # Valeur initiale de Vi en volts
L_initial = 1       # Inductance initiale en Henry
frequency_initial = 1 # Fréquence initiale en Hertz
period = 1 / frequency_initial
duration = 1.5       # Durée d'affichage

# Création de la figure et des axes
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 8))

# Génération du signal initial
t, command_signal = generate_square_wave(alpha_initial, period, duration)
voltage_signal = generate_voltage_signal(command_signal, Vi_initial)
alpha_Vi_signal = alpha_initial * Vi_initial * np.ones_like(t)
inductor_current, delta_i, Vo_R = generate_inductor_current(t, alpha_initial, Vi_initial, L_initial, period)

# Tracé du signal initial
ax1.plot(t, command_signal, 'b-', label='Signal de commande')
ax1.set_ylabel('Signal')
ax1.set_title('Signal de commande carré avec rapport cyclique de {:.2f}'.format(alpha_initial))
ax1.set_xlabel("Temps")
ax1.set_ylim(-0.5, 1.5)

ax2.plot(t, voltage_signal, 'r-', label='Tension Vd')
ax2.set_ylabel('Tension (V)')
ax2.set_title('Tension Vd en fonction du signal de commande')
ax2.set_xlabel("Temps")
ax2.set_ylim(-0.5 * Vi_initial, 1.5 * Vi_initial)

ax3.plot(t, alpha_Vi_signal, 'g-', label='Tension Vo')
ax3.set_xlabel('Temps')
ax3.set_ylabel('Tension Vo en sortie de filtre (V)')
ax3.set_title('Tension Vo en sortie du hacheur')
ax3.set_ylim(-0.5 * Vi_initial, 1.5 * Vi_initial)

ax4.plot(t, inductor_current, 'm-', label='Courant dans la bobine')
ax4.set_ylabel('Courant (A)')
ax4.set_title('Courant dans la bobine')
ax4.set_xlabel("Temps")
ax4.set_ylim(0, 40)  # Ajustez cette valeur en fonction de la variation attendue du courant
ax4.text(0.5, 0.5, 'ΔI = {:.2f} A\nVo/R = {:.2f} A'.format(delta_i, Vo_R), transform=ax4.transAxes,
         fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))

ax1.legend(loc='upper right')
ax2.legend(loc='upper right')
ax3.legend(loc='upper right')
ax4.legend(loc='upper right')

plt.tight_layout(pad=2.0, h_pad=2.5)  # Adjust spacing between plots

# Création de la fenêtre Tkinter
root = tk.Tk()
root.title("Simulation hacheur série avec paramètres : filtre idéal et conduction continue")

# Ajout du canevas matplotlib dans la fenêtre Tkinter
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

# Ajout du frame pour les sliders
slider_frame = tk.Frame(root)
slider_frame.pack(side=tk.RIGHT, fill=tk.Y)

# Ajout du slider pour régler le rapport cyclique
alpha_slider = tk.Scale(slider_frame, from_=0, to=1, resolution=0.01, orient=tk.HORIZONTAL, label="Rapport cyclique",
                         command=lambda event: update_plot())
alpha_slider.set(alpha_initial)
alpha_slider.pack(side=tk.TOP, fill=tk.X)

# Ajout du slider pour régler Vi
Vi_slider = tk.Scale(slider_frame, from_=0, to=100, resolution=1, orient=tk.HORIZONTAL, label="Vi (V)",
                     command=lambda event: update_plot())
Vi_slider.set(Vi_initial)
Vi_slider.pack(side=tk.TOP, fill=tk.X)

# Ajout du slider pour régler L
L_slider = tk.Scale(slider_frame, from_=1, to=10, resolution=0.1, orient=tk.HORIZONTAL, label="L (H)",
                     command=lambda event: update_plot())
L_slider.set(L_initial)
L_slider.pack(side=tk.TOP, fill=tk.X)

# Ajout du slider pour régler la fréquence
frequency_slider = tk.Scale(slider_frame, from_=1, to=10, resolution=0.1, orient=tk.HORIZONTAL, label="Fréquence (Hz)",
                            command=lambda event: update_plot())
frequency_slider.set(frequency_initial)
frequency_slider.pack(side=tk.TOP, fill=tk.X)

tk.mainloop()
