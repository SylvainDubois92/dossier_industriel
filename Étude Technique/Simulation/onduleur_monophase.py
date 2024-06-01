import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import UnivariateSpline

#%% Données extraites de la datasheet avec WebPlotDigitizer
datasheet_points = [
    (-0.08728721338208342, 0.7330465682891363),
    (-0.04274318585263914, 0.7321047377223366),
    (0.002408184359962817, 0.7320090780396828),
    (0.048164329210550805, 0.7310646794277113),
    (0.0927057886948231, 0.728427939047448),
    (0.1381563361699532, 0.7257892726333057),
    (0.18360559962249734, 0.7223031513124318),
    (0.22875183374475544, 0.7188176720028507),
    (0.2745015584824136, 0.7136359988572201),
    (0.3199469698672, 0.7076075128161508),
    (0.3647843965575354, 0.7003091284375693),
    (0.41022595587456384, 0.6917382776763046),
    (0.4556655891577133, 0.6818962445549419),
    (0.5017080710449697, 0.6699342901441638),
    (0.5465371515884963, 0.6571274488718257),
    (0.5916654093945513, 0.6417776008679995),
    (0.6373977998272988, 0.6251552864814895),
    (0.6828233088620028, 0.605991249386077),
    (0.7282456078402417, 0.5847085750238349),
    (0.773664054750723, 0.5608835359413974),
    (0.8184751589780463, 0.5362123259748142),
    (0.8638871857755974, 0.5081500123587176),
    (0.9092953605053911, 0.4775453340224257),
    (0.9543966538371407, 0.4443989329772312),
    (0.9994940951011329, 0.4087101672118415),
    (1.0439816256367949, 0.3704803207488423)
]

# Conversion des données extraites
currents_extracted = np.array([point[0] for point in datasheet_points]) * 22.5
inductances_extracted = np.array([point[1] for point in datasheet_points]) * 45

# Filtrer les données pour ne garder que les valeurs de courant >= 0
valid_indices = currents_extracted >= 0
currents_extracted = currents_extracted[valid_indices]
inductances_extracted = inductances_extracted[valid_indices]

# Ajouter les points symétriques négatifs
currents_sym = np.concatenate([-currents_extracted[::-1], currents_extracted])
inductances_sym = np.concatenate([inductances_extracted[::-1], inductances_extracted])

# Interpolation spline des données symétrisées
spline_interp = UnivariateSpline(currents_sym, inductances_sym, s=0)

# Fonction d'inductance du filtre basée sur l'interpolation
def inductance_filtre(I):
    return np.maximum(spline_interp(I), 1) / 1e6  # Conversion en Henrys
    #return 33e-6

#%% Paramètres de simulation
Udc = 650.5       # Tension continue en entrée du circuit en Volts
Fd = 200e3      # Fréquence de découpage du hacheur en Hertz
Td = 1 / Fd     # Période de découpage du hacheur
Fs = 50         # Fréquence de sortie
Ts = 1 / Fs     # Période de sortie
Te = 1e-7       # Pas de temps de simulation en secondes
Fe = 1 / Te     # Fréquence d'échantillonnage
Th = 1          # Durée de la simulation en secondes
C_filt = 50e-6  # Capacité du filtre
R = 10          # Résistance en Ohms
L = 0.5 * 33e-6 # Inductance de la charge (utilise la valeur initiale pour l'inductance du filtre)

#%% Fonctions de simulation
def alpha_a(t):
    return (1/2)+(np.sqrt(2)*230/Udc)*np.sin(2 * np.pi * Fs * t) # 230 eff

def Vii(alpha, t):
    pos = (t % Td) / Td
    return Udc / 2 if pos < alpha else -Udc / 2

def calc_Va_RL_serie(Vaanew, Vaold1, Vaold2, Ia_old1, L_filt):
    A_0 = 1 + L_filt / (R * Te + L) + L_filt * C_filt / (Te ** 2)
    A_1 = -2 * L_filt * C_filt / (Te ** 2)
    A_2 = L_filt * C_filt / (Te ** 2)
    B_1 = -L_filt / (Te + L / R)
    return (Vaanew - A_1 * Vaold1 - A_2 * Vaold2 - B_1 * Ia_old1) / A_0

def calc_Ia_RL_serie(Vanew, Ia_old1):
    return (Vanew + L * Ia_old1 / Te) / (R + L / Te)

def calc_Iaa(Vaanew, Vanew, Iaa_old1, L_filt):
    return Iaa_old1 + (Te / L_filt) * (Vaanew - Vanew)

def calc_Iacond(Vanew, Vaold1):
    return (C_filt / Te) * (Vanew - Vaold1)

def simulate():
    N = int(Th / Te) + 1
    T = np.linspace(0, Th, N)
    Vaa = np.zeros(N)
    Va = np.zeros(N)
    Ia = np.zeros(N)
    Iaa = np.zeros(N)
    Iacond = np.zeros(N)

    Va_old1 = Va_old2 = Ia_old1 = Iaa_old1 = 0

    for i in range(2, N):
        t_i = T[i]
        alpha_new = alpha_a(t_i)
        Vaa_new = Vii(alpha_new, t_i)
        
        L_filt = inductance_filtre(Iaa_old1)  # Calcul de l'inductance en fonction de Iaa
        
        Va_new = calc_Va_RL_serie(Vaa_new, Va_old1, Va_old2, Ia_old1, L_filt)
        Ia_new = calc_Ia_RL_serie(Va_new, Ia_old1)
        Iaa_new = calc_Iaa(Vaa_new, Va_new, Iaa_old1, L_filt)
        Iacond_new = calc_Iacond(Va_new, Va_old1)

        Vaa[i] = Vaa_new
        Va[i] = Va_new
        Ia[i] = Ia_new
        Iaa[i] = Iaa_new
        Iacond[i] = Iacond_new

        Va_old2, Va_old1 = Va_old1, Va_new
        Ia_old1 = Ia_new
        Iaa_old1 = Iaa_new

    return T, Vaa, Va, Ia, Iaa, Iacond

#%% Fonctions d'analyse spectrale
def analyse_spectrale(Ia):
    FFT_out = np.fft.rfft(Ia)
    spectrum_out_linear = np.abs(FFT_out / len(FFT_out))
    spectrum_out_dB = 10 * np.log10(spectrum_out_linear)
    frequency = np.linspace(0, Fe / 2, len(spectrum_out_linear))
    return frequency, spectrum_out_linear, spectrum_out_dB

def norme():
    normerangs = np.arange(2, 41)
    normeintensites = np.zeros(39)
    for n in normerangs:
        if n == 2:
            normeintensites[n - 2] = 1.08
        elif n == 3:
            normeintensites[n - 2] = 2.30
        elif n == 4:
            normeintensites[n - 2] = 0.43
        elif n == 5:
            normeintensites[n - 2] = 1.14
        elif n == 6:
            normeintensites[n - 2] = 0.30
        elif n == 7:
            normeintensites[n - 2] = 0.77
        elif n == 9:
            normeintensites[n - 2] = 0.40
        elif n == 11:
            normeintensites[n - 2] = 0.33
        elif n == 13:
            normeintensites[n - 2] = 0.21
        elif n % 2 == 0:
            normeintensites[n - 2] = 0.23 * 8 / n
        elif n % 2 == 1:
            normeintensites[n - 2] = 0.15 * 15 / n
        else:
            print("Problème dans la fonction norme")
    return normerangs, normeintensites

def analyse_expe_norme(rangs, freq, linearspec):
    I_rang = np.zeros(len(rangs))
    F_rang = np.zeros(len(rangs))
    for i in range(len(rangs)):
        F_n_i = int(rangs[i] * Fs / (Fe / 2) * len(linearspec))
        F_rang[i] = freq[F_n_i]
        I_rang[i] = linearspec[F_n_i]
    return F_rang, I_rang

#%% Fonction de tracé
def plot_spectrum_with_norms(frequency, spectrum_out_linear, simu_f_rangs, norme_intensites, Fs):
    plt.figure(dpi=200)  # Résolution de la figure
    plt.plot(frequency, spectrum_out_linear, 'g', label='$I_{a}$', linewidth=1)  # Trace le spectre en vert

    # Listes pour les croix rouges et bleues
    blue_x = []
    blue_y = []
    red_x = []
    red_y = []

    # Tracer les croix pour chaque point de la norme
    for i in range(len(simu_f_rangs)):
        idx = np.argmin(np.abs(frequency - simu_f_rangs[i]))
        norme_value = norme_intensites[i]
        spectrum_value = spectrum_out_linear[idx]

        if norme_value > spectrum_value:  # Vérifier si le point est au-dessus du spectre
            blue_x.append(simu_f_rangs[i])
            blue_y.append(norme_value)
        else:
            red_x.append(simu_f_rangs[i])
            red_y.append(norme_value)

    # Tracer les croix bleues et rouges
    plt.scatter(blue_x, blue_y, color='blue', marker='x', label='Norme respectée')
    plt.scatter(red_x, red_y, color='red', marker='x', label='Norme non respectée')

    plt.xlim(0, 40 * Fs)  # Limite de l'axe des x inchangée
    plt.xlabel('Fréquence (Hz)')
    plt.ylabel('Module en linéaire')
    plt.title("Spectre linéaire du courant absorbé par la charge avec les normes")
    plt.legend()
    plt.grid()
    plt.savefig("spectrum_linear.svg")  # Enregistrer la figure en SVG
    plt.show()

    # Conversion de la norme en dB
    norme_intensites_dB = 10 * np.log10(norme_intensites)

    # Création du graphe en dB
    spectrum_out_dB = 10 * np.log10(spectrum_out_linear)
    
    plt.figure(dpi=200)  # Résolution de la figure
    plt.plot(frequency, spectrum_out_dB, 'g', label='$I_{a}$', linewidth=1)  # Trace le spectre en vert

    # Listes pour les croix rouges et bleues en dB
    blue_x_dB = []
    blue_y_dB = []
    red_x_dB = []
    red_y_dB = []

    # Tracer les croix pour chaque point de la norme en dB
    for i in range(len(simu_f_rangs)):
        idx = np.argmin(np.abs(frequency - simu_f_rangs[i]))
        norme_value_dB = norme_intensites_dB[i]
        spectrum_value_dB = spectrum_out_dB[idx]

        if norme_value_dB > spectrum_value_dB:  # Vérifier si le point est au-dessus du spectre
            blue_x_dB.append(simu_f_rangs[i])
            blue_y_dB.append(norme_value_dB)
        else:
            red_x_dB.append(simu_f_rangs[i])
            red_y_dB.append(norme_value_dB)

    # Tracer les croix bleues et rouges en dB
    plt.scatter(blue_x_dB, blue_y_dB, color='blue', marker='x', label='Norme respectée')
    plt.scatter(red_x_dB, red_y_dB, color='red', marker='x', label='Norme non respectée')

    plt.xlim(0, 40 * Fs)  # Limite de l'axe des x inchangée
    plt.ylim(-50, 20)  # Limite de l'axe des y entre -50 dB et +20 dB
    plt.xlabel('Fréquence (Hz)')
    plt.ylabel('Module en dB')
    plt.title("Spectre en dB du courant absorbé par la charge avec les normes")
    plt.legend()
    plt.grid()
    plt.savefig("spectrum_dB.svg")  # Enregistrer la figure en SVG
    plt.show()

#%% Visualisation temporelle
def plot_last_two_periods(T, Vaa, Va, Ia, Iaa, Iacond, Ts):
    last_two_periods_idx = np.where((T >= T[-1] - 2 * Ts) & (T <= T[-1]))
    
    plt.figure(dpi=200)
    plt.plot(T[last_two_periods_idx], Vaa[last_two_periods_idx], label='Vaa')
    plt.xlabel('Temps (s)')
    plt.ylabel('Tension (V)')
    plt.title('Tensions bus DC découpée')
    plt.legend()
    plt.grid()
    plt.savefig("Vaa_last_two_periods.svg")
    plt.show()
    
    plt.figure(dpi=200)
    plt.plot(T[last_two_periods_idx], Va[last_two_periods_idx], label='Va')
    plt.ylim(-400, 400) 
    plt.xlabel('Temps (s)')
    plt.ylabel('Tension (V)')
    plt.title('Tension aux bornes de la charge')
    plt.legend()
    plt.grid()
    plt.savefig("Va_last_two_periods.svg")
    plt.show()
    
    plt.figure(dpi=200)
    plt.plot(T[last_two_periods_idx], Ia[last_two_periods_idx], label='Ia')
    plt.xlabel('Temps (s)')
    plt.ylabel('Courant (A)')
    plt.title('Intensité dans la charge')
    plt.legend()
    plt.grid()
    plt.savefig("Ia_last_two_periods.svg")
    plt.show()
    
    plt.figure(dpi=200)
    plt.plot(T[last_two_periods_idx], Iaa[last_two_periods_idx], label='Iaa')
    plt.xlabel('Temps (s)')
    plt.ylabel('Courant (A)')
    plt.title('Intensité dans Lfilt')
    plt.legend()
    plt.grid()
    plt.savefig("Iaa_last_two_periods.svg")
    plt.show()
    
    plt.figure(dpi=200)
    plt.plot(T[last_two_periods_idx], Iacond[last_two_periods_idx], label='Iacond')
    plt.xlabel('Temps (s)')
    plt.ylabel('Courant (A)')
    plt.title('Intensité dans le condensateur')
    plt.legend()
    plt.grid()
    plt.savefig("Iacond_last_two_periods.svg")
    plt.show()

#%% Fonction principale
def main():
    T, Vaa, Va, Ia, Iaa, Iacond = simulate()
    frequency, spectrum_out_linear, spectrum_out_dB = analyse_spectrale(Ia)
    norme_rangs, norme_intensites = norme()
    simu_f_rangs, simu_intensites = analyse_expe_norme(norme_rangs, frequency, spectrum_out_linear)
    
    # Affichage des spectres
    plt.figure(dpi=200)
    plt.plot(frequency/1000, spectrum_out_dB, 'r', label='$I_{a}$')
    plt.xlim(0, 1.1 * Fd/1000)  # Ajuster la limite de l'axe des x pour aller un peu au-delà de Fd
    plt.ylim(-80, 20) 
    plt.xlabel('Fréquence (kHz)')
    plt.ylabel('Module en dB')
    plt.title("Spectre en dB du courant absorbé par la charge")
    plt.legend()
    plt.grid()
    plt.savefig("spectrum_dB_main.svg")
    plt.show()
    
    plt.figure(dpi=200)
    plt.plot(frequency/1000, spectrum_out_linear, 'r', label='$I_{a}$')
    plt.xlim(0, 1.1 * Fd/1000)  # Ajuster la limite de l'axe des x pour aller un peu au-delà de Fd
    plt.xlabel('Fréquence (kHz)')
    plt.ylabel('Module en linéaire')
    plt.title("Spectre linéaire du courant absorbé par la charge")
    plt.legend()
    plt.grid()
    plt.savefig("spectrum_linear_main.svg")
    plt.show()
    
    plot_spectrum_with_norms(frequency, spectrum_out_linear, simu_f_rangs, norme_intensites, Fs)
    
    # Tableau de conformité à la norme
    result = pd.concat([
        pd.DataFrame(norme_rangs, columns=["Rangs"]),
        pd.DataFrame(simu_f_rangs, columns=["Fréquences simu"]),
        pd.DataFrame(simu_intensites, columns=["Simu Intensites"]),
        pd.DataFrame(norme_intensites, columns=["Norme Intensites"]),
        pd.DataFrame(norme_intensites - simu_intensites, columns=["Marge: N-S"])
    ], axis=1)
    
    print(result.to_string(index=False))
    
    # Affichage des deux dernières périodes des signaux temporels
    plot_last_two_periods(T, Vaa, Va, Ia, Iaa, Iacond, Ts)

main()
