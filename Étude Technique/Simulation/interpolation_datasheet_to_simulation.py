import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# Données extraites de la datasheet avec WebPlotDigitizer 
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
    return np.maximum(spline_interp(I), 1)  # Assurer que l'inductance reste positive

# Générer des points pour la courbe interpolée
currents_interp = np.linspace(-40, 40, 1000)  # Étendre de -40 A à 40 A
inductances_interp = inductance_filtre(currents_interp)

# Tracé des courbes
plt.figure(dpi=200)

# Tracer la courbe des données extraites avec des croix bleues
plt.plot(currents_extracted, inductances_extracted, marker='x', color='blue', linestyle='', label='Données extraites datasheet')

# Tracer la courbe interpolée symétrisée
plt.plot(currents_interp, inductances_interp, color='green', linestyle='-', label='Interpolation spline symétrisée')

# Ajouter une marque et une annotation à 40 A
inductance_40A = inductance_filtre(40)
plt.scatter(40, inductance_40A, color='red', zorder=5)
plt.text(40, inductance_40A, f'{inductance_40A:.2f} µH', color='red', fontsize=10, ha='right')

# Ajouter une marque et une annotation à -40 A
inductance_minus_40A = inductance_filtre(-40)
plt.scatter(-40, inductance_minus_40A, color='red', zorder=5)
plt.text(-40, inductance_minus_40A, f'{inductance_minus_40A:.2f} µH', color='red', fontsize=10, ha='left')

# Ajout des titres et labels
plt.title('Inductance en fonction du courant moyen DC')
plt.xlabel('Courant moyen DC (A)')
plt.ylabel('Inductance (µH)')
plt.grid(True)

# Ajout de la légende
plt.legend()

# Enregistrer le plot au format SVG
plt.savefig('inductance_interpolation_symmetrized.svg', format='svg')

# Affichage
plt.show()
