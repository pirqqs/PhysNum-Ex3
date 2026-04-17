import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
import glob
import math

# ============================================================
# PARAMÈTRES UTILISATEUR
# ============================================================

dossier = "Scan_dt_epsilon_1e-06"

disposition_graphes = {
    "vmax_vs_dt": True,
    "error_hmin_vs_dt": True,
    "error_vmax_vs_dt": True,
    "hmin_vs_epsilon": True,
    "error_hmin_vs_epsilon": True,
    "nsteps_vs_epsilon": True,
    "dt_and_distance_vs_time": True,
    "cost_comparison_fixed_vs_adaptive": True,
    "trajectories": True,
    "energy": True,
    "momentum": True,
    "distance_probe_earth": True,
    "altitude_probe": True,
    "perigee": True,
    "dt_used": True,
    "acceleration_probe": True,
    "drag_power": True,
    "convergence_perigee": True,
    "convergence_altitude_min": True
}

# ============================================================
# DOSSIER DE SORTIE
# ============================================================

dossier_figures = os.path.join(dossier, "figures")
os.makedirs(dossier_figures, exist_ok=True)

# ============================================================
# FONCTIONS UTILES
# ============================================================

def get_probe_speed(data, colonnes, indice_sonde):
    vx = data[:, colonnes[f"vx{indice_sonde}"]]
    vy = data[:, colonnes[f"vy{indice_sonde}"]]
    return np.sqrt(vx**2 + vy**2)

def get_hmin(data, colonnes):
    return np.min(data[:, colonnes["altitude_probe"]])

def get_vmax(data, colonnes, indice_sonde):
    v = get_probe_speed(data, colonnes, indice_sonde)
    return np.max(v)

def get_nsteps(data, colonnes):
    return len(data[:, colonnes["t"]]) - 1

def lire_entete_et_donnees(chemin_fichier):
    entete = None
    with open(chemin_fichier, "r") as f:
        for ligne in f:
            if ligne.startswith("#"):
                entete = ligne.strip()
                break

    if entete is None:
        raise RuntimeError(f"Aucune ligne d'entête commençant par '#' trouvée dans {chemin_fichier}")

    noms_colonnes = entete[1:].strip().split()
    dictionnaire_colonnes = {nom: i for i, nom in enumerate(noms_colonnes)}

    data = np.loadtxt(chemin_fichier, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] != len(noms_colonnes):
        raise RuntimeError(
            f"Incohérence entête/données dans {chemin_fichier} : "
            f"{len(noms_colonnes)} colonnes dans l'entête, {data.shape[1]} dans les données"
        )

    return noms_colonnes, dictionnaire_colonnes, data


def inferer_indices_corps(noms_colonnes):
    indices_corps = []
    i = 0
    while f"x{i}" in noms_colonnes and f"y{i}" in noms_colonnes and f"vx{i}" in noms_colonnes and f"vy{i}" in noms_colonnes:
        indices_corps.append(i)
        i += 1
    return indices_corps


def ligne_coloree(x, y, t, ax, vmin=None, vmax=None):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap="viridis")
    lc.set_array(t)

    if vmin is not None and vmax is not None:
        lc.set_clim(vmin, vmax)

    lc.set_linewidth(2)
    ax.add_collection(lc)
    ax.autoscale()

    return lc


def get_axes(cle_graphe, titre, nombre_jeux_donnees):
    if disposition_graphes[cle_graphe]:
        fig, ax = plt.subplots()
        axes = [ax] * nombre_jeux_donnees
    else:
        ncols = min(3, nombre_jeux_donnees)
        nrows = math.ceil(nombre_jeux_donnees / 3)

        fig, axarr = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
        axes = np.array(axarr).reshape(-1)

        for j in range(nombre_jeux_donnees, len(axes)):
            fig.delaxes(axes[j])

        axes = axes[:nombre_jeux_donnees]

    fig.suptitle(titre)
    return fig, axes


def tracer(ax, x, y, etiquette=None):
    ax.plot(x, y, label=etiquette)

# ============================================================
# LECTURE DES FICHIERS
# ============================================================

fichiers = sorted(glob.glob(os.path.join(dossier, "*.txt")))

jeux_donnees = []
valeurs_parametre = []
nom_parametre = None
colonnes_reference = None
dico_colonnes_reference = None
indices_corps = None

for fichier in fichiers:
    nom = os.path.basename(fichier)[:-4]
    morceaux = nom.split("_")

    try:
        nom_parametre_courant = morceaux[-2]
        valeur = float(morceaux[-1])
    except ValueError:
        print(f"Fichier ignoré car format inattendu : {fichier}")
        continue

    noms_colonnes, dico_colonnes, data = lire_entete_et_donnees(fichier)

    if colonnes_reference is None:
        colonnes_reference = noms_colonnes
        dico_colonnes_reference = dico_colonnes
        indices_corps = inferer_indices_corps(noms_colonnes)
    else:
        if noms_colonnes != colonnes_reference:
            raise RuntimeError(f"Colonnes incohérentes entre fichiers. Problème avec : {fichier}")

    jeux_donnees.append(data)
    valeurs_parametre.append(valeur)
    nom_parametre = nom_parametre_courant

print(f"{len(jeux_donnees)} jeux de données trouvés.")

if len(jeux_donnees) == 0:
    raise RuntimeError("Aucun jeu de données trouvé. Vérifie le nom du dossier et les fichiers de sortie.")

# Tri des jeux de données
ordre = np.argsort(valeurs_parametre)
valeurs_parametre = np.array(valeurs_parametre)[ordre]
jeux_donnees = [jeux_donnees[i] for i in ordre]

noms_colonnes = colonnes_reference
col = dico_colonnes_reference
Nbody = len(indices_corps)

# Interprétation par défaut des indices
indice_terre = 0
indice_lune = 1 if Nbody > 1 else None
indice_sonde = 2 if Nbody > 2 else (Nbody - 1)

colonnes_requises = ["t", f"x{indice_terre}", f"y{indice_terre}", f"x{indice_sonde}", f"y{indice_sonde}"]
for nom in colonnes_requises:
    if nom not in col:
        raise RuntimeError(f"Colonne requise '{nom}' absente de l'entête.")

# ============================================================
# PLAGE DE TEMPS COMMUNE
# ============================================================

tmin = min(data[:, col["t"]].min() for data in jeux_donnees)
tmax = max(data[:, col["t"]].max() for data in jeux_donnees)

cmap = plt.get_cmap("tab10")

# ============================================================
# GRAPHE : vitesse max en fonction de dt
# ============================================================

if disposition_graphes["vmax_vs_dt"]:
    fig, ax = plt.subplots()

    valeurs_vmax = np.array([get_vmax(data, col, indice_sonde) for data in jeux_donnees])

    ax.plot(valeurs_parametre, valeurs_vmax, marker="o", label="numérique")

    vmax_analytique = None
    if vmax_analytique is not None:
        ax.axhline(vmax_analytique, linestyle="--", label="analytique")

    ax.set_xlabel(nom_parametre)
    ax.set_ylabel(r"$v_{\max}$ [m/s]")
    ax.set_title("Vitesse maximale en fonction du pas de temps")
    ax.grid()

    if np.all(valeurs_parametre > 0):
        ax.set_xscale("log")

    ax.legend()
    fig.savefig(os.path.join(dossier_figures, "vmax_vs_dt.png"), dpi=300)

# ============================================================
# GRAPHE : erreur sur hmin en fonction de dt
# ============================================================

if disposition_graphes["error_hmin_vs_dt"]:
    fig, ax = plt.subplots()

    h_analytique = 10000.0
    valeurs_hmin = np.array([get_hmin(data, col) for data in jeux_donnees])
    erreur_hmin = np.abs(valeurs_hmin - h_analytique)

    ax.plot(valeurs_parametre, erreur_hmin, marker="o")
    ax.set_xlabel(nom_parametre)
    ax.set_ylabel(r"$|h_{\min}^{num} - h_{\min}^{ana}|$ [m]")
    ax.set_title("Erreur sur l'altitude minimale en fonction du pas de temps")
    ax.grid()

    if np.all(valeurs_parametre > 0):
        ax.set_xscale("log")
    if np.all(erreur_hmin > 0):
        ax.set_yscale("log")

    fig.savefig(os.path.join(dossier_figures, "error_hmin_vs_dt.png"), dpi=300)

# ============================================================
# GRAPHE : erreur sur vmax en fonction de dt
# ============================================================

if disposition_graphes["error_vmax_vs_dt"]:
    vmax_analytique = None

    if vmax_analytique is not None:
        fig, ax = plt.subplots()

        valeurs_vmax = np.array([get_vmax(data, col, indice_sonde) for data in jeux_donnees])
        erreur_vmax = np.abs(valeurs_vmax - vmax_analytique)

        ax.plot(valeurs_parametre, erreur_vmax, marker="o")
        ax.set_xlabel(nom_parametre)
        ax.set_ylabel(r"$|v_{\max}^{num} - v_{\max}^{ana}|$ [m/s]")
        ax.set_title("Erreur sur la vitesse maximale en fonction du pas de temps")
        ax.grid()

        if np.all(valeurs_parametre > 0):
            ax.set_xscale("log")
        if np.all(erreur_vmax > 0):
            ax.set_yscale("log")

        fig.savefig(os.path.join(dossier_figures, "error_vmax_vs_dt.png"), dpi=300)

# ============================================================
# GRAPHE : hmin en fonction de epsilon
# ============================================================

if disposition_graphes["hmin_vs_epsilon"] and nom_parametre == "epsilon":
    fig, ax = plt.subplots()

    valeurs_hmin = np.array([get_hmin(data, col) for data in jeux_donnees])

    ax.plot(valeurs_parametre, valeurs_hmin, marker="o", label="numérique")
    ax.axhline(10000.0, linestyle="--", label="analytique")

    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel(r"$h_{\min}$ [m]")
    ax.set_title("Altitude minimale en fonction de la tolérance adaptative")
    ax.grid()
    ax.set_xscale("log")
    ax.legend()

    fig.savefig(os.path.join(dossier_figures, "hmin_vs_epsilon.png"), dpi=300)

# ============================================================
# GRAPHE : erreur sur hmin en fonction de epsilon
# ============================================================

if disposition_graphes["error_hmin_vs_epsilon"] and nom_parametre == "epsilon":
    fig, ax = plt.subplots()

    h_analytique = 10000.0
    valeurs_hmin = np.array([get_hmin(data, col) for data in jeux_donnees])
    erreur_hmin = np.abs(valeurs_hmin - h_analytique)

    ax.plot(valeurs_parametre, erreur_hmin, marker="o")
    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel(r"$|h_{\min}^{num} - h_{\min}^{ana}|$ [m]")
    ax.set_title("Erreur sur l'altitude minimale en fonction de la tolérance adaptative")
    ax.grid()
    ax.set_xscale("log")

    if np.all(erreur_hmin > 0):
        ax.set_yscale("log")

    fig.savefig(os.path.join(dossier_figures, "error_hmin_vs_epsilon.png"), dpi=300)

# ============================================================
# GRAPHE : nombre de pas en fonction de epsilon
# ============================================================

if disposition_graphes["nsteps_vs_epsilon"] and nom_parametre == "epsilon":
    fig, ax = plt.subplots()

    valeurs_nsteps = np.array([get_nsteps(data, col) for data in jeux_donnees])

    ax.plot(valeurs_parametre, valeurs_nsteps, marker="o")
    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel("nombre de pas de temps")
    ax.set_title("Coût numérique en fonction de la tolérance adaptative")
    ax.grid()
    ax.set_xscale("log")

    fig.savefig(os.path.join(dossier_figures, "nsteps_vs_epsilon.png"), dpi=300)

# ============================================================
# GRAPHE : dt utilisé et distance en fonction du temps
# ============================================================

if disposition_graphes["dt_and_distance_vs_time"] and "dt_used" in col and "r_probe_earth" in col:
    fig, ax1 = plt.subplots()

    data = jeux_donnees[-1]
    t = data[:, col["t"]]
    dt_used = data[:, col["dt_used"]]
    r = data[:, col["r_probe_earth"]]

    ax1.plot(t, dt_used, label=r"$\Delta t$ utilisé")
    ax1.set_xlabel("t [s]")
    ax1.set_ylabel(r"$\Delta t$ utilisé [s]")
    ax1.grid()

    ax2 = ax1.twinx()
    ax2.plot(t, r, linestyle="--", label="distance sonde-Terre")
    ax2.set_ylabel("distance sonde-Terre [m]")

    ax1.set_title("Pas de temps adaptatif et distance sonde-Terre")

    fig.savefig(os.path.join(dossier_figures, "dt_and_distance_vs_time.png"), dpi=300)

# ============================================================
# GRAPHE : trajectoires dans l'espace réel
# ============================================================

if any(disposition_graphes.values()):
    fig, axes = get_axes("trajectories", "Trajectoires dans l'espace réel", len(jeux_donnees))
    lc_ref = None

    for i, data in enumerate(jeux_donnees):
        t = data[:, col["t"]]

        xT = data[:, col[f"x{indice_terre}"]]
        yT = data[:, col[f"y{indice_terre}"]]

        xS = data[:, col[f"x{indice_sonde}"]]
        yS = data[:, col[f"y{indice_sonde}"]]

        xS_rel = xS - xT
        yS_rel = yS - yT

        lc = ligne_coloree(xS_rel, yS_rel, t, axes[i], tmin, tmax)

        if i == 0:
            axes[i].scatter([0], [0], marker="o", s=40, label="Terre")

        axes[i].set_xlabel("x relatif à la Terre [m]")
        axes[i].set_ylabel("y relatif à la Terre [m]")
        axes[i].grid()
        axes[i].axis("equal")

        if not disposition_graphes["trajectories"]:
            axes[i].set_title(f"{nom_parametre} = {valeurs_parametre[i]}")
        else:
            axes[i].legend()

        if lc_ref is None:
            lc_ref = lc

    cbar = fig.colorbar(lc_ref, ax=axes)
    cbar.set_label("temps [s]")

    fig.savefig(os.path.join(dossier_figures, "trajectoires_sonde_relative_terre.png"), dpi=300)

# ============================================================
# GRAPHE : énergie mécanique totale
# ============================================================

if "E" in col:
    fig, axes = get_axes("energy", "Énergie mécanique totale", len(jeux_donnees))

    for i, data in enumerate(jeux_donnees):
        t = data[:, col["t"]]
        E = data[:, col["E"]]

        tracer(axes[i], t, E, etiquette=f"{nom_parametre}={valeurs_parametre[i]}")
        axes[i].set_xlabel("t [s]")
        axes[i].set_ylabel("E [J]")
        axes[i].grid()

        if not disposition_graphes["energy"]:
            axes[i].set_title(f"{nom_parametre} = {valeurs_parametre[i]}")

    if disposition_graphes["energy"]:
        axes[0].legend()

    fig.savefig(os.path.join(dossier_figures, "energie_all.png"), dpi=300)

# ============================================================
# GRAPHE : quantité de mouvement totale
# ============================================================

if "Px" in col and "Py" in col:
    fig, axes = get_axes("momentum", "Quantité de mouvement totale", len(jeux_donnees))

    for i, data in enumerate(jeux_donnees):
        t = data[:, col["t"]]
        Px = data[:, col["Px"]]
        Py = data[:, col["Py"]]
        norme_P = np.sqrt(Px**2 + Py**2)

        tracer(axes[i], t, norme_P, etiquette=f"{nom_parametre}={valeurs_parametre[i]}")
        axes[i].set_xlabel("t [s]")
        axes[i].set_ylabel(r"$|P|$ [kg m/s]")
        axes[i].grid()

        if not disposition_graphes["momentum"]:
            axes[i].set_title(f"{nom_parametre} = {valeurs_parametre[i]}")

    if disposition_graphes["momentum"]:
        axes[0].legend()

    fig.savefig(os.path.join(dossier_figures, "quantite_mouvement_all.png"), dpi=300)

# ============================================================
# GRAPHE : distance sonde-Terre
# ============================================================

if "r_probe_earth" in col:
    fig, axes = get_axes("distance_probe_earth", "Distance sonde-Terre", len(jeux_donnees))

    for i, data in enumerate(jeux_donnees):
        t = data[:, col["t"]]
        r = data[:, col["r_probe_earth"]]

        tracer(axes[i], t, r, etiquette=f"{nom_parametre}={valeurs_parametre[i]}")
        axes[i].set_xlabel("t [s]")
        axes[i].set_ylabel("distance sonde-Terre [m]")
        axes[i].grid()

        if not disposition_graphes["distance_probe_earth"]:
            axes[i].set_title(f"{nom_parametre} = {valeurs_parametre[i]}")

    if disposition_graphes["distance_probe_earth"]:
        axes[0].legend()

    fig.savefig(os.path.join(dossier_figures, "distance_sonde_terre_all.png"), dpi=300)

# ============================================================
# GRAPHE : altitude de la sonde
# ============================================================

if "altitude_probe" in col:
    fig, axes = get_axes("altitude_probe", "Altitude de la sonde", len(jeux_donnees))

    for i, data in enumerate(jeux_donnees):
        t = data[:, col["t"]]
        altitude = data[:, col["altitude_probe"]]

        tracer(axes[i], t, altitude, etiquette=f"{nom_parametre}={valeurs_parametre[i]}")
        axes[i].set_xlabel("t [s]")
        axes[i].set_ylabel("altitude [m]")
        axes[i].grid()

        if not disposition_graphes["altitude_probe"]:
            axes[i].set_title(f"{nom_parametre} = {valeurs_parametre[i]}")

    if disposition_graphes["altitude_probe"]:
        axes[0].legend()

    fig.savefig(os.path.join(dossier_figures, "altitude_sonde_all.png"), dpi=300)

# ============================================================
# GRAPHE : périgée courant
# ============================================================

if "perigee" in col:
    fig, axes = get_axes("perigee", "Périgée courant", len(jeux_donnees))

    for i, data in enumerate(jeux_donnees):
        t = data[:, col["t"]]
        perigee = data[:, col["perigee"]]

        tracer(axes[i], t, perigee, etiquette=f"{nom_parametre}={valeurs_parametre[i]}")
        axes[i].set_xlabel("t [s]")
        axes[i].set_ylabel("périgée [m]")
        axes[i].grid()

        if not disposition_graphes["perigee"]:
            axes[i].set_title(f"{nom_parametre} = {valeurs_parametre[i]}")

    if disposition_graphes["perigee"]:
        axes[0].legend()

    fig.savefig(os.path.join(dossier_figures, "perigee_all.png"), dpi=300)

# ============================================================
# GRAPHE : pas de temps utilisé
# ============================================================

if "dt_used" in col:
    fig, axes = get_axes("dt_used", "Pas de temps utilisé", len(jeux_donnees))

    for i, data in enumerate(jeux_donnees):
        t = data[:, col["t"]]
        dt_used = data[:, col["dt_used"]]

        tracer(axes[i], t, dt_used, etiquette=f"{nom_parametre}={valeurs_parametre[i]}")
        axes[i].set_xlabel("t [s]")
        axes[i].set_ylabel(r"$\Delta t$ utilisé [s]")
        axes[i].grid()

        if not disposition_graphes["dt_used"]:
            axes[i].set_title(f"{nom_parametre} = {valeurs_parametre[i]}")

    if disposition_graphes["dt_used"]:
        axes[0].legend()

    fig.savefig(os.path.join(dossier_figures, "dt_utilise_all.png"), dpi=300)

# ============================================================
# GRAPHE : accélération de la sonde
# ============================================================

if "a_probe" in col:
    fig, axes = get_axes("acceleration_probe", "Accélération de la sonde", len(jeux_donnees))

    for i, data in enumerate(jeux_donnees):
        t = data[:, col["t"]]
        a = data[:, col["a_probe"]]

        tracer(axes[i], t, a, etiquette=f"{nom_parametre}={valeurs_parametre[i]}")
        axes[i].set_xlabel("t [s]")
        axes[i].set_ylabel("accélération de la sonde [m/s²]")
        axes[i].grid()

        if not disposition_graphes["acceleration_probe"]:
            axes[i].set_title(f"{nom_parametre} = {valeurs_parametre[i]}")

    if disposition_graphes["acceleration_probe"]:
        axes[0].legend()

    fig.savefig(os.path.join(dossier_figures, "acceleration_sonde_all.png"), dpi=300)

# ============================================================
# GRAPHE : puissance de traînée
# ============================================================

if "P_drag" in col:
    fig, axes = get_axes("drag_power", "Puissance de traînée", len(jeux_donnees))

    for i, data in enumerate(jeux_donnees):
        t = data[:, col["t"]]
        Pdrag = data[:, col["P_drag"]]

        tracer(axes[i], t, Pdrag, etiquette=f"{nom_parametre}={valeurs_parametre[i]}")
        axes[i].set_xlabel("t [s]")
        axes[i].set_ylabel("puissance de traînée [W]")
        axes[i].grid()

        if not disposition_graphes["drag_power"]:
            axes[i].set_title(f"{nom_parametre} = {valeurs_parametre[i]}")

    if disposition_graphes["drag_power"]:
        axes[0].legend()

    fig.savefig(os.path.join(dossier_figures, "puissance_trainee_all.png"), dpi=300)

# ============================================================
# GRAPHE : convergence du périgée final
# ============================================================

if disposition_graphes["convergence_perigee"] and "perigee" in col:
    fig, ax = plt.subplots()

    perigee_final = np.array([data[-1, col["perigee"]] for data in jeux_donnees])

    ax.plot(valeurs_parametre, perigee_final, marker="o")
    ax.set_xlabel(nom_parametre)
    ax.set_ylabel("périgée final [m]")
    ax.set_title("Convergence du périgée final")
    ax.grid()

    if np.all(valeurs_parametre > 0):
        ax.set_xscale("log")

    fig.savefig(os.path.join(dossier_figures, "convergence_perigee.png"), dpi=300)

# ============================================================
# GRAPHE : convergence de l'altitude minimale
# ============================================================

if disposition_graphes["convergence_altitude_min"] and "altitude_probe" in col:
    fig, ax = plt.subplots()

    altitudes_min = np.array([np.min(data[:, col["altitude_probe"]]) for data in jeux_donnees])

    ax.plot(valeurs_parametre, altitudes_min, marker="o")
    ax.set_xlabel(nom_parametre)
    ax.set_ylabel("altitude minimale [m]")
    ax.set_title("Convergence de l'altitude minimale")
    ax.grid()

    if np.all(valeurs_parametre > 0):
        ax.set_xscale("log")

    fig.savefig(os.path.join(dossier_figures, "convergence_altitude_min.png"), dpi=300)

plt.show()