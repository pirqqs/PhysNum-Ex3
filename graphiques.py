import matplotlib.pyplot as plt
import numpy as np


def lire_donnees(nom_fichier="trajectory.dat"):
    with open(nom_fichier, "r", encoding="utf-8") as fichier:
        entete = fichier.readline().strip()

    colonnes = entete.lstrip("#").split()
    donnees = np.loadtxt(nom_fichier, comments="#")
    return colonnes, donnees


def indice(colonnes, nom):
    return colonnes.index(nom)


def tracer_trajectoires(colonnes, donnees):
    plt.figure(figsize=(8, 8))
    plt.plot(donnees[:, indice(colonnes, "x0")], donnees[:, indice(colonnes, "y0")], label="Terre")
    plt.plot(donnees[:, indice(colonnes, "x1")], donnees[:, indice(colonnes, "y1")], label="Lune")
    plt.plot(donnees[:, indice(colonnes, "x2")], donnees[:, indice(colonnes, "y2")], label="Artemis II")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Trajectoires dans le plan")
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("trajectoires.png", dpi=200)


def tracer_energie(colonnes, donnees):
    temps_jours = donnees[:, indice(colonnes, "t")] / 86400.0
    energie = donnees[:, indice(colonnes, "E")]

    plt.figure(figsize=(8, 5))
    plt.plot(temps_jours, energie)
    plt.xlabel("Temps [jours]")
    plt.ylabel("Energie mecanique totale [J]")
    plt.title("Evolution de l'energie")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("energie.png", dpi=200)


def tracer_distance(colonnes, donnees):
    temps_jours = donnees[:, indice(colonnes, "t")] / 86400.0
    distance = donnees[:, indice(colonnes, "r_probe_earth")] / 1000.0
    perigee = donnees[:, indice(colonnes, "perigee")] / 1000.0

    plt.figure(figsize=(8, 5))
    plt.plot(temps_jours, distance, label="Distance sonde-Terre")
    plt.plot(temps_jours, perigee, label="Perigee minimal")
    plt.xlabel("Temps [jours]")
    plt.ylabel("Distance [km]")
    plt.title("Distance de la sonde a la Terre")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("distance_sonde_terre.png", dpi=200)


def main():
    colonnes, donnees = lire_donnees("trajectory.dat")
    tracer_trajectoires(colonnes, donnees)
    tracer_energie(colonnes, donnees)
    tracer_distance(colonnes, donnees)
    plt.show()


if __name__ == "__main__":
    main()
