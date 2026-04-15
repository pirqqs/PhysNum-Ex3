#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <valarray>
#include <vector>

using namespace std;
using Etat = valarray<double>;

namespace
{
constexpr double G = 6.674e-11;
constexpr double PI = 3.14159265358979323846;

class Exercice4
{
private:
    size_t N_;
    vector<double> masses_;
    vector<double> rayons_;
    vector<int> trainee_active_;
    double temps_;
    double temps_fin_;
    double dt_;
    double dt_min_;
    double dt_max_;
    double precision_;
    double intervalle_sortie_;
    bool pas_adaptatif_;
    size_t indice_terre_;
    size_t indice_atmosphere_;
    size_t indice_sonde_;
    double rho0_;
    double lambda_;
    double Cx_;
    double section_;
    double perigee_cible_;
    Etat etat_;
    ofstream fichier_sortie_;
    double prochain_temps_sortie_;
    double distance_min_terre_;

public:
    Exercice4()
        : N_(3)
        , masses_({5.972e24, 7.349e22, 8500.0})
        , rayons_({6.371e6, 1.7374e6, 2.51})
        , trainee_active_({0, 0, 1})
        , temps_(0.0)
        , temps_fin_(1036800.0)
        , dt_(20.0)
        , dt_min_(1.0e-3)
        , dt_max_(300.0)
        , precision_(1.0e-8)
        , intervalle_sortie_(60.0)
        , pas_adaptatif_(true)
        , indice_terre_(0)
        , indice_atmosphere_(0)
        , indice_sonde_(2)
        , rho0_(1.2)
        , lambda_(7238.2)
        , Cx_(0.3)
        , section_(0.25 * PI * pow(5.02, 2))
        , perigee_cible_(10000.0)
        , etat_(4 * N_)
        , prochain_temps_sortie_(temps_)
        , distance_min_terre_(1.0e300)
    {
        const vector<double> x_initial = {0.0, 3.844e8, 3.14159e8};
        const vector<double> y_initial = {0.0, 0.0, 0.0};
        const vector<double> vx_initial = {0.0, 0.0, 0.0};
        const vector<double> vy_initial = {0.0, 1022.0, 1200.0};

        for (size_t i = 0; i < N_; ++i) {
            etat_[ix(i)] = x_initial[i];
            etat_[iy(i)] = y_initial[i];
            etat_[ivx(i)] = vx_initial[i];
            etat_[ivy(i)] = vy_initial[i];
        }

        distance_min_terre_ = distanceSondeTerre(etat_);
        fichier_sortie_.open("trajectory.dat");
        if (!fichier_sortie_) {
            throw runtime_error("Impossible d'ouvrir le fichier de sortie.");
        }
        fichier_sortie_ << scientific << setprecision(10);
    }

    void executer()
    {
        fichier_sortie_ << "# t";
        for (size_t i = 0; i < N_; ++i) {
            fichier_sortie_ << " x" << i << " y" << i << " vx" << i << " vy" << i;
        }
        fichier_sortie_ << " E Px Py r_probe_earth perigee\n";

        {
            double quantite_x = 0.0;
            double quantite_y = 0.0;
            quantiteMouvementTotale(quantite_x, quantite_y);

            fichier_sortie_ << temps_;
            for (size_t i = 0; i < N_; ++i) {
                fichier_sortie_ << ' ' << etat_[ix(i)] << ' ' << etat_[iy(i)] << ' ' << etat_[ivx(i)] << ' ' << etat_[ivy(i)];
            }
            fichier_sortie_ << ' ' << energieTotale()
                            << ' ' << quantite_x
                            << ' ' << quantite_y
                            << ' ' << distanceSondeTerre(etat_)
                            << ' ' << distance_min_terre_ - rayons_[indice_terre_]
                            << '\n';
        }

        while (temps_ < temps_fin_) {
            double dt_utilise = min(dt_, temps_fin_ - temps_);
            if (pas_adaptatif_) {
                dt_utilise = pasAdaptatif(dt_utilise);
            } else {
                etat_ = rk4(etat_, dt_utilise);
            }

            temps_ += dt_utilise;
            distance_min_terre_ = min(distance_min_terre_, distanceSondeTerre(etat_));

            if (collisionDetectee()) {
                double quantite_x = 0.0;
                double quantite_y = 0.0;
                quantiteMouvementTotale(quantite_x, quantite_y);

                fichier_sortie_ << temps_;
                for (size_t i = 0; i < N_; ++i) {
                    fichier_sortie_ << ' ' << etat_[ix(i)] << ' ' << etat_[iy(i)] << ' ' << etat_[ivx(i)] << ' ' << etat_[ivy(i)];
                }
                fichier_sortie_ << ' ' << energieTotale()
                                << ' ' << quantite_x
                                << ' ' << quantite_y
                                << ' ' << distanceSondeTerre(etat_)
                                << ' ' << distance_min_terre_ - rayons_[indice_terre_]
                                << '\n';
                break;
            }

            if (temps_ + 1e-12 >= prochain_temps_sortie_) {
                double quantite_x = 0.0;
                double quantite_y = 0.0;
                quantiteMouvementTotale(quantite_x, quantite_y);

                fichier_sortie_ << temps_;
                for (size_t i = 0; i < N_; ++i) {
                    fichier_sortie_ << ' ' << etat_[ix(i)] << ' ' << etat_[iy(i)] << ' ' << etat_[ivx(i)] << ' ' << etat_[ivy(i)];
                }
                fichier_sortie_ << ' ' << energieTotale()
                                << ' ' << quantite_x
                                << ' ' << quantite_y
                                << ' ' << distanceSondeTerre(etat_)
                                << ' ' << distance_min_terre_ - rayons_[indice_terre_]
                                << '\n';
                prochain_temps_sortie_ += intervalle_sortie_;
            }
        }

        cout << setprecision(10);
        cout << "Simulation terminee a t = " << temps_ << " s\n";
        cout << "Perigee minimal = " << distance_min_terre_ - rayons_[indice_terre_] << " m\n";
        cout << "Objectif = " << perigee_cible_ << " m\n";
    }

private:
    size_t ix(size_t i) const { return 2 * i; }
    size_t iy(size_t i) const { return 2 * i + 1; }
    size_t ivx(size_t i) const { return 2 * N_ + 2 * i; }
    size_t ivy(size_t i) const { return 2 * N_ + 2 * i + 1; }

    Etat f(const Etat& etat) const
    {
        Etat derivee(0.0, etat.size());

        for (size_t i = 0; i < N_; ++i) {
            derivee[ix(i)] = etat[ivx(i)];
            derivee[iy(i)] = etat[ivy(i)];
        }

        for (size_t i = 0; i < N_; ++i) {
            double acceleration_x = 0.0;
            double acceleration_y = 0.0;

            for (size_t j = 0; j < N_; ++j) {
                if (i == j) {
                    continue;
                }
                const double dx = etat[ix(j)] - etat[ix(i)];
                const double dy = etat[iy(j)] - etat[iy(i)];
                const double distance2 = dx * dx + dy * dy;
                const double distance3 = distance2 * sqrt(distance2) + 1e-30;
                acceleration_x += G * masses_[j] * dx / distance3;
                acceleration_y += G * masses_[j] * dy / distance3;
            }

            if (trainee_active_[i] != 0) {
                ajouterTrainee(etat, i, acceleration_x, acceleration_y);
            }

            derivee[ivx(i)] = acceleration_x;
            derivee[ivy(i)] = acceleration_y;
        }

        return derivee;
    }

    void ajouterTrainee(const Etat& etat, size_t i, double& acceleration_x, double& acceleration_y) const
    {
        const double dx = etat[ix(i)] - etat[ix(indice_atmosphere_)];
        const double dy = etat[iy(i)] - etat[iy(indice_atmosphere_)];
        const double distance = sqrt(dx * dx + dy * dy);
        const double altitude = distance - rayons_[indice_atmosphere_];
        if (altitude <= 0.0) {
            return;
        }

        const double rho = rho0_ * exp(-altitude / lambda_);
        if (rho < 1e-15) {
            return;
        }

        const double vx_relatif = etat[ivx(i)] - etat[ivx(indice_atmosphere_)];
        const double vy_relatif = etat[ivy(i)] - etat[ivy(indice_atmosphere_)];
        const double vitesse_relative = sqrt(vx_relatif * vx_relatif + vy_relatif * vy_relatif);
        const double factor = -0.5 * rho * section_ * Cx_ / masses_[i];
        acceleration_x += factor * vitesse_relative * vx_relatif;
        acceleration_y += factor * vitesse_relative * vy_relatif;
    }

    Etat rk4(const Etat& etat, double dt) const
    {
        const Etat k1 = f(etat);
        const Etat k2 = f(etat + 0.5 * dt * k1);
        const Etat k3 = f(etat + 0.5 * dt * k2);
        const Etat k4 = f(etat + dt * k3);
        return etat + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }

    double pasAdaptatif(double dt_essai)
    {
        double h = clamp(dt_essai, dt_min_, dt_max_);

        while (true) {
            const Etat etat_1_pas = rk4(etat_, h);
            const Etat etat_2_demi_pas = rk4(rk4(etat_, 0.5 * h), 0.5 * h);

            double erreur = 0.0;
            for (size_t k = 0; k < etat_.size(); ++k) {
                erreur = max(erreur, abs(etat_2_demi_pas[k] - etat_1_pas[k]) / max(1.0, abs(etat_2_demi_pas[k])));
            }

            if (erreur <= precision_ || h <= dt_min_) {
                etat_ = etat_2_demi_pas;
                const double facteur = (erreur > 0.0) ? 0.9 * pow(precision_ / erreur, 0.2) : 2.0;
                dt_ = clamp(h * clamp(facteur, 0.2, 2.0), dt_min_, dt_max_);
                return h;
            }

            h = max(dt_min_, 0.5 * h);
        }
    }

    double distanceSondeTerre(const Etat& etat) const
    {
        const double dx = etat[ix(indice_sonde_)] - etat[ix(indice_terre_)];
        const double dy = etat[iy(indice_sonde_)] - etat[iy(indice_terre_)];
        return sqrt(dx * dx + dy * dy);
    }

    bool collisionDetectee() const
    {
        for (size_t i = 0; i < N_; ++i) {
            for (size_t j = i + 1; j < N_; ++j) {
                const double dx = etat_[ix(j)] - etat_[ix(i)];
                const double dy = etat_[iy(j)] - etat_[iy(i)];
                const double distance = sqrt(dx * dx + dy * dy);
                if (distance <= rayons_[i] + rayons_[j]) {
                    cout << "Collision detectee entre les corps " << i << " et " << j << '\n';
                    return true;
                }
            }
        }
        return false;
    }

    double energieTotale() const
    {
        double energie = 0.0;

        for (size_t i = 0; i < N_; ++i) {
            const double vitesse2 = etat_[ivx(i)] * etat_[ivx(i)] + etat_[ivy(i)] * etat_[ivy(i)];
            energie += 0.5 * masses_[i] * vitesse2;
        }

        for (size_t i = 0; i < N_; ++i) {
            for (size_t j = i + 1; j < N_; ++j) {
                const double dx = etat_[ix(j)] - etat_[ix(i)];
                const double dy = etat_[iy(j)] - etat_[iy(i)];
                const double distance = sqrt(dx * dx + dy * dy);
                energie -= G * masses_[i] * masses_[j] / distance;
            }
        }

        return energie;
    }

    void quantiteMouvementTotale(double& quantite_x, double& quantite_y) const
    {
        quantite_x = 0.0;
        quantite_y = 0.0;
        for (size_t i = 0; i < N_; ++i) {
            quantite_x += masses_[i] * etat_[ivx(i)];
            quantite_y += masses_[i] * etat_[ivy(i)];
        }
    }

};
}

int main()
{
    try {
        Exercice4 exercice;
        exercice.executer();
    } catch (const exception& e) {
        cerr << "Erreur : " << e.what() << '\n';
        return 1;
    }

    return 0;
}
