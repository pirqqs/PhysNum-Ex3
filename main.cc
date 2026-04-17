#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <valarray>
#include <vector>
#include <string>
#include <algorithm>
#include "./ConfigFile.h"

using namespace std;
using Etat = valarray<double>;

class Engine
{
private:
    static constexpr double PI = 3.1415926535897932384626433832795028841971;

    // Paramètres physiques
    double G;
    double mT, mL, mA;
    double rT, rL, rA;
    double rho0, lambda, Cx, section;

    // Paramètres numériques
    double t;
    double tf;
    double dt;
    double dt_min;
    double dt_max;
    double epsilon;

    bool adaptive;
    unsigned int sampling;
    unsigned int last;

    // Structure du système
    size_t Nbody;
    Etat y_;
    vector<double> masses_;
    vector<double> rayons_;
    vector<int> trainee_active_;

    // Indices utiles
    size_t indice_terre;
    size_t indice_lune;
    size_t indice_sonde;
    size_t indice_atmosphere;

    // Conditions initiales
    double xT0, yT0, vxT0, vyT0;
    double xL0, yL0, vxL0, vyL0;
    double xA0, yA0, vxA0, vyA0;

    // Paramètres d'analyse
    double perigee_cible;
    double distance_min_terre;
    double dt_last_used;

    // Sortie
    ofstream* outputFile;

private:
    // Fonctions d'indices
    size_t ix(size_t i) const { return 2 * i; }
    size_t iy(size_t i) const { return 2 * i + 1; }
    size_t ivx(size_t i) const { return 2 * Nbody + 2 * i; }
    size_t ivy(size_t i) const { return 2 * Nbody + 2 * i + 1; }

    // Densité atmosphérique
    double densiteAtmosphere(double altitude) const
    {
        if (altitude <= 0.0) return rho0;
        return rho0 * exp(-altitude / lambda);
    }

    // Accélération gravitationnelle
    // + traînée éventuelle
    void accelerationCorps(const Etat& y, size_t i, double& ax, double& ay) const
    {
        ax = 0.0;
        ay = 0.0;

        // Gravitation due aux autres corps
        for (size_t j = 0; j < Nbody; ++j)
        {
            if (j == i) continue;

            const double dx = y[ix(j)] - y[ix(i)];
            const double dy = y[iy(j)] - y[iy(i)];
            const double r2 = dx * dx + dy * dy;
            const double r = sqrt(r2);
            const double r3 = r2 * r + 1e-30;

            ax += G * masses_[j] * dx / r3;
            ay += G * masses_[j] * dy / r3;
        }

        // Traînée atmosphérique éventuelle
        if (i < trainee_active_.size() && trainee_active_[i] != 0)
        {
            ajouterTrainee(y, i, ax, ay);
        }
    }

    void ajouterTrainee(const Etat& y, size_t i, double& ax, double& ay) const
    {
        const double dx = y[ix(i)] - y[ix(indice_atmosphere)];
        const double dy = y[iy(i)] - y[iy(indice_atmosphere)];
        const double r = sqrt(dx * dx + dy * dy);
        const double altitude = r - rayons_[indice_atmosphere];

        if (altitude < 0.0) return;

        const double rho = densiteAtmosphere(altitude);
        if (rho < 1e-15) return;

        const double vx_rel = y[ivx(i)] - y[ivx(indice_atmosphere)];
        const double vy_rel = y[ivy(i)] - y[ivy(indice_atmosphere)];
        const double v_rel = sqrt(vx_rel * vx_rel + vy_rel * vy_rel);

        if (v_rel < 1e-30) return;

        const double facteur = -0.5 * rho * section * Cx / masses_[i];

        ax += facteur * v_rel * vx_rel;
        ay += facteur * v_rel * vy_rel;
    }

    // =========================
    // Fonction f(y) = dy/dt
    // =========================
    Etat f(const Etat& y) const
    {
        Etat dydt(0.0, y.size());

        // Dérivées des positions = vitesses
        for (size_t i = 0; i < Nbody; ++i)
        {
            dydt[ix(i)] = y[ivx(i)];
            dydt[iy(i)] = y[ivy(i)];
        }

        // Dérivées des vitesses = accélérations
        for (size_t i = 0; i < Nbody; ++i)
        {
            double ax = 0.0;
            double ay = 0.0;
            accelerationCorps(y, i, ax, ay);

            dydt[ivx(i)] = ax;
            dydt[ivy(i)] = ay;
        }

        return dydt;
    }

    // =========================
    // RK4
    // =========================
    Etat rk4Step(const Etat& y, double h) const
    {
        const Etat k1 = f(y);
        const Etat k2 = f(y + 0.5 * h * k1);
        const Etat k3 = f(y + 0.5 * h * k2);
        const Etat k4 = f(y + h * k3);

        return y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }

    // =========================
    // Pas adaptatif
    // =========================
    double adaptiveStep(double h_try)
    {
        double h = std::clamp(h_try, dt_min, dt_max);

        while (true)
        {
            const Etat y_big = rk4Step(y_, h);
            const Etat y_half = rk4Step(rk4Step(y_, 0.5 * h), 0.5 * h);

            double err = 0.0;
            for (size_t k = 0; k < y_.size(); ++k)
            {
                const double scale = max(1.0, abs(y_half[k]));
                err = max(err, abs(y_half[k] - y_big[k]) / scale);
            }

            if (err <= epsilon || h <= dt_min)
            {
                y_ = y_half;

                double facteur = 2.0;
                if (err > 1e-30)
                {
                    facteur = 0.9 * pow(epsilon / err, 0.2);
                }
                facteur = clamp(facteur, 0.2, 2.0);

                dt = clamp(h * facteur, dt_min, dt_max);
                return h;
            }

            h = max(dt_min, 0.5 * h);
        }
    }

    // =========================
    // Un pas de temps
    // =========================
    void step()
    {
        if (adaptive)
        {
            dt_last_used = adaptiveStep(dt);
            t += dt_last_used;
        }
        else
        {
            y_ = rk4Step(y_, dt);
            t += dt;
            dt_last_used = dt;
        }

        distance_min_terre = min(distance_min_terre, distanceSondeTerre(y_));
    }

    // =========================
    // Diagnostics
    // =========================
    double energieTotale(const Etat& y) const
    {
        double E = 0.0;

        // énergie cinétique
        for (size_t i = 0; i < Nbody; ++i)
        {
            const double v2 = y[ivx(i)] * y[ivx(i)] + y[ivy(i)] * y[ivy(i)];
            E += 0.5 * masses_[i] * v2;
        }

        // énergie potentielle gravitationnelle
        for (size_t i = 0; i < Nbody; ++i)
        {
            for (size_t j = i + 1; j < Nbody; ++j)
            {
                const double dx = y[ix(j)] - y[ix(i)];
                const double dy = y[iy(j)] - y[iy(i)];
                const double r = sqrt(dx * dx + dy * dy) + 1e-30;

                E -= G * masses_[i] * masses_[j] / r;
            }
        }

        return E;
    }

    void quantiteMouvementTotale(const Etat& y, double& px, double& py) const
    {
        px = 0.0;
        py = 0.0;

        for (size_t i = 0; i < Nbody; ++i)
        {
            px += masses_[i] * y[ivx(i)];
            py += masses_[i] * y[ivy(i)];
        }
    }

    double distanceSondeTerre(const Etat& y) const
    {
        const double dx = y[ix(indice_sonde)] - y[ix(indice_terre)];
        const double dy = y[iy(indice_sonde)] - y[iy(indice_terre)];
        return sqrt(dx * dx + dy * dy);
    }

    double altitudeSonde(const Etat& y) const
    {
        return distanceSondeTerre(y) - rayons_[indice_terre];
    }

    double vitesseSonde(const Etat& y) const
    {
        const double vx = y[ivx(indice_sonde)];
        const double vy = y[ivy(indice_sonde)];
        return sqrt(vx * vx + vy * vy);
    }

    void accelerationSonde(const Etat& y, double& ax, double& ay) const
    {
        accelerationCorps(y, indice_sonde, ax, ay);
    }

    double accelerationNormeSonde(const Etat& y) const
    {
        double ax = 0.0, ay = 0.0;
        accelerationSonde(y, ax, ay);
        return sqrt(ax * ax + ay * ay);
    }

    double puissanceTraineeSonde(const Etat& y) const
    {
        if (!(indice_sonde < trainee_active_.size()) || trainee_active_[indice_sonde] == 0)
            return 0.0;

        const double dx = y[ix(indice_sonde)] - y[ix(indice_atmosphere)];
        const double dy = y[iy(indice_sonde)] - y[iy(indice_atmosphere)];
        const double r = sqrt(dx * dx + dy * dy);
        const double altitude = r - rayons_[indice_atmosphere];

        if (altitude < 0.0) return 0.0;

        const double rho = densiteAtmosphere(altitude);
        if (rho < 1e-15) return 0.0;

        const double vx_rel = y[ivx(indice_sonde)] - y[ivx(indice_atmosphere)];
        const double vy_rel = y[ivy(indice_sonde)] - y[ivy(indice_atmosphere)];
        const double v_rel = sqrt(vx_rel * vx_rel + vy_rel * vy_rel);

        if (v_rel < 1e-30) return 0.0;

        // Force de traînée
        const double facteur = -0.5 * rho * section * Cx;
        const double Fx = facteur * v_rel * vx_rel;
        const double Fy = facteur * v_rel * vy_rel;

        // Puissance = F · v
        const double vx = y[ivx(indice_sonde)];
        const double vy = y[ivy(indice_sonde)];

        return Fx * vx + Fy * vy;
    }

    bool collisionDetectee(const Etat& y) const
    {
        for (size_t i = 0; i < Nbody; ++i)
        {
            for (size_t j = i + 1; j < Nbody; ++j)
            {
                const double dx = y[ix(j)] - y[ix(i)];
                const double dy = y[iy(j)] - y[iy(i)];
                const double d = sqrt(dx * dx + dy * dy);

                if (d <= rayons_[i] + rayons_[j])
                {
                    cout << "Collision detectee entre les corps " << i
                         << " et " << j << endl;
                    return true;
                }
            }
        }
        return false;
    }

    // =========================
    // Sortie
    // =========================
    void printHeader()
    {
        *outputFile << "# t";
        for (size_t i = 0; i < Nbody; ++i)
        {
            *outputFile << " x" << i << " y" << i
                        << " vx" << i << " vy" << i;
        }

        *outputFile << " E Px Py";
        *outputFile << " r_probe_earth altitude_probe";
        *outputFile << " perigee dt_used";
        *outputFile << " a_probe P_drag";
        *outputFile << endl;
    }

    void printOut(bool write)
    {
        if ((!write && last >= sampling) || (write && last != 1))
        {
            double px = 0.0, py = 0.0;
            quantiteMouvementTotale(y_, px, py);

            const double r_probe_earth = distanceSondeTerre(y_);
            const double altitude_probe = altitudeSonde(y_);
            const double perigee = distance_min_terre - rayons_[indice_terre];
            const double a_probe = accelerationNormeSonde(y_);
            const double P_drag = puissanceTraineeSonde(y_);

            *outputFile << t;
            for (size_t i = 0; i < Nbody; ++i)
            {
                *outputFile << " " << y_[ix(i)]
                            << " " << y_[iy(i)]
                            << " " << y_[ivx(i)]
                            << " " << y_[ivy(i)];
            }

            *outputFile << " " << energieTotale(y_)
                        << " " << px
                        << " " << py
                        << " " << r_probe_earth
                        << " " << altitude_probe
                        << " " << perigee
                        << " " << dt_last_used
                        << " " << a_probe
                        << " " << P_drag
                        << endl;

            last = 1;
        }
        else
        {
            last++;
        }
    }

public:
    Engine(ConfigFile configFile)
    {
        // -------------------------
        // Valeurs par défaut
        // -------------------------
        G = 6.674e-11;

        mT = 5.972e24;
        mL = 7.3477e22;
        mA = 8500.0;

        rT = 6378.1e3;
        rL = 1737.4e3;
        rA = 0.5 * 5.02;

        rho0 = 0.0;
        lambda = 7238.2;
        Cx = 0.3;
        section = 0.25 * PI * 5.02 * 5.02;

        t = 0.0;
        tf = 2.0 * 24.0 * 3600.0;
        dt = 10.0;
        dt_min = 1e-3;
        dt_max = 300.0;
        epsilon = 1e-6;
        adaptive = false;
        sampling = 1;
        last = 0;

        Nbody = 3;

        indice_terre = 0;
        indice_lune = 1;
        indice_sonde = 2;
        indice_atmosphere = 0;

        xT0 = 0.0;
        yT0 = 0.0;
        vxT0 = 0.0;
        vyT0 = 0.0;

        xL0 = 384748e3;
        yL0 = 0.0;
        vxL0 = 0.0;
        vyL0 = 1022.0;

        xA0 = 314159e3;
        yA0 = 0.0;
        vxA0 = 0.0;
        vyA0 = 1200.0;

        perigee_cible = 10000.0;
        distance_min_terre = 1e300;
        dt_last_used = dt;

        // -------------------------
        // Lecture configuration
        // -------------------------
        G        = configFile.get<double>("G", G);

        mT       = configFile.get<double>("mT", mT);
        mL       = configFile.get<double>("mL", mL);
        mA       = configFile.get<double>("mA", mA);

        rT       = configFile.get<double>("rT", rT);
        rL       = configFile.get<double>("rL", rL);
        rA       = configFile.get<double>("rA", rA);

        rho0     = configFile.get<double>("rho0", rho0);
        lambda   = configFile.get<double>("lambda", lambda);
        Cx       = configFile.get<double>("Cx", Cx);
        section  = configFile.get<double>("S", section);

        tf       = configFile.get<double>("tfin", tf);
        dt       = configFile.get<double>("dt", dt);
        dt_min   = configFile.get<double>("dt_min", dt_min);
        dt_max   = configFile.get<double>("dt_max", dt_max);
        epsilon  = configFile.get<double>("epsilon", epsilon);
        adaptive = (configFile.get<int>("adaptive", adaptive ? 1 : 0) != 0);

        sampling = configFile.get<unsigned int>("nsampling", sampling);

        Nbody    = static_cast<size_t>(configFile.get<int>("Nbody", static_cast<int>(Nbody)));

        indice_terre      = static_cast<size_t>(configFile.get<int>("indice_terre", static_cast<int>(indice_terre)));
        indice_lune       = static_cast<size_t>(configFile.get<int>("indice_lune", static_cast<int>(indice_lune)));
        indice_sonde      = static_cast<size_t>(configFile.get<int>("indice_sonde", static_cast<int>(indice_sonde)));
        indice_atmosphere = static_cast<size_t>(configFile.get<int>("indice_atmosphere", static_cast<int>(indice_atmosphere)));

        xT0 = configFile.get<double>("xT0", xT0);
        yT0 = configFile.get<double>("yT0", yT0);
        vxT0 = configFile.get<double>("vxT0", vxT0);
        vyT0 = configFile.get<double>("vyT0", vyT0);

        xL0 = configFile.get<double>("xL0", xL0);
        yL0 = configFile.get<double>("yL0", yL0);
        vxL0 = configFile.get<double>("vxL0", vxL0);
        vyL0 = configFile.get<double>("vyL0", vyL0);

        xA0 = configFile.get<double>("xA0", xA0);
        yA0 = configFile.get<double>("yA0", yA0);
        vxA0 = configFile.get<double>("vxA0", vxA0);
        vyA0 = configFile.get<double>("vyA0", vyA0);

        perigee_cible = configFile.get<double>("h", perigee_cible);

        const int dragT = configFile.get<int>("drag_T", 0);
        const int dragL = configFile.get<int>("drag_L", 0);
        const int dragA = configFile.get<int>("drag_A", 1);

        trainee_active_ = {dragT, dragL, dragA};

        masses_ = {mT, mL, mA};
        rayons_ = {rT, rL, rA};

        if (Nbody < 1 || Nbody > 3)
            throw runtime_error("Nbody doit etre entre 1 et 3.");

        if (masses_.size() < Nbody || rayons_.size() < Nbody)
            throw runtime_error("Erreur de dimension des masses/rayons.");

        y_ = Etat(0.0, 4 * Nbody);

        vector<double> x_init = {xT0, xL0, xA0};
        vector<double> y_init = {yT0, yL0, yA0};
        vector<double> vx_init = {vxT0, vxL0, vxA0};
        vector<double> vy_init = {vyT0, vyL0, vyA0};

        for (size_t i = 0; i < Nbody; ++i)
        {
            y_[ix(i)] = x_init[i];
            y_[iy(i)] = y_init[i];
            y_[ivx(i)] = vx_init[i];
            y_[ivy(i)] = vy_init[i];
        }

        distance_min_terre = distanceSondeTerre(y_);

        const string outputName = configFile.get<string>("output", "trajectory.dat");
        outputFile = new ofstream(outputName.c_str());
        outputFile->precision(15);

        if (!(*outputFile))
            throw runtime_error("Impossible d'ouvrir le fichier de sortie.");
    }

    virtual ~Engine()
    {
        outputFile->close();
        delete outputFile;
    }

    void run()
    {
        t = 0.0;
        last = 0;
        dt_last_used = dt;

        printHeader();
        printOut(true);

        while (t < tf)
        {
            double t_before = t;
            step();

            if (collisionDetectee(y_))
            {
                printOut(true);
                break;
            }

            printOut(false);

            if (t <= t_before)
                throw runtime_error("Le temps n'avance plus. Verifie dt/dt_min.");
        }

        printOut(true);

        cout << setprecision(10);
        cout << "Fin de la simulation." << endl;
        cout << "Temps final atteint : " << t << " s" << endl;
        cout << "Perigee minimal : " << distance_min_terre - rayons_[indice_terre] << " m" << endl;
        cout << "Objectif demande : " << perigee_cible << " m" << endl;
    }
};

int main(int argc, char* argv[])
{
    ConfigFile configFile;

    for (int i = 1; i < argc; ++i)
        configFile.process(argv[i]);

    Engine* engine = nullptr;

    try
    {
        engine = new Engine(configFile);
        engine->run();
        delete engine;
        cout << "Simulation terminee." << endl;
        return 0;
    }
    catch (const exception& e)
    {
        delete engine;
        cerr << "Erreur : " << e.what() << endl;
        return 1;
    }
}
