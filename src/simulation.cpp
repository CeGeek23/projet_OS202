#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <array>
#include <vector>
#include <numeric>

#include "model.hpp"
#include "display.hpp"

#include <mpi.h>

using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType
{
    double length{1.};
    unsigned discretization{20u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};
};

void analyze_arg(int nargs, char* args[], ParamsType& params)
{
    if (nargs == 0) return;
    std::string key(args[0]);
    if (key == "-l"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour la longueur du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs - 2, &args[2], params);
        return;
    }
    auto pos = key.find("--longueur=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos + 11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs - 1, &args[1], params);
        return;
    }

    if (key == "-n"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour le nombre de cases par direction pour la discrétisation du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs - 2, &args[2], params);
        return;
    }
    pos = key.find("--number_of_cases=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos + 18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs - 1, &args[1], params);
        return;
    }

    if (key == "-w"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la direction du vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values = std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos + 1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs - 2, &args[2], params);
        return;
    }
    pos = key.find("--wind=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos + 7);
        params.wind[0] = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos + 1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs - 1, &args[1], params);
        return;
    }

    if (key == "-s"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la position du foyer initial !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values = std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la position du foyer initial" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos + 1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs - 2, &args[2], params);
        return;
    }
    pos = key.find("--start=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos + 8);
        params.start.column = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos + 1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs - 1, &args[1], params);
        return;
    }
}

ParamsType parse_arguments(int nargs, char* args[])
{
    if (nargs == 0) return {};
    if ((std::string(args[0]) == "--help"s) || (std::string(args[0]) == "-h"))
    {
        std::cout <<
            R"RAW(Usage : simulation [option(s)]
  Lance la simulation d'incendie en prenant en compte les [option(s)].
  Les options sont :
    -l, --longueur=LONGUEUR     Définit la taille LONGUEUR (réel en km) du carré représentant la carte de la végétation.
    -n, --number_of_cases=N     Nombre n de cases par direction pour la discrétisation
    -w, --wind=VX,VY            Définit le vecteur vitesse du vent (pas de vent par défaut).
    -s, --start=COL,ROW         Définit les indices I,J de la case où commence l'incendie (milieu de la carte par défaut)
)RAW";
        exit(EXIT_SUCCESS);
    }
    ParamsType params;
    analyze_arg(nargs, args, params);
    return params;
}

bool check_params(ParamsType& params)
{
    bool flag = true;
    if (params.length <= 0)
    {
        std::cerr << "[ERREUR FATALE] La longueur du terrain doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if (params.discretization <= 0)
    {
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit être positif et non nul !" << std::endl;
        flag = false;
    }

    if ((params.start.row >= params.discretization) || (params.start.column >= params.discretization))
    {
        std::cerr << "[ERREUR FATALE] Mauvais indices pour la position initiale du foyer" << std::endl;
        flag = false;
    }

    return flag;
}

void display_params(ParamsType const& params)
{
    std::cout << "Parametres définis pour la simulation : \n"
        << "\tTaille du terrain : " << params.length << std::endl
        << "\tNombre de cellules par direction : " << params.discretization << std::endl
        << "\tVecteur vitesse : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
        << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl;
}

void compute_load_distribution(int discretization, int globSize, std::vector<int>& recvcounts, std::vector<int>& displs)
{
    int rows_per_proc = discretization / (globSize - 1);
    int extra_rows = discretization % (globSize - 1);

    for (int i = 0; i < globSize - 1; ++i) {
        int num_rows = rows_per_proc + (i < extra_rows ? 1 : 0);
        recvcounts[i] = num_rows * discretization;
        displs[i] = (i > 0 ? displs[i - 1] + recvcounts[i - 1] : 0);
    }
}

int main(int nargs, char* args[])
{
    // On initialise MPI via MPI_Init
    MPI_Init(&nargs, &args);
    int globRank, globSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &globRank);
    MPI_Comm_size(MPI_COMM_WORLD, &globSize);

    if (globSize < 2)
    {
        std::cerr << "Cette simulation nécessite au moins 2 processus MPI" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // Abort permet de terminer MPI proprement en cas d'erreur
        return EXIT_FAILURE;
    }

    auto params = parse_arguments(nargs - 1, &args[1]);

    if (!check_params(params)) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // on termine MPI si les paramètres sont incorrects
        return EXIT_FAILURE;
    }

    // Taille totale de l'espace de simulation
    const int gridSize = params.discretization * params.discretization;

    std::vector<int> recvcounts(globSize - 1);
    std::vector<int> displs(globSize - 1);
    compute_load_distribution(params.discretization, globSize, recvcounts, displs);

    if (globRank == 0)
    {
        display_params(params);
        // On initialise l'affichage via SDL grâce à Displayer::init_instance
        auto displayer = Displayer::init_instance(params.discretization, params.discretization);

        // On alloue les buffers pour l'état de la simulation
        std::vector<std::uint8_t> vegetation(gridSize, 0);
        std::vector<std::uint8_t> fire(gridSize, 0);

        // Réception non bloquantes depuis les processus de calcul
        std::vector<MPI_Request> recvReqVeget(globSize - 1);
        std::vector<MPI_Request> recvReqFire(globSize - 1);

        for (int i = 1; i < globSize; ++i) {
            MPI_Irecv(vegetation.data() + displs[i - 1], recvcounts[i - 1], MPI_UINT8_T, i, 100, MPI_COMM_WORLD, &recvReqVeget[i - 1]);
            MPI_Irecv(fire.data() + displs[i - 1], recvcounts[i - 1], MPI_UINT8_T, i, 101, MPI_COMM_WORLD, &recvReqFire[i - 1]);
        }

        bool running = true;

        // variables pour mesurer les performances
        auto total_start_time = std::chrono::high_resolution_clock::now();
        auto total_end_time = total_start_time;
        auto display_start_time = total_start_time;
        auto display_end_time = display_start_time;

        while (running)
        {
            int flagVeget = 0, flagFire = 0;
            MPI_Status status;
            for (int i = 1; i < globSize; ++i) {
                MPI_Test(&recvReqVeget[i - 1], &flagVeget, &status);
                if (flagVeget)
                {
                    // Repost pour la prochaine mise à jour
                    MPI_Irecv(vegetation.data() + displs[i - 1], recvcounts[i - 1], MPI_UINT8_T, i, 100, MPI_COMM_WORLD, &recvReqVeget[i - 1]);
                }

                MPI_Test(&recvReqFire[i - 1], &flagFire, &status);
                if (flagFire)
                {
                    MPI_Irecv(fire.data() + displs[i - 1], recvcounts[i - 1], MPI_UINT8_T, i, 101, MPI_COMM_WORLD, &recvReqFire[i - 1]);
                }
            }

            // On mesure le temps d'affichage
            display_start_time = std::chrono::high_resolution_clock::now();
            // On met l'affichage à jour avec l'état reçu
            displayer->update(vegetation, fire);
            display_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> display_duration = display_end_time - display_start_time;
            std::cout << "Temps d'affichage : " << display_duration.count() << " seconds" << std::endl;

            SDL_Event event;
            while (SDL_PollEvent(&event))
            {
                if (event.type == SDL_QUIT)
                {
                    running = false;
                    int term = -1;
                    for (int i = 1; i < globSize; ++i) {
                        MPI_Send(&term, 1, MPI_INT, i, 200, MPI_COMM_WORLD); // On envoie un message de terminaison à tous les processus
                    }
                    break;
                }
            }
            //std::this_thread::sleep_for(0.1s);
        }

        total_end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> total_duration = total_end_time - total_start_time;
        std::cout << "Temps total de simulation : " << total_duration.count() << " seconds" << std::endl;
    }

    else
    {
        auto simu = Model(params.length, params.discretization, params.wind,
            params.start);
        double local_time = 0.0;
        int local_iter = 0;
        bool computing = true;

        // On se prépare à une réception non bloquante pour détecter un signal de terminaison envoyé par le processus 0
        int term_signal = 0;
        MPI_Request termReq;
        MPI_Irecv(&term_signal, 1, MPI_INT, 0, 200, MPI_COMM_WORLD, &termReq);

        // Diviser la grille entre les processus
        int start_row = (globRank - 1) * recvcounts[globRank - 1] / params.discretization;
        int num_rows = recvcounts[globRank - 1] / params.discretization;
        int local_grid_size = num_rows * params.discretization;

        // Allocation des buffers locaux pour les cartes de végétation et de feu
        std::vector<std::uint8_t> local_vegetation(local_grid_size, 0);
        std::vector<std::uint8_t> local_fire(local_grid_size, 0);

        while (simu.update() & computing)
        {
            auto iter_start = std::chrono::steady_clock::now();
            // Calcul de la prochaine itération
            simu.update();
            auto iter_end = std::chrono::steady_clock::now();
            double dt = std::chrono::duration<double>(iter_end - iter_start).count();
            local_time += dt;
            local_iter++;

            if ((simu.time_step() & 31) == 0) {
                std::cout << "Time step " << simu.time_step() << "\n===============" << std::endl;
            }

            // Récupération des cartes de végétation et de feu
            auto vegetation = simu.vegetal_map(); // Assuming vegetal_map() returns a vector
            auto fire = simu.fire_map(); // Assuming fire_map() returns a vector

            // Copie des données locales dans les buffers locaux
            std::copy(vegetation.begin() + start_row * params.discretization, vegetation.begin() + (start_row + num_rows) * params.discretization, local_vegetation.begin());
            std::copy(fire.begin() + start_row * params.discretization, fire.begin() + (start_row + num_rows) * params.discretization, local_fire.begin());

            // vérification non bloquante du signal de terminaison
            int flag = 0;
            MPI_Status status;
            MPI_Test(&termReq, &flag, &status);
            if (flag && term_signal == -1)
            {
                computing = false;
                break;
            }

            // Envoi non bloquant de l'état de simulation vers l'affichage
            MPI_Request sendReq;
            MPI_Isend(local_vegetation.data(), local_grid_size, MPI_UINT8_T, 0, 100, MPI_COMM_WORLD, &sendReq);
            MPI_Request_free(&sendReq); // On libère la requête
            MPI_Request sendReqFire;
            MPI_Isend(local_fire.data(), local_grid_size, MPI_UINT8_T, 0, 101, MPI_COMM_WORLD, &sendReqFire);
            MPI_Request_free(&sendReqFire); // On libère la requête

            //std::this_thread::sleep_for(0.1s);
        }

        // Réduction des temps de calculs locaux pour obtenir une statistique globale
        double global_total_time = 0.0;
        int global_total_iter = 0;
        MPI_Reduce(&local_time, &global_total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_iter, &global_total_iter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (globRank == 1 && global_total_iter > 0)
        {
            std::cout << "Temps de calcul moyen par itération : "
                << global_total_time / global_total_iter << "s" << std::endl;
        }

    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}