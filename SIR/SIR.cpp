#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>

//paramètres de l'infection
double beta = 1;           //taux de transmission
double gam = .1;           //taux de récupération
double N = 10000000;       //nombre d'individus dans la population

void exportToCSV(const std::vector<double>& v1,
                const std::vector<double>& v2,
                const std::vector<double>& v3,
                const std::vector<double>& v4,
                const std::vector<double>& v5,
                const std::string& filename) {
    /*
    Cette fonction prend cinq vecteurs et les exporte dans un fichier CSV.
    Ce fichier est utilisé pour produire les graphiques dans Python.
    */
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Erreur dans l'ouverture du fichier : " << filename << std::endl;
        return;
    }
    for (int i = 0; i < v1.size(); i++) {
        file << v1[i] << "," << v2[i] << "," << v3[i] << "," << v4[i] << "," << v5[i] << std::endl;
    }
    file.close();
}

void print_vector(const std::vector<double>& v) {
    /*
    Cette fonction affiche les vecteurs. Elle permet de simplifier le
    code plus bas afin de le rendre plus lisible. 
    */
    std::cout << "(";
    for (int i = 0; i < v.size(); i++) {
        if (i != v.size() - 1) {
            std::cout << v[i] << ",";
        } else{
            std::cout << v[i];
        }
    }
    std::cout << ")" << std::endl;
}

std::vector<double> scalar_vector_mult(const std::vector<double>& v, double s) {
    /*
    Cette fonction multiplie un vecteur par un scalaire et retourne
    le résultat. Il y a sûrement une façon plus élégante de le faire.
    */
    std::vector<double> result;
    for (int i = 0; i < v.size(); i++){
        result.push_back(v[i]*s);
    }
    return result;
}

std::vector<double> vector_add(const std::vector<double>& v1,const std::vector<double>& v2) {
    /*
    Cette fonction additionne deux vecteurs et retourne le résultat.
    Il y a sûrement une façon plus élégante de le faire, mais la
    lisibilité du code est améliorée pour le moment. 
    */
    std::vector<double> result;
    for (int i = 0; i < v1.size(); i++){
        result.push_back(v1[i] + v2[i]);
    }
    return result;
}

std::vector<double> f(const std::vector<double>& v, double t) {
    /*
    Y' = f(Y,t)     Le système d'EDOs est représenté par cette fonction.
    */
    double f0 = (-beta*v[0]*v[1]/N);
    double f1 = (beta*v[0]*v[1]/N) - gam*v[1];
    double f2 = gam*v[1];
    std::vector<double> result = {f0,f1,f2};
    return result;
}

int main() {
    //conditions initiales
    double I0 = 1;           //nombre de personnes infectées initialement
    double S0 = N - I0;      //nombre de personnes susceptibles initialement
    double R0 = 0;           //nombre de personnes ayant initialement récupérées

    std::cout << "Nombre de personnes dans la population : " << N << std::endl;
    std::cout << "Nombre initial d'infectés : " << I0 << std::endl;
    std::cout << "Nombre initial de susceptibles : " << S0 << std::endl;
    std::cout << "Nombre initial de personnes ayant récupérées de l'infection " << R0 << std::endl;

    //vecteurs vides des compartiments
    std::vector<double> S;               //nombre de personnes susceptibles
    std::vector<double> I;               //nombre de personnes infectées
    std::vector<double> R;               //nombre de personnes ayant récupérées
    std::vector<double> Y = {S0,I0,R0};  //vecteur d'état
    
    //ajout des conditions initiales
    S.push_back(S0);
    I.push_back(I0);
    R.push_back(R0);
    std::cout << "Vecteur initial (S,I,R) = ";
    print_vector(Y);

    //paramètres pour Runge-Kutta
    std::vector<double> k1(3);
    std::vector<double> k2(3);
    std::vector<double> k3(3);
    std::vector<double> k4(3);

    //paramètres temporels
    double dt = .1;                      //incréments de temps
    int Tmin = 0;                       //début de la simulation
    int Tmax = 100;                     //fin de la simulation
    int nb_iter = (Tmax - Tmin)/dt;     //nombre d'itérations
    std::vector<double> T;               //vecteurs des instants mesurés
    T.push_back(Tmin);                  //initialisation au premier instant

    //vecteur test (on vérifie que la population reste constante)
    std::vector<double> population_size;
    population_size.push_back(N);

    for (int n = 1; n <= nb_iter; n++) {
        //calcule des pentes k1, k2, k3 et k4 de la méthode RK4
        k1 = scalar_vector_mult(f(Y,T.back()), dt);
        std::vector<double> rk_incr1 = vector_add(Y,scalar_vector_mult(k1,.5));
        k2 = scalar_vector_mult(f(rk_incr1,T.back()+.5*dt),dt);
        std::vector<double> rk_incr2 = vector_add(Y,scalar_vector_mult(k2,.5));
        k3 = scalar_vector_mult(f(rk_incr2,T.back()+.5*dt),dt);
        std::vector<double> rk_incr3 = vector_add(Y,k3);
        k4 = scalar_vector_mult(f(rk_incr3,T.back()+dt),dt);
        std::vector<double> dY = vector_add(k1,scalar_vector_mult(k2,2));
        dY = vector_add(dY,scalar_vector_mult(k3,2));
        //calcule du nouveau vecteur d'état
        dY = vector_add(dY,k4);
        dY = scalar_vector_mult(dY,1.0/6.0);
        Y = vector_add(Y,dY);

        //ajout des valeurs calculées dans l'historique des compartiments
        S.push_back(Y[0]);
        I.push_back(Y[1]);
        R.push_back(Y[2]);
        T.push_back(T.back() + dt);
        population_size.push_back(Y[0]+Y[1]+Y[2]);
    }
    std::cout << "Vecteur final (S,I,R) = ";
    print_vector(Y);

    exportToCSV(T,S,I,R,population_size,"SIR.csv");
    int result = system("py plot_SIR_data.py SIR.csv");
    return 0;
}
