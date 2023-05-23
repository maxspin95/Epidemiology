#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>

//paramètres de l'infection
double beta = 1;           //taux de transmission
double gamma_R = .2;       //inverse de la durée moyenne d'une infection non sévère
double gamma_H = .1;       //inverse de la durée moyenne d'une hospitalisation pour une infection sévère
double eta = .2;           //inverse de la durée moyenne entre l'acquisition de l'infection sévère et l'hospitalisation
double N = 10000000;       //nombre d'individus dans la population
double x = .2;             //proportion des infections qui vont nécessiter une hospitalisation

void exportToCSV(const std::vector<double>& v1,
                const std::vector<double>& v2,
                const std::vector<double>& v3,
                const std::vector<double>& v4,
                const std::vector<double>& v5,
                const std::vector<double>& v6,
                const std::vector<double>& v7,
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
        file << v1[i] << "," << v2[i] << "," << v3[i] << "," << v4[i] << "," << v5[i] << "," << v6[i] << "," << v7[i] << std::endl;
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
    le résultat. Il y a sûrement une faÃ§on plus élégante de le faire.
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
    double lambda = beta*v[0]*(v[1] + v[2])/N;
    double f0 = -lambda;
    double f1 = (1-x)*lambda - gamma_R*v[1];
    double f2 = x*lambda - eta*v[2];
    double f3 = eta*v[2] - gamma_H*v[3];
    double f4 = gamma_R*v[1] + gamma_H*v[3];
    std::vector<double> result = {f0,f1,f2,f3,f4};
    return result;
}

int main() {
    //conditions initiales
    double I_R0 = 1;                     //nombre initial de personnes avec l'infection non sévère
    double I_H0 = 0;                     //nombre initial de personnes avec l'infection sévère
    double S0 = N - (I_R0 + I_H0);       //nombre initial de personnes susceptibles
    double H0 = 0;                       //nombre initial de personnes hospitalisées
    double R0 = 0;                       //nombre initial de personnes ayant récupérées

    std::cout << "Nombre de personnes dans la population : " << N << std::endl;
    std::cout << "Nombre initial de personnes avec l'infection non sévère : " << I_R0 << std::endl;
    std::cout << "Nombre initial de personnes susceptibles : " << S0 << std::endl;
    std::cout << "Nombre initial de personnes avec l'infection sévère : " << I_H0 << std::endl;
    std::cout << "Nombre initial de personnes hospitalisées : " << H0 << std::endl;
    std::cout << "Nombre initial de personnes ayant recuperees : " << R0 << std::endl;

    //vecteurs vides des compartiments
    std::vector<double> S;                           //nombre de personnes susceptibles
    std::vector<double> I_R;                         //nombre de personnes avec l'infection non sévère
    std::vector<double> I_H;                         //nombre de personnes avec l'infection sévère
    std::vector<double> H;                           //nombre de personne hospitalisée
    std::vector<double> R;                           //nombre de personnes ayant récupérées
    std::vector<double> Y = {S0,I_R0,I_H0,H0,R0};    //vecteur d'état

    //ajout des conditions initiales
    S.push_back(S0);
    I_R.push_back(I_R0);
    I_H.push_back(I_H0);
    H.push_back(H0);
    R.push_back(R0);
    std::cout << "Vecteur initial (S,I_R,I_H,H,R) = ";
    print_vector(Y);

    //paramètres pour Runge-Kutta
    std::vector<double> k1(5);
    std::vector<double> k2(5);
    std::vector<double> k3(5);
    std::vector<double> k4(5);

    //paramètres temporels
    double dt = .5;                      //incréments de temps
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
        //calcule du nouveau vecteur d'Ã©tat
        dY = vector_add(dY,k4);
        dY = scalar_vector_mult(dY,1.0/6.0);
        Y = vector_add(Y,dY);

        //ajout des valeurs calculées dans l'historique des compartiments
        S.push_back(Y[0]);
        I_R.push_back(Y[1]);
        I_H.push_back(Y[2]);
        H.push_back(Y[3]);
        R.push_back(Y[4]);
        T.push_back(T.back() + dt);
        population_size.push_back(Y[0]+Y[1]+Y[2]+Y[3]+Y[4]);
    }
    std::cout << "Vecteur final (S,I_R,I_H,H,R) = ";
    print_vector(Y);

    exportToCSV(T,S,I_R,I_H,H,R,population_size,"SIHR.csv");
    int result = system("py plot_SIHR_data.py SIHR.csv");
    return 0;

}
