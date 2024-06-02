#include <iomanip>
#include <iostream>
#include <clocale>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <locale.h>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>
#include <fixed/fixed_nonlinear_solver.h>
#include <cmath>
#include <gtest/gtest.h>   

using namespace std;
using namespace pde_solvers;

struct pipe
{
    double L; // протяженность участка
    double V; //скорость нефти
    double D_vnesh; //диаметр внешний
    double b; //толщина стенки
    double abc; //абсолютная шероховатость в м
    double z_0, z_L; //высотные отметки начала,конца
    double ro; //плотность
    double u; // кинематическая вязкость Стоксы в м2/с
    double Q; //расход м3/c
    double p_0;// давление в начале участка
    double p_L;// давление в начале участка
    double lambda;//коэфф.гидравл.сопротивления
    double t_w;//касательное напряжение трения
    double h; //шаг по координате расчетной сетки, в метрах
    double Re;
    double N; //КОЛ-ВО ТОЧЕК

    double get_inner_diameter() const {
        return D_vnesh - 2 * b;
    }
    double get_relative_roughness() const {
        return abc / get_inner_diameter();
    }

    double get_inner_area() const {
        double D = get_inner_diameter();
        double S = M_PI * D * D / 4;
        return V * S;
    }

    double get_V() const {
        double D = get_inner_diameter();
        return 4 * Q / (3.1415 * pow(D, 2));
    }

    double get_Re() const {
        double D = get_inner_diameter();
        return get_V() * D / u;
    }
    //касательное напряжение трения
    double get_t_w() const {
        return lambda / 8 * ro * pow(get_V(), 2);
    }

    double get_dx() const {
        return L / (N - 1);
    }

    double get_dt() const {
        return get_dx() / V;
    }

    double get_n() const {
        return (L / V) / get_dt();
    }
};



struct data_iterazii {
    vector<vector<double>> buffer;
    vector<double> pressure;
};



void party_layer(pipe myPipe, double parametr, vector<double>& current_layer, vector<double>& previous_layer) {

    //смещение предыдущего слоя и запись граничного условия
    current_layer[0] = parametr;
    for (size_t i = 1; i < myPipe.N; i++)
    {
        current_layer[i] = previous_layer[i - 1];
    }
}





vector<double>  euler(pipe myPipe, ring_buffer_t<vector<vector<double>>>& buffer, int i) {
    double p_0 = myPipe.p_0;
    double p_begin = myPipe.p_0;
    vector<vector<double>>& current_layer = buffer.current();
    vector<vector<double>>& previous_layer = buffer.previous();
    data_iterazii z;
    z.buffer = buffer.previous();
    z.pressure.clear(); // Очищаем вектор pressure перед заполнением
    if (i == 0) {
        size_t b = current_layer[0].size();
        for (size_t j = 0; j < b; j++) {
            double Re = myPipe.V * myPipe.get_inner_diameter() / previous_layer[1][j];
            double lambda = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            previous_layer[2][j] = p_0;
            double p_rachet = p_0 - myPipe.get_dx() * (lambda / myPipe.get_inner_diameter() * previous_layer[0][j] * pow(myPipe.V, 2) / 2 - M_G * previous_layer[0][j] * (myPipe.z_L - myPipe.z_0) / ((myPipe.N - 1) * myPipe.get_dx()));
            p_0 = p_rachet;
            z.pressure.push_back(previous_layer[2][j]); // Исправлено добавление в вектор
            z.buffer = buffer.previous();

        }
    }

    if (i != 0) {
        for (size_t j = 0; j < current_layer[0].size(); j++) {
            double Re = myPipe.V * myPipe.get_inner_diameter() / current_layer[1][j];
            double lambda = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            previous_layer[2][j] = p_0;
             double p_rachet = p_0 - myPipe.get_dx() * (lambda / myPipe.get_inner_diameter() * previous_layer[0][j] * pow(myPipe.V, 2) / 2 - M_G * previous_layer[0][j] * (myPipe.z_L - myPipe.z_0) / ((myPipe.N - 1) * myPipe.get_dx()));
            p_0 = p_rachet;
            z.buffer = buffer.previous();
        }
    }

    return z.pressure;
}



void excel(const string& filename, pipe myPipe, vector<vector<double>>& previous_layer, int i, vector <double>& data) {

    for (int j = 0; j < myPipe.N; j++) {
        if (i == 0 && j == 0) {
            ofstream outFile(filename, ios::out);
            setlocale(LC_ALL, "ru");
            outFile << "время,координата,плотность, вязкость, давление, разность давления" << "\n";
            outFile << i * myPipe.get_dt() << "," << (j)*myPipe.get_dx() << "," << previous_layer[0][j] << "," << previous_layer[1][j] << "," << previous_layer[2][j] << "," << previous_layer[2][j] - data[j] << "\n";
            outFile.close();
        }
        else {
            ofstream outFile(filename, ios::app); // Используйте ios::app для добавления данных
            outFile << i * myPipe.get_dt() << "," << (j)*myPipe.get_dx() << "," << previous_layer[0][j] << "," << previous_layer[1][j] << "," << previous_layer[2][j] << "," << previous_layer[2][j] - data[j] << "\n";
            outFile.close();

        }

    }
}
void excel_last_pressure(const string& filename, pipe myPipe, vector<vector<double>>& previous_layer, int i) {


    if (i == 0) {
        ofstream outFile(filename, ios::out);
        setlocale(LC_ALL, "ru");
        outFile << "нули,время,давление" << "\n";
        outFile << 0 << "," << i * myPipe.get_dt() << "," << previous_layer[2][myPipe.N - 1] << "\n";
        outFile.close();
    }
    else {
        ofstream outFile(filename, ios::app); // Используйте ios::app для добавления данных
        outFile << 0 << "," << i * myPipe.get_dt() << "," << previous_layer[2][myPipe.N - 1] << "\n";
        outFile.close();

    }


}

/// @brief Труба из заданий блока 3
pipe block3_pipe()
{
    pipe myPipe;
    myPipe.L = 100000;
    myPipe.p_0 = 6e6;
    myPipe.L = 100e3;
    myPipe.D_vnesh = 720e-3;
    myPipe.b = 10e-3;
    myPipe.z_0 = 100;
    myPipe.z_L = 50;
    myPipe.V = 0.5;
    myPipe.N = 101;
    myPipe.abc = 15e-6;
    return myPipe;
}

/// @brief Ступенчатая партия, вытеснение
TEST(BLOCK3, TASK1) 
{

    pipe myPipe = block3_pipe();

    vector<double> ro;

    ro = { 800 }; //значения плотности на входе трубы

    vector<double> u;

    u = { 10e-6 }; //значения кинемат. вязкости на входе трубы

    vector<double> pressure;
    pressure = { 0 };


    vector<double> ro_begin(myPipe.N, 900);
    vector<double> u_begin(myPipe.N, 15e-6);
    vector<double> pressure_begin(myPipe.N, 0);

    data_iterazii z;
    z.buffer;
    z.pressure = { 0 };

    ring_buffer_t<vector<vector<double>>> buffer(2, { ro_begin, u_begin,pressure_begin });

    // Объявляем переменную z вне цикла


    for (size_t i = 0; i < myPipe.get_n() + 2; i++) {
        if (i == 0) {
            euler(myPipe, buffer, i);

            excel_last_pressure("pressure_last_1.1.csv", myPipe, buffer.current(), i);
            party_layer(myPipe, ro[0], buffer.current()[0], buffer.previous()[0]);
            party_layer(myPipe, u[0], buffer.current()[1], buffer.previous()[1]);
            z.pressure = euler(myPipe, buffer, i);

            excel("3block_1.1.csv",myPipe, buffer.previous(), i, z.pressure);

            buffer.advance(1);
        }
        else {
            // Используем переменную z в блоке else
            euler(myPipe, buffer, i);
            excel_last_pressure("pressure_last_1.1.csv", myPipe, buffer.current(), i);
            party_layer(myPipe, ro[0], buffer.current()[0], buffer.previous()[0]);
            party_layer(myPipe, u[0], buffer.current()[1], buffer.previous()[1]);
            excel("3block_1.1.csv",myPipe, buffer.previous(), i, z.pressure);
            buffer.advance(1);
        }
    }


}


/// @brief Импульсная партия
TEST(BLOCK3, TASK2) {

    pipe myPipe = block3_pipe();


    double impulse_party_duration = 5*3600; // длительность захода импульсной партии, сек
    size_t impulse_party_layers = ceil(impulse_party_duration / myPipe.get_dt());
    // время моделирования - время вытеснения начальной партии плюс длительность импульсной партии
    size_t total_layers = myPipe.get_n() + impulse_party_layers + 2; 

    vector<double> ro;
    ro = {};
    for (int k = 0; k < impulse_party_layers; ++k) {
        ro.push_back(990);
    }

    // Заполняем оставшиеся 50 элементов значением 900
    for (int l = impulse_party_layers; l < total_layers; ++l) {
        ro.push_back(900);
    }


    vector<double> u;
    u = {  }; //значения кинемат. вязкости на входе трубы
    for (int k = 0; k < impulse_party_layers; ++k) {
        u.push_back(19e-6);
    }

    // Заполняем оставшиеся 50 элементов значением 900
    for (int l = impulse_party_layers; l < total_layers; ++l) {
        u.push_back(15e-6);
    }




    vector<double> pressure;
    pressure = { 0 };


    vector<double> ro_begin(myPipe.N, 900);
    vector<double> u_begin(myPipe.N, 15e-6);
    vector<double> pressure_begin(myPipe.N, 0);

    data_iterazii z;
    z.buffer;
    z.pressure = { 0 };

    ring_buffer_t<vector<vector<double>>> buffer(2, { ro_begin, u_begin,pressure_begin });

    // Объявляем переменную z вне цикла

    cout << impulse_party_layers;

    for (size_t i = 0; i < total_layers; i++) {
        if (i == 0) {
            euler(myPipe, buffer, i);
            excel_last_pressure("pressure_last_1.2.csv", myPipe, buffer.previous(), i);
            party_layer(myPipe, ro[i], buffer.current()[0], buffer.previous()[0]);
            party_layer(myPipe, u[i], buffer.current()[1], buffer.previous()[1]);
            z.pressure = euler(myPipe, buffer, i);
            excel("3block_1.2.csv", myPipe, buffer.previous(), i, z.pressure);
            buffer.advance(1);
        }
        else {
            // Используем переменную z в блоке else
            euler(myPipe, buffer, i);
            excel_last_pressure("pressure_last_1.2.csv", myPipe, buffer.previous(), i);
            party_layer(myPipe, ro[i], buffer.current()[0], buffer.previous()[0]);
            party_layer(myPipe, u[i], buffer.current()[1], buffer.previous()[1]);
            excel("3block_1.2.csv", myPipe, buffer.previous(), i, z.pressure);
            buffer.advance(1);
        }
    }







}


