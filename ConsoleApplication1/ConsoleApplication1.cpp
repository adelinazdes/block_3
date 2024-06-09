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


    double time;//время моделирования

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
    //шаг метода характеристик
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

/*vector<double>  euler(pipe myPipe, ring_buffer_t<vector<vector<double>>>& buffer, int i)*/ 
data_iterazii euler(pipe myPipe, ring_buffer_t<vector<vector<double>>>& buffer, int i) {
    double p_0 = myPipe.p_0;
    double p_begin = myPipe.p_0;
    vector<vector<double>>& current_layer = buffer.current();
    vector<vector<double>>& previous_layer = buffer.previous();
   
    data_iterazii z;
    
    vector<double> ro_profile(myPipe.N);
    vector<double> u_profile(myPipe.N);
    vector<double> Q_profile(myPipe.N);
    vector<double> pressure_begin_(myPipe.N);
    vector<double> pressure_last_(myPipe.N);

    std::vector<std::vector<std::vector<double>>> z_buffer(8, { ro_profile,  u_profile,Q_profile, pressure_begin_,pressure_last_ });
   /* z.buffer(8, { ro_profile,  u_profile,Q_profile, pressure_begin_,pressure_last_ };*/
    /*z.buffer[0] = { ro_profile };
    z.buffer[1] = { u_profile };
    z.buffer[2] = { Q_profile };
    z.buffer[3] = { pressure_begin_ };
    z.buffer[4] = { pressure_last_ };*/
    //vector<vector<double>>& 1_layer = z.buffer();
    //vector<vector<double>>& 2_layer = z.buffer();
    //vector<vector<double>> & 3_layer = z.buffer();
    //vector<vector<double>> & 4_layer = z.buffer();
    //vector<vector<double>> & 5_layer = z.buffer();
    //vector<vector<double>> & 6_layer = z.buffer();
    //vector<vector<double>> & 7_layer = z.buffer();
    //vector<vector<double>> & 8_layer = z.buffer();
    //вектор хранит значения давлений на первом слое
    z.pressure.clear(); // Очищаем вектор pressure перед заполнением
    if (i == 0) {
        size_t b = current_layer[0].size();
        for (size_t j = 0; j < b; j++) {
            double Re = myPipe.V * myPipe.get_inner_diameter() / previous_layer[1][j];
            double lambda = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            previous_layer[2][j] = p_0;
            double p_rachet = p_0 - myPipe.get_dx() * (lambda / myPipe.get_inner_diameter() * previous_layer[0][j] * pow(myPipe.V, 2) / 2 - M_G * previous_layer[0][j] * (myPipe.z_L - myPipe.z_0) / ((myPipe.N - 1) * myPipe.get_dx()));
            p_0 = p_rachet;
            z.pressure.push_back(previous_layer[2][j]); //добавление в вектор давления на первом слое
           
            ////БУФЕР для 3 задачи
            // z.buffer[0][j] = previous_layer[0][j];
            //z.buffer[1][j] = previous_layer[1][j];
            //z.buffer[2][j] = myPipe.Q;
            //z.buffer[3][j] = previous_layer[2][j];
            //z.buffer[4][j] = previous_layer[2][myPipe.N - 1];

        }
    }

    if (i != 0) {
        for (size_t j = 0; j < current_layer[0].size(); j++) {
            double Re = myPipe.V * myPipe.get_inner_diameter() / current_layer[1][j];
            double lambda = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            previous_layer[2][j] = p_0;
             double p_rachet = p_0 - myPipe.get_dx() * (lambda / myPipe.get_inner_diameter() * previous_layer[0][j] * pow(myPipe.V, 2) / 2 - M_G * previous_layer[0][j] * (myPipe.z_L - myPipe.z_0) / ((myPipe.N - 1) * myPipe.get_dx()));
            p_0 = p_rachet;
            ////БУФЕР для 3 задачи
            //z.buffer[0][j] = previous_layer[0][j];
            //z.buffer[1][j] = previous_layer[1][j];
            //z.buffer[2][j] = myPipe.Q;
            //z.buffer[3][j] = previous_layer[2][j];
            //z.buffer[4][j] = previous_layer[2][myPipe.N - 1];
        }
    }

    /*return z.pressure;*/
    return z;
}

void excel(const string& filename, pipe myPipe, vector<vector<double>>& previous_layer, int i, vector <double>& data) {

    for (int j = 0; j < myPipe.N; j++) {
        if (i == 0 && j == 0) {
            ofstream outFile(filename, ios::out);
            setlocale(LC_ALL, "ru");
            outFile << "время,координата,плотность, вязкость, давление, разность давления" << "\n";
            outFile << i * myPipe.time << "," << (j)*myPipe.get_dx() << "," << previous_layer[0][j] << "," << previous_layer[1][j] << "," << previous_layer[2][j] << "," << previous_layer[2][j] - data[j] << "\n";
            outFile.close();
        }
        else {
            ofstream outFile(filename, ios::app); // Используйте ios::app для добавления данных
            outFile << i * myPipe.time << "," << (j)*myPipe.get_dx() << "," << previous_layer[0][j] << "," << previous_layer[1][j] << "," << previous_layer[2][j] << "," << previous_layer[2][j] - data[j] << "\n";
            outFile.close();

        }

    }
}
void excel_last_pressure(const string& filename, pipe myPipe, vector<vector<double>>& previous_layer, int i) {

//давление в последней точке трубопровода на каждом слое
    if (i == 0) {
        ofstream outFile(filename, ios::out);
        setlocale(LC_ALL, "ru");
        outFile << "нули,время,давление" << "\n";
        outFile << 0 << "," << i * myPipe.time << "," << previous_layer[2][myPipe.N - 1] << "\n";
        outFile.close();
    }
    else {
        ofstream outFile(filename, ios::app); // Используйте ios::app для добавления данных
        outFile << 0 << "," << i * myPipe.time << "," << previous_layer[2][myPipe.N - 1] << "\n";
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
    myPipe.time = myPipe.get_dt();
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
            z= euler(myPipe, buffer, i);

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
    myPipe.time = myPipe.get_dt();

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
            z = euler(myPipe, buffer, i);
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


/// @brief Труба из заданий блока 3
pipe block3_task4_pipe()
{
    pipe myPipe;
    myPipe.L = 100000;
    myPipe.L = 100e3;
    myPipe.D_vnesh = 720e-3;
    myPipe.b = 10e-3;
    myPipe.z_0 = 100;
    myPipe.z_L = 50;
    myPipe.N = 101;
    myPipe.abc = 15e-6;
    return myPipe;
}


/////констурктор для буфера задачи 3
//int main() {
//    vector<double> ro_profile;
//    vector<double> u_profile;
//    vector<double> Q_profile;
//    vector<double> pressure_begin_;
//    vector<double> pressure_last_;
//
//    data_iterazii z;
//    z.buffer.resize(8);
//    z.buffer[0] = { ro_profile };
//    z.buffer[1] = { u_profile };
//    z.buffer[2] = { Q_profile };
//    z.buffer[3] = { pressure_begin_ };
//    z.buffer[4] = { pressure_last_ };
//
//    return 0;
//}





vector<double> diskretizatsia (pipe myPipe, vector<double>& input_Q) {
    vector<double> time_diskretizatsii; 
    time_diskretizatsii.push_back(0.0);
    for (size_t i = 0; i < input_Q.size(); i++) {
        myPipe.Q = input_Q[i];
        myPipe.V = myPipe.get_V();
        //шаг метода характеристик (зависимость от расхода)
        double step_characteristics_method = myPipe.get_dt();
        if (i == 0) {
            time_diskretizatsii.push_back(step_characteristics_method);
        }
        else {
            double time = time_diskretizatsii.back() + step_characteristics_method;
            time_diskretizatsii.push_back(time);
        }
    }
    return time_diskretizatsii;

}

vector<double> interpolyazia(vector<double> time_diskretizatsii, vector<double>& input) {
    vector<double> interpolyazia_value;
    double step_diskretizatsii = 2000;
    for (size_t i = 0; i < time_diskretizatsii.size(); i++) {
        double step = time_diskretizatsii[i] / step_diskretizatsii;
        double first_slice = floor(step);
        double second_slice = ceil(step);

        if (first_slice >= input.size()) {
            first_slice = input.size() - 1;
        }
        if (second_slice >= input.size()) {
            second_slice = input.size() - 1;
        }

        double x = time_diskretizatsii[i];
        double x0 = first_slice * step_diskretizatsii;
        double x1 = second_slice * step_diskretizatsii;
        double y0 = input[first_slice];
        double y1 = input[second_slice];

        double y;
        if (abs(x - x0) < 1e-4) {
            y = y0;
        }
        else {
            y = y0 + (x - x0) * (y0 - y1) / (x1 - x0);
        }
           //input[first_slice] + (time_diskretizatsii[i] - first_slice * step_diskretizatsii) * (input[first_slice] - input[second_slice]) / (second_slice * step_diskretizatsii - first_slice * step_diskretizatsii)
        interpolyazia_value.push_back(y);
    }
    return interpolyazia_value;
}








/// @brief задача 3
TEST(BLOCK3, TASK3) {

    pipe myPipe = block3_task4_pipe();
    
    double step_diskretizatsii = 2000;
    myPipe.time = step_diskretizatsii;
    //шаг метода характеристик (зависимость от расхода)
    double step_characteristics_method;

    vector<double> ro;
    ro = { 900,
        880,
        880,
        890,
        890,
        880,
        880,
        870 };
   
    vector<double> u;
    u = { 0.000015,
        0.000013,
        0.000013,
        0.000014,
        0.000014,
        0.000013,
        0.000013,
        0.000012 }; //значения кинемат. вязкости на входе трубы
   
    vector<double> pressure;
    pressure = { 6000000,
        5800000,
        5800000,
        5900000,
        5900000,
        5800000,
        5800000,
        5700000 };

    vector<double> Q;
    Q = { 0.192325,
        0.200000,
        0.210000,
        0.200000,
        0.180000,
        0.210000,
        0.210000,
        0.210000 };

    
    size_t  total_layers = Q.size();

    vector<double> ro_begin(myPipe.N, ro[0]);
    vector<double> u_begin(myPipe.N, u [0]);
    vector<double> pressure_begin(myPipe.N, pressure[0]);
    



    data_iterazii z;
    z.buffer;
    z.pressure = { 0 };

    vector<double> time_diskretizatsii= diskretizatsia (myPipe,Q);
    vector<double> ro_int= interpolyazia(time_diskretizatsii, ro);
    vector<double> u_int = interpolyazia(time_diskretizatsii, u);
    vector<double> pressure_int = interpolyazia(time_diskretizatsii, ro);
    vector<double> Q_int = interpolyazia(time_diskretizatsii, ro);


    //ring_buffer_t<vector<vector<double>>> buffer(2, { ro_begin, u_begin,pressure_begin });
    //    
    //for (size_t i = 0; i < total_layers; i++) {
    //    if (i == 0) {
    //        myPipe.p_0 = pressure[i];
    //        myPipe.Q = Q[i];
    //        myPipe.V = myPipe.get_V();
    //        //шаг метода характеристик (зависимость от расхода)
    //        step_characteristics_method = myPipe.get_dt();
    //        time_diskretizatsii.push_back(step_characteristics_method);
    //        euler(myPipe, buffer, i);
    //        excel_last_pressure("pressure_last_1.3.csv", myPipe, buffer.previous(), i);
    //        party_layer(myPipe, ro[i], buffer.current()[0], buffer.previous()[0]);
    //        party_layer(myPipe, u[i], buffer.current()[1], buffer.previous()[1]);
    //        z= euler(myPipe, buffer, i);
    //        excel("3block_1.3.csv", myPipe, buffer.previous(), i, z.pressure);
    //        buffer.advance(1);
    //        cout << step_characteristics_method;
    //        ring_buffer_t<vector<vector<double>>> new_buffer(total_layers, { i,pressure });
    //    }
    //    else {
    //        // Используем переменную z в блоке else
    //        myPipe.p_0 = pressure[i];
    //        myPipe.Q = Q[i];
    //        myPipe.V = myPipe.get_V();
    //        step_characteristics_method = time_diskretizatsii[i - 1] + myPipe.get_dt();
    //        time_diskretizatsii.push_back(step_characteristics_method);
    //        euler(myPipe, buffer, i);
    //        excel_last_pressure("pressure_last_1.3.csv", myPipe, buffer.previous(), i);
    //        party_layer(myPipe, ro[i], buffer.current()[0], buffer.previous()[0]);
    //        party_layer(myPipe, u[i], buffer.current()[1], buffer.previous()[1]);
    //        excel("3block_1.3.csv", myPipe, buffer.previous(), i, z.pressure);
    //        buffer.advance(1);

    //    }
    //}
        







}
