#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <math.h>
#include <iomanip>
using namespace std;

struct Rectangle
{
    int position; //подобласть 
    int x_left; 
    int x_right;  //координаты кравеого
    int y_down;
    int y_up;
};

int x_size = 0; //кол-во элементов по оси х
int y_size = 0; //кол-во элементов по оси у
int global_size = 0; //общее кол-во элементов
double h_x = 0; //шаг сетки по х
double h_y = 0; //шаг сетки по у
int L = 0; //кол-во подобластей
int cntKraev = 0; //кол-во краевых условий

vector<int> X_grid; //разбиение на сетку по х (координаты)
vector<int> Y_grid; //разбиение на сетку по у (координаты)
vector<vector <int> > nodes;
vector<vector<double> > G_l; //матрица жесткости
vector<vector<double> > M_l; //матрица масс
vector<double> b_l; //вектор правых частей
vector<vector<double> > localA; //матрица коэффициентов
vector<double> localB; //вектор коэффициентов
vector<double> B; //вектор глобальный
vector<double> b_prev; //предыдущая правая часть b^{j-1}
vector<Rectangle> rectangles;
vector<double> gridt; //временная сетка

//хранение матрицы в разреженном строчном формате
vector<double> di; //диагональ
vector<double> al;
vector<double> au;
vector<int> ig; //индексы строк ненулевых компонет
vector<int> jg; //индексы столбцов ненулевых компонент

//краевые условия
vector<vector<int> > kraevye_uslov; //кравеые условия
vector<vector<double> > val_A3; //значения 3 краевых условий для матрицы
vector<double> val_B; //значения 3 и 2 краевых условий для вектора b

//решение слау
vector<double> x0; //начальное приближение
vector<double> y;
vector<double> z;
vector<double> t;
vector<double> r;
vector<double> p;
vector<double> x;
vector<double> u_prev; //предыдущее решение u^n

int maxiter = 10000;
double e = 1e-18;

//параметры для Кранка-Николсона
double T = 5; //конечное время
double tau = 1/64.; //шаг по времени
int nt = static_cast<int>(T / tau) + 1; //кол-во временных шагов
double sigma = 1.0; //коэффициент sigma

double Func(int num, double x, double y, double t) 
{
    return cos(t);
}

double lambdaV (int num)
{
    switch(num){
        case 1: return 1.0;
        case 2: return 1.0;
        case 3: return 1.;
    }
}

double BettaV(int num)
{
    switch(num)
    {
        case 1: return 1.;
    }
}

double u_Betta(int num, double x, double y, double t)
{
    switch(num)
    {
        case 1: return sin(t);
    }
}

double Tetta(int num, double x, double y, double t)
{
    return 0;
    // switch(num)
    // {
    //     case 1: return -t;
    //     //case 2: return 0.0;
    // }
}

double UG(int num, double x, double y, double t)
{
    return sin(t);
    // switch(num)
    // {
    //     case 1: return x + 5 *t;
    //     case 2: return 1 + y * t;
    //     case 3: return 5 + y * t;
    // }
}

double u_g(int num, double x, double y, double t) {
    return sin(t);
}

/*получение глобального номера узла*/
int IndexOfUnknown(int i, int j) {
    int elementsPerRow = x_size - 1;
    int row = i / elementsPerRow;
    int col = i % elementsPerRow;
    int globalIndex;
    switch (j) {
        case 0: globalIndex = row * x_size + col; break;
        case 1: globalIndex = row * x_size + col + 1; break;
        case 2: globalIndex = (row + 1) * x_size + col; break;
        case 3: globalIndex = (row + 1) * x_size + col + 1; break;
        default: cout << "Invalid local index for element" << endl; return -1;
    }
    return globalIndex;
}

/*построение сетки и временной сетки*/
void GridBuilder()
{
    ifstream input("./test/sxod/matrix.txt");
    input >> x_size >> y_size; 
    X_grid.resize(x_size);
    Y_grid.resize(y_size);
    for (int i = 0; i < x_size; i++)
        input >> X_grid[i];
    for (int i = 0; i < y_size; i++)
        input >> Y_grid[i];
    global_size = x_size * y_size;

    //построение подобластей    
    input >> L;
    rectangles.resize(L);
    for (int i = 0; i < L; i++) 
        input >> rectangles[i].position >> rectangles[i].x_left >> rectangles[i].x_right >> rectangles[i].y_down >> rectangles[i].y_up;
    input.close();

    nodes.resize(global_size);
    int i = 0;
    for (double y_c : Y_grid) {
        for (double x_c : X_grid) {
            nodes[i].push_back(x_c);
            nodes[i].push_back(y_c);
            i++;
        }
    }

    //инициализация временной сетки
    gridt.resize(nt);
    for (int i = 0; i < nt; i++)
        gridt[i] = i * tau;
}

/*построение локальной матрицы жесткости*/
vector<vector<double> > LocalG_matrix(int k, int num_L) 
{
    vector<vector<double> > G(4, vector<double>(4));
    double lambda = lambdaV(k);

    h_x = rectangles[num_L].x_right - rectangles[num_L].x_left;
    h_y = rectangles[num_L].y_up - rectangles[num_L].y_down;
    double a1 = (lambda/6)*(h_y/h_x);
    double a2 = lambda*h_x/(6*h_y);

    G[0][1] = -2 * a1 + a2;
    G[0][2] = a1 - 2 * a2;
    G[0][3] = -a1 - a2;
    G[1][2] = -a1 - a2;
    G[1][3] = a1 - 2 * a2;
    G[2][3] = -2 * a1 + a2;

    G[0][0] = 2 * a1 + 2 * a2;
    G[1][1] = 2 * a1 + 2 * a2;
    G[2][2] = 2 * a1 + 2 * a2;
    G[3][3] = 2 * a1 + 2 * a2;

    G[1][0] = -2 * a1 + a2;
    G[2][0] = a1 - 2 * a2;
    G[2][1] = -a1 - a2;
    G[3][0] = -a1 - a2;
    G[3][1] = a1 - 2 * a2;
    G[3][2] = -2 * a1 + a2;

    return G;
}

/*построение локальной матрицы масс*/
vector<vector<double> > LocalM_matrix(int k, int num_L) 
{
    vector<vector<double> > M(4, vector<double>(4));
    h_x = rectangles[num_L].x_right - rectangles[num_L].x_left;
    h_y = rectangles[num_L].y_up - rectangles[num_L].y_down;
    double a = (h_x * h_y) / 36;
    M[0][1] = 2 * a;
    M[0][2] = 2 * a;
    M[0][3] = a;
    M[1][2] = a;
    M[1][3] = 2 * a;
    M[2][3] = 2 * a;

    M[0][0] = 4 * a;
    M[1][1] = 4 * a;
    M[2][2] = 4 * a;
    M[3][3] = 4 * a;

    M[1][0] = 2 * a;
    M[2][0] = 2 * a;
    M[2][1] = a;
    M[3][0] = a;
    M[3][1] = 2 * a;
    M[3][2] = 2 * a;
    
    return M;
}

/*построение локального вектора правых частей*/
vector<double> LocalB_vector(int num_L, double t) 
{
    h_x = rectangles[num_L].x_right - rectangles[num_L].x_left;
    h_y = rectangles[num_L].y_up - rectangles[num_L].y_down;
    vector<double> b(4);
    double f1 = Func(rectangles[num_L].position, rectangles[num_L].x_left, rectangles[num_L].y_down, t);
    double f2 = Func(rectangles[num_L].position, rectangles[num_L].x_right, rectangles[num_L].y_down, t);
    double f3 = Func(rectangles[num_L].position, rectangles[num_L].x_left, rectangles[num_L].y_up, t);
    double f4 = Func(rectangles[num_L].position, rectangles[num_L].x_right, rectangles[num_L].y_up, t);
    double a = (h_x * h_y) / 36;
    b[0] = a * (4 * f1 + 2 * f2 + 2 * f3 + f4);
    b[1] = a * (2 * f1 + 4 * f2 + f3 + 2 * f4);
    b[2] = a * (2 * f1 + f2 + 4 * f3 + 2 * f4);
    b[3] = a * (f1 + 2 * f2 + 2 * f3 + 4 * f4);

    cout << "VECTOR number of elem: " << num_L + 1 << " at t = " << t << endl;
    for (int i = 0; i < 4; i++) cout << b[i] << " ";
    cout << endl;
    return b;
}

/*вычисление локальной матрицы для Кранка-Николсона*/
vector<vector<double> > LocalMatrix(int num_L, int k) 
{
    vector<vector<double> > A(4, vector<double>(4));
    G_l = LocalG_matrix(num_L, k);
    M_l = LocalM_matrix(num_L, k);
    //матрица для левой части: (sigma/tau) * M + (1/2) * G
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            A[i][j] = (sigma / tau) * M_l[i][j] + (1.0 / 2.0) * G_l[i][j];
    }
    return A;
}

/*вычисление правой части для Кранка-Николсона*/
vector<double> LocalRHS(int num_L, int k, const vector<double>& u_prev, double t_n, double t_np1) 
{
    vector<double> rhs(4, 0.0);
    G_l = LocalG_matrix(num_L, k);
    M_l = LocalM_matrix(num_L, k);
    vector<double> b_n = LocalB_vector(num_L, t_n);     // b^{j-1} на t_n
    vector<double> b_np1 = LocalB_vector(num_L, t_np1); // b^j на t_{n+1}
    
    // Осреднение: (b^j + b^{j-1})/2
    vector<double> b_avg(4);
    for (int i = 0; i < 4; i++) {
        b_avg[i] = (b_np1[i] + b_n[i]) / 2.0;
    }
    
    //глобальные индексы узлов элемента
    vector<int> global_indices(4);
    for (int p = 0; p < 4; p++)
        global_indices[p] = IndexOfUnknown(k, p);
    
    //вычисление ((sigma/tau) * M - (1/2) * G) * u^n
    vector<double> temp(4, 0.0);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double mat_coeff = (sigma / tau) * M_l[i][j] - (1.0 / 2.0) * G_l[i][j];
            temp[i] += mat_coeff * u_prev[global_indices[j]];
        }
    }
    
    //правая часть: ((sigma/tau) * M - (1/2) * G) * u^n + (b^j + b^{j-1})/2
    for (int i = 0; i < 4; i++) {
        rhs[i] = temp[i] + b_avg[i];
    }
    
    return rhs;
}

/*построение портрета глобальной матрицы*/
void Portret_builders()
{
    vector<set<int> > list;
    list.resize(global_size);
    int m = 0;
    for (; m < L; m++) {
        vector<int> tmp;
        for (int p = 3; p >= 0; p--)
            tmp.push_back(IndexOfUnknown(m, p));
        for (int i = 0; i < 4; i++) {
            int ind1 = tmp[i];
            for (int j = i + 1; j < 4; j++) {
                int ind2 = tmp[j];
                list[ind1].insert(ind2);
            }
        }
    }
    
    ig.resize(global_size + 1);
    ig[0] = 0;
    ig[1] = 0;
    for (int i = 0; i < global_size; i++) {
        ig[i + 1] = ig[i] + list[i].size();
    }

    jg.resize(ig[global_size]);
    int k = 0;
    jg[k] = 0;
    for (int i = 0; i < global_size; i++) {
        for (int elem : list[i]) {
            jg[k] = elem;
            k++;
        }
    }

    di.resize(global_size);
    al.resize(ig[global_size]);
    au.resize(ig[global_size]);
}

/*построение глобальной матрицы*/
void GlobalMatrix(vector<vector<double> > &localA, int k)
{
    vector<int> tmp;
    for (int p = 0; p < 4; p++) {
        tmp.push_back(IndexOfUnknown(k, p));
        di[tmp[p]] += localA[p][p];
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) continue;
            int igg = tmp[i];
            int jgg = tmp[j];
            int row = max(igg, jgg);
            int col = min(igg, jgg);
            int index = ig[row];
            while (index < ig[row + 1] && jg[index] != col) index++;
            if (index == ig[row + 1]) {
                cout << "Error: Index not found in GlobalMatrix for row " << row << ", col " << col << endl;
                continue;
            }
            if (igg > jgg) {
                al[index] += localA[i][j];
            } else {
                au[index] += localA[i][j];
            }
        }
    }
}


/*вывод глобальной матрицы*/
void PrintGlobalMatrix(int global_size) {
    // Создаем временную двумерную матрицу для удобного отображения
    vector<vector<double> > fullMatrix;
    fullMatrix.resize(global_size, vector<double>(global_size, 0));

    // Заполнение диагональных элементов
    for (int i = 0; i < global_size; i++) {
        fullMatrix[i][i] = di[i];
    }

    // Заполнение недиагональных элементов
    for (int i = 0; i < global_size; i++) {
        for (int j = ig[i]; j < ig[i + 1]; j++) {
            int col = jg[j];        // Столбец из разреженной структуры
            fullMatrix[i][col] = al[j]; // Значение из массива al
            fullMatrix[col][i] = au[j]; // Симметричность матрицы
        }
    }

    // Вывод матрицы в консоль
    cout << "Global Matrix:" << endl;
    for (int i = 0; i < global_size; i++) {
        for (int j = 0; j < global_size; j++) {
            cout << setw(10) << fullMatrix[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Global Vector:" << endl;
    for (int i = 0; i < global_size; i++) {
        cout << setw(10) << B[i] << " ";
    }
    cout << endl;
}

/*сборка глобального вектора*/
vector<double> BuildGlobalVector(vector<double> &l_B, int k) 
{
    vector<double> tmp(4);
    for (int i = 0; i < 4; i++) tmp[i] = IndexOfUnknown(k, i);
    for (int i = 0; i < 4; i++) {
        int ind = tmp[i];
        B[ind] += l_B[i];
    }
    return B;
}

/*построение краевых условий*/
void Kraevye()
{
    ifstream input("./test/sxod/kraevie.txt");
    input >> cntKraev;
    kraevye_uslov.resize(cntKraev);
    for (int i = 0; i < cntKraev; i++) {
        kraevye_uslov[i].resize(4);
        for (int j = 0; j < 4; j++)
            input >> kraevye_uslov[i][j];
    }
    input.close();
}

/*учёт 3х краевых условий*/
void Third_K(const vector<int> &kraevye_uslov, double t) {
    vector<int> tmp;
    int i = kraevye_uslov[2];
    int j = kraevye_uslov[3];
    double betta = BettaV(kraevye_uslov[1]);
    double u_betta1 = u_Betta(kraevye_uslov[1], nodes[i][0], nodes[i][1], t);
    double u_betta2 = u_Betta(kraevye_uslov[1], nodes[j][0], nodes[j][1], t);
    tmp.push_back(i);
    tmp.push_back(j);
    double h;
    double x1 = nodes[i][0];
    double x2 = nodes[j][0];
    double y1 = nodes[i][1];
    double y2 = nodes[j][1];
    h = (x1 == x2) ? abs(y2 - y1) : abs(x2 - x1);
    double a = (betta * h) / 6.0;
    val_A3.resize(2, vector<double>(2));
    val_A3[0][0] = 2 * a;
    val_A3[0][1] = a;
    val_A3[1][0] = a;
    val_A3[1][1] = 2 * a;
    val_B.resize(2);
    val_B[0] = a * (2 * u_betta1 + u_betta2);
    val_B[1] = a * (u_betta1 + 2 * u_betta2);
    B[i] += val_B[0];
    B[j] += val_B[1];
    di[i] += val_A3[0][0];
    di[j] += val_A3[1][1];
    int ia = j, ja = i;
    int index = ig[ia];
    int flag = 1;
    for (; index < ig[ia + 1] && flag; index++)
        if (jg[index] == ja) flag = 0;
    index--;
    al[index] += val_A3[1][0];
    au[index] += val_A3[0][1];
}

/*учёт 2х краевых условий*/
void Second_K(const vector<int> &kraevye_uslov, double t)
{
    vector<int> tmp;
    double h;
    int i = kraevye_uslov[2];
    int j = kraevye_uslov[3];
    double tetta1 = Tetta(kraevye_uslov[1], nodes[i][0], nodes[i][1], t);
    double tetta2 = Tetta(kraevye_uslov[1], nodes[j][0], nodes[j][1], t);
    tmp.push_back(i);
    tmp.push_back(j);
    double x1 = nodes[i][0];
    double x2 = nodes[j][0];
    double y1 = nodes[i][1];
    double y2 = nodes[j][1];
    h = (x1 == x2) ? abs(y2 - y1) : abs(x2 - x1);
    double a = h / 6.0;
    val_B.resize(2);
    val_B[0] = a * (2 * tetta1 + tetta2);
    val_B[1] = a * (tetta1 + 2 * tetta2);
    B[i] += val_B[0];
    B[j] += val_B[1];
}

/*учёт 1х краевых условий*/
void First_K(const vector<int> &kraevye_uslov, double t) {
    vector<int> tmp;
    int p = kraevye_uslov[2];
    int s = kraevye_uslov[3];
    tmp.push_back(p);
    tmp.push_back(s);
    for (int globalNode : tmp) {
        double x_coord = nodes[globalNode][0];
        double y_coord = nodes[globalNode][1];
        double u_g = UG(kraevye_uslov[1], x_coord, y_coord, t);
        di[globalNode] = 1;
        B[globalNode] = u_g;
        for (int i = ig[globalNode]; i < ig[globalNode + 1]; i++) {
            al[i] = 0;
        }
        for (int i = ig[globalNode]; i < ig[global_size]; i++) {
            if (jg[i] == globalNode)
                au[i] = 0;
        }
    }
}

void AddKraevye(int num, vector<int> &kraevye_uslov, double t) 
{
    switch(num) {
        case 1: First_K(kraevye_uslov, t); break;
        case 2: Second_K(kraevye_uslov, t); break;
        case 3: Third_K(kraevye_uslov, t); break;
        default: cout << "Incorrect Boundary Condition Number" << endl; break;
    }
}

/*сборка глобальной матрицы и вектора*/
void BuildGlobalMatrixAndVector(const vector<double>& u_prev, double t_n, double t_np1) {
    //обнуление глобальных массивов
    fill(di.begin(), di.end(), 0.0);
    fill(al.begin(), al.end(), 0.0);
    fill(au.begin(), au.end(), 0.0);
    fill(B.begin(), B.end(), 0.0);
    
    //сборка матрицы и вектора
    for (int k = 0; k < L; k++) {
        localA = LocalMatrix(rectangles[k].position, k);
        GlobalMatrix(localA, k);
        localB = LocalRHS(k, k, u_prev, t_n, t_np1);
        B = BuildGlobalVector(localB, k);
    }

    //PrintGlobalMatrix(global_size);
    
    //учёт краевых условий на t_{n+1}
    for (int i = 0; i < cntKraev; i++)
        AddKraevye(kraevye_uslov[i][0], kraevye_uslov[i], t_np1);
}

double vector_multiplication(vector<double> &v1, vector<double> &v2)
{
    double sum = 0;
    for (int i = 0; i < global_size; i++)
        sum += v1[i] * v2[i];
    return sum;
}

vector<double> matrix_on_vector_multiplication(vector<double> &v1, vector<double> &v2) 
{
    for (int i = 0; i < global_size; i++)
        v2[i] = di[i] * v1[i];
    for (int i = 0; i < global_size; i++) {
        for (int j = ig[i]; j < ig[i+1]; j++) {
            v2[i] += al[j] * v1[jg[j]];
            v2[jg[j]] += au[j] * v1[i];
        }
    }
    return v2;
}

vector<double> vector_sum(vector<double> &v1, vector<double> &v2, vector<double> &v3, double a)
{
    for (int i = 0; i < global_size; i++)
        v3[i] = v1[i] + a * v2[i];
    return v3;
}

double Norma(vector<double> &v)
{
    double norma = 0;
    for (int i = 0; i < global_size; i++)
        norma += v[i] * v[i];
    return sqrt(norma);
}

void Calc_Zk(double a, vector<double> &v1)
{
    for (int i = 0; i < global_size; i++)
        z[i] = v1[i] + a * z[i];
}

void Calc_Pk(double a, vector<double> &v1)
{
    for (int i = 0; i < global_size; i++)
        p[i] = v1[i] + a * p[i];
}

void Calc_Xk_Rk_L(double alphaK)
{
    for (int i = 0; i < global_size; i++) {
        x[i] = x[i] + alphaK * z[i];
        r[i] = r[i] - alphaK * p[i];
    }
}

double calcResidual() {
    double normb = 0;
    for (int i = 0; i < global_size; i++) {
        normb += B[i] * B[i];
        r[i] = B[i] - di[i] * x[i];
        int m = ig[i + 1];
        for (int k = ig[i]; k < m; k++) {
            int j = jg[k];
            r[i] -= al[k] * x[j];
            r[j] -= au[k] * x[i];
        }
    }
    double mul = vector_multiplication(r, r);
    return sqrt(mul/normb);
}

/*ЛОС с диагональным предобуславливанием*/
void LOS_Dd()
{
    double alpha, betta;
    double normPr = Norma(B);
    matrix_on_vector_multiplication(x0, y);
    for (int i = 0; i < global_size; i++) r[i] = B[i] - y[i];
    for (int i = 0; i < global_size; i++) z[i] = r[i];
    for (int i = 0; i < global_size; i++) x[i] = x0[i];
    matrix_on_vector_multiplication(z, p);
    int k = 0;
    double discrepancy = calcResidual();
    for (; k < maxiter && discrepancy > e; k++) {
        double scal_pk_rk = vector_multiplication(r, p);
        double scal_pk = vector_multiplication(p, p);
        alpha = scal_pk_rk / scal_pk;
        Calc_Xk_Rk_L(alpha);
        matrix_on_vector_multiplication(r, y);
        double Ar_p = vector_multiplication(y, p);
        betta = -Ar_p / scal_pk;
        Calc_Zk(betta, r);
        Calc_Pk(betta, y);
        discrepancy = calcResidual();
    }
    cout << "Iteration: " << k << " RelDiscrepancy: " << discrepancy << endl;
}

void Output(ofstream& out, int time_step)
{
    out << fixed << setprecision(15);
    out << "Time step: " << time_step << " t = " << gridt[time_step] << endl;
    for (int i = 0; i < global_size; i++) {
        out << x[i] << endl;
    }
}

void SLAU()
{
    x0.resize(global_size);
    x.resize(global_size);
    y.resize(global_size);
    z.resize(global_size);
    t.resize(global_size);
    r.resize(global_size);
    p.resize(global_size);
    LOS_Dd();
}


int main() {
    GridBuilder();
    B.resize(global_size);
    u_prev.resize(global_size);
    for (int i = 0; i < global_size; i++) {
        u_prev[i] = 0.;
        //u_prev[i] = u_g(rectangles[0].position, nodes[i][0], nodes[i][1], 0.);
        //u_prev[i] = UG(rectangles[0].position, nodes[i][0], nodes[i][1], 0.0); // u(t=0) по u_g
    }
    b_prev.resize(global_size, 0.0); //инициализация b^{j-1}
    Kraevye();
    Portret_builders();
    
    ofstream out("out.txt");
    
    
    //вычисление b^{j-1} на t = 0
    for (int k = 0; k < L; k++) {
        vector<double> b_tmp = LocalB_vector(k, gridt[0]);
        B = BuildGlobalVector(b_tmp, k);
    }
    b_prev = B; //сохранение b^{j-1}
    
    //временной цикл
    for (int n = 0; n < nt - 1; n++) {
        double t_n = gridt[n];
        double t_np1 = gridt[n + 1];
        cout << "Time step: " << n + 1 << " t = " << t_np1 << endl;
        BuildGlobalMatrixAndVector(u_prev, t_n, t_np1);
        //PrintGlobalMatrix(global_size);
        SLAU();
        Output(out, n + 1);
        
        //обновление b^{j-1} для следующего шага
        fill(B.begin(), B.end(), 0.0);
        for (int k = 0; k < L; k++) {
            vector<double> b_tmp = LocalB_vector(k, t_np1);
            B = BuildGlobalVector(b_tmp, k);
        }
        b_prev = B;
        u_prev = x; //обновление решения для следующего шага
    }
    
    out.close();
    return 0;
}