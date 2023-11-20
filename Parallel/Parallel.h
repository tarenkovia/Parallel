#define _USE_MATH_DEFINES
#include <locale.h>
#include <time.h>
#include <ppl.h>
#include <vector>
#include <deque>
#include <list>
#include <algorithm>
#include <iostream>
#include "Task1.h"
#include "Task1_1.h"
//#define NNN 100

//int main()
//{
//    setlocale(LC_ALL, ".ACP");
//    double V[NNN];
//    double Time = clock();
//    for (int k = 1; k <= NNN; k++) {
//        V[k-1] = TaskOne(100*cos(k));
//    }
//    Time = (clock() - Time) / CLOCKS_PER_SEC;
//    std::cout << "Подзадачи завершены" << std::endl
//        << "Время посл. обр.: " << Time << " сек." << std::endl;
//    Time = clock();
//    Concurrency::parallel_for(0, NNN, [&V](size_t k) {V[k] = TaskOne(k); });
//    Time = (clock() - Time) / CLOCKS_PER_SEC;
//    std::cout << "Подзадачи завершены" << std::endl
//        << "Время парал. обр.: " << Time << " сек." << std::endl;
//}

//int main(int argc, char* argv[])
//{
//    setlocale(LC_ALL, ".ACP");
//    //Обработка векторов
//    std::cout << "Обработка векторов" << std::endl;
//    std::vector<double> V(NNN);
//    for (int k = 1; k <= NNN; k++)
//        V[k-1] = 100*cos(k);
//    std::vector<double> VP(V);//Конструируем копию
//    double Time = clock();
//    std::for_each(V.begin(), V.end(), [](double& x) {x = Task1_1(x); });
//    Time = (clock() - Time) / CLOCKS_PER_SEC;
//    std::cout << "Подзадачи завершены" << std::endl
//        << "Время посл. обр: " << Time << " сек." << std::endl;
//    Time = clock();
//    concurrency::parallel_for_each(VP.begin(), VP.end(), [](double& x) {x = Task1_1(x); });
//    Time = (clock() - Time) / CLOCKS_PER_SEC;
//    std::cout << "Подзадачи завершены" << std::endl
//        << "Время парал. обр.: " << Time << " сек." << std::endl;
//    //Обработка деков
//    std::cout << "Обработка деков" << std::endl;
//    std::deque<double> Dq;
//    for (int k = 0; k < NNN; k++)
//        Dq.push_back(k + 0.7);
//    std::deque<double> DqP(Dq);
//    Time = clock();
//    std::for_each(Dq.begin(), Dq.end(), [](double& x) {x = Task1_1(x); });
//    Time = (clock() - Time) / CLOCKS_PER_SEC;
//    std::cout << "Подзадачи завершены" << std::endl
//        << "Время посл. обр: " << Time << " сек." << std::endl;
//    Time = clock();
//    concurrency::parallel_for_each(DqP.begin(), DqP.end(), [](double& x) {x = Task1_1(x); });
//    Time = (clock() - Time) / CLOCKS_PER_SEC;
//    std::cout << "Подзадачи завершены" << std::endl
//        << "Время парал. обр.: " << Time << " сек." << std::endl;
//    //Обработка списков
//    std::cout << "Обработка списков" << std::endl;
//    std::list<double> Lst;
//    for (int k = 0; k < NNN; k++)
//        Lst.push_back(k + 0.7);
//    std::list<double> LstP(Lst);
//    Time = clock();
//    std::for_each(Lst.begin(), Lst.end(), [](double& x) {x = Task1_1(x); });
//    Time = (clock() - Time) / CLOCKS_PER_SEC;
//    std::cout << "Подзадачи завершены" << std::endl
//        << "Время посл. обр: " << Time << " сек." << std::endl;
//    Time = clock();
//    concurrency::parallel_for_each(LstP.begin(), LstP.end(), [](double& x) {x = Task1_1(x); });
//    Time = (clock() - Time) / CLOCKS_PER_SEC;
//    std::cout << "Подзадачи завершены" << std::endl
//        << "Время парал. обр.: " << Time << " сек." << std::endl;
//}

//#include <cmath>
//#define N1 100

//double Task2(double x)
//{
//	double Tmp = 0;
//	for (int k = 1; k <= std::max(20, 20 * (int)abs(x)); k++)
//		for (int j = 1; j <= std::max(20, 20 * (int)abs(x)); j++)
//			Tmp += ((x * x - x) * (k - 2 * j)) / (x * x + pow(k, 3) + pow(j, 3)) * cos((k + 2 * j) * x);
//	return Tmp;
//}
//
//int main(int argc, char* argv[])
//{
//	setlocale(LC_ALL, ".ACP");
//	std::vector<double> vX0(N1);
//	std::deque<double> dX0(N1);
//	std::list<double> lX0(N1);
//	for (int k = 1; k <= N1; k++)
//	{
//		vX0[k - 1] = 100 * cos(k);
//		dX0.push_back(vX0[k - 1]);
//		lX0.push_back(vX0[k - 1]);
//	}
//	std::vector<double> vY0(N1);
//	std::deque<double> dY0(N1);
//	std::list<double> lY0(N1);
//	for (int k = 1; k <= N1; k++)
//	{
//		vY0[k - 1] = 100 * cos(k);
//		dY0.push_back(vY0[k - 1]);
//		lY0.push_back(vY0[k - 1]);
//	}
//
//	double Time = clock();
//	std::vector<double>::iterator vIt = std::transform(vX0.begin(), vX0.end(), vY0.begin(), [](double x) -> double {return Task2(x); });
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	std::cout << "Преобразование завершено" << std::endl
//		<< "Время посл. обр: " << Time << " сек." << std::endl;
//	Time = clock();
//	std::deque<double>::iterator dIt = std::transform(dX0.begin(), dX0.end(), dY0.begin(), [](double x) -> double {return Task2(x); });
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	std::cout << "Преобразование завершено" << std::endl
//		<< "Время посл. обр: " << Time << " сек." << std::endl;
//	Time = clock();
//	std::list<double>::iterator lIt = std::transform(lX0.begin(), lX0.end(), lY0.begin(), [](double x) -> double {return Task2(x); });
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	std::cout << "Преобразование завершено" << std::endl
//		<< "Время посл. обр: " << Time << " сек." << std::endl;
//
//	Time = clock();
//	vIt = concurrency::parallel_transform(vX0.begin(), vX0.end(), vY0.begin(), [](double x) -> double {return Task2(x); });
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	std::cout << "Преобразование завершено" << std::endl
//		<< "Время парал. обр.: " << Time << " сек." << std::endl;
//	Time = clock();
//	dIt = concurrency::parallel_transform(dX0.begin(), dX0.end(), dY0.begin(), [](double x) -> double {return Task2(x); });
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	std::cout << "Преобразование завершено" << std::endl
//		<< "Время парал. обр.: " << Time << " сек." << std::endl;
//	Time = clock();
//	lIt = concurrency::parallel_transform(lX0.begin(), lX0.end(), lY0.begin(), [](double x) -> double {return Task2(x); });
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	std::cout << "Преобразование завершено" << std::endl
//		<< "Время парал. обр.: " << Time << " сек." << std::endl;
//}


#include <complex>
//#define N 50000000
//using namespace std;
//
//int main(int argc, char* argv[])
//{
//	setlocale(LC_ALL, ".ACP");
//	vector<complex<double>> Z0(N);
//	for (int k = 0; k < N; k++)
//		Z0[k] = exp(complex<double>(0, sqrt(k) / (1 + sqrt(k)) * sin(k)));
//	vector<complex<double>> Z(Z0);//Конструируем копию
//
//	double Time = clock();
//	sort(Z.begin(), Z.end(),
//		[](const complex<double>& left, const complex<double>& right) {
//			return left.real() < right.real(); });
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	cout << "Сортировка завершена" << endl
//		<< "Время посл. сортировки: " << Time << " сек." << endl;
//	Z = Z0;
//	Time = clock();
//	concurrency::parallel_sort(Z.begin(), Z.end(),
//		[](const complex<double>& left, const complex<double>& right) {
//			return left.real() < right.real(); });
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	cout << "Сортировка завершена" << endl
//		<< "Время парал. сортировки: " << Time << " сек." << endl;
//	Z = Z0;
//	Time = clock();
//	concurrency::parallel_buffered_sort(Z.begin(), Z.end(),
//		[](const complex<double>& left, const complex<double>& right) {
//			return left.real() < right.real(); });
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	cout << "Сортировка завершена" << endl
//		<< "Время парал. сортировки: " << Time << " сек." << endl;
//}


//#define N 50000000
//using namespace std;
//int main(int argc, char* argv[])
//{
//	setlocale(LC_ALL, ".ACP");
//	vector<size_t> X0(N);
//	for (int k = 0; k < N; k++)
//		X0[k] = size_t(4.0e+9 * sqrt(k) * sin(k)/(1+sqrt(k)));
//	vector<size_t> X(X0);//Конструируем копию
//
//	double Time = clock();
//	sort(X.begin(), X.end());
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	cout << "Сортировка завершена" << endl
//		<< "Время посл. сортировки: " << Time << " сек." << endl;
//	X = X0;
//	Time = clock();
//	concurrency::parallel_sort(X.begin(), X.end());
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	cout << "Сортировка завершена" << endl
//		<< "Время parallel_sort: " << Time << " сек." << endl;
//	X = X0;
//	Time = clock();
//	concurrency::parallel_buffered_sort(X.begin(), X.end());
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	cout << "Сортировка завершена" << endl
//		<< "Время parallel_buffered_sort: " << Time << " сек." << endl;
//	X = X0;
//	Time = clock();
//	concurrency::parallel_radixsort(X.begin(), X.end());
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	cout << "Сортировка завершена" << endl
//		<< "Время parallel_radix_sort: " << Time << " сек." << endl;
//
//}


//using namespace std;
//int main(int argc, char* argv[])
//{
//	setlocale(LC_ALL, ".ACP");
//	double V[NNN];
//	double Time = clock();
//	for (int k = 0; k < NNN; k++) V[k] = TaskOne(k);
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	cout << "Подзадачи завершены. Время: " << Time << " сек." << endl;
//	Concurrency::task_group TGr;
//	Time = clock();
//	for (int k = 0; k < NNN; k++) TGr.run([&V, k]() {V[k] = TaskOne(k); });
//	TGr.wait();
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	cout << "Подзадачи завершены. Время: " << Time << " сек." << endl;
//}


// TASK 6 - numb 5

//#include <concurrent_vector.h>
//#include <concurrent_queue.h>
//#include <queue>
//#define N1 400 //Число узлов по каждому направлению
//#define N2 50
//
//double Func(double x, double y)
//{double Tmp = 0;
//for (int k = 1; k <= N2; k++)
//    for (int j = 1; j <= N2; j++)
//        Tmp += cos(k*x)*sin(j*y) / ((1+ k + j)*sqrt(1+k*k*k*k+j*j*j*j));
//return Tmp;
//}
//
//struct Point
//{ double x, y, f; 
//Point(double _x, double _y, double _f)
//{
//	x = _x; y = _y; f = _f;
//}
//};
//
//int main(int argc, char *argv[])
//{setlocale(LC_ALL, ".ACP");
// std::vector<Point> Pts;
// double h = 2.0 * M_PI / N1;
// double Time = clock();
// for (int k = 0; k < N1; k++)
//	 for (int j = 0; j < N1; j++){
//	 double x = h*k;
//	 double y = h*j;
//	 double f = Func(x, y);
//	 if (f <= 0)
//		 Pts.push_back(Point(x, y, f));
//	 }
// Time = (clock() - Time) / CLOCKS_PER_SEC;
// std::cout << "Вычисления завершены. Размер контейнера: "<<Pts.size() << std::endl
//	 << "Время : " << Time << " сек." << std::endl;
// Pts.clear();
// concurrency::concurrent_vector<Point> CPts;
// Time = clock();
// concurrency::parallel_for(0, N1,
//	 [&CPts,h](int k){
//	 for (int j = 0; j < N1; j++){
//		 double x = h*k;
//		 double y = h*j;
//		 double f = Func(x, y);
//		 if (f >= 0)
//			 CPts.push_back(Point(x, y, f));}}
//     );
// Time = (clock() - Time) / CLOCKS_PER_SEC;
// std::cout << "Вычисления завершены. Размер контейнера: " << CPts.size() << std::endl
//	 << "Время : " << Time << " сек." << std::endl;
// CPts.clear();
//}

//int main(int argc, char* argv[])
//{
//	setlocale(LC_ALL, ".ACP");
//	std::queue<Point> Pts;
//	double h = 2.0 * M_PI / N1;
//	double Time = clock();
//	for (int k = 0; k < N1; k++)
//		for (int j = 0; j < N1; j++) {
//			double x = h * k;
//			double y = h * j;
//			double f = Func(x, y);
//			if (f >= 0)
//				Pts.push(Point(x, y, f));
//		}
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	std::cout << "Вычисления завершены. Размер контейнера: " << Pts.size() << std::endl
//		<< "Время : " << Time << " сек." << std::endl;
//
//
//	concurrency::concurrent_queue<Point> CPts;
//	Time = clock();
//	concurrency::parallel_for(0, N1,
//		[&CPts, h](int k) {
//			for (int j = 0; j < N1; j++) {
//				double x = h * k;
//				double y = h * j;
//				double f = Func(x, y);
//				if (f >= 0)
//					CPts.push(Point(x, y, f));
//			}}
//	);
//	Time = (clock() - Time) / CLOCKS_PER_SEC;
//	std::cout << "Вычисления завершены. Размер контейнера: " << CPts.unsafe_size() << std::endl
//		<< "Время : " << Time << " сек." << std::endl;
//	CPts.clear();
//}

// TASK 7 - numb 5
//#pragma once //Concurrent_Simps.h
//#include <ppl.h>
//#define _USE_MATH_DEFINES //Simps_Tst.cpp
//#include <iostream>
//#include <cmath>
//#include <time.h>
//#include <locale.h>
//using namespace std;
//using namespace concurrency;
//typedef double (*Double_Func_Double)(double);
//double Simps(double a, double b, int N, Double_Func_Double);
//double Concurrent_Simps(double a, double b, int N, Double_Func_Double);
//namespace MethodCall {
//	static void* ObjAddr = nullptr; //без __declspec(thread)
//	template <class Ty>
//	double Sub_Int_Func(double x)
//	{
//		return (*((Ty*)ObjAddr))(x);
//	}
//	template <class Ty>
//	double Simpson(double a, double b, int N, Ty const& Obj)
//	{
//		ObjAddr = (void*)&Obj;
//		return Simps(a, b, N, Sub_Int_Func<Ty>);
//	}
//	template <class Ty>
//	double Concurrent_Simpson(double a, double b, int N, Ty const& Obj)
//	{
//		ObjAddr = (void*)&Obj;
//		return Concurrent_Simps(a, b, N, Sub_Int_Func<Ty>);
//	}
//};
//
//double Simps(double a, double b, int N, Double_Func_Double Func)
//{
//	double h = (b - a) / (2 * N);
//	double S1 = 0, S2 = 0;
//	for (int k = 1; k < N; k++) {
//		double Tmp = a + (2 * k - 1) * h; S1 += Func(Tmp); S2 += Func(Tmp + h);
//	}
//	S1 += Func(b - h); return h * (Func(a) + Func(b) + 4.0 * S1 + 2.0 * S2) / 3.0;
//}
//double Concurrent_Simps(double a, double b, int N, Double_Func_Double Func)
//{
//	double h = (b - a) / (2 * N);
//	combinable<double> CS1([]() {return 0.0; }), CS2([]() { return 0.0; });
//	parallel_for(1, N, [a, h, Func, &CS1, &CS2](int k)
//		{double Tmp = a + (2 * k - 1) * h;
//	CS1.local() += Func(Tmp); CS2.local() += Func(Tmp + h); });
//	double S1 = CS1.combine([](double x, double y) {return x + y; });
//	double S2 = CS2.combine([](double x, double y) {return x + y; });
//	S1 += Func(b - h); return h * (Func(a) + Func(b) + 4.0 * S1 + 2.0 * S2) / 3.0;
//}
//
//class My_Class {
//private:
//	double y;
//	const int N = 200000; //Число разбиений отрезка интегрир.
//	const double _Inf = 5000.0; //Факт. верхний предел интегрир.
//public:
//	My_Class(double _y = 0.0) { y = _y; }
//	double Sub_Int_Func(double x)
//	{
//		double Tmp = 15 * log(10.0) + log(abs(x)) - sqrt(x);
//		int N = Tmp > 0 ? ceil(Tmp * Tmp + 1) : 1;
//		double P = exp(-abs(x));
//		double Tmp2 = pow(cos(y), 2);
//		for (int k = 0; k <= N; k++)
//			P *= cos(x / (Tmp2 + exp(sqrt((double)k))));
//		return P;
//	}
//	double Quad()
//	{
//		return MethodCall::Simpson(0.0, _Inf, N,
//			[this](double x) {return Sub_Int_Func(x); });
//	}
//	double Concurrent_Quad()
//	{
//		return MethodCall::Concurrent_Simpson(0.0, _Inf, N,
//			[this](double x) {return Sub_Int_Func(x); });
//	}
//};
//int main(void)
//{
//	setlocale(LC_ALL, ".ACP");
//	double y; cout << "y="; cin >> y;
//	My_Class TObj(y);
//	double Tms = clock();
//	double F = TObj.Quad();
//	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
//	cout.precision(8);
//	cout << "F=" << F << endl << "Время=" << Tms << " с" << endl;
//	Tms = clock();
//	F = TObj.Concurrent_Quad();
//	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
//	cout << "F=" << F << endl << "Время=" << Tms << " с" << endl;
//}

#include <agents.h> // conveyor-transformer.cpp
#include <cmath>
#include <iostream>
#include <locale.h>
#include <time.h>
using namespace concurrency;
using namespace std;
int const N1 = 100; //Размерность рядов

int const Chunk = 100; //Порция обрабатываемых элементов
double Func1(double x)//Первая ступень конвейера
{
	double Tmp = 0, Tmpa = x * x;
	for (int n = 0; n <= N1; n++) {
		double Tmpxn = x * x + n, Tmpn3 = n * n * n;
		for (int k = 0; k <= N1; k++) {
			double Tmpk = k * k;
			for (int j = 0; j <= N1; j++)
				Tmp += Tmpxn / (1 + Tmpa + Tmpn3 + k * k + j * j);
		}
	}
	return Tmp;
}

double Func2(double x)//Вторая ступень конвейера
{
	double Tmp = 0, Tmpa = x * x;
	for (int n = 0; n <= N1; n++) {
		double Tmpn = n * n;
		for (int k = 0; k <= N1; k++) {
			double Tmpk = k * k * k;
			for (int j = 0; j <= N1; j++)
				Tmp += 1.0 / (1 + Tmpa + pow((Tmpn + j * j), 3 / 2) + Tmpk);
		}
	}
	return x * x * Tmp;
}

int main()
{
	int const NNN = 1000; //Длина массива на входе конвейера
	setlocale(LC_ALL, ".ACP");
	vector<double> X(NNN), Y(NNN), Z(NNN);
	for (int k = 0; k < X.size(); k++)//Заполнение входного вектора
		X[k] = 10 * sin((double)k);
	double Tms = clock();
	for (int k = 0; k < X.size(); k++) {//Последоват. вычисления
		Y[k] = Func1(X[k]); Y[k] = Func2(Y[k]);
	}

	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout << "Время последовательного алгоритма: " << Tms << " c." << endl;
	transformer<double, double> t1(Func1), //Ступени конвейера
		t2(Func2);
	t1.link_target(&t2); // t1 -> t2
	Tms = clock();
	for (int k = 0; k < Chunk; k++) //Параллельные вычисления
		send(t1, X[k]);
	for (int k = Chunk; k < X.size(); k++) {
		Z[k - Chunk] = receive(t2); send(t1, X[k]);
	}
	for (int k = X.size() - Chunk; k < X.size(); k++) {
		Z[k] = receive(t2);
	}
	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout << "Время параллельного алгоритма: " << Tms << " c."
		<< endl;
}

