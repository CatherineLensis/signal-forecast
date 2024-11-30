#define _USE_MATH_DEFINES
//Лазько Екатерина 0522Б1ИСкс1
#include "MyForm.h"
#include "math.h"
#include <chrono>
#include "FFT.h"
#include <cmath>
#include <cstdlib> 
#include <ctime>

using namespace System;
using namespace System::Windows::Forms;

[STAThreadAttribute]

int main(array <String^>^ args)
{
	Application::SetCompatibleTextRenderingDefault(false);
	Application::EnableVisualStyles();

	START::MyForm form;
	Application::Run(% form);
}

System::Void START::MyForm::button1_Click(System::Object^ sender, System::EventArgs^ e)
{
	{
		this->chart1->Series[0]->Points->Clear();
		this->chart1->Series[1]->Points->Clear();
		this->chart1->Series[2]->Points->Clear();
		this->chart2->Series[0]->Points->Clear();
		this->chart2->Series[1]->Points->Clear();
		this->chart2->Series[2]->Points->Clear();
		this->chart2->Series[3]->Points->Clear();

	}

	double f = 0;
	double f1 = 0, f2 = 0, f3 = 0;

	double fd = 0, t1 = 0, fi_1 = 0, t2 = 0, fi_2 = 0, fi_3 = 0;



	f1 = System::Convert::ToDouble(textBox2->Text);
	f2 = System::Convert::ToDouble(textBox5->Text);
	f3 = System::Convert::ToDouble(textBox8->Text);
	t1 = System::Convert::ToDouble(textBox11->Text);
	t2 = System::Convert::ToDouble(textBox12->Text);
	fd = System::Convert::ToDouble(textBox10->Text);
	double K = System::Convert::ToDouble(textBox1->Text);
	double TH = System::Convert::ToDouble(textBox3->Text);

	int N = static_cast<int>(System::Convert::ToDouble(textBox14->Text));

	std::vector<double> X(N);


	for (int i = 0; i < N; i++) {
		double t = (1 / fd) * i;
		if (t < t1) f = f1;
		else if (t >= t1 && t < t2) f = f2;
		else f = f3; //(t >= t2 && t <= N / fd)
		X[i] = sin(2 * M_PI * f * t);
		//this->chart1->Series[1]->Points->AddXY(t, X[i]);//выводит график с разрывом
	};

	///////////////Избавление от разрыва////////////////////////


	double u = 0;
	double n1 = t1 * fd, n2 = t2 * fd;
	for (int i = 0; i < N; i++) {
		double t = (1.0 / fd) * i;
		if (i < n1) {
			u += 2 * M_PI * f1 * (1.0 / fd);
		}
		else if (i >= n1 && i < n2) {
			u += 2 * M_PI * f2 * (1.0 / fd);
		}
		else {
			u += 2 * M_PI * f3 * (1.0 / fd);
		}
		X[i] = sin(u);
		this->chart1->Series[2]->Points->AddXY(t, X[i]);  // График без разрывов
	}

	////////////ПРОГНОЗ///////////////////////
	double dt = 1 / fd;
	std::vector<double> Xpr(N);
	std::vector <double> delta(N); // Отклонение прогноза от исходного сигнала
	for (int i = 2; i < N; i++) {
		Xpr[i] = 2 * cos(2 * M_PI * f2 * dt) * X[i - 1] - X[i - 2];
		delta[i] = ((X[i] - Xpr[i])* (X[i] - Xpr[i]));
	}

	for (int i = 2; i < N; i++) {
		double t = (1 / fd) * i;
		t = Math::Round(t, 2); // Округление
		this->chart2->Series[1]->Points->AddXY(t, pow(X[i] - Xpr[i], 2));
	}
	//Детектирование сигнала второй частоты
	std::vector <double> delta_W(N, 0); // Cвёртка дельты и ф-ии окна
	for (int i = 0; i < N; ++i) // Продвижение по отсчётам
	{
		delta_W[i] = 0;
		for (int j = 0; j < K; ++j) // Продвижение по окну
		{
			int index = i + j - K / 2;
			if (index < 0 || index >= N) {
				continue; // За пределами массива пропускаем
			}
			double weight = 0.5 - 0.5 * cos(2 * M_PI * j / (K - 1)); // Окно Ханна
			delta_W[i] += delta[index] * weight;
		}
		delta_W[i] /= K; // Нормализация
		double t = i * dt;
		t = Math::Round(t, 2); // Округление
		this->chart2->Series[2]->Points->AddXY(t, delta_W[i]); // Отрисовка
		this->chart2->Series[3]->Points->AddXY(t, TH); // Отрисовка порога
	}

	//for (int i = 0; i < N; ++i) // Продвижение по отсчётам
	//{
	//	
	//	for (int j = 0; j < K; ++j) // Продвижение по окну
	//	{
	//		
	//		if (i + j - (K / 2) < 0 || i + j - (K / 2) >= N)
	//			delta_W[i] += 0;
	//		else
	//			delta_W[i] += delta[i + j - (K / 2)] / K;

	//		 // Среднее значение в окне
	//	}
	//	double t = i * dt;
	//	t = Math::Round(t, 2); // Округление
	//	this->chart2->Series[2]->Points->AddXY(t, delta_W[i]); // Отрисовка
	//	this->chart2->Series[3]->Points->AddXY(t, TH); // Отрисовка порога
	//}

	// Оценка длительности сигнала второй частоты
	double t_1 = -1;
	double t_2 = -1;
	bool in_second_frequency = false; // Находится ли сигнал в области второй частоты? (по умолчанию нет)
	for (int i = 0; i < N; ++i)
	{
		// Проверяем, пересекает ли сигнал порог
		if (delta_W[i] >= TH && !in_second_frequency)
		{
			t_2 = i * dt;
			in_second_frequency = true;
		}
		else if (delta_W[i] < TH && in_second_frequency)
		{
			t_1 = i * dt;
			in_second_frequency = false;
		}
	}
	t_1 = Math::Round(t_1, 2); // Округление
	t_2 = Math::Round(t_2, 2); // Округление
	textBox4->Text = Convert::ToString(t_1);
	textBox6->Text = Convert::ToString(t_2);


	///////////////////СГЛАЖИВАНИЕ////////////////////
	 //k = Convert::ToDouble(K->Text); //высота окна чем выше окно тем больше берём высоких частот
	 //th = Convert::ToDouble(TH->Text); //порог


	return System::Void();
}
