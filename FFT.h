#include <complex>
#include <vector>
#include <algorithm>
#include <random>

typedef std::complex<double> Complex;

class cmplx {
public:
	double real;
	double image;

	cmplx() { real = image = 0.; }
	cmplx(double x, double y) { real = x; image = y; }
	cmplx& operator = (cmplx&);
	friend cmplx operator * (cmplx& x, cmplx& y);
	friend cmplx operator / (cmplx& x, cmplx& y);
	friend cmplx operator / (cmplx& x, double y);
	friend cmplx operator + (cmplx& x, cmplx& y);
	friend cmplx operator - (cmplx& x, cmplx& y);
};

//----------------------------------------------------------------
cmplx operator * (cmplx& x, cmplx& y)
{
	cmplx z;
	z.real = x.real * y.real - x.image * y.image;
	z.image = x.real * y.image + y.real * x.image;
	return z;
}

//----------------------------------------------------------------
cmplx operator / (cmplx& x, cmplx& y)
{
	cmplx z;
	double y2 = y.real * y.real + y.image * y.image;
	if (y2 < 10e-40)  return z;
	z.real = (x.real * y.real + x.image * y.image) / y2;
	z.image = (y.real * x.image - x.real * y.image) / y2;
	return z;
}

cmplx operator / (cmplx& x, double y)
{
	cmplx z;
	if (y < 10e-40)  return z;
	z.real = x.real / y;
	z.image = x.image / y;
	return z;
}

cmplx operator + (cmplx& x, cmplx& y)
{
	cmplx z;
	z.real = x.real + y.real;
	z.image = x.image + y.image;
	return z;
}

cmplx operator - (cmplx& x, cmplx& y)
{
	cmplx z;
	z.real = x.real - y.real;
	z.image = x.image - y.image;
	return z;
}

cmplx& cmplx::operator = (cmplx& c)
{
	real = c.real;
	image = c.image;
	return *this;
}

// комплексное сопряжение
cmplx conjg(cmplx c) { return cmplx(c.real, -c.image); }
cmplx conjg(double real, double image) { return cmplx(real, -image); }

void Быстрое_Фурье_преобразование(cmplx* F, long v_size, int is) // БДВПФ
{
	cmplx  temp, w, c;
	long i, i1, j, j1, istep;
	long m, mmax;
	long n = v_size;
	double fn, r1, theta;
	fn = (double)n;
	double r = M_PI * is;
	j = 1;
	for (i = 1; i <= n; i++)
	{
		i1 = i - 1;
		if (i < j)
		{
			j1 = j - 1;
			temp = F[j1];
			F[j1] = F[i1];
			F[i1] = temp;
		}
		m = n / 2;
		while (j > m) { j -= m;	m = (m + 1) / 2; }
		j += m;
	}
	mmax = 1;
	while (mmax < n)
	{
		istep = 2 * mmax;
		r1 = r / (double)mmax;
		for (m = 1; m <= mmax; m++)
		{
			theta = r1 * (double)(m - 1);
			w = cmplx(cos(theta), sin(theta));
			for (i = m - 1; i < n; i += istep)
			{
				j = i + mmax;
				c = F[j];
				temp = w * c;
				F[j] = F[i] - temp;
				F[i] = F[i] + temp;
			}
		}
		mmax = istep;
	}
	if (is > 0)  for (i = 0; i < n; i++) { F[i].real /= fn;  F[i].image /= fn; }
	return;
}

double Sum12() {
    static std::mt19937 rng(static_cast<unsigned int>(std::time(0)));
    std::uniform_int_distribution<int> dist(-32768, 32767);
    double sum = 0.0;
    for (int i = 0; i < 12; i++) {
        sum += static_cast<double>(dist(rng));
    }
    return sum;
}
