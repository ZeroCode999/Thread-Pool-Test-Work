#include <iostream>
#include <thread>
#include <chrono>
#include <vector>

#include "MyThreadPool.h"

using namespace std;
using namespace std::chrono;

// standard matrix multiplication
void mult_v0(int M, int N, int K, const double * A, const double * B, double * C)
{
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			C[i * N + j] = 0;
			for (int k = 0; k < K; ++k)
				C[i * N + j] += A[i * K + k] * B[k * N + j];
		}
	}
}

void mult_v1(int M, int N, int K, const double * A, const double * B, double * C)
{
	for (int i = 0; i < M; ++i)
	{
		// remove the constant part from the inner loop
		double * c = C + i * N;
		for (int j = 0; j < N; ++j)
			c[j] = 0;
		// sequential line-by-line traversal for all three matrices
		for (int k = 0; k < K; ++k)
		{
			const double * b = B + k * N;
			double  a = A[i * K + k];
			for (int j = 0; j < N; ++j)
				c[j] += a * b[j];
		}
	}
}

void mult_v1_sub_iter4(int shift, int M, int N, int K, const double* A, const double* B, double* C);

void mult_v1_threaded4(int M, int N, int K, const double* A, const double* B, double* C)
{
	thread threadCalc0(mult_v1_sub_iter4, 0, M, N, K, A, B, C);
	thread threadCalc1(mult_v1_sub_iter4, 1, M, N, K, A, B, C);
	thread threadCalc2(mult_v1_sub_iter4, 2, M, N, K, A, B, C);
	thread threadCalc3(mult_v1_sub_iter4, 3, M, N, K, A, B, C);
	threadCalc0.join();
	threadCalc1.join();
	threadCalc2.join();
	threadCalc3.join();

}

void mult_v1_sub_iter4(int shift, int M, int N, int K, const double* A, const double* B, double* C)
{
	for (int i = shift; i < M; i += 4)
	{
		double* c = C + i * N;
		for (int j = 0; j < N; ++j)
			c[j] = 0;
		for (int k = 0; k < K; ++k)
		{
			const double* b = B + k * N;
			double  a = A[i * K + k];
			for (int j = 0; j < N; ++j)
				c[j] += a * b[j];
		}
	}
}

void mult_v1_sub_iter8(int shift, int M, int N, int K, const double* A, const double* B, double* C);

void mult_v1_threaded8(int M, int N, int K, const double* A, const double* B, double* C)
{
	vector<thread> threads;
	for (int i = 0; i < 8; i++)
	{
		threads.push_back(thread(mult_v1_sub_iter8, i, M, N, K, A, B, C));
	}
	for (int i = 0; i < 8; i++)
	{
		threads[i].join();
	}
}


void mult_v1_sub_iter8(int shift, int M, int N, int K, const double* A, const double* B, double* C)
{
	for (int i = shift; i < M; i += 8)
	{
		double* c = C + i * N;
		for (int j = 0; j < N; ++j)
			c[j] = 0;
		for (int k = 0; k < K; ++k)
		{
			const double* b = B + k * N;
			double  a = A[i * K + k];
			for (int j = 0; j < N; ++j)
				c[j] += a * b[j];
		}
	}
}

void mult_v1_sub_iter_pool4(int shift, int M, int N, int K, const double* A, const double* B, double* C);

void mult_v1_threadpool4(int M, int N, int K, const double* A, const double* B, double* C)
{
	MyThreadPool pool(4);

	for (int i = 0; i < 4; i++)
	{
		pool.add_task(mult_v1_sub_iter_pool4, i, M, N, K, A, B, C);
	}

}

void mult_v1_sub_iter_pool4(int shift, int M, int N, int K, const double* A, const double* B, double* C)
{
	for (int i = shift; i < M; i += 4)
	{
		double* c = C + i * N;
		for (int j = 0; j < N; ++j)
			c[j] = 0;
		for (int k = 0; k < K; ++k)
		{
			const double* b = B + k * N;
			double  a = A[i * K + k];
			for (int j = 0; j < N; ++j)
				c[j] += a * b[j];
		}
	}
}

void mult_v1_sub_iter_pool8(int shift, int M, int N, int K, const double* A, const double* B, double* C);

void mult_v1_threadpool8(int M, int N, int K, const double* A, const double* B, double* C)
{
	MyThreadPool pool(8);

	for (int i = 0; i < 8; i++)
	{
		pool.add_task(mult_v1_sub_iter_pool8, i, M, N, K, A, B, C);
	}

}

void mult_v1_sub_iter_pool8(int shift, int M, int N, int K, const double* A, const double* B, double* C)
{
	for (int i = shift; i < M; i += 8)
	{
		double* c = C + i * N;
		for (int j = 0; j < N; ++j)
			c[j] = 0;
		for (int k = 0; k < K; ++k)
		{
			const double* b = B + k * N;
			double  a = A[i * K + k];
			for (int j = 0; j < N; ++j)
				c[j] += a * b[j];
		}
	}
}

void solve()
{
	const int N = 1000, M = 1000, K = 1000;

	double * A = new double[M*K];
	for (int i = 0; i < M*K; i++)
	{
		A[i] = rand() % 10000;
	}

	double * B = new double[K*N];
	for (int i = 0; i < K*N; i++)
	{
		B[i] = rand() % 10000;
	}

	cout << "mult_v0" << endl;
	double * C = new double[M*N];
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	mult_v0(M, N, K, A, B, C);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	cout << duration_cast<microseconds>(t2 - t1).count() << endl << endl;

	cout << endl;
	cout << endl;
	cout << "mult_v1" << endl;
	t1 = high_resolution_clock::now();
	mult_v1(M, N, K, A, B, C);
	t2 = high_resolution_clock::now();
	cout << duration_cast<microseconds>(t2 - t1).count() << endl << endl;

	cout << endl;
	cout << endl;
	cout << "mult_v1_threaded4" << endl;
	t1 = high_resolution_clock::now();
	mult_v1_threaded4(M, N, K, A, B, C);
	t2 = high_resolution_clock::now();
	cout << duration_cast<microseconds>(t2 - t1).count() << endl << endl;

	cout << endl;
	cout << endl;
	cout << "mult_v1_threadpool4" << endl;
	t1 = high_resolution_clock::now();
	mult_v1_threadpool8(M, N, K, A, B, C);
	t2 = high_resolution_clock::now();
	cout << duration_cast<microseconds>(t2 - t1).count() << endl << endl;

	cout << endl;
	cout << endl;
	cout << "mult_v1_threaded8" << endl;
	t1 = high_resolution_clock::now();
	mult_v1_threaded4(M, N, K, A, B, C);
	t2 = high_resolution_clock::now();
	cout << duration_cast<microseconds>(t2 - t1).count() << endl << endl;

	cout << endl;
	cout << endl;
	cout << "mult_v1_threadpool8" << endl;
	t1 = high_resolution_clock::now();
	mult_v1_threadpool8(M, N, K, A, B, C);
	t2 = high_resolution_clock::now();
	cout << duration_cast<microseconds>(t2 - t1).count() << endl << endl;

	delete[] A, B, C;
}

int main()
{
	solve();
}