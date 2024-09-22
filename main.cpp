#include <iostream>
#include <iomanip>
#include <thread>
#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <windows.h>

using namespace boost::multiprecision;

// Define the precision
using BigFloat = number<cpp_dec_float<100000>>;

struct GaussLegendreParams {
    BigFloat a;
    BigFloat b;
    BigFloat t;
    BigFloat p;
    unsigned int iterations;
    BigFloat pi;
};

// Function to perform a chunk of the Gauss-Legendre algorithm
void gauss_legendre_chunk(GaussLegendreParams& params) {
    BigFloat a = params.a;
    BigFloat b = params.b;
    BigFloat t = params.t;
    BigFloat p = params.p;

    for (unsigned int i = 0; i < params.iterations; ++i) {
        BigFloat a_next = (a + b) / BigFloat(2);
        BigFloat b_next = sqrt(a * b);
        BigFloat t_next = t - p * (a - a_next) * (a - a_next);
        a = a_next;
        b = b_next;
        t = t_next;
        p *= BigFloat(2);
    }

    // Update the shared pi value
    params.pi = (a + b) * (a + b) / (BigFloat(4) * t);
}

int main() {
    // Initialize the parameters
    BigFloat a = BigFloat(1);
    BigFloat b = BigFloat(1) / sqrt(BigFloat(2));
    BigFloat t = BigFloat(1) / BigFloat(4);
    BigFloat p = BigFloat(1);
    unsigned int total_iterations = 500; // Increase for better precision
    unsigned int num_threads = std::thread::hardware_concurrency(); // Use available threads
    unsigned int iterations_per_thread = total_iterations / num_threads;

    // Vector to hold thread parameters
    std::vector<GaussLegendreParams> params(num_threads);
    std::vector<std::thread> threads;

    // Start measuring time
    LARGE_INTEGER start, end, frequency;
    QueryPerformanceFrequency(&frequency); // Get the frequency of the high-resolution performance counter
    QueryPerformanceCounter(&start);

    // Create threads
    for (unsigned int i = 0; i < num_threads; ++i) {
        params[i] = { a, b, t, p, iterations_per_thread, BigFloat(0) };
        threads.emplace_back(gauss_legendre_chunk, std::ref(params[i]));
    }

    std::cout << "Number of created threads: " << num_threads << std::endl;
    std::cout << "Calculating PI, please wait..." << std::endl;

    // Join threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Combine results
    BigFloat pi = BigFloat(0);
    for (const auto& param : params) {
        pi += param.pi;
    }
    pi /= num_threads; // Average the results

    // Stop measuring time
    QueryPerformanceCounter(&end);

    // Calculate elapsed time
    double elapsedTime = static_cast<double>(end.QuadPart - start.QuadPart) /
        static_cast<double>(frequency.QuadPart);

    // Print the result
    std::cout << std::setprecision(100000) << std::fixed; // Set precision for output
    std::cout << "Calculated value of PI: " << pi << std::endl;

    std::cout << std::setprecision(10) << std::fixed;
    std::cout << "Execution time: " << elapsedTime << " seconds" << std::endl;

    return 0;
}
