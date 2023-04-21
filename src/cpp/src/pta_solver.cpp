#include "pta_solver.h"

#include <chrono>
#include <condition_variable>
#include <ctime>
#include <functional> // for _1, _2
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>

//
// Return a timestamp std::string to we can see how long things take
//
std::string timestamp(void)
{
    auto now = std::chrono::system_clock::now();
    auto seconds = std::chrono::time_point_cast<std::chrono::seconds>(now);
    auto mseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - seconds);
    auto date = std::chrono::system_clock::to_time_t(now);

    struct tm local_time;
    localtime_r(&date, &local_time);

    char buffer[128];
    buffer[0] = '\0';
    auto buffer_size = sizeof(buffer) - 1;
    auto out = strftime(buffer, buffer_size, "%H:%M:%S", &local_time);
    out += snprintf(buffer + out, buffer_size - out, ".%03d ", (int)mseconds.count());

    return std::string(buffer);
}

template <typename ret, typename T, typename... Rest>
using fn = std::function<ret(T, Rest...)>;

template <typename ret, typename T, typename... Rest>
ret wrap_my_slow_function(fn<ret, T, Rest...> f, T t, Rest... rest)
{
    std::mutex my_mutex;
    std::condition_variable my_condition_var;
    ret result = 0;

    std::unique_lock<std::mutex> my_lock(my_mutex);

    //
    // Spawn a thread to call my_function_that_might_block().
    // Pass in the condition variables and result by reference.
    //
    std::thread my_thread([&]()
                          {
    result = f(t, rest...);
    // Unblocks one of the threads currently waiting for this condition.
    my_condition_var.notify_one(); });

    //
    // Detaches the thread represented by the object from the calling
    // thread, allowing them to execute independently from each other. B
    //
    my_thread.detach();

    if (my_condition_var.wait_for(my_lock, std::chrono::seconds(1)) == std::cv_status::timeout)
    {
        //
        // Throw an exception so the caller knows we failed
        //
        // Timed out at       :
        throw std::runtime_error("Timeout");
    }

    return result;
}

PTASolver::PTASolver(double bulk_fugacity, double temperature, call_mono_eos eos)
{
    this->EquationOfState = eos;
    this->BulkFugacity = bulk_fugacity;
    this->Temperature = temperature;
}

PTASolver::PTASolver(double bulk_fugacity, double temperature, call_mix_eos eos)
{
    this->MixEquationOfState = eos;
    this->BulkFugacity = bulk_fugacity;
    this->Temperature = temperature;
}

double PTASolver::equilibrium(double P,
                              double f_eps)
{
    if (std::isinf(P) == 1 || std::isnan(P) == 1)
    {
        return 10000.;
    }
    struct mono_eos adsorbed = this->EquationOfState(P, this->Temperature);
    return this->BulkFugacity * f_eps - adsorbed.fug;
}

double PTASolver::findOptimizedPressure(double InitialPressureEstimate, double f_eps)
{

    auto equilibrium = [&](double x)
    {
        return this->equilibrium(x, f_eps);
    };

    auto fnn = [&](double P)
    {
        return brent_zeroin(equilibrium, P, this->DEFAULT_TOL);
    };

    try
    {
        auto f1 = fn<double, double>(fnn);
        return wrap_my_slow_function(f1, InitialPressureEstimate);

        // Success, no timeout
    }
    catch (std::runtime_error &e)
    {
        //
        // Do whatever you need here upon timeout failure
        //
        // throw
        std::cout << "Timeout... re-running with increase of 10% on initial Pz" << std::endl;
    }

    return PTASolver::findOptimizedPressure(InitialPressureEstimate * 1.1, f_eps);
}