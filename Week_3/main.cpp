/**
 * 
 * THRESHOLD VALUE IN A BINARY DECISION
 * 
 * Requierements:
 *  - gnuplot ($ brew install gnuplot) for plotting the results.
 *  - gnuplot-cpp from https://code.google.com/archive/p/gnuplot-cpp/downloads
 *  - Gauss class "gauss.h" and "gauss.cpp" (in the same repository)
 * 
 * Two stimulus are assumed to produce spikes at different firing rates. The firing
 * rates are represented as a Gaussian distribution. To decide wheter the stimulus 1 
 * "s1" or stimulus 2 "s2" is occurring (binary decision), a threshold on the firing
 * rate response must be set. This program will calculate this threshold of decision.
 * 
 * The class "Gauss" represents such gaussian distributions plus the prior probability
 * and the cost of a bad prediction.
 * 
 * This is developed for the Computational Neuroscience Course at:
 * https://www.coursera.org/learn/computational-neuroscience
 * 
 * Usage (with g++ compiler):
 * g++ main.cpp -std=c++11 (-I<PATH_TO_gnuplot-cpp_FOLDER>) 
 * 
 * Franklin Y. Alvarez - 05/05/2020
 * 
 **/

#include <iostream>
#include "gauss.cpp"

// Remember to include the path to this header file.
// Download from https://code.google.com/archive/p/gnuplot-cpp/downloads
#include "gnuplot_i.hpp"

int main()
{
    float err = 0.00001;
    float thr1, thr2, u1, u2, sd1, sd2, cost1, cost2, prior1, prior2;

    // TODO: Validation of the input
    // Gauss distribution under first stimulus "s1"
    std::cout << "Please indicate the parameters of the gaussian distribution of the stimulus s1:" << std::endl;
    std::cout << " - Mean of s1: ";
    std::cin >> u1;
    std::cout << " - Standard deviation of s1: ";
    std::cin >> sd1;
    std::cout << " - Prior probability of s1[0 - 1]: ";
    std::cin >> prior1;
    std::cout << " - Cost of mistakingly guessing s1: ";
    std::cin >> cost1;
    Gauss s1(u1, sd1, prior1, cost1);

    // Gaussian distribution under second stimulus "s2"
    std::cout << "Please indicate the parameters of the gaussian distribution of the stimulus s2:" << std::endl;
    std::cout << " - Mean of s2: ";
    std::cin >> u2;
    std::cout << " - Standard deviation of s2: ";
    std::cin >> sd2;
    std::cout << " - Prior probability of s2[0 - 1]: ";
    std::cin >> prior2;
    std::cout << " - Cost of mistakingly guessing s2: ";
    std::cin >> cost2;
    Gauss s2(u2, sd2, prior2, cost2);

    // Print the gaussians in the terminal
    s1.dump();
    s2.dump();

    // To calculate the threshold using the gaussian distribution equation
    // f(x) = {k*exp[-(1/2)*((x-mean)/sd)^2]}/[sd*sqrt(2 * pi)]
    // k = prior / cost
    float k1 = s1.get_prior() / s1.get_losscost();
    float k2 = s2.get_prior() / s2.get_losscost();

    if (s1.get_stddev() == s2.get_stddev())
    // If the standard deviation is the same in both distributions there is
    // only 1 solution.
    {
        thr1 = ( 2*pow(s1.get_stddev(), 2)*log(k1/k2) + 
                 pow(s2.get_mean(), 2) - pow(s1.get_mean(), 2) ) /
               ( 2 * (s2.get_mean() - s1.get_mean()) );
        thr2 = thr1;
        std::cout << "The only threshold is located at " << thr1 << std::endl;
    }
    else
    // By solving the equation f1(x) = f2(x) a quadratic function is obtained
    // x1 = (- b + sqrt(b^2 - 4ac)) / (2 * a)
    // x2 = (- b - sqrt(b^2 - 4ac)) / (2 * a)
    {
        // a = sd2^2 - sd1^2
        float a = pow(s2.get_stddev(), 2) - pow(s1.get_stddev(), 2);

        // b = 2 * (sd1^2 * mean2 - sd2^2 * mean1)
        float b = 2 * ( pow(s1.get_stddev(), 2) * s2.get_mean() - 
                        pow(s2.get_stddev(), 2) * s1.get_mean() );
        
        // c = (sd2*mean1)^2 - (sd1*mean2)^2 - 2*ln(k1*sd2/(k2*sd1))*(sd1*sd2)^2
        float c = pow(s2.get_stddev()*s1.get_mean(),2) -
                  pow(s1.get_stddev()*s2.get_mean(),2) -
                  pow(s2.get_stddev(), 2)*pow(s1.get_stddev(), 2) *
                  2 * log(k1*s2.get_stddev()/(k2*s1.get_stddev()));
        
        float root = pow(b, 2) - 4*a*c;
        if (root < 0)
        // There is no solution
        {
            std::cout << "There is no possible threshold";
        }
        else if (root == 0)
        // Only one solution
        {
            thr1 = -b / (2 * a);
            thr2 = thr1;
            std::cout << "The only threshold is located at " << thr1 << std::endl;
        }
        else
        // There are two solutions
        {
            thr1 = (-b + sqrt(root)) / (2 * a);
            thr2 = (-b - sqrt(root)) / (2 * a);
            std::cout << "Threshold are located at:" << std::endl;
            std::cout << " - thr1 = " << thr1 << std::endl;
            std::cout << " - thr2 = " << thr2 << std::endl;
        }
        
    }

    // Plot both gaussians distributions
    Gnuplot::set_terminal_std("qt");
    Gnuplot g1("lines");
    g1.set_title("Binary decision between two stimulus");
    g1.set_xlabel("Firing rate");
    g1.set_ylabel("Probability");
    s1.gnu_plot(g1);
    s2.gnu_plot(g1);

    // To not terminate the program before enjoying the plot!
    std::cout << "PRESS ENTER TO EXIT" << std::endl;
    getchar();
    getchar();
    
    return 0;
}