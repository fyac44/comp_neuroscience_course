/**
 * 
 * GAUSSIAN DISTRIBUTION CLASS
 * 
 * Requierements:
 *  - gnuplot ($ brew install gnuplot) for plotting the results.
 *  - gnuplot-cpp from https://code.google.com/archive/p/gnuplot-cpp/downloads
 * 
 * The class is a representation of the Gaussian distribution funtion of the
 * firing rate in a neuron given a specific stimulus. The function is weighted
 * by the prior probability of the stimulus ocurrence and the cost of guessing
 * mistakingly such stimulus.
 * 
 * This is developed for the Computational Neuroscience Course at:
 * https://www.coursera.org/learn/computational-neuroscience
 * 
 * Franklin Y. Alvarez - 05/05/2020
 * 
 **/

#ifndef GAUSS_FYAC_H_
#define GAUSS_FYAC_H_

#define DEFAULT_MEAN 0.0
#define DEFAULT_SD 1.0
#define DEFAULT_PRIOR 0.5
#define DEFAULT_LOSS_COST 1.0
#define PI 3.14159265358979323846
#define PLOT_SAMPLES 5000
#define PLOT_SD 5.0

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "gnuplot_i.hpp"

class Gauss
{
private:

    static int class_count;
    int id_;
	float mean_;
    float sd_;
    float prior_;
    float loss_cost_;


    // Functions
    void identify();


public:
    // Constructors
	Gauss();
	Gauss(const float, const float);
    Gauss(const float, const float, const float, const float);
	Gauss(const Gauss &other);

    // Destructor
    virtual ~Gauss();

    // Overload
	Gauss& operator=(const Gauss &other);

    // Setters
	void set_mean(const float);
	void set_stddev(const float);
    void set_id(const int);
    void set_prior(const float);
    void set_losscost(const float);
    
    // Getters
    float get_mean() const;
	float get_stddev() const;
    int get_id() const;
    float get_prior() const;
    float get_losscost() const;

    // Functions
    void dump() const;
    float evaluate(const float) const;
    void gnu_plot(Gnuplot &) const;
};

#endif