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

#include "gauss.h"

// ***************** Initialize Static ID *********************
int Gauss::class_count = 1;

// ***************** Constructors *****************************
Gauss::Gauss()
{
	set_mean(DEFAULT_MEAN);
	set_stddev(DEFAULT_SD);
    set_prior(DEFAULT_PRIOR);
    set_losscost(DEFAULT_LOSS_COST);
    identify();
}

Gauss::Gauss(const float mean, const float sd)
{
	set_mean(mean);
	set_stddev(sd);
    set_prior(DEFAULT_PRIOR);
    set_losscost(DEFAULT_LOSS_COST);
    identify();
}

Gauss::Gauss(const float mean,
             const float sd,
             const float prior,
             const float cost) 
{
	set_mean(mean);
	set_stddev(sd);
    set_prior(prior);
    set_losscost(cost);
    identify();
}

Gauss::Gauss(const Gauss &other)
{
	*this = other;
}

// ***************** Destructors ******************************
Gauss::~Gauss() { }

// ***************** Overload *********************************
Gauss& Gauss::operator=(const Gauss &other)
{
	set_mean(other.get_mean());
	set_stddev(other.get_stddev());
    identify();
	return *this ;
}

// ***************** Setters **********************************
void Gauss::set_mean(const float mean) { mean_ = mean; }
void Gauss::set_stddev(const float sd) { sd_ = sd; }
void Gauss::set_id(const int id)       { id_ = id; }
void Gauss::set_prior(const float prior)
{
    if (prior >= 0.0 && prior <= 1.0) { prior_ = prior; }
    else { prior_ = DEFAULT_PRIOR; }
}
void Gauss::set_losscost(const float cost) { loss_cost_ = cost; }

// ***************** Getters **********************************
float Gauss::get_mean() const   { return mean_; }
float Gauss::get_stddev() const { return sd_; }
int Gauss::get_id() const { return id_; }
float Gauss::get_losscost() const { return loss_cost_; }
float Gauss::get_prior() const { return prior_; }

// ***************** Functions ********************************

void Gauss::identify()
/*----------------------------------------------------------------------------
Private function to enumarate the instances of the class Gauss
-----------------------------------------------------------------------------*/
{
    set_id(Gauss::class_count);
    Gauss::class_count++;
}

void Gauss::dump() const
/*----------------------------------------------------------------------------
Print in terminal all the properties of the instance
-----------------------------------------------------------------------------*/
{
	std::cout << "Gaussian #: " << get_id() << std::endl;
	std::cout << "  - Mean: " << get_mean() << std::endl;
	std::cout << "  - Standard deviation: " << get_stddev() << std::endl;
    std::cout << "  - Prior probabillity: " << get_prior() << std::endl;
    std::cout << "  - Loss cost: " << get_losscost() << std::endl;
}

void Gauss::gnu_plot(Gnuplot & g) const 
/*----------------------------------------------------------------------------
Uses an intance of gnuplot-cpp class to plot the gaussian distribution.

The distribution is represented for 5 standard deviations (PLOT_SD) in 5000
samples (PLOT_SAMPLES)
-----------------------------------------------------------------------------*/
{
    // Define x and y axis
    float min_val = get_mean() - PLOT_SD * get_stddev();
    float max_val = get_mean() + PLOT_SD * get_stddev();
    float step = ( max_val - min_val ) / PLOT_SAMPLES;
    std::vector<double> x(PLOT_SAMPLES), y(PLOT_SAMPLES);
    for (int i = 0; i < PLOT_SAMPLES; i++)
    {
        x[i] = min_val + (float) i * step;
        y[i] = evaluate(x[i]);
    }
    std::string name = "Gaussian distribution of s" + std::to_string(get_id());
    g.plot_xy(x, y, name);

}

float Gauss::evaluate(const float x) const
/*----------------------------------------------------------------------------
Evaluation of the Gaussian distribution function "f_x" in a given value "x".
-----------------------------------------------------------------------------*/
{
    float f_x;
    float mean = get_mean();
    float sd = get_stddev();

    f_x = exp(- 0.5 * pow( (x - mean) / sd , 2.0) ) / (sd * sqrt( 2.0 * PI));
    f_x = get_prior() * f_x / get_losscost();

    return f_x;
}



