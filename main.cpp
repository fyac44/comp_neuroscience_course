/**
 * 
 * SPIKE-TRIGGERED AVERAGE CALCULATION
 * 
 * Requierements:
 *  - C++11 standard
 *  - gnuplot ($ brew install gnuplot) for plotting the results.
 *  - gnuplot-cpp from https://code.google.com/archive/p/gnuplot-cpp/downloads
 * 
 * Quiz 2 from the Computational Neuroscience Course in Coursera developed in C++.
 * 
 * The data provided in the quiz was dumped into the text file "data.txt"
 * 
 * Usage (with g++ compiler):
 * g++ main.cpp -std=c++11 (-I<PATH_TO_gnuplot-cpp_FOLDER>) 
 * 
 * Franklin Y. Alvarez - 29/04/2020
 * 
 **/

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Remember to include the path to this header file.
// Download from https://code.google.com/archive/p/gnuplot-cpp/downloads
#include "gnuplot_i.hpp"

using std::cout;
using std::endl;

void read_data(std::vector<double> &,
               std::vector<bool> &,
               const std::string);

void compute_sta(const std::vector<double> &,
                 const std::vector<bool> &,
                 double *,
                 const short);

void plot_sta(double *,
              short,
              short);

int main()
{
    std::vector<bool> rho;              // Spike response of the neuron
    std::vector<double> stim;           // Stimulation signal

    // Variables
    const short sampling_perdiod = 2;   // [us]
    const short w_len = 300;            // window lenght [us]
    short num_samples = w_len / sampling_perdiod;
    
    //  Initialize STA
    double * p_sta = new double[num_samples];
    for (int p = 0; p < num_samples; p++) p_sta[p] = 0.;
    
    // Read data from data.txt
    read_data(stim, rho, "data.txt");

    // Call the spike-triggered average function
    compute_sta(stim, rho, p_sta, num_samples);

    // Plot the results
    plot_sta(p_sta, num_samples, sampling_perdiod);

    delete[] p_sta;

    return 0;
}

void read_data(std::vector<double> & stim,
               std::vector<bool> & rho,
               const std::string filename)
/*----------------------------------------------------------------------------
Function to fill the stimulus (stim) and the spikes (rho) generate in a single
neuron from (filename).

Each line in filename represents a time step. They should be formatted as:

<data_rho data_stim>

data_rho is a boolean that represents wether a spike occurred (1) or not (0)
data_stim is a double with the stimulus value at the corresponding time step
-----------------------------------------------------------------------------*/
{
    std::ifstream data_file (filename);
    std::string line;
    bool data_rho;
    double data_stim;
    std::stringstream data_linestream;
    while (getline(data_file, line)) 
    {
        std::stringstream data_linestream(line);
        data_linestream >> data_rho;
        rho.push_back(data_rho);
        data_linestream >> data_stim;
        stim.push_back(data_stim);
    }
}

void compute_sta(const std::vector<double> & stim, 
                 const std::vector<bool> & rho,
                 double * sta,
                 const short num_samples)
/*----------------------------------------------------------------------------
Function to compute the spike-triggered average (sta) from the stimulus (stim) 
and spikes (rho).

The size of "sta" is given by the number of samples (num_samples) that will be
taken into account for the average.

Kindly reminder: The "sta" is the average stimulus previous to a spike.
-----------------------------------------------------------------------------*/
{
    // Generate a vector with the location of each spikes (indexes).
    // Start looking for spikes when the minimum length of the previous stimulus
    // is granted.
    std::vector<int> spike_indexes;
    int signal_size = stim.size();
    int index = num_samples;
    while (index < signal_size)
    {
        auto itr_index = std::find(rho.begin() + index, rho.end(), true);
        index = std::distance(rho.begin(), itr_index);
        spike_indexes.push_back(index);
        index++;
    }

    // Average the previous stimulus for each spike
    int num_spikes = spike_indexes.size();
    cout << "Spikes found = " << num_spikes << endl;    
    for (int i : spike_indexes)
    {
        for (int k = 0; k < num_samples; k++)
        {
            sta[k] = sta[k] + stim[i-num_samples+k]/num_spikes;
        }
    }
}

void plot_sta(double * sta, short num_samples, short sampling_period)
/*----------------------------------------------------------------------------
Function that uses the gnuplot tools to plot the STA through the gnuplot-cpp
functions.

The x-axis are negative numbers representing the time before a spike.
The y-axis is the average stimulus level at that moment for each spike.
-----------------------------------------------------------------------------*/
{
    // Define x and y axis
    std::vector<double> x(num_samples), y(num_samples);
    for (int i = 0; i < num_samples; i++)
    {
        x[i] = sampling_period * i - sampling_period * num_samples;
        y[i] = sta[i];
    }

    // Usage of the gnuplot-cpp functions
    Gnuplot::set_terminal_std("qt");
    Gnuplot g1("lines");
    g1.set_title("Spike-Triggered Average");
    g1.set_xlabel("Time [us]");
    g1.set_ylabel("Stimulus level");
    g1.plot_xy(x, y, "STA");

    // To not terminate the program before enjoying the plot!
    cout << "Press ENTER to exit: ";
    getchar();

}
