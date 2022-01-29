#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <thread>
#include <signal.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Wrapped C functions cannot be stopped in python with CTRL-C.
// Use signal(SIGINT, signalHandler) in functions that have long execution time.
void signalHandler(int signum) {std::cerr << "Killed by the user.\n"; exit(signum);}


double interpolate(double *yData, size_t size, double x) {
       
    int i = int(std::floor(x));
    if (x>=size-2)  i = size-2;
    double xL = i, yL = yData[i];
    double xR = i+1, yR = yData[i+1];
    double dydx = (yR-yL) / (xR-xL);
    return yL + dydx * (x-xL);
}
    
extern "C"
{

float *grid_otf_C(float *data, int nspectra, int nchan,                  // Data has dimensions (nchan, nspectra)
                   double *weights, double* x_pix, double *y_pix,        // These have dimension (nspectra)
                   int nx, int ny,                                       // X and Y size of datacube 
                   const char* kernel,
                   double r_support_pix_sqrd, double pre_delta_pix_sqrd,
                   double *pre_conv_fn, int pre_conv_size,
                   double cap_dist_sqrd_pix, double max_conv_fn,
                   double cutoff_conv_fn, int threads) 
{
    
    signal(SIGINT, signalHandler);
    
    std::string kern = kernel;
    float *data_cube = new float[nchan*ny*nx];
    for (auto i=0; i<nchan*ny*nx; i++) data_cube[i] = 0;
    
    // Starting main loop over all spectra
    std::cout << "Processing " << nspectra << " spectra ...    " << std::flush;
        
    int count = 0;
#pragma omp parallel for num_threads(threads) shared(count) schedule(dynamic)
    for (auto i=0; i<nx*ny; i++) {
        
        // Updating loop progress
        count++;
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        if(tid==0) {
            for(int i=0; i<4; i++) std::cout << '\b';
            int per = lround(float(count)/(float(nx*ny))*100);
            std::cout << std::setw(3) << per << "%" << std::flush;
        }
        
        std::vector<size_t> keep;
        std::vector<double> combined_weight;
        double coverage=0, wsum=0;

        int x = i % nx;
        int y = int(i/nx);
        
        for (auto k=0; k<nspectra; k++) {
            double xdist = x_pix[k]-x;
            double ydist = y_pix[k]-y;
                
            if (kern=="nearest") {
                if ((xdist>=-0.5) && (xdist<0.5) && (ydist>=-0.5) && (ydist<0.5)) {
                    keep.push_back(k);
                    coverage += 1;
                    combined_weight.push_back(weights[k]);
                    wsum += weights[k];
                }
            }
            else {
                double pdist_sqrd = xdist*xdist + ydist*ydist;
                if (pdist_sqrd<=r_support_pix_sqrd) {
                    keep.push_back(k);

                    double c = max_conv_fn;
                    if (pdist_sqrd>=cap_dist_sqrd_pix) {
                        // INTERPOLATE TO GET THE CONV. FN. FOR EACH DATA POINT
                        double pre_grid_x = pdist_sqrd / pre_delta_pix_sqrd;
                        c = interpolate(pre_conv_fn,pre_conv_size,pre_grid_x);
                    }
                    coverage += c;
                    combined_weight.push_back(c*weights[k]);
                    wsum += c*weights[k];
                }
            }
        }

        if (keep.size()==0) continue;

        // PLACE A MINIMUM THRESHOLD NEEDED TO CONSIDER A GRID POINT
        if (coverage>cutoff_conv_fn && wsum>0) {
            for (auto k=0; k<keep.size(); k++) {
                for (auto z=0; z<nchan; z++) {
                    double dataval = data[z+keep[k]*nchan];
                    if (dataval==dataval) 
                        data_cube[i+z*nx*ny] += dataval*combined_weight[k]/wsum;
                }
            }
        }
    }
    
    for(int i=0; i<4; i++) std::cout << '\b';
    std::cout << " done!" << std::endl;
    
    return data_cube;
}

}
