#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <pthread.h> 
#include <unistd.h>
#include <sched.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0

int num_threads;
int num_processors;
pthread_t* tid; // array of thread ids


double* band_power;


struct thread_inputs {
  signal* sig;
  int thread;
  double bandwidth;
  int filter_order;
  double* filter_coeffs;
  double* band_power;
  int start_band;
  int blocksize;
};

void* worker(void* arg);

void usage() {
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads num_processors\n");
}

double avg_power(double* data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double* data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double* data, int num) {

  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double* data, int num) {

  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc;
  }
}


int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub) {
  double Fc        = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data,sig->num_samples);

  double signal_power = avg_power(sig->data,sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart,THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  double filter_coeffs[filter_order + 1];
  //double band_power[num_bands];

  //int parallel = 1;
    tid         = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);
    if (!tid) {
      fprintf(stderr, "cannot allocate memory\n");
      exit(-1);
    }


    band_power = (double *)malloc(num_bands * sizeof(double));
    
    if (num_bands >= num_threads) {
      for (int i = 0; i < num_threads; i++ ) {
        // input block
        struct thread_inputs* cur_thread_inputs = (struct thread_inputs*) malloc(sizeof(struct thread_inputs));
        cur_thread_inputs->sig = sig;
        cur_thread_inputs->thread = i;
        cur_thread_inputs->bandwidth = bandwidth;
        cur_thread_inputs->filter_order = filter_order;
        cur_thread_inputs->filter_coeffs = filter_coeffs;
        cur_thread_inputs->blocksize = num_bands / num_threads;
        cur_thread_inputs->start_band = cur_thread_inputs->blocksize * i;
        if (i == num_threads - 1) {
          cur_thread_inputs->blocksize += num_bands % num_threads;
        }
        
        
        int returncode = pthread_create(&(tid[i]),  // thread id gets put here
                                        NULL, // use default attributes
                                        worker, // thread will begin in this function
                                        (void*) cur_thread_inputs // we'll give it i as the argument
                                        );
        if (returncode != 0) {
          perror("Failed to start thread");
          exit(-1);
        }
      }

      for (int i = 0; i < num_threads; i++) {
        int returncode = pthread_join(tid[i], NULL);
        if (returncode != 0) {
          perror("join failed");
          exit(-1);
        }
      }
    }
  else {
    for (int i = 0; i < num_bands; i++) {
      struct thread_inputs* cur_thread_inputs = (struct thread_inputs*) malloc(sizeof(struct thread_inputs));
      cur_thread_inputs->sig = sig;
      cur_thread_inputs->thread = i;
      cur_thread_inputs->bandwidth = bandwidth;
      cur_thread_inputs->filter_order = filter_order;
      cur_thread_inputs->filter_coeffs = filter_coeffs;
      cur_thread_inputs->blocksize = 1;
      cur_thread_inputs->start_band = i;

      int returncode = pthread_create(&(tid[i]),  // thread id gets put here
                                        NULL, // use default attributes
                                        worker, // thread will begin in this function
                                        (void*) cur_thread_inputs // we'll give it i as the argument
                                        );
      if (returncode != 0) {
        perror("Failed to start thread");
        exit(-1);
      }
    }

    for (int i = 0; i < num_bands; i++) {
      int returncode = pthread_join(tid[i], NULL);
      if (returncode != 0) {
        perror("join failed");
        exit(-1);
      }
    }
  }


  
  


  unsigned long long tend = get_cycle_count();
  double end = get_seconds();

  resources rend;
  get_resources(&rend,THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

  for (int band = 0; band < num_bands; band++) {
    double band_low  = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n\
User time        %lf seconds\n\
System time      %lf seconds\n\
Page faults      %ld\n\
Page swaps       %ld\n\
Blocks of I/O    %ld\n\
Signals caught   %ld\n\
Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  return wow;
}

int main(int argc, char* argv[]) {

  if (argc != 8) {
    usage();
    return -1;
  }

  char sig_type    = toupper(argv[1][0]);
  char* sig_file   = argv[2];
  double Fs        = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  int num_bands    = atoi(argv[5]);
  num_threads = atoi(argv[6]);
  num_processors = atoi(argv[7]);

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n\
file:     %s\n\
Fs:       %lf Hz\n\
order:    %d\n\
bands:    %d\n",
         sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file,
         Fs,
         filter_order,
         num_bands);

  printf("Load or map file\n");

  signal* sig;
  switch (sig_type) {
    case 'T':
      sig = load_text_format_signal(sig_file);
      break;

    case 'B':
      sig = load_binary_format_signal(sig_file);
      break;

    case 'M':
      sig = map_binary_format_signal(sig_file);
      break;

    default:
      printf("Unknown signal type\n");
      return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end   = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}

void* worker(void* arg) {
  struct thread_inputs* args = (struct thread_inputs* ) arg;
  int blocksize = args->blocksize;
  int start = args->start_band;
  int end = start + blocksize;  //int blocksize = arg->num_bands / num_threads; // note: floor
  signal* sig = args->sig;
  int thread = args->thread;
  double bandwidth = args-> bandwidth;
  int filter_order = args->filter_order;
  double filter_coeffs[filter_order + 1];
  //double* band_power = args->band_power;



  

  // put ourselves on the desired processor
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(args->thread % num_processors, &set);
  if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
    perror("Can't setaffinity"); // hopefully doesn't fail
    exit(-1);
  }

  // This figures out the chunk of the vector I should
  // work on based on my id
  //int mystart =  * blocksize;



  for (int band = start; band < end; band++) {
    // Make the filter
    generate_band_pass(sig->Fs,
                       band * bandwidth + 0.0001, // keep within limits
                       (band + 1) * bandwidth - 0.0001,
                       filter_order,
                       filter_coeffs);
    hamming_window(filter_order,filter_coeffs);

    // Convolve
    // convolve_and_compute_power(sig->num_samples,
    //                            sig->data,
    //                            filter_order,
    //                            filter_coeffs,
    //                            &(band_power[band]));
    convolve_and_compute_power(sig->num_samples,
                               sig->data,
                               filter_order,
                               filter_coeffs,
                               &(band_power[band]));

 


  }

  // Now I sum that chunk and put the result in partial_sum
  // partial_sum[myid] = 0.0;
  // for (int i = mystart; i < myend; i++) {
  //   partial_sum[myid] += vector[i];
  // }

  // Done.  The master thread will sum up the partial sums
  pthread_exit(NULL);           // finish - no return value

}




