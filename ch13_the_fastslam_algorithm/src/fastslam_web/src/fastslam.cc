#include <iomanip>
#include "particle.h"
#include <uv.h>
#define NAPI_EXPERIMENTAL
#include <assert.h>
#include <node_api.h>

typedef struct {
  uv_thread_t thread;
  napi_threadsafe_function tsfn;
} ThreadData;

extern std::map <uint32_t, char> m;

void fastslam(uint n, ThreadData* data){
  std::string filename1 = "aa3_lsr2.h5";
  std::string filename2 = "aa3_dr.h5";
  std::string filename3 = "measurement.h5";
  std::string filename4 = "timeLsr.bin";
  std::string include_dir = "/app/data/";

  std::string filename5 = include_dir + "poses.dat";
  std::ofstream ostrm(filename5, std::ios::binary);

  std::string filename6 = include_dir + "ldmarks.dat";
  std::ofstream ostrm1(filename6, std::ios::binary);

  arma::mat speed, steering, time,  measurmt;
  uvec timelsr;
  speed.load(arma::hdf5_name(include_dir + filename2, "speed"));
  steering.load(arma::hdf5_name(include_dir + filename2, "steering"));
  time.load(arma::hdf5_name(include_dir + filename2, "time"));
  timelsr.load(include_dir + filename4);
  measurmt.load(arma::hdf5_name(include_dir + filename3, "zz"));


  static float const dt = 25e-3;
  Particle particles(n);

  uint stindex = 0;
  uvec jvec = find(timelsr >= time(stindex),1,"first");
  uint llsr = timelsr.n_elem;
  uint j = jvec(0);

  //mat::fixed<5000,2> poses = zeros<mat>(5000,2);
  //uint max_threads = mkl_get_max_threads();
  //mkl_set_num_threads(max_threads);

  /*std::chrono::steady_clock::time_point time_;
  std::chrono::steady_clock::time_point start_time_;
  std::chrono::steady_clock::time_point start_iter_;
  float deltatime_(0.0f);
  float deltatime_iter_(0.0f);
  start_time_ = std::chrono::steady_clock::now();*/

  for (size_t i = stindex; i < time.n_elem; i++) {
    //start_iter_ = std::chrono::steady_clock::now();

    mat z = zeros<mat>(2,0);
    while (j < llsr && timelsr(j) < time(i)+dt) {
      z = join_horiz(z,measurmt(uvec({2*j,2*j+1}),find(measurmt(2*j,span::all) != 0)));
      ++j;
    }

    particles.motion(speed(i),steering(i),dt,z);

    /*time_ = std::chrono::steady_clock::now();
    deltatime_ = std::chrono::duration_cast<std::chrono::duration<float>>(time_-start_time_).count();
    deltatime_iter_ = std::chrono::duration_cast<std::chrono::duration<float>>(time_-start_iter_).count();
    std::cout << "iter: " << std::setw(10)  << i
          << std::setw(15) << "iter time: " << std::setw(10) << deltatime_iter_
          << std::setw(15) << "avg: " << std::setw(10) << deltatime_/(i+1) <<
          '\r' ;*/

    //ostrm << particles << '\n';

    // the call into javascript is queued to be handled by the main event thread.
    float* poses = (float*) malloc(2 * n * sizeof(float));
    particles.cpy(poses);
    assert(napi_call_threadsafe_function(data->tsfn,
                                               (void*)poses,
                                               napi_tsfn_blocking) == napi_ok);

    /*if (i%1000==0) {
      //ostrm1.seekp(0,std::ios::beg);
      //particles.print_tree(ostrm1);
    }*/

    if (i%200 == 0) {
      //this has to be made thread safe.
      {
        if (m[(uint32_t) data->thread] == 0) {
          m.erase((uint32_t) data->thread);
          return;
        }
      }
    }


  }
    //ostrm1.seekp(0,std::ios::beg);
    //particles.print_tree(ostrm1);


  std::cout << '\n';
}
