#include "seif.h"

void seif(){
  std::string filename1 = "aa3_lsr2.h5";
  std::string filename2 = "aa3_dr.h5";
  std::string filename3 = "measurement.h5";
  std::string filename4 = "timeLsr.bin";
  std::string include_dir = "./data/";

  std::string filename5 = include_dir + "poses.dat";
  std::ofstream ostrm(filename5, std::ios::binary);

  arma::mat speed, steering, time,  measurmt;
  uvec timelsr;
  speed.load(arma::hdf5_name(include_dir + filename2, "speed"));
  steering.load(arma::hdf5_name(include_dir + filename2, "steering"));
  time.load(arma::hdf5_name(include_dir + filename2, "time"));
  timelsr.load(include_dir + filename4);
  measurmt.load(arma::hdf5_name(include_dir + filename3, "zz"));

  vec μ = zeros<vec>(3);
  vec ξ = zeros<vec>(3);
  mat Ω = 10e4*eye<mat>(3,3);
  std::vector<uint> m0;
  SpMat<ushort> Λ;

  static uint const n(20);
  static float const dt = 25e-3;

  uint stindex = 0;
  uvec jvec = find(timelsr >= time(stindex),1,"first");
  uint llsr = timelsr.n_elem;
  uint j = jvec(0);

  //mat::fixed<5000,2> poses = zeros<mat>(5000,2);

  //uint max_threads = mkl_get_max_threads();
  //mkl_set_num_threads(max_threads);

  std::chrono::steady_clock::time_point time_;
  std::chrono::steady_clock::time_point start_time_;
  std::chrono::steady_clock::time_point start_iter_;
  float deltatime_(0.0f);
  float deltatime_iter_(0.0f);
  start_time_ = std::chrono::steady_clock::now();


  for (size_t i = stindex; i < time.n_elem; i++) {
    start_iter_ = std::chrono::steady_clock::now();

    if (i > stindex && i%5000==0) {
      /*std::string filename = include_dir + "iter_" + std::to_string(i) + ".h5";
      poses.save(hdf5_name(filename,"poses"));
      μ.save(hdf5_name(filename, "μ"));
      ξ.save(hdf5_name(filename, "ξ"));
      Ω.save(hdf5_name(filename, "Ω"));
      Col<uint>(m0).save(hdf5_name(filename, "m0"));
      Λ.save(include_dir + "Λ_iter_" + std::to_string(i) + ".bin");*/

      /*std::string filename = include_dir + "iter_" + std::to_string(i) + ".dat";
      poses.save(include_dir + "poses_" + std::to_string(i) + ".dat");
      μ.save(include_dir + "μ_" + std::to_string(i) + ".dat");
      ξ.save(include_dir + "ξ_" + std::to_string(i) + ".dat");
      Ω.save(include_dir + "Ω_" + std::to_string(i) + ".dat");
      Col<uint>(m0).save(include_dir + "m0_" + std::to_string(i) + ".dat");
      Λ.save(include_dir + "Λ_" + std::to_string(i) + ".dat");*/
    }

    motion(speed(i),steering(i),dt,m0,μ,ξ,Ω,Λ);

    estimate(m0,μ,ξ,Ω);

    std::vector<uint> m3;
    std::vector<uint> dest;
    while (j <= llsr && timelsr(j) < time(i)+dt) {
      mat z1 = measurmt(uvec({2*j,2*j+1}),find(measurmt(2*j,span::all) != 0));
      if (!z1.empty()) {
        std::vector<uint> c;
        c = correspondence(z1,m0,μ,ξ,Ω,Λ);
        measurement(z1,c,μ,ξ,Ω,Λ);
        std::sort(c.begin(),c.end());
        std::set_union(m3.begin(),m3.end(),c.begin(),c.end(),std::back_inserter(dest));
        std::swap(m3,dest);
      }

      j++;
    }

    if (!m3.empty()) {

      std::vector<uint> m1, m4, m5, m6, m2(m0);

      std::sort(m2.begin(),m2.end());
      std::set_difference(m2.begin(),m2.end(),
            m3.begin(),m3.end(),
            std::inserter(m4,m4.begin()));
      auto last = std::unique(m3.begin(),m3.end());
      m4.insert(m4.end(),m3.begin(),last);
      m4.erase(m4.begin(),m4.begin()+std::max(0,int(m4.size()-n)));

      std::vector<uint> m8(m4);
      std::sort(m8.begin(),m8.end());

      std::set_difference(m2.begin(),m2.end(),
            m8.begin(),m8.end(),
            std::inserter(m1,m1.begin()));

      std::set_difference(m3.begin(),last,
            m8.begin(),m8.end(),
            std::inserter(m5,m5.begin()));

      std::set_union(m1.begin(),m1.end(),m5.begin(),m5.end(),std::back_inserter(m6));
      m0 = std::move(m4);
      if (!m6.empty()) {
        sparsification(m0,m6,μ,ξ,Ω,Λ);
      }
    }

    time_ = std::chrono::steady_clock::now();
    deltatime_ = std::chrono::duration_cast<std::chrono::duration<float>>(time_-start_time_).count();
    deltatime_iter_ = std::chrono::duration_cast<std::chrono::duration<float>>(time_-start_iter_).count();
    std::cout << "iter: " << std::setw(10)  << i <<
          "      iter time: " << std::setw(10) << deltatime_iter_ <<
          "      avg: " << std::setw(10) << deltatime_/(i+1) <<
          '\r' ;

    ostrm  << std::setw(15)  << μ(0) <<
              std::setw(15)  << μ(1) << '\n';

  }
  std::cout << '\n';
}
