{
  'targets': [
    {
      'target_name': 'addon',
      'sources': [
      '../tmp/linear_model.cc',
      '../tmp/particle.cc',
      '../tmp/fastslam.cc',
      'addon.cc'
      ],
      'include_dirs': [
                    'src',
                    '.',
                    '/opt/intel/compilers_and_libraries/linux/mkl/include/',
                    '/usr/include/hdf5/serial'
      ],
      "link_settings": {
        "libraries": [
          "-lmkl_core",
          "-lmkl_def",
          "-lmkl_intel_thread",
          "-lmkl_core",
          "-liomp5",
          "-larmadillo",
          "-lm",
          "-ldl",
          "-lhdf5_cpp",
          "-lhdf5",
        #  "-lfslam"
      ],
        "ldflags": [
            "-L.",
            "-L/app",
            "-L/usr/local/hdf5/lib",
            "-L./build/Release"
            "-L./build/Release/obj.target",
            "-L/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin",
            "-L/opt/intel/compilers_and_libraries/linux/lib/intel64",
            "-Wl,-rpath=./build/Release",
            "-Wl,-rpath=/app",
            "-Wl,-rpath,@loader_path",
            "-Wl,-rpath=/usr/local/lib64",
            "-Wl,-rpath=/usr/local/lib",
            "-Wl,-rpath=/usr/lib",
            "-Wl,-rpath=/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin",
            "-Wl,-rpath=/opt/intel/compilers_and_libraries/linux/lib/intel64"
        ]
      },
    'configurations': {
      'Release': {
        'cflags_cc': [
          '-O3',
          '-std=c++17',
          '-m64'
        ],
        'cflags_cc!': [
          '-fno-rtti',
          '-pthread',
          '-fno-omit-frame-pointer',
          '-fno-exceptions',
          '-std=gnu++1y',
         ]
      }
    }
    }
  ]
}
