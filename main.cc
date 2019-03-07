#include <string>
void fastslam(uint);

int main(int argc, char** argv){
  uint i=100;
  if (argc > 1) {
    i = (uint) std::stoi(argv[1]);
  }

  fastslam(i);

  return 0;
}
