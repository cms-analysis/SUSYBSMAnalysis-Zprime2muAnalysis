#include "EMuBackgroundStudy.h"

int main(int argc, char** argv) {
  if (argc < 1) {
    printf("usage: %s ntuple.root [max # of events to run over]\n", argv[0]);
    return 1;
  }
    
  int max_entries = argc > 2 ? atoi(argv[1]) : -1;
  EMuBackgroundStudy bs(argv[1], max_entries);

  bs.Loop();

  return 0;
}
