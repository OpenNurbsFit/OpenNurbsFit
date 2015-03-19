

#include "TomGine/tgTomGineThread.h"

TomGine::tgTomGineThread viewer(800,600, "Example FitSurfaceDepth");

int main(int argc, char *argv[])
{


  viewer.Update();
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Escape);
  return 0;
}
