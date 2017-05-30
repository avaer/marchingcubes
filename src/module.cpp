#include <node.h>

#include "MarchingCubes.h"
#include "LookUpTable.h"

//_____________________________________________________________________________
// run the MC algorithm
void run()
//-----------------------------------------------------------------------------
{
	strcpy(formula, fun_def[9]);

  if( strlen(formula) <= 0 ) return false ;
  //if( export_iso && strlen( out_filename->get_text() ) <= 0 ) export_iso = 0 ;

  // Init data
  mc.set_resolution( size_x, size_y, size_z ) ;
  mc.init_all() ;

  // Parse formula
  FunctionParser fparser ;
  fparser.Parse( (const char*)formula, "x,y,z,c,i" ) ;
  if( fparser.EvalError() )
  {
    printf( "parse error\n" ) ;
    return false ;
  }

  // Fills data structure
  int i,j,k ;
  float val[5], w ;
  float rx = (xmax-xmin) / (size_x - 1) ;
  float ry = (ymax-ymin) / (size_y - 1) ;
  float rz = (zmax-zmin) / (size_z - 1) ;
  unsigned char buf[sizeof(float)] ;
  for( i = 0 ; i < size_x ; i++ )
  {
    val[X] = (float)i * rx  + xmin ;
    for( j = 0 ; j < size_y ; j++ )
    {
      val[Y] = (float)j * ry  + ymin ;
      for( k = 0 ; k < size_z ; k++ )
      {
        val[Z] = (float)k * rz  + zmin ;

        if( csg_root )
        {
          val[3] = csg_root->eval( val[X],val[Y],val[Z] ) ;
        }
        if( isofile  )
        {
          fread (buf, sizeof(float), 1, isofile);
          val[4] = * (float*) buf ;
        }

        w = fparser.Eval(val) - isoval ;
        mc.set_data( w, i,j,k ) ;
      }
    }
  }
  //if( export_iso ) mc.writeISO( out_filename->get_text() ) ;

  mc.run();
}

void Init(v8::Local<v8::Object> exports) {
  NODE_SET_METHOD(exports, "add", Add);
}

NODE_MODULE(addon, Init)
