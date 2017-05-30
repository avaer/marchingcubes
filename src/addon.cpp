#include <node.h>

#include "MarchingCubes.h"
#include "LookUpTable.h"
#include "csg.h"

using std::pow;

float fn(int x, int y, int z) {
  return ( ( pow((8*x),2) + pow((8*y-2),2) + pow((8*z),2) + 16 - 1.85*1.85 ) * ( pow((8*x),2) + pow((8*y-2),2)+ pow((8*z),2) + 16 - 1.85*1.85 ) - 64 * ( pow((8*x),2) + pow((8*y-2),2) ) ) * ( ( pow((8*x),2) + ((8*y-2)+4)*((8*y-2)+4) + pow((8*z),2) + 16 - 1.85*1.85 ) * ( pow((8*x),2) + ((8*y-2)+4)*((8*y-2)+4) + pow((8*z),2) + 16 - 1.85*1.85 ) - 64 * ( ((8*y-2)+4)*((8*y-2)+4) + pow((8*z),2) ) ) + 1025;
}

// isovalue defining the isosurface
float isoval = 0.0f ;

// original/topological MC switch
int   originalMC = 0 ;

// grid extension
float xmin=-1.0f, xmax=1.0f,  ymin=-1.0f, ymax=1.0f,  zmin=-1.0f, zmax=1.0f ;

//_____________________________________________________________________________
// run the MC algorithm
void run()
//-----------------------------------------------------------------------------
{
  // grid size control
  int   size_x=50, size_y=50, size_z=50 ;

  // Init data
  MarchingCubes mc;
  mc.set_resolution( size_x, size_y, size_z ) ;
  mc.init_all() ;

  // Fills data structure
  int i,j,k ;
  float w ;
  float rx = (xmax-xmin) / (size_x - 1) ;
  float ry = (ymax-ymin) / (size_y - 1) ;
  float rz = (zmax-zmin) / (size_z - 1) ;
  unsigned char buf[sizeof(float)] ;
  for( i = 0 ; i < size_x ; i++ )
  {
    for( j = 0 ; j < size_y ; j++ )
    {
      for( k = 0 ; k < size_z ; k++ )
      {

        w = fn(i, j, k) - isoval ;
        mc.set_data( w, i,j,k ) ;
      }
    }
  }
  //if( export_iso ) mc.writeISO( out_filename->get_text() ) ;

  mc.run();
}

void March(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();

  // Check the number of arguments passed.
  if (args.Length() < 1) {
    // Throw an Error that is passed back to JavaScript
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong number of arguments")));
    return;
  }
  // Check the argument types
  if (!args[0]->IsObject()) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong arguments")));
    return;
  }

  v8::Local<v8::Object> opts = args[0]->ToObject();
  v8::Local<v8::Value> width = opts->Get(v8::String::NewFromUtf8(isolate, "width"));
  v8::Local<v8::Value> height = opts->Get(v8::String::NewFromUtf8(isolate, "height"));
  v8::Local<v8::Value> depth = opts->Get(v8::String::NewFromUtf8(isolate, "depth"));
  if (!(width->IsNumber() && height->IsNumber() && depth->IsNumber())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  /* MyObject* obj1 = node::ObjectWrap::Unwrap<MyObject>(
      args[0]->ToObject());
  MyObject* obj2 = node::ObjectWrap::Unwrap<MyObject>(
      args[1]->ToObject());

  double sum = obj1->value() + obj2->value(); */
  args.GetReturnValue().Set(v8::Number::New(isolate, 1));
}

void Init(v8::Local<v8::Object> exports) {
  NODE_SET_METHOD(exports, "march", March);
}

NODE_MODULE(addon, Init)
