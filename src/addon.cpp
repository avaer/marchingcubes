#include <iostream>
#include <node.h>

#include "MarchingCubes.h"
#include "LookUpTable.h"
#include "csg.h"

using std::pow;
using std::abs;

float fn(float x, float y, float z) {
  //return pow((2*pow((1.3*x),2)+pow((1.3*y),2)+pow((1.3*z),2)-1),3)-(1/10)*pow((1.3*x),2)*pow((1.3*z),3)-pow((1.3*y),2)*pow((1.3*z),3);
  return pow(2*x,2)+pow(y,2)+pow(z,2);
}

// grid extension
float xmin=-1.0f, xmax=1.0f,  ymin=-1.0f, ymax=1.0f,  zmin=-1.0f, zmax=1.0f ;

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

  // begin mc
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
  for( i = 0 ; i < size_x ; i++ )
  {
    for( j = 0 ; j < size_y ; j++ )
    {
      for( k = 0 ; k < size_z ; k++ )
      {
        w = fn(
          (((float)(i) / (float)(size_x)) - ((float)(size_x) / 2)) * 20,
          (((float)(j) / (float)(size_y)) - ((float)(size_y) / 2)) * 20,
          (((float)(k) / (float)(size_z)) - ((float)(size_z) / 2)) * 20
        );
        w = (abs(i - (size_x / 2)) <= 4 && abs(j - (size_y / 2)) <= 4 && abs(k - (size_z / 2)) <= 4) ? -0.9 : 0.9;
        mc.set_data( w, i, j, k ) ;
      }
    }
  }

  // mc.set_method(true);
  mc.run();
  // end mc

  unsigned int numTrigs = mc.ntrigs();
  unsigned int numVerts = numTrigs * 3;
  unsigned int numPositions = numVerts * 3;
  std::cerr << "num positions: " << numTrigs << ":" << numPositions << "\n";
  v8::Local<v8::Float32Array> positions = v8::Float32Array::New(v8::ArrayBuffer::New(isolate, numPositions * 4), 0, numPositions);
  v8::Local<v8::Float32Array> normals = v8::Float32Array::New(v8::ArrayBuffer::New(isolate, numPositions * 4), 0, numPositions);

  auto triangles = mc.triangles();
  auto vertices = mc.vertices();
  for (unsigned int i = 0; i < numTrigs; i++) {
    const Triangle &triangle = triangles[i];
    auto a = vertices[triangle.v1];
    auto b = vertices[triangle.v2];
    auto c = vertices[triangle.v3];

    unsigned int baseIndex = i * 3 * 3;
    positions->Set(baseIndex + 0, v8::Number::New(isolate, a.x));
    positions->Set(baseIndex + 1, v8::Number::New(isolate, a.y));
    positions->Set(baseIndex + 2, v8::Number::New(isolate, a.z));
    positions->Set(baseIndex + 3, v8::Number::New(isolate, b.x));
    positions->Set(baseIndex + 4, v8::Number::New(isolate, b.x));
    positions->Set(baseIndex + 5, v8::Number::New(isolate, b.x));
    positions->Set(baseIndex + 6, v8::Number::New(isolate, c.x));
    positions->Set(baseIndex + 7, v8::Number::New(isolate, c.x));
    positions->Set(baseIndex + 8, v8::Number::New(isolate, c.x));

    normals->Set(baseIndex + 0, v8::Number::New(isolate, a.x));
    normals->Set(baseIndex + 1, v8::Number::New(isolate, a.y));
    normals->Set(baseIndex + 2, v8::Number::New(isolate, a.z));
    normals->Set(baseIndex + 3, v8::Number::New(isolate, b.x));
    normals->Set(baseIndex + 4, v8::Number::New(isolate, b.x));
    normals->Set(baseIndex + 5, v8::Number::New(isolate, b.x));
    normals->Set(baseIndex + 6, v8::Number::New(isolate, c.x));
    normals->Set(baseIndex + 7, v8::Number::New(isolate, c.x));
    normals->Set(baseIndex + 8, v8::Number::New(isolate, c.x));
  }

  args.GetReturnValue().Set(positions);
}

void Init(v8::Local<v8::Object> exports) {
  NODE_SET_METHOD(exports, "march", March);
}

NODE_MODULE(addon, Init)
