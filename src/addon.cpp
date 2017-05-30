#include <node.h>

#include "MarchingCubes.h"
#include "LookUpTable.h"
#include "csg.h"

using std::pow;
using std::abs;
using std::sqrt;

float fn(float x, float y, float z) {
  //return pow((2*pow((1.3*x),2)+pow((1.3*y),2)+pow((1.3*z),2)-1),3)-(1/10)*pow((1.3*x),2)*pow((1.3*z),3)-pow((1.3*y),2)*pow((1.3*z),3);
  return sqrt(pow(x*2,2)+pow(y,2)+pow(z,2)) - 0.5f;
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
  v8::Local<v8::Value> data = opts->Get(v8::String::NewFromUtf8(isolate, "data"));
  if (!(width->IsNumber() && height->IsNumber() && depth->IsNumber() && data->IsFloat32Array())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  v8::Local<v8::Float32Array> dataArray = data.As<v8::Float32Array>();

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
  unsigned int index;
  float rx = (xmax-xmin) / (size_x - 1) ;
  float ry = (ymax-ymin) / (size_y - 1) ;
  float rz = (zmax-zmin) / (size_z - 1) ;
  for( i = 0 ; i < size_x ; i++ )
  {
    for( j = 0 ; j < size_y ; j++ )
    {
      for( k = 0 ; k < size_z ; k++ )
      {
        index = i + (j * size_x) + (k * size_x * size_y);
        w = (float)dataArray->Get(index)->NumberValue();
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
    positions->Set(baseIndex + 4, v8::Number::New(isolate, b.y));
    positions->Set(baseIndex + 5, v8::Number::New(isolate, b.z));
    positions->Set(baseIndex + 6, v8::Number::New(isolate, c.x));
    positions->Set(baseIndex + 7, v8::Number::New(isolate, c.y));
    positions->Set(baseIndex + 8, v8::Number::New(isolate, c.z));

    normals->Set(baseIndex + 0, v8::Number::New(isolate, a.nx));
    normals->Set(baseIndex + 1, v8::Number::New(isolate, a.ny));
    normals->Set(baseIndex + 2, v8::Number::New(isolate, a.nz));
    normals->Set(baseIndex + 3, v8::Number::New(isolate, b.nx));
    normals->Set(baseIndex + 4, v8::Number::New(isolate, b.ny));
    normals->Set(baseIndex + 5, v8::Number::New(isolate, b.nz));
    normals->Set(baseIndex + 6, v8::Number::New(isolate, c.nx));
    normals->Set(baseIndex + 7, v8::Number::New(isolate, c.ny));
    normals->Set(baseIndex + 8, v8::Number::New(isolate, c.nz));
  }

  v8::Local<v8::Object> result = v8::Object::New(isolate);
  result->Set(v8::String::NewFromUtf8(isolate, "positions"), positions);
  result->Set(v8::String::NewFromUtf8(isolate, "normals"), normals);

  args.GetReturnValue().Set(result);
}

void Init(v8::Local<v8::Object> exports) {
  NODE_SET_METHOD(exports, "march", March);
}

NODE_MODULE(addon, Init)
