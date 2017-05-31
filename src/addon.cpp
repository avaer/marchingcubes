#include <node.h>

#include "MarchingCubes.h"
#include "LookUpTable.h"
#include "csg.h"

// grid extension
float xmin=-1.0f, xmax=1.0f,  ymin=-1.0f, ymax=1.0f,  zmin=-1.0f, zmax=1.0f ;

v8::Local<v8::Value> DoMarchCubes(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto normalsKey = v8::String::NewFromUtf8(isolate, "normals");

  // Check the number of arguments passed.
  if (args.Length() < 1) {
    // Throw an Error that is passed back to JavaScript
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong number of arguments")));
    return v8::Null(isolate);
  }
  // Check the argument types
  if (!args[0]->IsObject()) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong arguments")));
    return v8::Null(isolate);
  }

  v8::Local<v8::Object> opts = args[0]->ToObject();
  v8::Local<v8::Value> width = opts->Get(v8::String::NewFromUtf8(isolate, "width"));
  v8::Local<v8::Value> height = opts->Get(v8::String::NewFromUtf8(isolate, "height"));
  v8::Local<v8::Value> depth = opts->Get(v8::String::NewFromUtf8(isolate, "depth"));
  v8::Local<v8::Value> data = opts->Get(v8::String::NewFromUtf8(isolate, "data"));
  if (!(width->IsNumber() && height->IsNumber() && depth->IsNumber() && data->IsFloat32Array())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return v8::Null(isolate);
  }

  v8::Local<v8::Float32Array> dataArray = data.As<v8::Float32Array>();

  // begin mc
  // grid size control
  int size_x = width->ToInt32()->Value(), size_y = height->ToInt32()->Value(), size_z = depth->ToInt32()->Value();
  int half_size_x = size_x / 2, half_size_y = size_y / 2, half_size_z = size_z / 2;

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
  unsigned int numNormals = numVerts * 3;
  v8::Local<v8::Float32Array> positions = v8::Float32Array::New(v8::ArrayBuffer::New(isolate, numPositions * 4), 0, numPositions);
  v8::Local<v8::Float32Array> normals = v8::Float32Array::New(v8::ArrayBuffer::New(isolate, numNormals * 4), 0, numNormals);

  auto triangles = mc.triangles();
  auto vertices = mc.vertices();
  for (unsigned int i = 0; i < numTrigs; i++) {
    const Triangle &triangle = triangles[i];
    const Vertex &a = vertices[triangle.v1];
    const Vertex &b = vertices[triangle.v2];
    const Vertex &c = vertices[triangle.v3];

    unsigned int baseIndex = i * 3 * 3;
    positions->Set(baseIndex + 0, v8::Number::New(isolate, a.x - half_size_x));
    positions->Set(baseIndex + 1, v8::Number::New(isolate, a.y - half_size_y));
    positions->Set(baseIndex + 2, v8::Number::New(isolate, a.z - half_size_z));
    positions->Set(baseIndex + 3, v8::Number::New(isolate, b.x - half_size_x));
    positions->Set(baseIndex + 4, v8::Number::New(isolate, b.y - half_size_y));
    positions->Set(baseIndex + 5, v8::Number::New(isolate, b.z - half_size_z));
    positions->Set(baseIndex + 6, v8::Number::New(isolate, c.x - half_size_x));
    positions->Set(baseIndex + 7, v8::Number::New(isolate, c.y - half_size_y));
    positions->Set(baseIndex + 8, v8::Number::New(isolate, c.z - half_size_z));

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
  result->Set(positionsKey, positions);
  result->Set(normalsKey, normals);

  return result;
}

void MarchCubes(const v8::FunctionCallbackInfo<v8::Value>& args) {
  args.GetReturnValue().Set(DoMarchCubes(args));
}
void MarchCubesPlanet(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto normalsKey = v8::String::NewFromUtf8(isolate, "normals");

  v8::Local<v8::Object> marchingCubes = DoMarchCubes(args).As<v8::Object>();
  v8::Local<v8::Value> positions = marchingCubes->Get(positionsKey).As<v8::Float32Array>();
  v8::Local<v8::Value> normals = marchingCubes->Get(normalsKey).As<v8::Float32Array>();

  args.GetReturnValue().Set(marchingCubes);
}

void Init(v8::Local<v8::Object> exports) {
  NODE_SET_METHOD(exports, "marchCubes", MarchCubes);
  NODE_SET_METHOD(exports, "marchCubesPlanet", MarchCubesPlanet);
}

NODE_MODULE(addon, Init)
