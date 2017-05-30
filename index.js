const marchingcubes = require('./build/Release/addon');

console.log('got marching cubes', marchingcubes.march({
  width: 50,
  height: 50,
  depth: 50,
}));
