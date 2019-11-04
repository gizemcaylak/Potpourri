current_path=cd;
cd([current_path, '/src/cpp/SPADIS']);
mex('-g', '-O', '-largeArrayDims', '-lut', '-output', [current_path, '/binaries/spadis'], '*.cpp');
cd([current_path]);

current_path=cd;
cd([current_path, '/src/cpp/ET']);
mex('-g', '-O', '-largeArrayDims', '-lut', '-output', [current_path, '/binaries/et'], '*.cpp');
cd([current_path]);





