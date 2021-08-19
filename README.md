Navigate to the folder containing this file and clone morphologica:

    git clone https://github.com/ABRG-Models/morphologica.git

Then:

    mkdir build
    cd build
    cmake ..
    make
    cd ..
    ./build/model config.json logs 1
    
This should display the dynamics and produce logs/data.h5, the contents of which you can display in python by running:

    python analysis.py

