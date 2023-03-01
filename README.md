All files in the ```Fortran/``` directory are written by Dr. Kevin Patrick Furlong. Compile the file ```tqtec.f90```  to ```tqtec``` and ```readtqtec.f90``` to ```tqtec``` in the ```Fortran/``` directory to use the scripts in ```Model/``` to generate batches of tectonic-thermal models. To generate models, navigate to ```DataGeneration/Python/``` and run the command 

```
python DataGenerator.py [-h] [-n Number of Samples] [-c Chunk Size < N] [-p Orogen Parameters] [-o Output Dir] [-t Temporary Dir]
``` 

