# IoT-SmartBins

# packages to install - you can use pip install:
1. matplotlib
2. pygame
3. numpy
4. rand

# instructions:
If you want to see the animation, run the file "main_ANIMATE.py"
If you want to run the simulation and see some outputs, run the file "main_SIMUL.py"

If you are using Windows:
- make sure that you have python3 installed already. If you don't have it, you will have to install it.
- open a prompt command window (cmd) and type "py main_ANIMATE.py" or "py main_SIMUL.py".

If you are using Linux:
- make sure you have python3 installed already. If you don't have it, install python3 package: "sudo apt install python3"
- in a terminal, type "python3 main_ANIMATE.py" or "python3 main_SIMUL.py".

# animation:

# simulation:
You can change some parameters from simulation very easily.

At "main_SIMUL.py" file:
1. TIME_STEP: time of each iteraction. It's the time between everything is updated in the simulation.
2. TIME_OF_SIMULATION: total time of a simulation.
3. NUMBER_OF_SIMULATION: number of simulations using a single map; and the same commertial points and the same bin's initial position.

At "consts.py" file:
1. TIME_TO_EMPTY_BINS: time for passing the trucks, the period.
2. TIME_TO_CLEAN_STREETS: time for street sweepers working, the period.
3. MIN_PERCENTAGE_BIN: minimal filling of the bin for which the trucks will empty it.
4. MAX_TIME_IGNORING_A_BIN: maximal time that the trucks can ignore a certain bin.
5. COST_TRUCK_PER_KM: cost for each km traveled by the trucks.
6. COST_TRUCK_PER_DAY: cost for each time the trucks run through the city - it doesn't depend on the distance traveled.
7. COST_OF_CLEANING_A_STREET: cost for removing all the trash from a street.
8. EPSILON: if the distance between 2 points is less than EPSILON, the program will see them as the same.
