# time consts
MINUTE = 60
HOUR = 3600
DAY = 3600*24
MONTH = 30*DAY
YEAR = 365*DAY

# bins/streets are emptied/cleaned once each:
PERIOD_EMPTY_BINS = 3*DAY
PERIOD_CLEAN_STREETS = 3*DAY

# minimal percentage of a bin to get emptied by a truck:
MIN_PERCENTAGE_BIN = 0.50

# maximal days a bin can be ignored:
MAX_TIME_IGNORING_A_BIN = 7*DAY

# costs
COST_TRUCK_PER_KM = 15
COST_TRUCK_PER_DAY = 10
COST_OF_CLEANING_A_STREET = 10

# epsilons
EPSILON_DIST = 1e-3
EPSILON_POT = 1e-20
EPSILON_CHARGE = 1e-10 
EPSILON_TEMPTR = 1e-10

# time to update potentials
PERIOD_UPDATE_POTS = 1*DAY


### FOR POTENTIALS ###
# division of streets for each order of potential (in meters)
SIZE_OF_DIVISION = [20, 2, 0.2]

# maximal range for each order of potential (in meters)
MAXIMAL_RANGE = [1e3, 1e2, 10]

# maximal order of potential
MAXIMAL_ORDER = 3

# weight of street sweeper opinion on reports
STR_SWEEPER_WEIGHT = 10

### FOR TEMPERATURE ###
# division of streets (in meters)
SIZE_OF_DIV_TEMPTR = 20

# distance to be considered in the extremity (in meters)
MAX_DIST_EXTREMITY = 3
if MAX_DIST_EXTREMITY > SIZE_OF_DIV_TEMPTR:
    raise ValueError("MAX_DIST_EXTREMITY > SIZE_OF_DIV_TEMPTR not allowed for now.")

# rho/A for streets
RHO_A = 10