# time consts
MINUTE = 60
HOUR = 3600
DAY = 3600*24
MONTH = 30*DAY
YEAR = 365*DAY

# bins/streets are emptied/cleaned once each:
TIME_TO_EMPTY_BINS = 3*DAY
TIME_TO_CLEAN_STREETS = 3*DAY

# minimal percentage of a bin to get emptied by a truck:
MIN_PERCENTAGE_BIN = 0.50

# maximal days a bin can be ignored:
MAX_TIME_IGNORING_A_BIN = 7*DAY

# costs
COST_TRUCK_PER_KM = 15
COST_TRUCK_PER_DAY = 10
COST_OF_CLEANING_A_STREET = 10