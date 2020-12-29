""" Author: AbuBakr Ghaznavi
    Description: Finds hazards using GPS data
    Data: 11/29/2020
"""

# Import necessary files
import sys
from statistics import median
import numpy as np
import pandas as pd
import datetime
import re
# Define a window size to take the average of for turn classification
WINDOW_SIZE = 20
# Create regex to find multiples
duplicate_finder = re.compile("GP.{3}")

# Class for giving generic data access for the command sentences
class GPCMD():
	# Get the longitude of the data
	def long(self):
		coeff = 1
		if self.horizontal == "W":
			coeff = -1
		return coeff * longFrac(self.longitude)
	# Get the latitude of the data
	def lat(self):
		coeff = 1
		if self.vertical == "S":
			coeff = -1	
		return coeff * latFrac(self.latitude)
	# Get the speed of the data
	def spd(self):
		return float(self.speed)
	# Only using the timestamp as a point of reference
	def ts(self):
		dt = datetime.datetime.strptime(self.timestamp, "%H%M%S.%f")
		return dt
	# Get the altitude of the data
	def alt(self):
		return float(self.antenna_height_geoid)
	# Get the angle of the data
	def angle(self):
		return float(self.true_course)
	# Get the number of the satellites
	def num_sat(self):
		return int(self.num_satellites)

# Get the GPGGA data
class GPGGA(GPCMD):
	def __init__(self, line):
		pieces = line.split(",")
		self.code = pieces[0].strip()
		self.timestamp = pieces[1].strip()
		self.latitude = pieces[2].strip()
		self.vertical = pieces[3].strip()
		self.longitude = pieces[4].strip()
		self.horizontal = pieces[5].strip()
		self.quality_indicator = pieces[6].strip()
		self.num_satellites = pieces[7].strip()
		self.hori_dilution = pieces[8].strip()
		self.antenna_height_geoid = pieces[9].strip()
		self.antenna_height_meters = pieces[10].strip()
		self.geoidal_seperation = pieces[11].strip()
		self.unit_seperation = pieces[12].strip()
		self.update_age = pieces[13].strip()
		self.diff_refr_id = pieces[14].strip()

	def __repr__(self):
		gps_map = {"Code": self.code,
			   "Timestamp": self.timestamp,
 			   "Latitude": self.latitude,
			   "Vertical": self.vertical,
			   "Longitude": self.longitude,
			   "Horizontal": self.horizontal,
			   "Number of satellites": self.num_satellites,
			   "Horizontal dillution": self.hori_dilution,
			   "Antenna height geoid": self.antenna_height_geoid,
			   "Antenna height meters": self.antenna_height_meters,
			   "Geoidal seperation unit":self.geoidal_seperation,
			   "Meters seperation": self.unit_seperation,
			   "Update age":self.update_age,
			   "Difference reference ID":self.diff_refr_id}
		string_format = ""
		for key in sorted(gps_map.keys()):
			string_format += key + ":" +  gps_map[key]
			string_format += "\n"	
		return string_format

# Get the Recommended minimum requirements
class GPRMC(GPCMD):
	def __init__(self, line):
		pieces = line.split(",")
		self.code = pieces[0].strip()
		self.timestamp = pieces[1].strip()
		self.validity = pieces[2].strip()
		self.latitude = pieces[3].strip()
		self.vertical = pieces[4].strip()
		self.longitude = pieces[5].strip()
		self.horizontal = pieces[6].strip()
		# Speed in knots
		self.speed = pieces[7].strip()
		self.true_course = pieces[8].strip()
		self.datestamp = pieces[9].strip()
		self.variation = None
		self.variation_horizontal = None
		if len(pieces) >  10:
			self.variation = pieces[10].strip()
		if len(pieces) > 11:
			self.variation_horizontal = pieces[11].strip()
	def __repr__(self):
		gps_map = {"Timestamp": self.timestamp,
			   "Validity": self.validity,	
			   "Latitude": self.latitude,	
			   "Vertical": self.vertical,	
			   "Longitude": self.longitude,	
			   "Horizontal": self.horizontal,	
			   "Speed": self.speed,	
			   "True course": self.true_course,	
			   "Datestamp": self.true_course,
			   "Variation": self.variation,
			   "Variation horizontal": self.variation_horizontal,
			   "Code": self.code}
		          
		string_format = ""
		for key in sorted(gps_map.keys()):
			val = gps_map[key]
			if val == None:
				str_val = "None"
			else:
				str_val = str(val)
			string_format += key + ":" + str_val
			string_format += "\n"	
		return string_format

# Define a method for reading the GPS recorder
class GPSReader():
	sentence_codes = ["$GPRMC","$GPGGA"]
	# Split file header to read fields
	def __init__(self, filename):
		self.filehandle = open(filename, "r")
		self.version = self.filehandle.readline().replace("Vers","").strip()
		self.use_serial_feedback = self.filehandle.readline().split("=")[1].strip().lower() == "true"
		self.development_mode = self.filehandle.readline().split("=")[1].strip().lower() == "true"
		self.rmc_only = self.filehandle.readline().split("=")[1].strip().lower() == "true"
		# Read the blank line to align with the rest of the file
		self.filehandle.readline()
	
	# Get the code from the line
	@staticmethod
	def getCode(line):
		comma_index = line.find(",")
		return line[0:comma_index]
	
	# Get the command sentence
	def getGPCmd(self):
		code = ""
		while code not in GPSReader.sentence_codes:
			
			line = self.filehandle.readline()
			# Skip over lines with double sentences
			if len(re.findall(duplicate_finder, line)) > 1:
				continue
			# If we reached EOF
			if len(line) == 0:
				return None
			code = GPSReader.getCode(line).strip()
				
		# Throw error if code is not GPS code
		if not(code in GPSReader.sentence_codes):
			raise Exception("Invalid GP code")
		if code == "$GPGGA":
			return GPGGA(line)
		elif code == "$GPRMC":
			return GPRMC(line)

	# Define an iterator to reading the lines as commands
	def iterCmds(self):	
		while True:
			GPCmd = self.getGPCmd()
			if GPCmd == None:
				self.filehandle.close()
				break
			yield GPCmd	

# Function to convert to fractional degrees
def longFrac(long_str):
	# Get first three digits for degree part
	deg_part = long_str[0:3]
	minute_part = long_str[3:]
	return float(deg_part) + float(minute_part)/60
	
# Convert the latitude to fractional	
def latFrac(long_str):
	# Get first three digits for degree part
	deg_part = long_str[0:2]
	minute_part = long_str[2:]
	return float(deg_part) + float(minute_part)/60
	
# Define a placemarker class to write the data	
class KMLPlacemarker():
	# Define headers and formats to use
	KML_HEADER = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
 """

	KML_TRAILER = """</Document>
</kml>"""
	KML_PlacemarkFMT = """
  <Placemark>
    <name>%s</name>
    <description>%s</description>
    <Style id="%s">
      <IconStyle>
        <color> %s </color>
        <Icon>
           <href>http://maps.google.com/mapfiles/kml/paddle/1.png</href>
        </Icon>
      </IconStyle>
    </Style>
    <Point>
      <coordinates>%s</coordinates>
    </Point>
  </Placemark>
	"""
	def __init__(self, filename):
		self.filename = filename
		self.filehandle = open(self.filename, "w")
		# Write the header to the file
		self.filehandle.write(KMLPlacemarker.KML_HEADER + "\n")
		
	# Write a placemark to the data
	def writePlacemark(self, color, values):
		# Define the colors and styles for each hazard
		if color.lower() == "red":
			name = "Stop pin"
			description = "Red PIN for a stop"
			hex_color = "ff0000ff"
			style_id = "redStyle"
		elif color.lower() == "cyan":
			name = "Right turn"
			description = "Cyan PIN for right turn"
			hex_color = "ffffff00"
			style_id = "cyanStyle"
		elif color.lower() == "yellow":
			name = "Left turn"
			description = "Yellow PIN for left turn"
			hex_color = "ff00ffff"
			style_id = "yellowStyle"
		else:
			name = "Placemark"
			hex_color = "ff00A5ff"
			description = "Generic placemark"
			style_id = "genericStyle"
		# Make the coordinate string that will be used for the placemark
		values = [str(coordinate) for coordinate in values]
		coordinate_string = ",".join(values)
		placemark_string = KMLPlacemarker.KML_PlacemarkFMT % (name, description, style_id, hex_color, coordinate_string)
		self.filehandle.write(placemark_string + "\n") 	
	# Write the trailer	
	def endFile(self):
		self.filehandle.write(KMLPlacemarker.KML_TRAILER)	
	# Close the file handle
	def close(self):
		self.filehandle.close()

# Read in the arguments and define arguments
if len(sys.argv) > 2:
	TEST_FILENAMES = sys.argv[1:-1]
	print(TEST_FILENAMES)
	OUT_FILENAME = sys.argv[-1]
else:
	print("Usage GPStoKML.py <gps_file1> <gps_file2> ...  <outfile>.kml")
	sys.exit(1)

SPEED_TOLERANCE = 0.04
ANGLE_SPEED_TOLERANCE = 0.01
# Have max speed tolerance to exclude lane changes
MAX_SPEED_TOLERANCE = 40


def angle_dist(a, b):
	d = a - b
	d = ( d + 180 ) % 360 - 180
	return d

def calcAngleDiffs(df):
	# Keep an initial zero for the sake of indexing
	diffs = [0]
	num_rows = len(df.index)
	for i in range(num_rows - 1):
		row = df.iloc[i]
		next_row = df.iloc[i + 1]
		first_angle = row["angle"]
		second_angle = next_row["angle"]
		# Calculate the angle delta
		angle_difference = angle_dist(second_angle, first_angle)
		# Apply the speed as a factor (average between them)
		# Get the average speed
		avg_speed = (row["speed"] + next_row["speed"]) / 2.0
		# Exclude angles when there is no speed factor
		if avg_speed < ANGLE_SPEED_TOLERANCE:
			angle_difference = 0
		if avg_speed > MAX_SPEED_TOLERANCE:
			angle_difference = 0		
		diffs.append(angle_difference)
	return np.array(diffs)

# Find the range from the end for which there is no movement at the end
def find_stop_range(df):
	# Go to the end of the data frame and check the last time the speed is under the threshold
	last_index = len(df.index) - 1
	cur_idx = last_index
	for index, row in df[::-1].iterrows():
		current_speed = row["speed"]
		if current_speed > SPEED_TOLERANCE:
			return (last_index, cur_idx)
		cur_idx -= 1
	return (start_index, cur_idx + 1)

# Find the range for which there is no motion at the start
def find_start_range(df):
	start_index = 0
	cur_idx = start_index
	for index, row in df[::-1].iterrows():
		current_speed = row["speed"]
		if current_speed < SPEED_TOLERANCE:
			return (start_index, cur_idx)
		cur_idx += 1	
	# Leave one value in the range to represent
	return (start_index, cur_idx - 1)
	
def agglo_turns(df):
	# Use on the cross index of the dataframe and a sliding window approach for smoothing
	# Adjustable parameter for this function
	WIN_TOLERANCE = 2
	in_turn = False 
	window_size = WINDOW_SIZE
	middle_indices_arr = []
	num_rows = len(df.index)
	medians = []
	for idx in range(num_rows - window_size):
		middle_idx = int((idx + window_size) // 2)
		window = []
		for window_cell in range(window_size):
			cell_cross = df.iloc[idx + window_cell]["cross"]
			window.append(cell_cross)
		window_mean = np.mean(np.array(window))
		# Mean high enough to enter turn
		high_mean = window_mean > WIN_TOLERANCE
		# If mean was comprised of 1 high value surrouned by negatives then it will not count
		count_pos = len([num for num in window if num > 0])
		count_prop = count_pos / window_size
		if count_prop <= 0.5:
			high_mean = False	
		if not(in_turn) and high_mean:
			middle_indices_arr.append(middle_idx)
			in_turn = True
		if not(in_turn) and not(high_mean):
			# Do nothing if we aren't in high mean and not gonna be
			continue
		if in_turn and high_mean:
			# Append to our middle_indices
			middle_indices_arr.append(middle_idx)
		# Tough state where we handle transition to not in_turn
		# Save the median of our indices
		if in_turn and not(high_mean):
			in_turn = False
			medians.append(int(median(middle_indices_arr)))
			middle_indices_arr = []
	return medians		
		

def r_agglo_turns(df):
	# Use on the cross index of the dataframe and a sliding window approach for smoothing
	# Adjustable parameter for this function
	WIN_TOLERANCE = -2
	in_turn = False 
	window_size = WINDOW_SIZE
	middle_indices_arr = []
	num_rows = len(df.index)
	medians = []
	for idx in range(num_rows - window_size):
		middle_idx = int((idx + window_size) // 2)
		window = []
		for window_cell in range(window_size):
			cell_cross = df.iloc[idx + window_cell]["cross"]
			window.append(cell_cross)
		window_mean = np.mean(np.array(window))
		# Mean high enough to enter turn
		high_mean = window_mean < WIN_TOLERANCE
		# If mean was comprised of 1 high value surrouned by negatives then it will not count
		count_neg = len([num for num in window if num < 0])
		count_prop = count_neg / window_size
		if count_prop <= 0.5:
			high_mean = False	
		if not(in_turn) and high_mean:
			middle_indices_arr.append(middle_idx)
			in_turn = True
		if not(in_turn) and not(high_mean):
			# Do nothing if we aren't in high mean and not gonna be
			continue
		if in_turn and high_mean:
			# Append to our middle_indices
			middle_indices_arr.append(middle_idx)
		# Tough state where we handle transition to not in_turn
		# Save the median of our indices
		if in_turn and not(high_mean):
			in_turn = False
			medians.append(int(median(middle_indices_arr)))
			middle_indices_arr = []
	return medians		

def detect_stops(df):
	# Use on the cross index of the dataframe and a sliding window approach for smoothing
	# Adjustable parameter for this function
	STOP_TOL = 12
	stopped = False
	SUFFICIENT_STOP_COUNT = 40
	SUFFFICIENT_GO_COUNT = 4
	TOO_MANY_COUNT = 400
	stopped_idx = []
	current_stopped_idx = None
	stop_count = 0 
	go_count = 0
	# Need to see 20 slows in a row before entering stopped state
	num_idx = 0
	for index, row in df.iterrows():
		speed = row["speed"]
		if stopped == False:
			if speed < STOP_TOL:
				stop_count += 1
			else:
				stop_count = 0
			if stop_count == SUFFICIENT_STOP_COUNT:
				stopped = True
				current_stopped_idx = num_idx
		if stop_count == TOO_MANY_COUNT:
			current_stopped_idx = None
		if stopped == True:
			if speed >= STOP_TOL:
				go_count += 1
			else:
				stop_count += 1
				go_count = 0
			if go_count == SUFFFICIENT_GO_COUNT:
				stopped = False
				go_count = 0
				stop_count = 0
				if current_stopped_idx != None:
					stopped_idx.append(current_stopped_idx)
		num_idx += 1
	return stopped_idx	

def process_gps(gps_file):
	reader = GPSReader(gps_file)
	data_mat = []
	# Get all the latitudes and longitudes
	cmds = [item for item in reader.iterCmds()]

	# Pair up sentences of seperate codes
	for i in range(len(cmds) - 1):
		if i % 2 == 0:
			continue
		# Means that we have a pair of the same code type
		first_code = cmds[i].code
		second_code = cmds[i + 1].code
		if first_code == second_code:
			continue
		

		first_code = cmds[i].code
		second_code = cmds[i + 1].code

		

		if first_code == "$GPGGA":
			gpgga = cmds[i]
			gprmc = cmds[i + 1]
		elif first_code == "$GPRMC":
			gpgga = cmds[i + 1]
			gprmc = cmds[i]

		# Check for low fixquality or invalid 
		if gprmc.validity == "V" or gpgga.quality_indicator == "0":
			continue


		# Flatten out the relevant data fields to be used for pandas array
		# Most relevant fields will come from the GPRMC data
		data_row = [gprmc.ts(), gpgga.long(), gpgga.lat(), gpgga.alt(), gprmc.spd(), gprmc.angle(), gpgga.num_sat()]
		data_mat.append(data_row)

	# Make the data frame from the data
	column_list = ["time","longitude","latitude","altitude","speed","angle","num_sat"]
	gps_data = pd.DataFrame(data_mat, columns=column_list)	
	angle_dif = calcAngleDiffs(gps_data)
	gps_data["cross"] = 5 * (1 + np.square(angle_dif)) * np.sign(angle_dif) 
	#with pd.option_context('display.max_rows', None, 'display.max_columns', None):
	#	print(gps_data[["speed","angle","cross"]])

	stop_range = find_stop_range(gps_data)
	start_range = find_start_range(gps_data)
	print(stop_range)
	print(start_range)
	# Drop non moving rows on both ends of the track
	gps_data.drop(gps_data.index[start_range[0]:start_range[1]], inplace = True)
	gps_data.drop(gps_data.index[stop_range[1]:stop_range[0]], inplace = True)
	right_turns = r_agglo_turns(gps_data)
	left_turns = agglo_turns(gps_data)
	print(right_turns)
	print(left_turns)
	stopped_idx = detect_stops(gps_data)
	gps_data["zero"] = np.zeros(len(gps_data.index))


	# Write in each of the hazards
	for turn_idx in right_turns:
		row = gps_data.iloc[turn_idx]
		row_values = [row["longitude"], row["latitude"], row["speed"]]
		placemarker.writePlacemark("cyan", row_values)
	for turn_idx in left_turns:
		row = gps_data.iloc[turn_idx]
		row_values = [row["longitude"], row["latitude"], row["speed"]]
		placemarker.writePlacemark("yellow", row_values)
	for stop_idx in stopped_idx:
		row = gps_data.iloc[stop_idx]
		row_values = [row["longitude"], row["latitude"], row["speed"]]
		placemarker.writePlacemark("red", row_values)

placemarker = KMLPlacemarker(OUT_FILENAME)
# Iterate through files and process them
for filename in TEST_FILENAMES:
	print("Working on processing file:", filename)
	process_gps(filename)

placemarker.endFile()
placemarker.close()
